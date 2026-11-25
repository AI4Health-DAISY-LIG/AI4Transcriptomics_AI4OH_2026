"""Step 2: quality control and exploration of the raw counts.

Outputs:
  results/step2_qc/
    library_sizes.png                   total counts per sample, coloured by condition
    read_counts_summary.csv             summary statistics of sequencing depth by group
    read_counts_by_<variable>.png       read count scatter and boxplots by batch, site, and condition
    gene_filter.txt                     genes before and after filtering, with the rule used

  results/step2_qc/{RLOG2_PCA, VST_PCA, INLIERS_VST_PCA}/
    pca_variance.csv                    variance explained by the first 10 components
    pca_scores.csv                      PC1 to PC10 per sample
    pca_by_<variable>.png               PC1 vs PC2 coloured by condition and metadata variables
    PCA_top_variables.csv               R² / η² variance explained and FDR per variable on PC1 & PC2
    PC1_PC2_top_variables.png           horizontal bar chart of top variables driving PC1 & PC2

  data/processed/
    counts_filtered.csv                 filtered raw counts (genes x samples), input of step 3
    metadata.csv                        one row per sample column of the counts, with condition

Decisions marked TODO(student) are yours. Record them in worksheets/step2_qc.md.

Run:  uv run python scripts/step2_qc_exploration.py
"""

from __future__ import annotations

import matplotlib
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

from infer_sex import verify_sample_sex
from helpers.statistics import benjamini_hochberg

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402
from step2_view_genes_filtering_effects import (
    display_genes_and_sample_size_for_biomarkers,
)
from norm_before_pca import (
    python_pydeseq2_vst,
    deseq2_rlog_exact_python,
)
from adjustText import adjust_text

import config  # noqa: E402

OUT = config.RESULTS / "step2_qc"
AUDIT = config.RESULTS / "step1_audit"

OUT_RLOG2_PCA = config.RESULTS / "step2_qc" / "RLOG2_PCA"
OUT_VST_PCA = config.RESULTS / "step2_qc" / "VST_PCA"
OUT_INLIERS_VST_PCA = config.RESULTS / "step2_qc" / "INLIERS_VST_PCA"


for path in [OUT_RLOG2_PCA, OUT_VST_PCA, OUT_INLIERS_VST_PCA]:
    path.mkdir(parents=True, exist_ok=True)


def load_counts(with_mapping=True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Raw counts (genes x samples) and metadata aligned on the sample columns."""
    counts = pd.read_csv(config.RAW / config.COUNTS_FILE, index_col=0)

    if with_mapping:
        samples = pd.read_csv(AUDIT / "sample_table.csv")
    else:
        samples = pd.read_csv(AUDIT / "sample_table.csv", index_col=0)

    # this function is used in the step1 to infer gender so we skip the rest if it's already aligned
    if "char:infer_gender" in samples.columns:
        return counts, samples

    # TODO(student): read results/step1_audit/counts_header_check.txt first. This join keeps
    # only the count columns whose name equals a GEO sample_id. If a column has no match,
    # find out which biopsy it is and decide: rename it, or leave it out and say why.

    # rename sample_id 01-0022 to 01-0022-1 to match the geo record
    samples.loc[samples["gsm"] == "GSM3186680", "char:sample_id"] = "01-0022-1"

    key = "char:sample_id" if "char:sample_id" in samples.columns else "title"
    meta = samples.set_index(key).reindex(counts.columns)

    missing = meta["condition"].isna()
    if missing.any():
        print(
            f"[warn] {missing.sum()} count column(s) dropped, no GEO record under that name: "
            f"{list(counts.columns[missing])}"
        )
    return counts, meta


variables = (
    ["condition"]
    + config.TECHNICAL_VARS
    + config.BIOLOGICAL_VARS
    + config.PATHOLOGY_VARS
    + config.INFER_VARS
)


def _sig_stars(fdr: float) -> str:
    if pd.isna(fdr) or fdr >= 0.05:
        return ""
    if fdr < 0.001:
        return "***"
    if fdr < 0.01:
        return "**"
    return "*"


def plot_top_pc_drivers(
    pcs: pd.DataFrame,
    metadata: pd.DataFrame,
    ratio,
    path,
    variables: list | set | None = None,
    ntop: int = 20,
    title: str = "Top Variables Separating Samples on PC1 & PC2",
) -> pd.DataFrame:
    """
    Ranks variables in `metadata` by how well they separate samples along PC1 and PC2
    (using R² for numeric and η² for categorical variables), plots the top `ntop`,
    and exports the full results table (with p-values and FDR) to a CSV file.
    """
    records = []

    for col in metadata.columns:
        if variables is not None and col not in variables:
            continue

        s = metadata[col]
        valid = s.notna()
        if valid.sum() < 3:
            continue

        s_valid = s[valid]
        numeric = pd.to_numeric(s_valid, errors="coerce")
        is_numeric = numeric.notna().all() and numeric.nunique() > 8

        row = {
            "variable": col.replace("char:", ""),
            "original_column": col,
            "type": "numeric" if is_numeric else "categorical",
            "metric": "R2" if is_numeric else "eta2",
            "n_samples": int(valid.sum()),
        }

        for pc in ["PC1", "PC2"]:
            y = pcs.loc[valid, pc].to_numpy(dtype=float)
            ss_total = np.sum((y - y.mean()) ** 2)

            if ss_total == 0:
                row[pc] = 0.0
                row[f"{pc}_pvalue"] = np.nan
                continue

            if is_numeric:
                x = numeric.to_numpy(dtype=float)
                if np.std(x) == 0:
                    row[pc] = 0.0
                    row[f"{pc}_pvalue"] = np.nan
                else:
                    r, p_val = stats.pearsonr(x, y)
                    row[pc] = float(r**2)
                    row[f"{pc}_pvalue"] = float(p_val)
            else:
                groups = s_valid.astype(str)
                if groups.nunique() < 2 or groups.nunique() >= len(groups):
                    row[pc] = 0.0
                    row[f"{pc}_pvalue"] = np.nan
                else:
                    df_tmp = pd.DataFrame({"y": y, "g": groups})
                    group_means = df_tmp.groupby("g")["y"].transform("mean")
                    ss_between = np.sum((group_means - y.mean()) ** 2)
                    row[pc] = float(ss_between / ss_total)

                    # One-way ANOVA p-value corresponding to eta-squared
                    grouped_y = [grp["y"].to_numpy() for _, grp in df_tmp.groupby("g")]
                    _, p_val = stats.f_oneway(*grouped_y)
                    row[f"{pc}_pvalue"] = float(p_val)

        records.append(row)

    scores = pd.DataFrame(records)

    if not scores.empty:
        # Add Benjamini-Hochberg FDR-adjusted p-values
        scores["PC1_fdr"] = benjamini_hochberg(scores["PC1_pvalue"])
        scores["PC2_fdr"] = benjamini_hochberg(scores["PC2_pvalue"])

        # Add combined variance explained across PC1 & PC2 for overall ranking
        scores["combined_score"] = scores["PC1"] + scores["PC2"]
        scores = scores.sort_values("combined_score", ascending=False).reset_index(
            drop=True
        )

    # Export to CSV
    scores.to_csv(path / "PCA_top_variables.csv", index=False)

    # --- Plotting ---
    fig, axes = plt.subplots(1, 2, figsize=(11, max(4, ntop * 0.4)), sharey=False)
    palette = {"numeric": "#2b5c8f", "categorical": "#d95f02"}

    for ax, pc, var_ratio in zip(axes, ["PC1", "PC2"], ratio[:2]):
        top_df = scores.nlargest(ntop, pc).iloc[::-1]
        colors = [palette[t] for t in top_df["type"]]

        bars = ax.barh(
            top_df["variable"], top_df[pc], color=colors, edgecolor="none", height=0.7
        )

        ax.set_xlabel("Variance Explained ($R^2$ / $\\eta^2$)", fontsize=9)
        ax.set_title(f"{pc} ({var_ratio:.0%})", fontsize=10, fontweight="bold")
        ax.tick_params(axis="y", labelsize=8)
        ax.grid(axis="x", linestyle="--", alpha=0.4)
        ax.spines[["top", "right"]].set_visible(False)

        ax.set_xlim(0, 1.15)

        # Add value + FDR significance labels at the end of each bar
        for bar, fdr in zip(bars, top_df[f"{pc}_fdr"]):
            w = bar.get_width()
            stars = _sig_stars(fdr)
            ax.text(
                w + 0.015,
                bar.get_y() + bar.get_height() / 2,
                f"{w:.2f}{stars}",
                va="center",
                ha="left",
                fontsize=7.5,
                fontweight="bold" if stars else "normal",
                color="#111111" if stars else "#666666",
            )

    handles = [
        plt.Rectangle((0, 0), 1, 1, color=palette["numeric"], label="Numeric ($R^2$)"),
        plt.Rectangle(
            (0, 0), 1, 1, color=palette["categorical"], label="Categorical ($\\eta^2$)"
        ),
    ]

    fig.legend(
        handles=handles,
        loc="upper right",
        frameon=False,
        fontsize=8,
        bbox_to_anchor=(0.98, 0.98),
    )

    fig.suptitle(title, fontsize=11, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path / "PC1_PC2_top_variables.png", dpi=150)
    plt.close(fig)

    return scores


def scatter_by(
    pcs: pd.DataFrame,
    values: pd.Series,
    var: str,
    ratio,
    title: str,
    path,
    labels: pd.Series | None = None,
) -> None:
    """PC1 against PC2, coloured by one variable with non-overlapping point labels."""
    fig, ax = plt.subplots(figsize=(6, 5))
    numeric = pd.to_numeric(values, errors="coerce")

    if numeric.notna().any() and numeric.nunique() > 8:
        na = numeric.isna()
        ax.scatter(
            pcs.loc[na, "PC1"],
            pcs.loc[na, "PC2"],
            c="lightgrey",
            s=35,
            label="no value",
        )
        sc = ax.scatter(
            pcs.loc[~na, "PC1"],
            pcs.loc[~na, "PC2"],
            c=numeric[~na],
            cmap="viridis",
            s=35,
        )
        fig.colorbar(sc, ax=ax, label=var.replace("char:", ""))
        if na.any():
            ax.legend(fontsize=7)
    else:
        levels = values.astype(object).where(values.notna(), "no value").astype(str)
        for level in sorted(levels.unique()):
            idx = levels == level
            ax.scatter(
                pcs.loc[idx, "PC1"], pcs.loc[idx, "PC2"], label=level, s=35, alpha=0.8
            )
        ax.legend(fontsize=7, title=var.replace("char:", ""))

    # Create text objects for each point (defaults to pcs.index if labels is None)
    point_labels = pcs.index.astype(str) if labels is None else labels.astype(str)
    texts = [
        ax.text(x, y, lbl, fontsize=7, alpha=0.85)
        for x, y, lbl in zip(pcs["PC1"], pcs["PC2"], point_labels)
    ]

    ax.set_xlabel(f"PC1 ({ratio[0]:.0%})")
    ax.set_ylabel(f"PC2 ({ratio[1]:.0%})")
    ax.set_title(title, fontsize=9)
    fig.tight_layout()

    # Repel labels from each other and from the scatter points
    adjust_text(
        texts,
        x=pcs["PC1"].to_numpy(),
        y=pcs["PC2"].to_numpy(),
        ax=ax,
        arrowprops=dict(arrowstyle="-", color="black", lw=0.9, alpha=0.6),
        force_text=(0.3, 0.5),
    )

    fig.savefig(path, dpi=150)
    plt.close(fig)


def compare_read_counts(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    group_variables: list[str],
    output_dir,
) -> pd.DataFrame:
    """
    Compare raw sequencing read counts per sample across batches, sites,
    and biological groups.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Calculate total reads per sample
    library_sizes = counts.sum(axis=0)

    # Align metadata with count columns
    df = metadata.reindex(counts.columns).copy()
    df["sample"] = counts.columns
    df["read_count"] = library_sizes.reindex(counts.columns).values
    df["read_count_millions"] = df["read_count"] / 1e6

    summary_records = []

    for variable in group_variables:

        if variable not in df.columns:
            print(f"[warn] {variable} not found in metadata. Skipping.")
            continue

        plot_df = df.dropna(subset=[variable]).copy()

        if plot_df.empty:
            continue

        # Convert categories to strings for plotting
        plot_df[variable] = plot_df[variable].astype(str)

        # Summary statistics
        summary = (
            plot_df.groupby(variable)["read_count"]
            .agg(
                n_samples="count",
                mean="mean",
                median="median",
                std="std",
                min="min",
                max="max",
            )
            .reset_index()
            .rename(columns={variable: "group"})
        )

        summary["variable"] = variable
        summary_records.append(summary)

        # Order groups by median read count
        order = summary.sort_values("median")["group"].tolist()

        # Create figure
        fig, axes = plt.subplots(
            1,
            2,
            figsize=(14, 5),
            gridspec_kw={"width_ratios": [1.5, 1]},
        )

        # Individual sample read counts
        sample_order = plot_df.sort_values(variable)["sample"]

        colors = plt.cm.tab20(np.linspace(0, 1, plot_df[variable].nunique()))

        color_map = dict(zip(order, colors))

        for category in order:
            subset = plot_df[plot_df[variable] == category]

            axes[0].scatter(
                subset["sample"],
                subset["read_count_millions"],
                label=category,
                color=color_map[category],
                s=45,
                alpha=0.85,
            )

        axes[0].set_xticks(range(len(sample_order)))
        axes[0].set_xticklabels(
            sample_order,
            rotation=90,
            fontsize=7,
        )
        axes[0].set_ylabel("Total reads (millions)")
        axes[0].set_title(f"Read counts per sample by {variable}")
        axes[0].legend(
            title=variable,
            fontsize=8,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
        )
        axes[0].grid(axis="y", alpha=0.3)

        # Boxplot comparing group distributions
        box_data = [
            plot_df.loc[
                plot_df[variable] == category,
                "read_count_millions",
            ].values
            for category in order
        ]

        axes[1].boxplot(
            box_data,
            tick_labels=order,
            showmeans=True,
            showfliers=True,
        )

        axes[1].set_ylabel("Total reads (millions)")
        axes[1].set_title(f"Read count distribution by {variable}")
        axes[1].tick_params(axis="x", rotation=45)
        axes[1].grid(axis="y", alpha=0.3)

        fig.suptitle(
            f"Sequencing depth comparison: {variable}",
            fontsize=12,
        )

        fig.tight_layout()

        safe_name = variable.replace("char:", "").replace(" ", "_").replace(".", "_")

        fig.savefig(
            output_dir / f"read_counts_by_{safe_name}.png",
            dpi=150,
            bbox_inches="tight",
        )

        plt.close(fig)

    # Save combined summary
    if summary_records:
        summary_df = pd.concat(summary_records, ignore_index=True)

        summary_df.to_csv(
            output_dir / "read_counts_summary.csv",
            index=False,
        )

        return summary_df

    return pd.DataFrame()


def runPCA(
    data: pd.DataFrame,
    meta_inliers: pd.DataFrame,
    output_folder_path: Path,
    method_name: str,
    ntop: int = 1500,
):
    print(output_folder_path)
    top_genes = data.var(axis=1).sort_values(ascending=False).index[:ntop]
    x = data.loc[top_genes].T
    x = x - x.mean(axis=0)
    pca = PCA(n_components=10)
    names = [f"PC{i + 1}" for i in range(10)]
    pcs = pd.DataFrame(pca.fit_transform(x), index=x.index, columns=names)
    ratio = pca.explained_variance_ratio_

    pd.DataFrame({"component": names, "variance_explained": ratio}).to_csv(
        output_folder_path / "pca_variance.csv", index=False
    )
    pcs.to_csv(output_folder_path / "pca_scores.csv")

    for var in variables:
        if var not in meta_inliers.columns:
            continue
        safe = var.replace("char:", "").replace(" ", "_").replace(".", "_")
        scatter_by(
            pcs,
            meta_inliers[var],
            var,
            ratio,
            f"PCA of {method_name}, top {ntop} variable genes, coloured by {safe}",
            output_folder_path / f"pca_by_{safe}.png",
        )
    plot_top_pc_drivers(
        pcs=pcs,
        metadata=meta_inliers,
        ratio=ratio,
        path=output_folder_path,
        variables=variables,
        ntop=10,
        title=f"Top Variables Separating Samples on PC1 & PC2 | {method_name}",
    )
    print("figures written to", output_folder_path)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    config.PROCESSED.mkdir(parents=True, exist_ok=True)
    counts, meta = load_counts(with_mapping=False)
    counts = counts.loc[:, ~meta["condition"].isna()]
    meta = meta.loc[counts.columns]
    print(f"counts: {counts.shape[0]} genes x {counts.shape[1]} samples")

    group_variables = [
        "char:batch",
        "char:site",
        "condition",
    ]

    # Compare sequencing depth across technical and biological variables
    read_count_summary = compare_read_counts(
        counts=counts,
        metadata=meta,
        group_variables=group_variables,
        output_dir=OUT,
    )

    print("\nRead count summary:")
    print(read_count_summary)

    # Library sizes
    lib = counts.sum(axis=0)
    colors = meta["condition"].map({"Control": "tab:blue", "FSHD": "tab:red"})
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.bar(range(len(lib)), lib.values / 1e6, color=colors.values)
    ax.set_xticks(range(len(lib)))
    ax.set_xticklabels(lib.index, rotation=90, fontsize=6)
    ax.set_ylabel("counts (millions)")
    ax.set_title("Library size per sample (blue: control, red: FSHD)")
    fig.tight_layout()
    fig.savefig(OUT / "library_sizes.png", dpi=150)
    plt.close(fig)

    # Gene filtering.
    # TODO(student): choose and justify the rule. The default keeps a gene when at least
    # `min_samples` samples have `min_count` counts or more, with min_samples = smallest group.
    # Check what your rule does to the four biomarker genes of the paper (config.BIOMARKERS).
    min_count = 10
    min_samples = int(meta["condition"].value_counts().min())  # is 9
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    filtered = counts.loc[keep]

    (OUT / "gene_filter.txt").write_text(
        f"rule: >= {min_count} counts in >= {min_samples} samples\n"
        f"genes before: {len(counts)}\ngenes after: {len(filtered)}\n",
        encoding="utf-8",
    )
    print(f"genes kept: {len(filtered)} / {len(counts)}")
    filtered.to_csv(config.PROCESSED / "counts_filtered.csv")
    meta.to_csv(config.PROCESSED / "metadata.csv")

    # Read counts_filtered.csv and check if the 4 biomarkers remain after filtering, plus report how many samples are kept
    display_genes_and_sample_size_for_biomarkers()

    # rlog function coded in python
    rlog_py, _ = deseq2_rlog_exact_python(filtered)

    # Because rlog is not available in pydeseq2 I use vst to have a similar result.
    vst_df = python_pydeseq2_vst(filtered)

    # top gene selection count
    ntop = 1500

    # make the PCA for each different methods
    for name, data, path in zip(
        ["rlog2", "VST"], [rlog_py, vst_df], [OUT_RLOG2_PCA, OUT_VST_PCA]
    ):
        runPCA(
            data,
            meta_inliers=meta,
            output_folder_path=path,
            method_name=name,
            ntop=ntop,
        )

    # second analysis without outliers

    inliers = [s for s in filtered.columns if s not in config.OUTLIERS]
    meta_inliers = meta.loc[meta.index.isin(inliers)].copy()

    # Because rlog is not available in pydeseq2 I use vst to have a similar result.
    vst_df = python_pydeseq2_vst(filtered[inliers])

    runPCA(
        vst_df,
        meta_inliers=meta_inliers,
        output_folder_path=OUT_INLIERS_VST_PCA,
        method_name="VST Inliers only",
        ntop=ntop,
    )


if __name__ == "__main__":
    main()
