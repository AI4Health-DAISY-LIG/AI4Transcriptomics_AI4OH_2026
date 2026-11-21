"""Step 2: quality control and exploration of the raw counts.

Outputs (results/step2_qc/):
  library_sizes.png                 total counts per sample, coloured by condition
  gene_filter.txt                   genes before and after filtering, with the rule used
  pca_variance.csv                  variance explained by the first 10 components
  pca_scores.csv                    PC1 to PC10 per sample
  pca_by_<variable>.png             PC1 vs PC2 coloured by condition and by each variable
  data/processed/counts_filtered.csv   filtered raw counts (genes x samples), input of step 3
  data/processed/metadata.csv          one row per sample column of the counts, with condition

Decisions marked TODO(student) are yours. Record them in worksheets/step2_qc.md.

Run:  uv run python scripts/step2_qc_exploration.py
"""

from __future__ import annotations

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402
from step2_view_genes_filtering_effects import (
    display_genes_and_sample_size_for_biomarkers,
)
from pca import deseq2_rlog_exact_python, python_pydeseq2_vst
from adjustText import adjust_text

import config  # noqa: E402

OUT = config.RESULTS / "step2_qc"
AUDIT = config.RESULTS / "step1_audit"


def load_counts() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Raw counts (genes x samples) and metadata aligned on the sample columns."""
    counts = pd.read_csv(config.RAW / config.COUNTS_FILE, index_col=0)
    samples = pd.read_csv(AUDIT / "sample_table.csv")
    # TODO(student): read results/step1_audit/counts_header_check.txt first. This join keeps
    # only the count columns whose name equals a GEO sample_id. If a column has no match,
    # find out which biopsy it is and decide: rename it, or leave it out and say why.

    # # remove a list of samples sample that have
    # counts = counts.drop(columns=["01-0037", "32-0016"])

    # samples = samples.loc[~samples["title"].isin(["01-0037", "32-0016"])]

    # rename sample_id 01-0022 to 01-0022-1 to match the geo record
    samples.loc[samples["gsm"] == "GSM3186680", "char:sample_id"] = "01-0022-1"

    key = "char:sample_id" if "char:sample_id" in samples.columns else "title"
    meta = samples.set_index(key).reindex(counts.columns)
    meta.to_csv(config.RAW / "tmp.csv")
    missing = meta["condition"].isna()
    if missing.any():
        print(
            f"[warn] {missing.sum()} count column(s) dropped, no GEO record under that name: "
            f"{list(counts.columns[missing])}"
        )
    return counts, meta


def plot_top_pc_drivers(
    pcs: pd.DataFrame,
    metadata: pd.DataFrame,
    ratio,
    path,
    top_n: int = 20,
    title: str = "Top Variables Separating Samples on PC1 & PC2",
) -> pd.DataFrame:
    """
    Ranks variables in `metadata` by how well they separate samples along PC1 and PC2
    (using R² for numeric and η² for categorical variables) and plots the top `top_n`.
    """

    records = []

    # Iterate through each column in the metadata to evaluate its association with PCs
    for col in metadata.columns:
        # Filter to keep only the variable used in the PCA plot
        if not (col.startswith("char:")) and col != "condition":
            continue

        # Specific exclusion for a known non-informative or redundant variable
        if col == "char:subject status":
            continue

        s = metadata[col]

        # Data Quality Check: Ensure there are enough valid data points to compute statistics reliably
        valid = s.notna()
        if valid.sum() < 3:
            continue

        # Determine Variable Type (Numeric vs Categorical)
        # Attempt to convert to numeric : if successful for all non-NaN values and has sufficient variance, treat as numeric
        s_valid = s[valid]
        numeric = pd.to_numeric(s_valid, errors="coerce")

        # Heuristic for categorical data
        is_numeric = numeric.notna().all() and numeric.nunique() > 8

        row = {
            "variable": col.replace("char:", ""),  # Clean variable name for plotting
            "type": "numeric" if is_numeric else "categorical",
        }

        # Calculate association strength (Variance Explained) for PC1 and PC2
        for pc in ["PC1", "PC2"]:
            y = pcs.loc[valid, pc].to_numpy(
                dtype=float
            )  # The Principal Component scores (target variable)

            # Total Sum of Squares (SS_total): Variance in the PC scores themselves
            ss_total = np.sum((y - y.mean()) ** 2)

            # If PC has no variance (constant), association is zero
            if ss_total == 0:
                row[pc] = 0.0
                continue

            if is_numeric:
                x = numeric.to_numpy(dtype=float)

                # Check for constant numeric variable (no variance to explain)
                if np.std(x) == 0:
                    row[pc] = 0.0
                else:
                    # For Numeric Variables: Use Pearson Correlation Coefficient squared (R²)
                    # R² represents the proportion of variance in PC scores explained by the linear relationship with variable x
                    r = np.corrcoef(x, y)[0, 1]
                    row[pc] = float(r**2)
            else:
                # For Categorical Variables: Use Eta-squared (η²)
                # η² represents the proportion of variance in PC scores explained by group membership

                groups = s_valid.astype(str)

                # Check for valid categorical structure: Must have at least 2 groups to compare means and not be unique per sample (which would make it an ID, not a category)
                if groups.nunique() < 2 or groups.nunique() >= len(groups):
                    row[pc] = 0.0
                else:
                    df_tmp = pd.DataFrame({"y": y, "g": groups})

                    # Calculate Sum of Squares Between Groups
                    group_means = df_tmp.groupby("g")["y"].transform("mean")
                    ss_between = np.sum((group_means - y.mean()) ** 2)

                    # Eta-squared
                    row[pc] = float(ss_between / ss_total)

        records.append(row)

    # Convert list of dictionaries into a DataFrame for easy sorting and plotting
    scores = pd.DataFrame(records)

    fig, axes = plt.subplots(1, 2, figsize=(11, max(4, top_n * 0.4)), sharey=False)
    palette = {"numeric": "#2b5c8f", "categorical": "#d95f02"}

    # Iterate through axes to plot PC1 and PC2 results separately
    for ax, pc, var_ratio in zip(axes, ["PC1", "PC2"], ratio[:2]):

        # Sort scores by the current PC's score (descending) and take top N
        # .iloc[::-1] reverses order so highest value is at the top of the bar chart
        top_df = scores.nlargest(top_n, pc).iloc[::-1]
        colors = [palette[t] for t in top_df["type"]]

        # Create horizontal bar chart
        bars = ax.barh(
            top_df["variable"], top_df[pc], color=colors, edgecolor="none", height=0.7
        )

        # Configure axes and labels
        ax.set_xlim(0, 1.08)  # Add slight padding to the right for text labels
        ax.set_xlabel("Variance Explained ($R^2$ / $\\eta^2$)", fontsize=9)
        ax.set_title(f"{pc} ({var_ratio:.0%})", fontsize=10, fontweight="bold")

        # Style adjustments for readability
        ax.tick_params(axis="y", labelsize=8)
        ax.grid(axis="x", linestyle="--", alpha=0.4)
        ax.spines[["top", "right"]].set_visible(False)

        # Add value labels at the end of each bar
        for bar in bars:
            w = bar.get_width()
            ax.text(
                w + 0.015,
                bar.get_y() + bar.get_height() / 2,
                f"{w:.2f}",
                va="center",
                ha="left",
                fontsize=7.5,
                color="#333333",
            )

    # Create a shared legend for the entire figure indicating variable types
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=palette["numeric"], label="Numeric ($R^2$)"),
        plt.Rectangle(
            (0, 0), 1, 1, color=palette["categorical"], label="Categorical ($\\eta^2$)"
        ),
    ]

    # Place legend outside the plot area to avoid overlapping bars
    fig.legend(
        handles=handles,
        loc="upper right",
        frameon=False,
        fontsize=8,
        bbox_to_anchor=(0.98, 0.98),
    )

    # Finalize figure layout and save
    fig.suptitle(title, fontsize=11, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=150)
    plt.close(fig)


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


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    config.PROCESSED.mkdir(parents=True, exist_ok=True)
    counts, meta = load_counts()
    counts = counts.loc[:, ~meta["condition"].isna()]
    meta = meta.loc[counts.columns]
    print(f"counts: {counts.shape[0]} genes x {counts.shape[1]} samples")

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

    rlog_py, p_var = deseq2_rlog_exact_python(filtered)

    # Because rlog is not available in pydeseq2 I use vst to have a similar result.
    vst_df = python_pydeseq2_vst(filtered)

    # make the PCA for each different methods
    for name, data in zip(["rlog2", "VST"], [rlog_py, vst_df]):
        n_top = 1500
        top_genes = data.var(axis=1).sort_values(ascending=False).index[:n_top]
        x = data.loc[top_genes].T
        x = x - x.mean(axis=0)

        pca = PCA(n_components=10)
        names = [f"PC{i + 1}" for i in range(10)]
        pcs = pd.DataFrame(pca.fit_transform(x), index=x.index, columns=names)
        ratio = pca.explained_variance_ratio_

        pd.DataFrame({"component": names, "variance_explained": ratio}).to_csv(
            OUT / "pca_variance.csv", index=False
        )
        pcs.to_csv(OUT / "pca_scores.csv")

        variables = (
            ["condition"]
            + config.TECHNICAL_VARS
            + config.BIOLOGICAL_VARS
            + config.PATHOLOGY_VARS
        )
        for var in variables:
            if var not in meta.columns:
                continue
            safe = var.replace("char:", "").replace(" ", "_").replace(".", "_")
            scatter_by(
                pcs,
                meta[var],
                var,
                ratio,
                f"PCA of {name}, top {n_top} variable genes, coloured by {var}",
                OUT / f"pca_by_{safe}_{name}.png",
            )
        plot_top_pc_drivers(
            pcs,
            meta,
            ratio,
            OUT / f"pca_PC1_and_PC2_top_variables_{name}.png",
            10,
            "Top Variables Separating Samples on PC1 & PC2",
        )
        print("figures written to", OUT)
        print("\nNext: fill worksheets/step2_qc.md")


if __name__ == "__main__":
    main()
