"""Step 4: compare your analysis with what the paper reports.

Part A, the biomarker score. The paper scores each biopsy from four DUX4 target genes
and deposits that score and the resulting group on GEO. You write `biomarker_score`
from the paper's Methods; the script then compares your score with the deposited one
and with the MRI and pathology variables.

Part B, the claims. You fill paper/paper_claims.csv; the script joins each claimed gene
with your differential expression results.

Inputs:
  data/raw/<counts>, data/processed/metadata.csv, data/processed/gene_names.csv
  results/step3_de/de_*.csv, paper/paper_claims.csv
Outputs (results/step4_comparison/):
  biomarker_score.csv                 your score next to the authors' score and group
  score_mine_vs_deposited.png
  score_correlations.csv              Spearman of your score with the deposited variables
  claims_vs_results_<contrast>.csv    each claimed gene with your padj, log2FC, direction check
  summary.txt

Run:  uv run python scripts/step4_compare_with_paper.py
"""

from __future__ import annotations

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np
import config  # noqa: E402
from pydeseq2.dds import DeseqDataSet
import seaborn as sns
from scipy.stats import mannwhitneyu

OUT = config.RESULTS / "step4_comparison"
ALPHA = 0.05


def biomarker_score(
    counts: pd.DataFrame, meta: pd.DataFrame, gene_names: pd.DataFrame
) -> pd.Series | None:
    """Your reproduction of the paper's FSHD score, one value per sample of `counts`.

    counts      raw counts, genes (GENCODE ids) x samples
    meta        one row per sample, with `condition`
    gene_names  gene_id, gene_id_noversion, gene_name, gene_type

    TODO(student): read "RNA-Seq analysis" in the paper's Methods and implement the score
    for config.BIOMARKERS. The paper names a DESeq2 transformation that PyDESeq2 does not
    have: choose a substitute, and state it in worksheets/step4_comparison.md. Return a Series
    indexed like counts.columns.
    """
    all_samples = counts.columns
    outliers = ["01-0037", "32-0016"]
    inliers = [s for s in all_samples if s not in outliers]

    # Fit VST parameters ONLY on clean inlier samples
    meta_inliers = meta.loc[inliers].copy()
    counts_inliers_T = counts[inliers].T.round().astype(int)

    dds_fit = DeseqDataSet(
        counts=counts_inliers_T,
        metadata=meta_inliers,
        design="~1",
        quiet=True,
    )
    dds_fit.vst()  # Learns size factors & dispersion trend from inliers

    # Apply the learned transformation to ALL samples (including outliers)
    counts_all_T = counts[all_samples].T.round().astype(int)
    vst_all_T = dds_fit.vst_transform(counts_all_T)

    vst_df = pd.DataFrame(
        vst_all_T.T,
        index=counts.index,
        columns=all_samples,
    )

    # Filter to the 4 biomarker genes
    gene_filtered = gene_names.loc[gene_names["gene_name"].isin(config.BIOMARKERS)]
    stripped_index = vst_df.index.str.split(".").str[0]
    mask = vst_df.index.isin(gene_filtered["gene_id"]) | stripped_index.isin(
        gene_filtered["gene_id_noversion"]
    )
    vst_biomarkers = vst_df.loc[mask]

    # Compute control mean using ONLY clean control samples
    inlier_control_mask = meta_inliers["condition"].str.lower().str.contains("control")
    clean_controls = meta_inliers.index[inlier_control_mask]
    control_mean = vst_biomarkers.loc[:, clean_controls].mean(axis=1)

    # Score ALL samples against the clean control baseline
    log2_diff_from_ctrl = vst_biomarkers.sub(control_mean, axis=0)
    scores = log2_diff_from_ctrl.sum(axis=0)

    return scores.reindex(all_samples)


def plot_stir_comparison(
    scores: pd.Series,
    meta: pd.DataFrame,
    stir_col: str = "STIR",
    condition_col: str = "condition",
    save_path: str | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Compare FSHD biomarker scores between Normal MRI and STIR+ MRI (> 0) samples."""
    # Align scores and metadata
    df = (
        pd.DataFrame({"score": scores}).join(meta, how="inner").dropna(subset=["score"])
    )

    # Assign group labels
    def assign_group(row: pd.Series) -> str | None:
        stir_val = pd.to_numeric(row.get(stir_col), errors="coerce")
        if pd.isna(stir_val):
            return None
        return "STIR+ MRI (>0)" if stir_val > 0 else "Normal MRI"

    df["group"] = df.apply(assign_group, axis=1)
    df = df.dropna(subset=["group"])

    order = [g for g in ["STIR+ MRI (>0)", "Normal MRI"] if g in df["group"].unique()]

    # Statistical comparison: STIR- vs STIR+
    neg_scores = df.loc[df["group"] == "Normal MRI", "score"]
    pos_scores = df.loc[df["group"] == "STIR+ MRI (>0)", "score"]
    _, p_val = mannwhitneyu(pos_scores, neg_scores, alternative="two-sided")

    # Build boxplot + jittered individual points
    fig, ax = plt.subplots(figsize=(6.5, 5), dpi=150)
    palette = {
        "Normal MRI": "#66c2a5",
        "STIR+ MRI (>0)": "#fc8d62",
    }

    sns.boxplot(
        data=df,
        x="group",
        y="score",
        order=order,
        palette=palette,
        width=0.45,
        fliersize=0,
        boxprops={"alpha": 0.75},
        ax=ax,
    )
    sns.stripplot(
        data=df,
        x="group",
        y="score",
        order=order,
        color="#222222",
        size=6,
        jitter=0.15,
        alpha=0.85,
        ax=ax,
    )

    # Annotate sample counts (n) on x-axis labels & p-value bracket
    counts_by_group = df["group"].value_counts()
    ax.set_xticklabels([f"{g}\n(n={counts_by_group[g]})" for g in order])

    idx_neg = order.index("Normal MRI")
    idx_pos = order.index("STIR+ MRI (>0)")
    y_max = df["score"].max()
    y_span = df["score"].max() - df["score"].min()
    y_bar = y_max + 0.08 * y_span
    h = 0.03 * y_span

    ax.plot(
        [idx_neg, idx_neg, idx_pos, idx_pos],
        [y_bar, y_bar + h, y_bar + h, y_bar],
        lw=1.3,
        c="black",
    )
    ax.text(
        (idx_neg + idx_pos) / 2,
        y_bar + h + 0.01 * y_span,
        f"Mann-Whitney p = {p_val:.4f}",
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
    )

    ax.set_ylim(top=y_bar + 0.14 * y_span)
    ax.set_xlabel("")
    ax.set_ylabel("FSHD Biomarker Score (VST log2 diff)")
    ax.set_title("Biomarker Score by Muscle STIR Status")
    sns.despine(ax=ax)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, bbox_inches="tight")

    return fig, ax


def compare_score(score: pd.Series, meta: pd.DataFrame) -> list[str]:
    table = pd.DataFrame(
        {
            "my_score": score,
            "condition": meta["condition"],
            "deposited_score": pd.to_numeric(meta[config.SCORE_COL], errors="coerce"),
            "deposited_group": pd.to_numeric(meta[config.GROUP_COL], errors="coerce"),
        }
    )
    table.to_csv(OUT / "biomarker_score.csv")
    plot_stir_comparison(
        score,
        meta,
        stir_col="char:stir",
        save_path=OUT / "Stir_plot.png",
    )
    rows = []
    for var in [config.SCORE_COL] + [
        v for v in config.PATHOLOGY_VARS if v != config.SCORE_COL
    ]:
        other = pd.to_numeric(meta[var], errors="coerce")
        both = score.notna() & other.notna()
        rows.append(
            {
                "variable": var.replace("char:", ""),
                "n": int(both.sum()),
                "spearman": round(score[both].corr(other[both], method="spearman"), 3),
            }
        )
    corr = pd.DataFrame(rows)
    corr.to_csv(OUT / "score_correlations.csv", index=False)

    fig, ax = plt.subplots(figsize=(5, 5))
    for group, part in table.dropna(subset=["deposited_score"]).groupby(
        "deposited_group"
    ):
        ax.scatter(
            part["deposited_score"], part["my_score"], s=30, label=f"group {int(group)}"
        )

    ax.axline((0, 0), slope=1, color="blue", linestyle="--", label="y = x")
    ax.set_xlabel("score deposited on GEO (rnaseq.score)")
    ax.set_ylabel("my score")
    ax.legend(fontsize=7, title="authors' group")
    fig.tight_layout()
    fig.savefig(OUT / "score_mine_vs_deposited.png", dpi=150)
    plt.close(fig)

    by_group = table.groupby(table["deposited_group"].fillna(0).astype(int))["my_score"]
    return [
        "[score] Spearman of my score with the deposited variables:",
        corr.to_string(index=False),
        "[score] my score by authors' group (0 = control): min, median, max",
        by_group.agg(["count", "min", "median", "max"]).round(2).to_string(),
    ]


def compare_claims(
    res: pd.DataFrame, claims: pd.DataFrame, symbol_to_id: dict[str, str]
) -> pd.DataFrame:
    lookup = res.copy()
    lookup.index = lookup.index.to_series().astype(str).str.split(".").str[0]
    lookup = lookup[~lookup.index.duplicated()]
    rows = []
    for _, c in claims.iterrows():
        gene = str(c["gene"]).strip()
        g = symbol_to_id.get(gene, gene.split(".")[0])
        if g not in lookup.index:
            rows.append({**c.to_dict(), "found_in_my_results": False})
            continue
        h = lookup.loc[g]
        my_dir = "up" if h["log2FoldChange"] > 0 else "down"
        rows.append(
            {
                **c.to_dict(),
                "found_in_my_results": True,
                "my_log2FC": round(float(h["log2FoldChange"]), 3),
                "my_padj": h["padj"],
                "my_significant": (
                    bool(h["padj"] < ALPHA) if pd.notna(h["padj"]) else False
                ),
                "my_direction": my_dir,
                "direction_agrees": my_dir
                == str(c.get("direction_in_paper", "")).strip().lower(),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary: list[str] = []
    gene_names = pd.read_csv(config.GENE_NAMES)
    meta = pd.read_csv(config.PROCESSED / "metadata.csv", index_col=0)
    counts = pd.read_csv(config.RAW / config.COUNTS_FILE, index_col=0)[meta.index]

    # Part A
    score = biomarker_score(counts, meta, gene_names)
    if score is None:
        summary.append("[score] biomarker_score() is not written yet: part A skipped")
    else:
        summary += compare_score(score.reindex(meta.index), meta)

    # Part B
    claims = pd.read_csv(config.ROOT / "paper" / "paper_claims.csv", dtype=str)
    claims = claims[claims["gene"].notna() & (claims["gene"].str.strip() != "")]
    if claims.empty:
        summary.append(
            "[claims] paper/paper_claims.csv has no gene rows yet: part B skipped"
        )
    symbol_to_id = dict(zip(gene_names["gene_name"], gene_names["gene_id_noversion"]))
    for de_file in (
        sorted((config.RESULTS / "step3_de").glob("de_*.csv"))
        if not claims.empty
        else []
    ):
        tag = de_file.stem.removeprefix("de_")
        mine = claims[claims["my_contrast"].fillna("").str.strip().isin(["", tag])]
        if mine.empty:
            continue
        table = compare_claims(pd.read_csv(de_file, index_col=0), mine, symbol_to_id)
        table.to_csv(OUT / f"claims_vs_results_{tag}.csv", index=False)
        found = table[table["found_in_my_results"]]
        summary += [
            f"[{tag}] claimed genes: {len(table)}  found in my results: {len(found)}",
            f"[{tag}] significant in my results (padj < {ALPHA}): {int(found.get('my_significant', pd.Series(dtype=bool)).sum())}",
            f"[{tag}] same direction as the paper: {int(found.get('direction_agrees', pd.Series(dtype=bool)).sum())}",
        ]
    (OUT / "summary.txt").write_text("\n".join(summary), encoding="utf-8")
    print("\n".join(summary))
    print("\nNext: fill worksheets/step4_comparison.md, then write your report")


if __name__ == "__main__":
    main()
