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

import config  # noqa: E402

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
    return None


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
