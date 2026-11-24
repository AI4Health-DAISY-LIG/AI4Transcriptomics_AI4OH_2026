"""Step 4: compare your analysis with what the paper reports.

Part A, the biomarker score. The paper scores each biopsy from four DUX4 target genes
and deposits that score and the resulting group on GEO. You write `biomarker_score`
from the paper's Methods; the script then compares your score with the deposited one
and with the MRI and pathology variables.

Part B, the claims. You fill paper/paper_claims.csv; the script joins each claimed gene
with your differential expression results.

Inputs:
  data/raw/<counts>, data/processed/metadata.csv, data/processed/gene_names.csv
  data/raw/GSE115650_biopsy_FPKM.csv.gz (group medians of the figure genes)
  results/step3_de/track_a/de_*.csv.gz (the paper's DE, reproduced in R) and size_factors.csv
  results/step3_de/track_b/groups/de_*.csv.gz (independent DE, groups vs controls)
  paper/paper_claims.csv
Outputs (results/step4_comparison/):
  biomarker_score.csv                 your score next to the authors' score and group
  score_mine_vs_deposited.png
  score_correlations.csv              Spearman of your score with the deposited variables
  score_sensitivity.csv               PyDESeq2 vs deposited size factors; VST score next to the
                                      log2(normalised count + 0.5) score
  score_clusters.txt                  4 clusters of the four VST values vs the deposited groups
  claims_vs_results_<track>__<contrast>.csv   each claimed gene with your padj, log2FC, direction
  figure_genes_group_medians.csv      median log10(FPKM + 0.5) per group for the claimed genes
  summary.txt

Score (Part A): rlog is not in PyDESeq2. Substitute: blind VST (design ~1) fitted on the full
count matrix of the 43 biopsies, then the paper's formula on the four markers: value minus the
control mean of each gene, summed over the genes, minus the minimum over the 43 biopsies (the
GEO convention, authors' getRankSummary). The exact R reproduction is r/score_reproduction.R.

Run:  uv run python scripts/step4_compare_with_paper.py
"""

from __future__ import annotations

import matplotlib
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import config  # noqa: E402

OUT = config.RESULTS / "step4_comparison"
ALPHA = 0.05
DE = config.RESULTS / "step3_de"
NOTES: list[str] = []          # extra summary lines written by biomarker_score()


def size_factors_of(dds) -> pd.Series:
    """PyDESeq2 keeps the size factors in obs (0.5) or obsm (older)."""
    if "size_factors" in dds.obs.columns:
        return dds.obs["size_factors"]
    if "size_factors" in dds.obsm.keys():
        return pd.Series(np.asarray(dds.obsm["size_factors"]), index=dds.obs_names)
    raise SystemExit(f"size factors not found; obs: {list(dds.obs.columns)}, "
                     f"obsm: {list(dds.obsm.keys())}")


def paper_formula(values: pd.DataFrame, control: pd.Series) -> pd.Series:
    """Paper / authors' code: per gene, value minus the control mean; sum over the genes;
    minus the minimum over all biopsies."""
    relative = values.sub(values.loc[:, control].mean(axis=1), axis=0)
    total = relative.sum(axis=0)
    return total - total.min()


def marker_ids(counts: pd.DataFrame, gene_names: pd.DataFrame) -> dict[str, str]:
    ids = {}
    for symbol in config.BIOMARKERS:
        names = gene_names.loc[gene_names["gene_name"] == symbol]
        bare = set(names["gene_id_noversion"]) | set(names["gene_id"].str.split(".").str[0])
        hit = [g for g in counts.index if g.split(".")[0] in bare]
        if len(hit) != 1:
            raise SystemExit(f"{symbol}: {len(hit)} matching count rows")
        ids[symbol] = hit[0]
    return ids


def biomarker_score(counts: pd.DataFrame, meta: pd.DataFrame, gene_names: pd.DataFrame) -> pd.Series | None:
    """Your reproduction of the paper's FSHD score, one value per sample of `counts`.

    counts      raw counts, genes (GENCODE ids) x samples
    meta        one row per sample, with `condition`
    gene_names  gene_id, gene_id_noversion, gene_name, gene_type

    Substitute for rlog: blind PyDESeq2 VST on the full matrix (see the module docstring).
    Also writes score_sensitivity.csv and score_clusters.txt.
    """
    from pydeseq2.dds import DeseqDataSet

    control = meta["condition"] == config.REFERENCE_LEVEL
    ids = marker_ids(counts, gene_names)
    dds = DeseqDataSet(counts=counts.T.round().astype(int), metadata=meta[["condition"]],
                       design="~ condition", quiet=True)
    if not hasattr(dds, "vst"):
        raise SystemExit("this PyDESeq2 has no DeseqDataSet.vst(); check the installed version")
    dds.vst(use_design=False)                      # blind: intercept-only design
    vst = pd.DataFrame(np.asarray(dds.layers["vst_counts"]), index=dds.obs_names,
                       columns=dds.var_names).T
    score = paper_formula(vst.loc[list(ids.values())], control)

    sf = size_factors_of(dds)
    deposited = pd.read_csv(DE / "track_a" / "size_factors.csv",
                            index_col="sample")["size_factor_deposited"].reindex(sf.index)
    rel = ((sf - deposited).abs() / deposited).max()
    normalised = counts.loc[list(ids.values())].div(sf.reindex(counts.columns), axis=1)
    log2_score = paper_formula(np.log2(normalised + 0.5), control)
    deposited_score = pd.to_numeric(meta[config.SCORE_COL], errors="coerce")
    fshd = deposited_score.notna()
    pd.DataFrame({"size_factor_pydeseq2": sf, "size_factor_deposited": deposited,
                  "score_vst": score, "score_log2_norm_0.5": log2_score,
                  "score_deposited": deposited_score}).to_csv(OUT / "score_sensitivity.csv")

    group = pd.to_numeric(meta[config.GROUP_COL], errors="coerce").fillna(0).astype(int)
    clusters = pd.Series(fcluster(linkage(vst.loc[list(ids.values())].T.to_numpy(),
                                          method="complete", metric="euclidean"),
                                  4, criterion="maxclust"), index=vst.columns)
    (OUT / "score_clusters.txt").write_text(
        "4 clusters (complete linkage, Euclidean) of the four VST values (rows) vs the deposited "
        "group (columns, 0 = control); no relabelling applied\n"
        + pd.crosstab(clusters, group).to_string() + "\n", encoding="utf-8")

    def spearman(a: pd.Series) -> float:
        return round(a[fshd].corr(deposited_score[fshd], method="spearman"), 4)

    NOTES.extend([
        f"[score] markers: " + ", ".join(f"{k} {v}" for k, v in ids.items()),
        f"[score] PyDESeq2 vs deposited size factors: max relative difference {rel:.2e}",
        f"[score] Spearman with the deposited score (34 FSHD): VST {spearman(score)}, "
        f"log2(normalised + 0.5) {spearman(log2_score)}; max |VST - deposited| "
        f"{(score - deposited_score)[fshd].abs().max():.3f}",
        "[score] clusters of the VST values vs deposited groups: score_clusters.txt",
    ])
    return score


def compare_score(score: pd.Series, meta: pd.DataFrame) -> list[str]:
    table = pd.DataFrame({
        "my_score": score,
        "condition": meta["condition"],
        "deposited_score": pd.to_numeric(meta[config.SCORE_COL], errors="coerce"),
        "deposited_group": pd.to_numeric(meta[config.GROUP_COL], errors="coerce"),
    })
    table.to_csv(OUT / "biomarker_score.csv")
    rows = []
    for var in [config.SCORE_COL] + [v for v in config.PATHOLOGY_VARS if v != config.SCORE_COL]:
        other = pd.to_numeric(meta[var], errors="coerce")
        both = score.notna() & other.notna()
        rows.append({"variable": var.replace("char:", ""), "n": int(both.sum()),
                     "spearman": round(score[both].corr(other[both], method="spearman"), 3)})
    corr = pd.DataFrame(rows)
    corr.to_csv(OUT / "score_correlations.csv", index=False)

    fig, ax = plt.subplots(figsize=(5, 5))
    for group, part in table.dropna(subset=["deposited_score"]).groupby("deposited_group"):
        ax.scatter(part["deposited_score"], part["my_score"], s=30, label=f"group {int(group)}")
    ax.set_xlabel("score deposited on GEO (rnaseq.score)")
    ax.set_ylabel("my score")
    ax.legend(fontsize=7, title="authors' group")
    fig.tight_layout()
    fig.savefig(OUT / "score_mine_vs_deposited.png", dpi=150)
    plt.close(fig)
    by_group = table.groupby(table["deposited_group"].fillna(0).astype(int))["my_score"]
    return (["[score] Spearman of my score with the deposited variables:", corr.to_string(index=False),
             "[score] my score by authors' group (0 = control): min, median, max",
             by_group.agg(["count", "min", "median", "max"]).round(2).to_string()])


def compare_claims(res: pd.DataFrame, claims: pd.DataFrame, symbol_to_id: dict[str, str]) -> pd.DataFrame:
    lookup = res.copy()
    lookup.index = lookup.index.to_series().astype(str).str.split(".").str[0]
    lookup = lookup[~lookup.index.duplicated()]
    rows = []
    for _, c in claims.iterrows():
        gene = str(c["gene"]).strip()
        g = symbol_to_id.get(gene, gene.split(".")[0])
        if g not in lookup.index:
            rows.append({**c.to_dict(), "gene_id": g, "found_in_my_results": False})
            continue
        h = lookup.loc[g]
        my_dir = "up" if h["log2FoldChange"] > 0 else "down"
        rows.append({
            **c.to_dict(),
            "gene_id": g,
            "found_in_my_results": True,
            "my_log2FC": round(float(h["log2FoldChange"]), 3),
            "my_padj": h["padj"],
            "my_significant": bool(h["padj"] < ALPHA) if pd.notna(h["padj"]) else False,
            "my_direction": my_dir,
            "direction_agrees": my_dir == str(c.get("direction_in_paper", "")).strip().lower(),
        })
    return pd.DataFrame(rows)


def de_files() -> list[tuple[str, "Path"]]:
    """(tag, file) for the reproduced paper contrasts (Track A) and the Track B groups fit."""
    files = [(f"track_a/{f.name.removeprefix('de_').removesuffix('.csv.gz')}", f)
             for f in sorted((DE / "track_a").glob("de_*.csv.gz"))]
    files += [(f"track_b/groups/{f.name.removeprefix('de_').removesuffix('.csv.gz')}", f)
              for f in sorted((DE / "track_b" / "groups").glob("de_*.csv.gz"))]
    return files


def figure_gene_medians(meta: pd.DataFrame, claims: pd.DataFrame,
                        symbol_to_id: dict[str, str]) -> list[str]:
    """Median log10(FPKM + 0.5) per group, the scale of Figures 4 and 5, for the claimed genes.
    Group-vs-control tests alone do not show the 'progressive' order the captions describe."""
    path = config.RAW / "GSE115650_biopsy_FPKM.csv.gz"
    if not path.exists():
        return [f"[figures] {path.name} not found: group medians skipped"]
    fpkm = pd.read_csv(path, index_col=0)[meta.index]
    fpkm.index = fpkm.index.str.split(".").str[0]
    group = pd.to_numeric(meta[config.GROUP_COL], errors="coerce").fillna(0).astype(int)
    labels = group.map({0: "Control", 1: "G1", 2: "G2", 3: "G3", 4: "G4"})
    outliers = meta.index.isin(["01-0037", "32-0016"])
    rows = []
    for gene in sorted(claims["gene"].str.strip().unique()):
        gid = symbol_to_id.get(gene, "").split(".")[0]
        if gid not in fpkm.index:
            continue
        values = np.log10(fpkm.loc[gid] + 0.5)
        med = values.groupby(labels).median()
        med["G1_without_outliers"] = values[(labels == "G1") & ~outliers].median()
        med["G1_outliers"] = values[outliers].median()
        rows.append({"gene": gene, "gene_id": gid, **med.round(3).to_dict(),
                     "G2_le_G3_le_G4": bool(med["G2"] <= med["G3"] <= med["G4"]),
                     "G2_G3_G4_above_control": bool(min(med["G2"], med["G3"], med["G4"])
                                                    > med["Control"])})
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "figure_genes_group_medians.csv", index=False)
    return [f"[figures] {len(table)} genes; group medians G2 <= G3 <= G4: "
            f"{int(table['G2_le_G3_le_G4'].sum())}; G2, G3 and G4 all above controls: "
            f"{int(table['G2_G3_G4_above_control'].sum())}"]


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
        summary += NOTES + compare_score(score.reindex(meta.index), meta)

    # Part B
    claims = pd.read_csv(config.ROOT / "paper" / "paper_claims.csv", dtype=str)
    claims = claims[claims["gene"].notna() & (claims["gene"].str.strip() != "")]
    if claims.empty:
        summary.append("[claims] paper/paper_claims.csv has no gene rows yet: part B skipped")
    symbol_to_id = dict(zip(gene_names["gene_name"], gene_names["gene_id_noversion"]))
    keys = set(claims["my_contrast"].fillna("").str.strip()) - {""}
    unmatched = sorted(keys - {tag for tag, _ in de_files()})
    summary.append(f"[claims] contrast keys: {len(keys)}; matching no result file: "
                   + (", ".join(unmatched) if unmatched else "none"))
    for tag, de_file in de_files() if not claims.empty else []:
        mine = claims[claims["my_contrast"].fillna("").str.strip().isin(["", tag])]
        if mine.empty:
            continue
        table = compare_claims(pd.read_csv(de_file, index_col=0), mine, symbol_to_id)
        table.to_csv(OUT / f"claims_vs_results_{tag.replace('/', '__')}.csv", index=False)
        found = table[table["found_in_my_results"]]
        summary += [
            f"[{tag}] claimed genes: {len(table)}  found in my results: {len(found)}",
            f"[{tag}] significant in my results (padj < {ALPHA}): {int(found.get('my_significant', pd.Series(dtype=bool)).sum())}",
            f"[{tag}] same direction as the paper: {int(found.get('direction_agrees', pd.Series(dtype=bool)).sum())}",
        ]
    if not claims.empty:
        summary += figure_gene_medians(meta, claims, symbol_to_id)
    (OUT / "summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))
    print("\nNext: fill worksheets/step4_comparison.md, then write your report")


if __name__ == "__main__":
    main()
