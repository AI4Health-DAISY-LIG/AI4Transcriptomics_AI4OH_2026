"""Bonus: choose the genes for the Translator queries, by the rule in worksheets/bonus_translator.md.

View A (required gene query): primary FSHD vs control genes that are robust (planned five
sensitivities and Rochester), |log2FC| >= 1, and follow the DUX4-target score within FSHD.
View B (optional query): genes whose expression follows the DUX4-target score within FSHD,
robust in the score checks, |log2FC across the G1-G4 span| >= 1, not a Yao 2014 DUX4 target.
Both views exclude non-autosomal, non-coding, biomarker, immediate-early and composition genes.

Inputs:  results/step3_de/track_b/ (primary, score, groups, robustness tables)
         bonus_lists/ (reference lists, see bonus_lists/sources.tsv)
Outputs (results/bonus/translator/):
  list_matching.csv          symbols of each list found in our annotation
  yao2014_matching.csv       each Yao 2014 gene: matched by Ensembl ID, by symbol, or absent
  pool_<view>_<filter>.csv   top 10 of each pool, with flags and the reason for each skip
  influence_checks.csv       view A refits without the largest-Cook's biopsy
  known_regulators.csv       the known DUX4 regulators in both views
  top10_by_run_<view>.csv    exploration: each run's own top 10 after the gene filters, and
                             the rank of every such gene in every run (no effect on the choice)
  selection.txt              the chosen genes and their Translator question

Run:  uv run --locked python scripts/bonus_select_genes.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import config
import step3_differential_expression as de

LISTS = config.ROOT / "bonus_lists"
OUT = config.RESULTS / "bonus" / "translator"
TRACK_B = config.RESULTS / "step3_de" / "track_b"
ALPHA = 0.05
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
IEG_ALIASES = {"CCN1": "CYR61", "CCN2": "CTGF"}       # GENCODE v24 still uses the old names
COMPOSITION = {"HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ",
               "ADIPOQ", "LEP", "PLIN1", "PLIN4", "CIDEA", "CIDEC", "FABP4",
               "S100A8", "S100A9", "S100A12", "FCGR3B", "CXCR1", "CXCR2", "CSF3R", "ELANE",
               "MPO", "CAMP", "LTF"}
KNOWN_REGULATORS = ["SMCHD1", "DNMT3B", "LRIF1", "CHD4", "HDAC2", "MTA2", "MBD2", "RBBP4",
                    "CHAF1A", "CHAF1B", "MBD3L2", "MAPK14", "MAPK11"]
NOT_CHROMATIN = {"RNA modification", "Scaffold protein", "Protein modification"}
TOP = 10


def reference_lists() -> dict[str, set[str]]:
    ieg = {s.upper() for s in pd.read_csv(LISTS / "Tullai2007_NHLBI_ESBL_ieg.tsv",
                                          sep="\t")["source_symbol"]}
    ieg |= {IEG_ALIASES[s] for s in ieg if s in IEG_ALIASES}
    tf = set(pd.read_csv(LISTS / "TF_names_v_1.01.txt", header=None)[0])
    epi = pd.read_csv(LISTS / "EpiGenes_main.csv")
    functions = epi["Function"].fillna("#").str.split(r",\s*").map(set)
    # Drop only entries whose functions are all unrelated to chromatin: RNA modification (with
    # or without scaffold) or protein modification of non-histone proteins. 713 of 801 remain.
    keep = ~functions.map(lambda f: f <= NOT_CHROMATIN and f != {"Scaffold protein"})
    chromatin = set(epi.loc[keep, "HGNC_symbol"])
    regeneration = set(pd.read_csv(LISTS / "QuickGO_GO0043403_human.tsv",
                                   sep="\t")["source_symbol"].str.upper())
    yao = pd.read_csv(LISTS / "Yao2014_TableS1_DUX4_induced_genes.tsv", sep="\t")
    return {"immediate_early": ieg, "transcription_factor": tf, "chromatin": chromatin,
            "regeneration": regeneration, "yao2014_ensembl": set(yao["ensembl_gene_id"]),
            "yao2014_symbol": set(yao["source_gene_symbol"]),
            "composition": COMPOSITION, "known_regulator": set(KNOWN_REGULATORS)}


def flags(table: pd.DataFrame, lists: dict[str, set[str]]) -> pd.DataFrame:
    name = table["gene_name"].fillna("")
    out = pd.DataFrame(index=table.index)
    out["autosomal"] = table["chromosome"].isin(AUTOSOMES)
    out["protein_coding"] = table["gene_type"] == "protein_coding"
    out["biomarker"] = name.isin(config.BIOMARKERS)
    out["immediate_early"] = name.isin(lists["immediate_early"])
    out["composition"] = name.isin(lists["composition"]) | name.str.startswith("HLA-")
    # Yao 2014 by Ensembl ID; for the IDs absent from GENCODE v24 (retired or renamed since
    # 2014), by symbol, so that targets such as DUXA and the PRAMEF genes are not missed.
    out["yao2014_target"] = table.index.str.split(".").str[0].isin(lists["yao2014_ensembl"]) \
        | name.isin(lists["yao2014_symbol_fallback"])
    for key in ["transcription_factor", "chromatin", "regeneration", "known_regulator"]:
        out[key] = name.isin(lists[key])
    return out


def ranked(table: pd.DataFrame, padj: str, stat: str) -> pd.DataFrame:
    order = pd.DataFrame({"padj": table[padj].to_numpy(), "abs_stat": -table[stat].abs().to_numpy(),
                          "gene_id": table.index.to_numpy()})
    return table.iloc[order.sort_values(["padj", "abs_stat", "gene_id"]).index]


class Influence:
    """View A check: refit the primary model without one biopsy (cached by biopsy)."""

    def __init__(self) -> None:
        self.x, self.meta = de.sample_table()
        self.fits: dict[str, pd.Series] = {}

    def lfc_without(self, sample: str) -> pd.Series:
        if sample not in self.fits:
            spec = dict(de.FITS["primary"], exclude=de.FITS["primary"]["exclude"] + [sample])
            m = de.fit_samples(self.meta, spec)
            dds = de.DeseqDataSet(counts=self.x.loc[m.index],
                                  metadata=m[["batch", "sex", "condition"]],
                                  design=de.MAIN, refit_cooks=True, quiet=True)
            dds.deseq2()
            stats = de.DeseqStats(dds, contrast=["condition", "FSHD", de.REF], alpha=ALPHA,
                                  quiet=True)
            stats.summary()
            self.fits[sample] = stats.results_df["log2FoldChange"]
            print(f"  influence refit without {sample} done")
        return self.fits[sample]


RUNS_A = {"primary": "primary", "ladder_batch": "ladder_batch", "site": "site",
          "all43": "all43", "restore_32_0004": "restore_32_0004",
          "omit_32_0016": "omit_32_0016", "rochester": "rochester", "ladder_none": "ladder_none"}
RUNS_B = ["score", "score_nosex", "score_all34", "score_omit_32_0006", "score_omit_01_0025"]


def top10_by_run(runs: dict[str, pd.DataFrame], size_col: str, lists: dict[str, set[str]],
                 main_eligible: pd.Index, main_run: str) -> pd.DataFrame:
    """Each run's own top 10: gene filters + padj < 0.05 and |size| >= 1 in that run, ranked
    by that run's padj. Returns the union with every gene's rank in every run."""
    ranks = {}
    for run, table in runs.items():
        f = flags(table, lists)
        ok = (table["padj"] < ALPHA) & (table[size_col].abs() >= 1) & f["autosomal"] & \
            f["protein_coding"] & ~f["biomarker"] & ~f["immediate_early"] & ~f["composition"]
        ranks[run] = pd.Series(range(1, int(ok.sum()) + 1),
                               index=ranked(table[ok], "padj", "stat").index)
    top = sorted(set().union(*[r.index[:TOP] for r in ranks.values()]))
    main = runs[main_run]
    out = pd.DataFrame({"gene_name": main["gene_name"].reindex(top)}, index=pd.Index(top))
    for run, r in ranks.items():
        out[f"rank_{run}"] = r.reindex(top)
    out[f"{main_run}_{size_col}"] = main[size_col].reindex(top)
    out[f"{main_run}_padj"] = main["padj"].reindex(top)
    out["eligible_by_rule"] = out.index.isin(main_eligible)
    out["in_main_top10"] = out[f"rank_{main_run}"] <= TOP
    return out.sort_values(f"rank_{main_run}", na_position="last")


def pick(pools: list[tuple[str, pd.DataFrame]], check) -> tuple[str | None, list[dict]]:
    """First gene, pool by pool in rank order, that passes `check` (None if none)."""
    log = []
    for pool_name, pool in pools:
        for gene in pool.index:
            ok, note = check(gene)
            log.append({"pool": pool_name, "gene_id": gene, "gene_name": pool.loc[gene, "gene_name"],
                        "passes": ok, "note": note})
            if ok:
                return gene, log
    return None, log


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    lists = reference_lists()
    primary = pd.read_csv(TRACK_B / "primary" / "de_FSHD_vs_Control.csv.gz", index_col=0)
    ids = set(primary.index.str.split(".").str[0])
    yao = pd.read_csv(LISTS / "Yao2014_TableS1_DUX4_induced_genes.tsv", sep="\t")
    yao["matched_by"] = np.where(yao["ensembl_gene_id"].isin(ids), "ensembl_id",
                                 np.where(yao["source_gene_symbol"].isin(primary["gene_name"]),
                                          "symbol", "not in the counts"))
    yao.to_csv(OUT / "yao2014_matching.csv", index=False)
    lists["yao2014_symbol_fallback"] = set(yao.loc[yao["matched_by"] == "symbol",
                                                   "source_gene_symbol"])
    score = pd.read_csv(TRACK_B / "score" / "de_score_slope.csv.gz", index_col=0)
    g4g1 = pd.read_csv(TRACK_B / "groups" / "de_G4_vs_G1.csv.gz", index_col=0)
    robust_a = pd.read_csv(TRACK_B / "robustness.csv.gz", index_col=0)
    robust_b = pd.read_csv(TRACK_B / "score_robustness.csv.gz", index_col=0)

    names = set(primary["gene_name"].dropna())
    matching = [{"list": key, "symbols": len(values),
                 "found": len(values & ids) if key == "yao2014_ensembl" else len(values & names)}
                for key, values in lists.items() if key != "yao2014_symbol_fallback"]
    pd.DataFrame(matching).to_csv(OUT / "list_matching.csv", index=False)

    # View A.
    a = primary.join(flags(primary, lists))
    a["score_slope"] = score["log2FoldChange"]
    a["score_padj"] = score["padj"]
    a["robust_with_rochester"] = robust_a["robust_with_rochester"].reindex(a.index).fillna(False)
    a["follows_score"] = (a["score_padj"] < ALPHA) & \
        (np.sign(a["score_slope"]) == np.sign(a["log2FoldChange"]))
    eligible_a = a[(a["padj"] < ALPHA) & (a["log2FoldChange"].abs() >= 1)
                   & a["robust_with_rochester"].astype(bool) & a["follows_score"]
                   & a["autosomal"] & a["protein_coding"] & ~a["biomarker"]
                   & ~a["immediate_early"] & ~a["composition"]]
    eligible_a = ranked(eligible_a, "padj", "stat")
    pools_a = [("A_regeneration", eligible_a[eligible_a["regeneration"]]),
               ("A_tf_chromatin", eligible_a[eligible_a["transcription_factor"]
                                             | eligible_a["chromatin"]]),
               ("A_all", eligible_a)]

    influence = Influence()

    def check_a(gene: str) -> tuple[bool, str]:
        sample = a.loc[gene, "max_cooks_sample"]
        lfc = influence.lfc_without(sample).get(gene, np.nan)
        full = a.loc[gene, "log2FoldChange"]
        ok = bool(np.isfinite(lfc) and np.sign(lfc) == np.sign(full) and abs(lfc) >= 0.5 * abs(full))
        return ok, f"without {sample}: log2FC {lfc:.3f} vs {full:.3f}"

    gene_a, log_a = pick(pools_a, check_a)
    pd.DataFrame(log_a).to_csv(OUT / "influence_checks.csv", index=False)

    # View B (Yao 2014 targets stay in the tables, flagged, but cannot be picked).
    b = score.join(flags(score, lists))
    b["robust_b"] = robust_b["robust_b"].reindex(b.index).fillna(False).astype(bool)
    b["primary_log2FC"] = primary["log2FoldChange"]
    b["primary_padj"] = primary["padj"]
    b["g4_vs_g1_log2FC"] = g4g1["log2FoldChange"]
    eligible_b = b[(b["padj"] < ALPHA) & (b["lfc_span"].abs() >= 1) & b["robust_b"]
                   & b["autosomal"] & b["protein_coding"] & ~b["biomarker"]
                   & ~b["immediate_early"] & ~b["composition"]]
    eligible_b = ranked(eligible_b, "padj", "stat")
    pools_b = [("B_chromatin", eligible_b[eligible_b["chromatin"]]),
               ("B_tf", eligible_b[eligible_b["transcription_factor"]])]
    gene_b, log_b = pick(pools_b, lambda g: (not b.loc[g, "yao2014_target"],
                                             "Yao 2014 DUX4-induced gene"
                                             if b.loc[g, "yao2014_target"] else "eligible"))

    skipped = {(r["pool"], r["gene_id"]): r["note"] for r in log_a + log_b}
    cols_a = ["gene_name", "baseMean", "log2FoldChange", "padj", "stat", "score_slope",
              "score_padj", "max_cooks_sample", "transcription_factor", "chromatin",
              "regeneration", "yao2014_target"]
    cols_b = ["gene_name", "baseMean", "log2FoldChange", "lfc_span", "padj", "stat",
              "g4_vs_g1_log2FC", "primary_log2FC", "primary_padj", "transcription_factor",
              "chromatin", "regeneration", "yao2014_target"]
    for pool_name, pool in pools_a + pools_b:
        top = pool.head(TOP)[cols_a if pool_name.startswith("A") else cols_b].copy()
        top.insert(0, "rank", range(1, len(top) + 1))
        top["check"] = [skipped.get((pool_name, g), "not reached") for g in top.index]
        top.to_csv(OUT / f"pool_{pool_name}.csv", float_format="%.6g", index_label="gene_id")

    # Exploration only: does the top 10 move between runs? (no effect on the choice)
    runs_a = {k: pd.read_csv(TRACK_B / v / "de_FSHD_vs_Control.csv.gz", index_col=0)
              for k, v in RUNS_A.items()}
    runs_b = {k: pd.read_csv(TRACK_B / k / "de_score_slope.csv.gz", index_col=0) for k in RUNS_B}
    explore = {"A": top10_by_run(runs_a, "log2FoldChange", lists, eligible_a.index, "primary"),
               "B": top10_by_run(runs_b, "lfc_span", lists, eligible_b.index, "score")}
    for view, table in explore.items():
        table.to_csv(OUT / f"top10_by_run_{view}.csv", float_format="%.6g", index_label="gene_id")

    known = []
    for symbol in KNOWN_REGULATORS:
        ids = primary.index[primary["gene_name"] == symbol]
        for gene in ids:
            known.append({"gene_name": symbol, "gene_id": gene,
                          "primary_log2FC": primary.loc[gene, "log2FoldChange"],
                          "primary_padj": primary.loc[gene, "padj"],
                          "robust_with_rochester": bool(a.loc[gene, "robust_with_rochester"]),
                          "score_slope": score.loc[gene, "log2FoldChange"],
                          "score_lfc_span": score.loc[gene, "lfc_span"],
                          "score_padj": score.loc[gene, "padj"],
                          "g4_vs_g1_log2FC": g4g1.loc[gene, "log2FoldChange"],
                          "g4_vs_g1_padj": g4g1.loc[gene, "padj"]})
        if len(ids) == 0:
            known.append({"gene_name": symbol, "gene_id": "not in the counts"})
    pd.DataFrame(known).to_csv(OUT / "known_regulators.csv", float_format="%.6g", index=False)

    lines = ["list matching (symbols found in our annotation):"]
    lines += [f"  {m['list']}: {m['found']} / {m['symbols']}" for m in matching]
    lines += [f"view A eligible: {len(eligible_a)}; pools: "
              + ", ".join(f"{n} {len(p)}" for n, p in pools_a),
              f"view B eligible: {len(eligible_b)} ({int(eligible_b['yao2014_target'].sum())} "
              "Yao 2014 targets, not pickable); pools: "
              + ", ".join(f"{n} {len(p)}" for n, p in pools_b)]
    for view, table in explore.items():
        outside = table[~table["in_main_top10"]]
        lines.append(f"top 10 by run, view {view}: {len(table)} genes in the union; "
                     f"{len(outside)} not in the main top 10: "
                     + (", ".join(outside["gene_name"].fillna("?")) or "none"))
    if gene_a is not None:
        up = a.loc[gene_a, "log2FoldChange"] > 0
        lines.append(f"view A gene: {a.loc[gene_a, 'gene_name']} ({gene_a}), log2FC "
                     f"{a.loc[gene_a, 'log2FoldChange']:.3f}, padj {a.loc[gene_a, 'padj']:.3g}; "
                     f"Translator: 'What chemicals may {'decrease' if up else 'increase'} the "
                     "activity of a gene?'")
    else:
        lines.append("view A gene: none passed")
    if gene_b is not None:
        up = b.loc[gene_b, "log2FoldChange"] > 0
        lines.append(f"view B gene: {b.loc[gene_b, 'gene_name']} ({gene_b}), log2FC across the "
                     f"span {b.loc[gene_b, 'lfc_span']:.3f}, padj {b.loc[gene_b, 'padj']:.3g}; "
                     f"rises with the DUX4-target score: {up}; optional Translator question "
                     f"with '{'decrease' if up else 'increase'}'")
    else:
        lines.append("view B gene: none (rule not relaxed)")
    (OUT / "selection.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
