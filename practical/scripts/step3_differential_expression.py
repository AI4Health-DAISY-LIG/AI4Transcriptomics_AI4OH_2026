"""Step 3, Track B: independent differential expression with PyDESeq2, on raw counts.

Track A (r/de_reproduction.R) reproduces the authors' tables. This script is the independent analysis,
with the design and the sample rules fixed before any DE result (23 Nov):

  primary          39 biopsies   ~ batch + sex + condition         main result
  ladder_none      39            ~ condition                       worksheet: without batch
  ladder_batch     39            ~ batch + condition               worksheet: with batch; no-sex sensitivity
  site             39            ~ batch + sex + site + condition  site sensitivity
  all43            43            ~ batch + sex + condition         nothing excluded
  restore_32_0004  40            ~ batch + sex + condition         32-0004 put back
  omit_32_0016     38            ~ batch + sex + condition         planned outlier sensitivity
  groups           39            ~ batch + sex + group             G1-G4 vs controls, descriptive
  groups_none      39            ~ group                           worksheet: groups without batch
  rochester        24            ~ batch + sex + condition         Rochester biopsies only (same site)
  groups_g3_check  41            ~ batch + sex + group             32-0008 and 32-0009 put back

rochester and groups_g3_check were added on 23 Nov after the first run of the nine fits above,
before any gene-level result was examined: rochester removes the case/control site difference
(it cannot remove the banking of the controls); groups_g3_check tests whether the low-mapping
libraries 32-0008 and 32-0009 explain the gap with the published G3 count. The groups fit also
gives G4 vs G1 (DUX4-marker-high vs DUX4-marker-low FSHD biopsies), exploratory: the groups
were built from DUX4 targets.

Added on 24 Nov as further checks:
  rochester_norefit   24        as rochester, Cook's replacement off   checks the 15,409 refits
  score               30 FSHD   ~ batch + sex + score                  DUX4-target score slope
  score_nosex         30 FSHD   ~ batch + score
  score_all34         34 FSHD   ~ batch + sex + score
  score_omit_32_0006  29 FSHD   ~ batch + sex + score                  highest score
  score_omit_01_0025  29 FSHD   ~ batch + sex + score                  highest score leverage
The score models give "regulatory genes associated with the DUX4-target score within FSHD"
(bonus view B). Their slope is also expressed as a log2FC across the G1-to-G4 score span
(SCORE_SPAN, medians 0.35207 and 15.33630 in the primary set), on the scale of view A.

The 39 biopsies: all 43 minus config.TRACK_B_EXCLUDE (three libraries mapped below 90%, and
32-0004, whose recorded sex disagrees with its expression). Sex is the sex inferred from
expression (results/step2_qc/sex_check.csv); the GEO labels stay in the sample table. Site comes
from the ID prefix (config.SITE_BY_PREFIX). PyDESeq2 defaults: Wald test, Cook's outlier
refitting and filtering, independent filtering at alpha = 0.05.

Robustness (genes with padj < 0.05 in primary): same sign and at least 50% of the primary
|log2FC| in each planned adjusted sensitivity (ladder_batch, site, all43, restore_32_0004,
omit_32_0016) -> robust_planned; robust_with_rochester adds the rochester fit. A finite
log2FC is required. Agreement with ladder_none and rochester_norefit is reported, not required.
X/Y genes are flagged as not interpretable. Genes significant only without 32-0016 are listed
apart as sensitive to its inclusion, with log2FC, standard error and dispersion in both fits.
The groups_none vs groups difference is the combined effect of batch and sex adjustment.

A fit-free preflight runs first: sample IDs unique and aligned with the counts, covariates
complete, cohort sizes, and the rank of every design matrix.

Inputs:  data/processed/counts_filtered.csv, metadata.csv, gene_names.csv (steps 1 and 2)
         data/raw/gencode.v24.annotation.gtf.gz (chromosomes), results/step2_qc/sex_check.csv
Outputs (results/step3_de/track_b/):
  <fit>/de_<level>_vs_Control.csv.gz   all genes, with gene name, chromosome, type, dispersion
  <fit>/summary.txt                    samples, genes tested, significant genes
  <fit>/ma_*.png, volcano_*.png, pvalue_histogram_*.png   for primary and the two ladder fits
  robustness.csv.gz, sensitive_to_32_0016.csv, sensitivity_comparison.csv, preflight.txt,
  score_robustness.csv.gz (view B: same sign and >= 50% of the slope in score_nosex,
  score_all34 and both omissions; same sign in G4 vs G1), summary.txt
  batch_overlap.csv               genes at padj < 0.05 with and without batch (ladder_none vs
                                  ladder_batch): both, only one, and the totals

Run:  uv run --locked python scripts/step3_differential_expression.py              (all fits)
      uv run --locked python scripts/step3_differential_expression.py --fit primary (one fit)
      uv run --locked python scripts/step3_differential_expression.py --preflight   (no fit)
      uv run --locked python scripts/step3_differential_expression.py --no-fit      (summaries
                                                    and comparisons from the saved fits only)
"""

from __future__ import annotations

import argparse
import gzip
import re

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from pydeseq2.dds import DeseqDataSet  # noqa: E402
from pydeseq2.ds import DeseqStats  # noqa: E402

import config  # noqa: E402

OUT = config.RESULTS / "step3_de" / "track_b"
ALPHA = 0.05
GZIP = {"method": "gzip", "mtime": 0}  # byte-stable gzip output on reruns
REF = config.REFERENCE_LEVEL
MAIN = "~ batch + sex + condition"
SCORE = "~ batch + sex + score"
# Median DUX4-target score of G1 and G4 in the 30 FSHD biopsies of the primary set (GEO
# rnaseq.score): 0.35207 and 15.33630. Frozen; preflight recomputes it as a check.
SCORE_SPAN = 14.98423
FITS = {
    "primary":         {"exclude": config.TRACK_B_EXCLUDE, "design": MAIN},
    "ladder_none":     {"exclude": config.TRACK_B_EXCLUDE, "design": "~ condition"},
    "ladder_batch":    {"exclude": config.TRACK_B_EXCLUDE, "design": "~ batch + condition"},
    "site":            {"exclude": config.TRACK_B_EXCLUDE, "design": "~ batch + sex + site + condition"},
    "all43":           {"exclude": [], "design": MAIN},
    "restore_32_0004": {"exclude": config.SENSITIVITY_EXCLUDE, "design": MAIN},
    "omit_32_0016":    {"exclude": config.TRACK_B_EXCLUDE + ["32-0016"], "design": MAIN},
    "groups":          {"exclude": config.TRACK_B_EXCLUDE, "design": "~ batch + sex + group"},
    "groups_none":     {"exclude": config.TRACK_B_EXCLUDE, "design": "~ group"},
    "rochester":       {"exclude": config.TRACK_B_EXCLUDE, "design": MAIN, "site": "Rochester"},
    "groups_g3_check": {"exclude": ["01-0037", "32-0004"], "design": "~ batch + sex + group"},
    "rochester_norefit": {"exclude": config.TRACK_B_EXCLUDE, "design": MAIN, "site": "Rochester",
                          "refit_cooks": False},
    "score":              {"exclude": config.TRACK_B_EXCLUDE, "design": SCORE, "fshd_only": True},
    "score_nosex":        {"exclude": config.TRACK_B_EXCLUDE, "design": "~ batch + score",
                           "fshd_only": True},
    "score_all34":        {"exclude": [], "design": SCORE, "fshd_only": True},
    "score_omit_32_0006": {"exclude": config.TRACK_B_EXCLUDE + ["32-0006"], "design": SCORE,
                           "fshd_only": True},
    "score_omit_01_0025": {"exclude": config.TRACK_B_EXCLUDE + ["01-0025"], "design": SCORE,
                           "fshd_only": True},
}
EXTRA_CONTRASTS = {"groups": [("G4", "G1")]}
PLANNED_FITS = ["ladder_batch", "site", "all43", "restore_32_0004", "omit_32_0016"]
ROBUSTNESS_FITS = PLANNED_FITS + ["rochester"]
INFO_FITS = ["ladder_none", "rochester_norefit"]
SCORE_CHECKS = ["score_nosex", "score_all34", "score_omit_32_0006", "score_omit_01_0025"]
SCORE_TABLE = "de_score_slope.csv.gz"
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
PLOT_FITS = {"primary", "ladder_none", "ladder_batch"}
TABLE = "de_FSHD_vs_Control.csv.gz"


def gene_annotation() -> pd.DataFrame:
    """gene_name and gene_type from step 1, chromosome from the GENCODE v24 GTF."""
    annot = (pd.read_csv(config.GENE_NAMES).drop_duplicates("gene_id")
             .set_index("gene_id")[["gene_name", "gene_type"]])
    chrom = {}
    with gzip.open(config.RAW / config.ANNOTATION_FILE, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 9)
            if fields[2] == "gene":
                chrom[fields[8].split('gene_id "', 1)[1].split('"', 1)[0]] = fields[0]
    annot["chromosome"] = pd.Series(chrom)
    return annot


def sample_table() -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(config.PROCESSED / "counts_filtered.csv", index_col=0)
    meta = pd.read_csv(config.PROCESSED / "metadata.csv", index_col=0)
    sex = pd.read_csv(config.RESULTS / "step2_qc" / "sex_check.csv", index_col=0)["inferred_sex"]
    for label, ids in [("counts", counts.columns), ("metadata", meta.index),
                       ("sex_check", sex.index)]:
        if not ids.is_unique or set(ids) != set(counts.columns):
            raise SystemExit(f"{label}: sample IDs not unique or not the same 43 as the counts")
    meta = meta.loc[counts.columns]
    meta = pd.DataFrame({
        "condition": meta["condition"].astype(str),
        "batch": meta["char:batch"].astype(str),
        "sex": sex.reindex(meta.index),
        "site": meta.index.str[:2].map(config.SITE_BY_PREFIX),
        "group": pd.to_numeric(meta[config.GROUP_COL], errors="coerce"),
        "score": pd.to_numeric(meta[config.SCORE_COL], errors="coerce"),
    }, index=meta.index)
    meta["group"] = ("G" + meta["group"].astype("Int64").astype(str)).where(meta["group"].notna())
    meta.loc[meta["condition"] == REF, "group"] = REF
    missing = meta.drop(columns="score").isna().any(axis=1) | \
        (meta["score"].isna() & (meta["condition"] != REF))
    if missing.any():
        raise SystemExit(f"missing covariates:\n{meta[missing]}")
    return counts.T.round().astype(int), meta


def fit_samples(meta: pd.DataFrame, spec: dict) -> pd.DataFrame:
    """The biopsies of one fit: all minus its exclusions, restricted to one site if given."""
    m = meta.drop(index=spec["exclude"])
    if spec.get("fshd_only"):
        m = m[m["condition"] != REF]
    return m[m["site"] == spec["site"]] if "site" in spec else m


def design_matrix(m: pd.DataFrame, variables: list[str]) -> pd.DataFrame:
    """Intercept plus treatment-coded dummies (first level as reference), as the fit uses."""
    dummies = [m[v].astype(float) if pd.api.types.is_numeric_dtype(m[v]) else
               pd.get_dummies(m[v], prefix=v, drop_first=True, dtype=float) for v in variables]
    return pd.concat([pd.Series(1.0, index=m.index, name="intercept")] + dummies, axis=1)


def describe(m: pd.DataFrame, factor: str) -> str:
    if factor == "score":
        return f"FSHD={len(m)}, score {m['score'].min():.2f}-{m['score'].max():.2f}"
    return ", ".join(f"{k}={v}" for k, v in m[factor].value_counts().sort_index().items())


def preflight(x: pd.DataFrame, meta: pd.DataFrame) -> list[str]:
    """Cohort sizes and design ranks for every fit, without fitting anything."""
    fshd = fit_samples(meta, FITS["score"])
    medians = fshd.groupby("group")["score"].median()
    lines = [f"genes in the input: {x.shape[1]}; samples: {x.shape[0]}",
             f"score span G1-G4 medians: {medians['G4'] - medians['G1']:.5f} "
             f"(frozen SCORE_SPAN {SCORE_SPAN})"]
    for name, spec in FITS.items():
        m = fit_samples(meta, spec)
        variables = re.findall(r"[a-z_]+", spec["design"])
        mat = design_matrix(m, variables)
        rank = np.linalg.matrix_rank(mat.to_numpy())
        zero = int((x.loc[m.index].sum(axis=0) == 0).sum())
        sizes = describe(m, variables[-1])
        lines.append(f"{name}: n={len(m)} ({sizes}); design rank {rank}/{mat.shape[1]}; "
                     f"all-zero genes {zero}")
        if rank < mat.shape[1]:
            raise SystemExit("\n".join(lines) + f"\n{name}: design matrix not full rank")
    return lines


def plots(res: pd.DataFrame, sig: pd.DataFrame, tag: str, out) -> None:
    label = res["gene_name"].fillna(res.index.to_series())
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(np.log10(res["baseMean"] + 1), res["log2FoldChange"], s=4, c="lightgrey")
    ax.scatter(np.log10(sig["baseMean"] + 1), sig["log2FoldChange"], s=4, c="tab:red")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("log10(base mean + 1)")
    ax.set_ylabel(f"log2 fold change {tag}")
    ax.set_title(f"MA plot ({len(sig)} genes with padj < {ALPHA})")
    fig.tight_layout()
    fig.savefig(out / f"ma_{tag}.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(res["log2FoldChange"], -np.log10(res["pvalue"]), s=4, c="lightgrey")
    ax.scatter(sig["log2FoldChange"], -np.log10(sig["pvalue"]), s=4, c="tab:red")
    for gene in sig.index[:15]:
        xy = (sig.loc[gene, "log2FoldChange"], -np.log10(sig.loc[gene, "pvalue"]))
        ax.annotate(str(label.loc[gene]), xy, fontsize=6)
    ax.set_xlabel(f"log2 fold change {tag}")
    ax.set_ylabel("-log10 p-value")
    fig.tight_layout()
    fig.savefig(out / f"volcano_{tag}.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(res["pvalue"].dropna(), bins=50)
    ax.set_xlabel("p-value")
    ax.set_ylabel("genes")
    ax.set_title(tag)
    fig.tight_layout()
    fig.savefig(out / f"pvalue_histogram_{tag}.png", dpi=150)
    plt.close(fig)


def gene_vector(dds: DeseqDataSet, key: str, default=None) -> pd.Series:
    """A per-gene result of the fit: PyDESeq2 0.5 keeps these in dds.var, older ones in varm."""
    if key in dds.var.columns:
        return dds.var[key]
    if key in dds.varm.keys():
        return pd.Series(np.asarray(dds.varm[key]), index=dds.var_names)
    if default is not None:
        return pd.Series(default, index=dds.var_names)
    raise SystemExit(f"{key} not found; var: {list(dds.var.columns)}; varm: {list(dds.varm.keys())}")


def run_fit(name: str, x: pd.DataFrame, meta: pd.DataFrame, annot: pd.DataFrame) -> None:
    spec = FITS[name]
    out = OUT / name
    out.mkdir(parents=True, exist_ok=True)
    m = fit_samples(meta, spec)
    variables = re.findall(r"[a-z_]+", spec["design"])
    factor = variables[-1]
    dds = DeseqDataSet(counts=x.loc[m.index], metadata=m[variables], design=spec["design"],
                       refit_cooks=spec.get("refit_cooks", True), quiet=True)
    dds.deseq2()
    dispersion = gene_vector(dds, "dispersions")
    if "cooks" not in dds.layers.keys():
        raise SystemExit(f"no Cook's distances; layers: {list(dds.layers.keys())}")
    cooks = pd.DataFrame(np.asarray(dds.layers["cooks"]), index=dds.obs_names,
                         columns=dds.var_names)
    max_cooks_sample = cooks.fillna(-np.inf).idxmax(axis=0)
    refitted = int(gene_vector(dds, "refitted", default=False).astype(bool).sum())

    lines = [f"[{name}] design: {spec['design']}  samples: {len(m)}  excluded: "
             f"{', '.join(spec['exclude']) or 'none'}"
             + (f"; {spec['site']} only" if "site" in spec else "")
             + ("; FSHD only" if spec.get("fshd_only") else "")
             + ("; Cook's replacement off" if spec.get("refit_cooks") is False else ""),
             f"  {describe(m, factor)}; genes refitted after Cook's outlier replacement: {refitted}"]
    if factor == "score":
        columns = list(dds.obsm["design_matrix"].columns)
        if "score" not in columns:
            raise SystemExit(f"no score column in the design matrix: {columns}")
        vector = np.zeros(len(columns))
        vector[columns.index("score")] = 1.0
        contrasts = [("score_slope", vector)]
    else:
        pairs = [(level, REF) for level in sorted(set(m[factor]) - {REF})]
        contrasts = [(f"{level}_vs_{ref}", [factor, level, ref])
                     for level, ref in pairs + EXTRA_CONTRASTS.get(name, [])]
    for tag, contrast in contrasts:
        stats = DeseqStats(dds, contrast=contrast, alpha=ALPHA, quiet=True)
        stats.summary()
        res = annot.reindex(stats.results_df.index).join(stats.results_df)
        if factor == "score":
            res["lfc_span"] = res["log2FoldChange"] * SCORE_SPAN
        res["dispersion"] = dispersion
        res["max_cooks_sample"] = max_cooks_sample
        res = res.sort_values("padj")
        res.to_csv(out / f"de_{tag}.csv.gz", float_format="%.6g", index_label="gene_id")
        sig = res[res["padj"] < ALPHA]
        n_up = int((sig["log2FoldChange"] > 0).sum())
        cooks_na = int((res["pvalue"].isna() & (res["baseMean"] > 0)).sum())
        filtered = int((res["pvalue"].notna() & res["padj"].isna()).sum())
        lines.append(
            f"  {tag}: finite p-value {int(res['pvalue'].notna().sum())}, finite padj "
            f"{int(res['padj'].notna().sum())}, Cook's-filtered {cooks_na}, independently "
            f"filtered {filtered}; padj < {ALPHA}: {len(sig)} "
            f"(up {n_up}, down {len(sig) - n_up}); also |log2FC| >= 1: "
            f"{int((sig['log2FoldChange'].abs() >= 1).sum())}, >= 2: "
            f"{int((sig['log2FoldChange'].abs() >= 2).sum())}"
            + (f"; |log2FC across the G1-G4 span| >= 1: {int((sig['lfc_span'].abs() >= 1).sum())}"
               if factor == "score" else ""))
        if name in PLOT_FITS:
            plots(res, sig, tag, out)
    (out / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


def robustness() -> list[str]:
    needed = ["primary"] + ROBUSTNESS_FITS + INFO_FITS
    missing = [f for f in needed if not (OUT / f / TABLE).exists()]
    if missing:
        return [f"robustness: skipped, fits not run yet: {', '.join(missing)}"]
    tables = {f: pd.read_csv(OUT / f / TABLE, index_col=0) for f in needed}
    primary = tables["primary"]
    sig = primary[primary["padj"] < ALPHA]
    table = sig[["gene_name", "chromosome", "gene_type", "baseMean", "log2FoldChange", "padj"]].copy()
    table = table.rename(columns={"log2FoldChange": "lfc_primary", "padj": "padj_primary"})
    for f in ROBUSTNESS_FITS + INFO_FITS:
        lfc = tables[f]["log2FoldChange"].reindex(sig.index)
        table[f"lfc_{f}"] = lfc
        ok = lfc.notna() & (np.sign(lfc) == np.sign(sig["log2FoldChange"])) & \
             (lfc.abs() >= 0.5 * sig["log2FoldChange"].abs())
        table[f"holds_{f}"] = ok
    table["robust_planned"] = table[[f"holds_{f}" for f in PLANNED_FITS]].all(axis=1)
    table["robust_with_rochester"] = table["robust_planned"] & table["holds_rochester"]
    table["autosomal"] = table["chromosome"].isin(AUTOSOMES)
    table["biomarker"] = table["gene_name"].isin(config.BIOMARKERS)
    table.to_csv(OUT / "robustness.csv.gz", float_format="%.6g", index_label="gene_id",
                 compression=GZIP)

    omit = tables["omit_32_0016"]
    only = omit.index[(omit["padj"] < ALPHA) & ~(primary["padj"].reindex(omit.index) < ALPHA)]
    cols = ["log2FoldChange", "lfcSE", "dispersion", "padj"]
    masked = (omit.loc[only, ["gene_name", "chromosome"] + cols]
              .join(primary.loc[only, cols], lsuffix="_without_32_0016", rsuffix="_primary"))
    masked.sort_values("padj_without_32_0016").to_csv(OUT / "sensitive_to_32_0016.csv",
                                                      float_format="%.6g", index_label="gene_id")
    comparison = []
    for f in ROBUSTNESS_FITS + INFO_FITS:
        other = tables[f]
        in_primary = primary["padj"] < ALPHA
        in_other = other["padj"].reindex(primary.index) < ALPHA
        paired = sig.index
        comparison.append({
            "fit": f,
            "finite_padj": int(other["padj"].notna().sum()),
            "significant": int((other["padj"] < ALPHA).sum()),
            "both": int((in_primary & in_other).sum()),
            "primary_only": int((in_primary & ~in_other).sum()),
            "fit_only": int((~in_primary & in_other).sum()),
            "lfc_spearman_primary_genes": round(float(sig["log2FoldChange"].corr(
                other["log2FoldChange"].reindex(paired), method="spearman")), 4),
            "median_lfc_ratio": round(float((other["log2FoldChange"].reindex(paired)
                                            / sig["log2FoldChange"]).median()), 3),
            "median_se_ratio": round(float((other["lfcSE"].reindex(paired)
                                           / sig["lfcSE"]).median()), 3),
        })
    comparison = pd.DataFrame(comparison)
    comparison.to_csv(OUT / "sensitivity_comparison.csv", index=False)

    autosomal = table["autosomal"]
    return [f"robustness: {len(table)} genes with padj < {ALPHA} in primary; robust (planned five) "
            f"{int(table['robust_planned'].sum())} "
            f"({int((table['robust_planned'] & autosomal).sum())} on chr1-22); with Rochester "
            f"{int(table['robust_with_rochester'].sum())} "
            f"({int((table['robust_with_rochester'] & autosomal).sum())} on chr1-22); "
            f"holds without sex {int(table['holds_ladder_batch'].sum())}, site "
            f"{int(table['holds_site'].sum())}, all43 {int(table['holds_all43'].sum())}, 32-0004 "
            f"restored {int(table['holds_restore_32_0004'].sum())}, 32-0016 omitted "
            f"{int(table['holds_omit_32_0016'].sum())}, Rochester only "
            f"{int(table['holds_rochester'].sum())}; reported only: ~ condition "
            f"{int(table['holds_ladder_none'].sum())}, Rochester without Cook's replacement "
            f"{int(table['holds_rochester_norefit'].sum())}",
            f"sensitive to 32-0016 (padj < {ALPHA} only without it): {len(masked)}",
            "paired comparison with primary (genes significant in primary for ratios):",
            comparison.to_string(index=False)]


def batch_overlap() -> list[str]:
    """Worksheet check: genes at padj < 0.05 without batch (ladder_none, ~ condition) and with
    batch (ladder_batch, ~ batch + condition), on the same 39 biopsies."""
    paths = {f: OUT / f / TABLE for f in ["ladder_none", "ladder_batch"]}
    missing = [f for f, path in paths.items() if not path.exists()]
    if missing:
        return [f"batch overlap: skipped, not run yet: {', '.join(missing)}"]
    sig = {f: set(pd.read_csv(path, index_col=0).query("padj < @ALPHA").index)
           for f, path in paths.items()}
    none, batch = sig["ladder_none"], sig["ladder_batch"]
    row = {"without_batch": len(none), "with_batch": len(batch), "both": len(none & batch),
           "only_without_batch": len(none - batch), "only_with_batch": len(batch - none)}
    pd.DataFrame([row]).to_csv(OUT / "batch_overlap.csv", index=False)
    return ["batch overlap (padj < 0.05, ~ condition vs ~ batch + condition): "
            + ", ".join(f"{k} {v}" for k, v in row.items())]


def score_robustness() -> list[str]:
    """View B: genes with padj < 0.05 for the score slope, checked in the four score
    sensitivities (same sign, >= 50% of the slope) and against G4 vs G1 (same sign)."""
    paths = {f: OUT / f / SCORE_TABLE for f in ["score"] + SCORE_CHECKS}
    paths["g4_vs_g1"] = OUT / "groups" / "de_G4_vs_G1.csv.gz"
    missing = [f for f, path in paths.items() if not path.exists()]
    if missing:
        return [f"score robustness: skipped, not run yet: {', '.join(missing)}"]
    tables = {f: pd.read_csv(path, index_col=0) for f, path in paths.items()}
    slopes = tables["score"]
    sig = slopes[slopes["padj"] < ALPHA]
    table = sig[["gene_name", "chromosome", "gene_type", "baseMean", "log2FoldChange",
                 "lfc_span", "padj"]].rename(columns={"log2FoldChange": "slope"}).copy()
    for f in SCORE_CHECKS:
        slope = tables[f]["log2FoldChange"].reindex(sig.index)
        table[f"slope_{f}"] = slope
        table[f"holds_{f}"] = slope.notna() & (np.sign(slope) == np.sign(sig["log2FoldChange"])) \
            & (slope.abs() >= 0.5 * sig["log2FoldChange"].abs())
    g4 = tables["g4_vs_g1"]["log2FoldChange"].reindex(sig.index)
    table["lfc_g4_vs_g1"] = g4
    table["holds_g4_vs_g1"] = g4.notna() & (np.sign(g4) == np.sign(sig["log2FoldChange"]))
    table["robust_b"] = table[[c for c in table.columns if c.startswith("holds_")]].all(axis=1)
    table["autosomal"] = table["chromosome"].isin(AUTOSOMES)
    table.to_csv(OUT / "score_robustness.csv.gz", float_format="%.6g", index_label="gene_id",
                 compression=GZIP)
    big = table["lfc_span"].abs() >= 1
    return [f"score (view B): {len(table)} genes with padj < {ALPHA}; robust "
            f"{int(table['robust_b'].sum())}; robust with |log2FC across the span| >= 1 "
            f"{int((table['robust_b'] & big).sum())}; holds without sex "
            f"{int(table['holds_score_nosex'].sum())}, all 34 {int(table['holds_score_all34'].sum())}"
            f", without 32-0006 {int(table['holds_score_omit_32_0006'].sum())}, without 01-0025 "
            f"{int(table['holds_score_omit_01_0025'].sum())}, same sign in G4 vs G1 "
            f"{int(table['holds_g4_vs_g1'].sum())}"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--fit", action="append", choices=list(FITS),
                        help="run only this fit (repeatable); default: all fits")
    parser.add_argument("--preflight", action="store_true", help="checks only, no fit")
    parser.add_argument("--no-fit", action="store_true",
                        help="rebuild summaries and comparisons from the saved fits only")
    args = parser.parse_args()
    fits = args.fit or list(FITS)
    OUT.mkdir(parents=True, exist_ok=True)
    if not args.no_fit:  # --no-fit reads only the saved fit tables, not the counts
        x, meta = sample_table()
        checks = preflight(x, meta)
        (OUT / "preflight.txt").write_text("\n".join(checks) + "\n", encoding="utf-8")
        print("\n".join(checks))
        if args.preflight:
            return
        annot = gene_annotation()
        for name in fits:
            run_fit(name, x, meta, annot)
    lines = [(OUT / f / "summary.txt").read_text(encoding="utf-8").rstrip()
             for f in FITS if (OUT / f / "summary.txt").exists()]
    robustness_lines = robustness() + batch_overlap() + score_robustness()
    lines += robustness_lines
    (OUT / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(robustness_lines))


if __name__ == "__main__":
    main()
