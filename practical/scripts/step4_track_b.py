"""Step 4, Track B: independent analysis of the DUX4 score and of the Translator panel genes.

Not the reproduction (r/score_reproduction.R). Three parts, all descriptive:

1. Alternative DUX4 score. Blind DESeq2 VST of the full matrix (r/vst_track_b.R); per gene, the
   value minus its median in the 9 controls; the score is the mean over a panel:
     S1        Yao et al. 2014 Table S1 DUX4-induced genes, as mapped for the bonus
               (results/bonus/translator/yao2014_matching.csv)
     S3        the authors' own mapping of Yao Table S3 (data/raw/authors/mk.rda, via the R step)
     *_lfo     the same without every feature of LEUTX, KHDC1L, PRAMEF2 and TRIM43 (the four
               markers that define the published score), by deposited or source symbol;
               S1_lfo is the headline score
     *_z       sensitivity: floored z-scores, each control SD floored at the median of the positive
             control SDs of the panel
   Each score is compared with the deposited GEO score: Spearman with MRI and pathology on the
   endpoint-complete FSHD biopsies, and the paired bootstrap of the difference in rho.
2. Tissue state of the Translator panel: the genes of every committed top-10 list and pool of the
   bonus selection that are not known DUX4 targets, regulators or members of their families,
   minus those failing the robustness rule of their view, plus POSTN (chosen from the
   literature). Spearman with MRI, pathology and the GEO score; partial correlation given the
   GEO score (Pearson correlation of rank residuals). Groups, fixed from Step 3 before this
   analysis: "follows the DUX4 level" (score padj < 0.05, G4 vs G1 padj < 0.05, and the
   same sign as FSHD vs controls when that is significant) and "case-control only" (the
   other panel genes: selected through FSHD vs controls, outside that joint criterion; the
   label says nothing about expression being uniform across FSHD muscles).
3. RUNX1-POSTN and SOX11-POSTN expression correlations.

Bootstrap: 5,000 resamples of the FSHD biopsies with replacement, seed 20261124, the same draws
for every score and gene; pointwise 95% percentile intervals; resamples where a correlation is
undefined are dropped and counted. Nothing here is adjusted for multiple testing.

Run, from practical/:
  uv run --locked python scripts/step4_track_b.py --panel      (writes panel_genes.csv)
  pixi run --manifest-path r/pixi.toml Rscript r/vst_track_b.R
  uv run --locked python scripts/step4_track_b.py
Outputs (results/step4_comparison/track_b/):
  panel_genes.csv, scores.csv, panels.csv, score_correlations.csv, delta_rho.csv,
  tissue_state.csv, postn_correlations.csv, scores_by_batch.png, postn_pairs.png, summary.txt
"""

from __future__ import annotations

import argparse

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import config  # noqa: E402

OUT = config.RESULTS / "step4_comparison" / "track_b"
DE = config.RESULTS / "step3_de" / "track_b"
BONUS = config.RESULTS / "bonus" / "translator"
ALPHA = 0.05
N_BOOT = 5000
SEED = 20261124
ENDPOINTS = {"stir": "char:stir", "t1": "char:t1", "fat_fraction": "char:fat.fraction",
             "pathology": "char:path.score", "active": "char:active",
             "inflammation": "char:inflam"}
SCORE_ENDPOINTS = ["stir", "t1", "fat_fraction", "pathology"]   # the paper's four
LOW_MAPPING = config.SENSITIVITY_EXCLUDE
FAMILIES = ("PRAMEF", "HNRNPCL")
LITERATURE = {"POSTN": "periostin, tissue-remodelling indicator (literature, after selection)"}


# --- Panel -----------------------------------------------------------------------------------
def build_panel() -> pd.DataFrame:
    sources: dict[str, list[str]] = {}
    for name in ["top10_by_run_A", "top10_by_run_B", "pool_A_all", "pool_A_tf_chromatin",
                 "pool_A_regeneration", "pool_B_tf", "pool_B_chromatin"]:
        table = pd.read_csv(BONUS / f"{name}.csv")
        for gene in table["gene_id"]:
            sources.setdefault(gene, []).append(name)
    yao = set(pd.read_csv(BONUS / "yao2014_matching.csv")["feature"].dropna())
    known = set(pd.read_csv(BONUS / "known_regulators.csv")["gene_id"])
    primary = pd.read_csv(DE / "primary" / "de_FSHD_vs_Control.csv.gz", index_col=0)
    score = pd.read_csv(DE / "score" / "de_score_slope.csv.gz", index_col=0)
    g4g1 = pd.read_csv(DE / "groups" / "de_G4_vs_G1.csv.gz", index_col=0)
    robust_a = pd.read_csv(DE / "robustness.csv.gz", index_col=0)["robust_with_rochester"]
    robust_b = pd.read_csv(DE / "score_robustness.csv.gz", index_col=0)["robust_b"]
    for symbol, why in LITERATURE.items():
        sources.setdefault(primary.index[primary["gene_name"] == symbol][0], []).append(why)

    rows = []
    for gene, where in sources.items():
        name = primary.loc[gene, "gene_name"]
        p_sig = primary.loc[gene, "padj"] < ALPHA
        if gene in yao or gene in known or name.startswith(FAMILIES):
            status = "left out: known DUX4 target, regulator or family member"
        elif p_sig and not robust_a.get(gene, False):
            status = "left out: fails the View A robustness rule"
        elif not p_sig and gene in robust_b.index and not robust_b[gene]:
            status = "left out: fails the View B robustness rule"
        else:
            status = "panel"
        follows = (score.loc[gene, "padj"] < ALPHA and g4g1.loc[gene, "padj"] < ALPHA
                   and (not p_sig or np.sign(score.loc[gene, "log2FoldChange"])
                        == np.sign(primary.loc[gene, "log2FoldChange"])))
        rows.append({"gene_id": gene, "gene_name": name, "status": status,
                     "group": "follows the DUX4 level" if follows else "case-control only",
                     "low_count": primary.loc[gene, "baseMean"] < 20,
                     "primary_log2FC": primary.loc[gene, "log2FoldChange"],
                     "primary_padj": primary.loc[gene, "padj"],
                     "score_lfc_span": score.loc[gene, "lfc_span"],
                     "score_padj": score.loc[gene, "padj"],
                     "g4_vs_g1_padj": g4g1.loc[gene, "padj"],
                     "sources": "; ".join(where)})
    panel = pd.DataFrame(rows).sort_values(["status", "group", "gene_name"])
    OUT.mkdir(parents=True, exist_ok=True)
    panel.to_csv(OUT / "panel_genes.csv", index=False, float_format="%.6g")
    kept = panel[panel["status"] == "panel"]
    print(f"{len(panel)} genes in the lists; panel {len(kept)} "
          f"({(kept['group'] == 'follows the DUX4 level').sum()} follow the DUX4 level); "
          + "; ".join(f"{k}: {v}" for k, v in panel["status"].value_counts().items()))
    return panel


# --- Helpers ---------------------------------------------------------------------------------
def ranks(x: np.ndarray) -> np.ndarray:
    """Average ranks within each row (rows are resamples)."""
    return pd.DataFrame(x).rank(axis=1).to_numpy()


def corr_rows(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Pearson correlation row by row; NaN where a row is constant."""
    a = a - a.mean(axis=1, keepdims=True)
    b = b - b.mean(axis=1, keepdims=True)
    den = np.sqrt((a ** 2).sum(axis=1) * (b ** 2).sum(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, (a * b).sum(axis=1) / den, np.nan)


def residual_rows(y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Residuals of y on z with an intercept, row by row."""
    zc = z - z.mean(axis=1, keepdims=True)
    yc = y - y.mean(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        beta = (zc * yc).sum(axis=1, keepdims=True) / (zc ** 2).sum(axis=1, keepdims=True)
    return yc - beta * zc


def _spearman(xs: np.ndarray, ys: np.ndarray, gs: np.ndarray | None) -> np.ndarray:
    """Row-wise Spearman of 2-D arrays; partial given gs: ranks are taken within each row,
    then x and y ranks are residualised on the gs ranks."""
    rx, ry = ranks(xs), ranks(ys)
    if gs is not None:
        rg = ranks(gs)
        rx, ry = residual_rows(rx, rg), residual_rows(ry, rg)
    return corr_rows(rx, ry)


def spearman(x: np.ndarray, y: np.ndarray, given: np.ndarray | None = None) -> float:
    """Spearman (partial given `given`) on one sample."""
    return float(_spearman(x[None], y[None], None if given is None else given[None])[0])


def spearman_boot(x: np.ndarray, y: np.ndarray, idx: np.ndarray,
                  given: np.ndarray | None = None) -> np.ndarray:
    """The same on each resample: idx holds positions into x, y and given, one row per draw."""
    return _spearman(x[idx], y[idx], None if given is None else given[idx])


def interval(values: np.ndarray) -> tuple[float, float, int]:
    ok = values[np.isfinite(values)]
    return float(np.percentile(ok, 2.5)), float(np.percentile(ok, 97.5)), int(len(values) - len(ok))


def label(lo: float, hi: float) -> str:
    return "weaker" if hi < 0 else "stronger" if lo > 0 else "difference unresolved"


# --- Analysis --------------------------------------------------------------------------------
def load() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    vst = pd.read_csv(OUT / "vst_selected.csv.gz", index_col=0)
    features = pd.read_csv(OUT / "vst_features.csv", index_col=0)
    meta = pd.read_csv(config.PROCESSED / "metadata.csv", index_col=0).loc[vst.columns]
    sex = pd.read_csv(config.RESULTS / "step2_qc" / "sex_check.csv", index_col=0)
    meta["batch"] = meta["char:batch"]
    meta["inferred_sex"] = sex["inferred_sex"].reindex(meta.index)
    return vst, features, meta


def scores(vst: pd.DataFrame, features: pd.DataFrame, control: pd.Series) -> tuple[pd.DataFrame, list[dict]]:
    centred = vst.sub(vst.loc[:, control].median(axis=1), axis=0)
    sd = vst.loc[:, control].std(axis=1, ddof=1)
    # Marker features: the deposited symbol, or the symbol in Yao's source table (Yao lists
    # KHDC1L under two IDs; ENSG00000243501 carries the deposited symbol RP11-257K9.8).
    yao = pd.read_csv(BONUS / "yao2014_matching.csv")
    source_markers = set(yao.loc[yao["source_gene_symbol"].isin(config.BIOMARKERS), "feature"])
    markers = features["gene_name"].isin(config.BIOMARKERS) | features.index.isin(source_markers)
    out, info = {}, []
    for panel in ["S1", "S3"]:
        members = features.index[features[f"in_{panel}"]]
        floor = sd.loc[members][sd.loc[members] > 0].median()     # per panel
        z = centred.div(sd.clip(lower=floor), axis=0)
        for suffix, keep in [("", members), ("_lfo", members[~markers.loc[members]])]:
            out[f"{panel}{suffix}"] = centred.loc[keep].mean(axis=0)
            out[f"{panel}{suffix}_z"] = z.loc[keep].mean(axis=0)
            info.append({"score": f"{panel}{suffix}", "features": len(keep),
                         "markers_removed": int(len(members) - len(keep)),
                         "floored_sd": int((sd.loc[keep] < floor).sum()),
                         "sd_floor": round(float(floor), 5)})
    return pd.DataFrame(out), info


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--panel", action="store_true", help="write panel_genes.csv only")
    if parser.parse_args().panel:
        build_panel()
        return

    vst, features, meta = load()
    control = meta["condition"] == config.REFERENCE_LEVEL
    fshd = ~control
    geo = pd.to_numeric(meta[config.SCORE_COL], errors="coerce")
    table, info = scores(vst, features, control)
    table["geo"] = geo
    table.insert(0, "batch", meta["batch"])
    table.insert(1, "condition", meta["condition"])
    table.to_csv(OUT / "scores.csv", float_format="%.6g")
    names = [c for c in table.columns if c not in ("batch", "condition", "geo")]

    panels = pd.DataFrame(info)
    top_control = table.loc[control, names].max()
    panels["fshd_above_highest_control"] = [int((table.loc[fshd, s["score"]]
                                                 > top_control[s["score"]]).sum()) for s in info]
    panels["spearman_with_geo_fshd"] = [round(table.loc[fshd, s["score"]].corr(
        geo[fshd], method="spearman"), 4) for s in info]
    panels.to_csv(OUT / "panels.csv", index=False)

    rng = np.random.default_rng(SEED)
    n_fshd = int(fshd.sum())
    fshd_ids = meta.index[fshd]
    subsets = {"all": np.ones(n_fshd, bool),
               "without_low_mapping": ~fshd_ids.isin(LOW_MAPPING),
               "without_32_0004": fshd_ids != "32-0004"}

    geo_f = geo[fshd].to_numpy()

    # 1. Scores against MRI and pathology, and the change in rho (same draws for every score).
    corr_rows_out, delta_rows = [], []
    for endpoint in SCORE_ENDPOINTS:
        y_all = pd.to_numeric(meta.loc[fshd, ENDPOINTS[endpoint]], errors="coerce").to_numpy()
        for sub_name, sub in subsets.items():
            keep = sub & np.isfinite(y_all)
            pos = np.flatnonzero(keep)
            idx = pos[rng.integers(0, len(pos), size=(N_BOOT, len(pos)))]
            y = y_all
            base_full = spearman(geo_f[pos], y[pos])
            base_boot = spearman_boot(geo_f, y, idx)
            for s in ["geo"] + names:
                x = table.loc[fshd, s].to_numpy()
                full = spearman(x[pos], y[pos])
                boot = spearman_boot(x, y, idx)
                lo, hi, bad = interval(boot)
                corr_rows_out.append({"endpoint": endpoint, "subset": sub_name, "score": s,
                                      "n": len(pos), "rho": round(full, 4), "ci_low": round(lo, 4),
                                      "ci_high": round(hi, 4), "undefined_resamples": bad})
                if s != "geo":
                    d = boot - base_boot
                    dlo, dhi, dbad = interval(d)
                    delta_rows.append({"endpoint": endpoint, "subset": sub_name, "score": s,
                                       "n": len(pos), "delta_rho": round(full - base_full, 4),
                                       "ci_low": round(dlo, 4), "ci_high": round(dhi, 4),
                                       "label": label(dlo, dhi), "undefined_resamples": dbad})
    corr = pd.DataFrame(corr_rows_out)
    delta = pd.DataFrame(delta_rows)
    corr.to_csv(OUT / "score_correlations.csv", index=False)
    delta.to_csv(OUT / "delta_rho.csv", index=False)

    # 2. Tissue state of the panel genes.
    panel = pd.read_csv(OUT / "panel_genes.csv")
    panel = panel[panel["status"] == "panel"].set_index("gene_id")
    track_b = ~fshd_ids.isin(config.TRACK_B_EXCLUDE)
    tissue = []
    for set_name, sub in [("34 FSHD", np.ones(n_fshd, bool)), ("30 Track B", track_b)]:
        for endpoint, col in list(ENDPOINTS.items()) + [("dux4_score_geo", None)]:
            y_all = geo_f if col is None else \
                pd.to_numeric(meta.loc[fshd, col], errors="coerce").to_numpy()
            pos = np.flatnonzero(sub & np.isfinite(y_all))
            idx = pos[rng.integers(0, len(pos), size=(N_BOOT, len(pos)))]
            for gene, row in panel.iterrows():
                x = vst.loc[gene, fshd_ids].to_numpy()
                full = spearman(x[pos], y_all[pos])
                lo, hi, _ = interval(spearman_boot(x, y_all, idx))
                rec = {"set": set_name, "gene": row["gene_name"], "group": row["group"],
                       "low_count": row["low_count"], "variable": endpoint, "n": len(pos),
                       "rho": round(full, 4), "ci_low": round(lo, 4), "ci_high": round(hi, 4)}
                if col is not None:
                    pf = spearman(x[pos], y_all[pos], geo_f[pos])
                    plo, phi, _ = interval(spearman_boot(x, y_all, idx, geo_f))
                    rec.update({"partial_rho_given_geo": round(pf, 4),
                                "partial_ci_low": round(plo, 4), "partial_ci_high": round(phi, 4)})
                tissue.append(rec)
    tissue = pd.DataFrame(tissue)
    tissue.to_csv(OUT / "tissue_state.csv", index=False)

    # 3. POSTN pairs.
    by_name = {r["gene_name"]: g for g, r in panel.iterrows()}
    pairs = []
    pair_draws = {}                                  # the same draws for both pairs
    for set_name, sub in [("34 FSHD", np.ones(n_fshd, bool)), ("30 Track B", track_b)]:
        pos = np.flatnonzero(sub)
        pair_draws[set_name] = (pos, pos[rng.integers(0, len(pos), size=(N_BOOT, len(pos)))])
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for ax, other in zip(axes, ["RUNX1", "SOX11"]):
        x = vst.loc[by_name[other], fshd_ids].to_numpy()
        y = vst.loc[by_name["POSTN"], fshd_ids].to_numpy()
        for set_name, (pos, idx) in pair_draws.items():
            full = spearman(x[pos], y[pos])
            lo, hi, _ = interval(spearman_boot(x, y, idx))
            pairs.append({"pair": f"{other}-POSTN", "set": set_name, "n": len(pos),
                          "rho": round(full, 4), "ci_low": round(lo, 4), "ci_high": round(hi, 4)})
        for batch, part in meta.loc[fshd_ids].groupby("batch"):
            sel = fshd_ids.isin(part.index)
            ax.scatter(x[sel], y[sel], s=18, label=f"batch {batch}")
        ax.set_xlabel(f"{other} (VST)")
        ax.set_ylabel("POSTN (VST)")
        ax.set_title(f"{other} vs POSTN, 34 FSHD biopsies", fontsize=9)
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "postn_pairs.png", dpi=150)
    plt.close(fig)
    pd.DataFrame(pairs).to_csv(OUT / "postn_correlations.csv", index=False)

    fig, ax = plt.subplots(figsize=(6, 4))
    for i, (batch, part) in enumerate(table.groupby("batch")):
        for j, cond in enumerate(["Control", "FSHD"]):
            v = part.loc[part["condition"] == cond, "S1_lfo"]
            ax.scatter(np.full(len(v), i + 0.2 * j - 0.1), v, s=14,
                       marker="o" if cond == "FSHD" else "s",
                       color="tab:red" if cond == "FSHD" else "tab:blue",
                       label=cond if i == 0 else None)
    ax.set_xticks(range(table["batch"].nunique()), [f"batch {b}" for b in
                                                   sorted(table["batch"].unique())])
    ax.set_ylabel("S1 score without the four markers")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "scores_by_batch.png", dpi=150)
    plt.close(fig)

    head = delta[(delta["subset"] == "all") & (delta["score"].isin(["S1_lfo", "S3_lfo"]))]
    lines = [f"bootstrap: {N_BOOT} resamples, seed {SEED}; descriptive, not adjusted for "
             "multiple testing",
             "scores: " + "; ".join(f"{r['score']} {r['features']} features "
                                    f"(markers removed {r['markers_removed']}), "
                                    f"Spearman with GEO {r['spearman_with_geo_fshd']}, FSHD above "
                                    f"the highest control {r['fshd_above_highest_control']}/34"
                                    for _, r in panels.iterrows()),
             "headline, change in rho against the GEO score (all endpoint-complete FSHD):"]
    lines += [f"  {r['endpoint']:<12} {r['score']:<7} n {r['n']:2d}  delta {r['delta_rho']:+.3f} "
              f"[{r['ci_low']:+.3f}, {r['ci_high']:+.3f}]  {r['label']}"
              for _, r in head.iterrows()]
    geo_rho = corr[(corr["subset"] == "all") & (corr["score"] == "geo")]
    lines += ["GEO score rho: " + ", ".join(f"{r['endpoint']} {r['rho']:.3f}"
                                            for _, r in geo_rho.iterrows())]
    key = tissue[(tissue["set"] == "34 FSHD") & tissue["gene"].isin(["RUNX1", "SOX11", "POSTN"])]
    lines += ["tissue state, 34 FSHD (rho; partial rho given the GEO score):"]
    for gene, part in key.groupby("gene"):
        lines.append(f"  {gene}: " + ", ".join(
            f"{r['variable']} {r['rho']:+.2f}"
            + ("" if pd.isna(r.get("partial_rho_given_geo")) else
               f" ({r['partial_rho_given_geo']:+.2f})") for _, r in part.iterrows()))
    lines += ["POSTN pairs: " + "; ".join(f"{r['pair']} {r['set']} rho {r['rho']:+.3f} "
                                          f"[{r['ci_low']:+.3f}, {r['ci_high']:+.3f}]"
                                          for r in pairs)]
    (OUT / "summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
