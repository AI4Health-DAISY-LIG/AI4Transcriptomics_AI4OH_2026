"""Step 3: differential expression with PyDESeq2 on raw counts, with sensitivity runs.

Inputs:  data/processed/counts_filtered.csv and metadata.csv (from step 2)
Outputs (results/step3_de/):
  de_<level>_vs_Control.csv, plots, summary.txt   primary model
  sensitivity/<run>/                              the same files for each other run
  comparison.csv                                  every run compared with the primary one
  blood_by_condition.csv                          blood variables per site and condition

Record your decisions in worksheets/step3_de.md.

Run:  uv run python scripts/step3_differential_expression.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from formulaic import model_matrix  # noqa: E402
from pydeseq2.dds import DeseqDataSet  # noqa: E402
from pydeseq2.ds import DeseqStats  # noqa: E402

import config  # noqa: E402

OUT = config.RESULTS / "step3_de"
ALPHA = 0.05
BLOOD = ["hemoglobin_pct", "neutrophil_pct"]
FULL = ["batch", "sex", "site", *BLOOD]
REFERENCE_RUN = {"site01_no_blood": "site01"}


def without(*terms: str) -> list[str]:
    return [t for t in FULL if t not in terms]


# name: (covariates, sample selection). The first run is the primary model.
RUNS = {
    "primary": (FULL, "inliers"),
    "no_batch": (without("batch"), "inliers"),
    "no_site": (without("site"), "inliers"),
    "no_blood": (without(*BLOOD), "inliers"),
    "site01": (without("site"), "site01"),
    "site01_no_blood": (without("site", *BLOOD), "site01"),
    "batch_III_IV": (FULL, "batch_III_IV"),
    "without_32-0014": (FULL, "without_32-0014"),
    "with_outliers": (FULL, "all"),
}


def select(meta: pd.DataFrame, subset: str) -> pd.DataFrame:
    inlier = ~meta.index.isin(config.OUTLIERS)
    if subset == "all":
        return meta
    if subset == "inliers":
        return meta[inlier]
    if subset == "site01":
        return meta[inlier & (pd.to_numeric(meta["site"]) == 1)]
    if subset == "batch_III_IV":
        return meta[inlier & meta["batch"].isin(["III", "IV"])]
    if subset == "without_32-0014":
        if "32-0014" not in meta.index:
            raise ValueError("sample 32-0014 not found: check the sample names")
        return meta[inlier & (meta.index != "32-0014")]
    raise ValueError(f"unknown subset: {subset}")


def plots(res: pd.DataFrame, sig: pd.DataFrame, tag: str, out: Path) -> None:
    label = res["gene_name"] if "gene_name" in res.columns else res.index.to_series()
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


def run(
    x: pd.DataFrame,
    meta: pd.DataFrame,
    factor: str,
    covariates: list[str],
    names: pd.Series | None,
    out: Path,
) -> tuple[list[str], dict[str, tuple[pd.DataFrame, int, int]]]:
    """Fit one model for `factor` and test each of its levels against the controls."""
    meta = meta.dropna(subset=[factor, *covariates])
    x = x.loc[meta.index]
    design = "~ " + " + ".join([*covariates, factor])
    data = meta[[factor, *covariates]]

    # A model whose terms are fully redundant cannot be estimated.
    matrix = np.asarray(model_matrix(design, data), dtype=float)
    if np.linalg.matrix_rank(matrix) < matrix.shape[1]:
        return [f"design: {design}", "SKIPPED: terms are fully confounded"], {}

    dds = DeseqDataSet(counts=x, metadata=data, design=design, quiet=True)
    dds.deseq2()
    counts_per_level = meta[factor].value_counts().sort_index()
    lines = [
        f"design: {design}",
        "samples per level: "
        + ", ".join(f"{k}={v}" for k, v in counts_per_level.items()),
    ]
    results = {}
    for level in sorted(set(meta[factor]) - {config.REFERENCE_LEVEL}):
        tag = f"{level}_vs_{config.REFERENCE_LEVEL}"
        stats = DeseqStats(
            dds,
            contrast=[factor, level, config.REFERENCE_LEVEL],
            alpha=ALPHA,
            quiet=True,
        )
        stats.summary()
        res = stats.results_df.copy()
        if names is not None:
            res.insert(0, "gene_name", names.reindex(res.index))
        res = res.sort_values("padj")
        res.to_csv(out / f"de_{tag}.csv", float_format="%.5g")
        sig = res[res["padj"] < ALPHA]
        n_up = int((sig["log2FoldChange"] > 0).sum())
        lines.append(
            f"[{tag}] genes tested: {len(res)}  padj < {ALPHA}: {len(sig)}  "
            f"up: {n_up}  down: {len(sig) - n_up}"
        )
        plots(res, sig, tag, out)
        results[tag] = (
            res,
            int(counts_per_level[config.REFERENCE_LEVEL]),
            int(counts_per_level[level]),
        )
    return lines, results


def compare(res: pd.DataFrame, ref: pd.DataFrame) -> dict:
    """Compare one contrast of one run with the same contrast of a reference run."""
    sig = set(res.index[res["padj"] < ALPHA])
    ref_sig = set(ref.index[ref["padj"] < ALPHA])
    lfc = res["log2FoldChange"].reindex(ref.index)
    ref_lfc = ref["log2FoldChange"]
    same_sign = np.sign(lfc) == np.sign(ref_lfc)  # a missing value counts as False

    shared = list(sig & ref_sig)
    ref_genes = list(ref_sig)
    top = ref.sort_values("padj").head(100).index
    up = int((res.loc[list(sig), "log2FoldChange"] > 0).sum())

    return {
        "significant": len(sig),
        "up": up,
        "down": len(sig) - up,
        "ref_significant": len(ref_sig),
        "shared_same_direction": int(same_sign[shared].sum()),
        "shared_opposite_direction": int((~same_sign[shared]).sum()),
        "only_this_run": len(sig - ref_sig),
        "only_reference": len(ref_sig - sig),
        # rank agreement over all genes (says nothing about effect size)
        "lfc_spearman": round(lfc.rank().corr(ref_lfc.rank()), 3),
        # over the reference's significant genes: 1 = same size, 0.5 = halved
        "median_effect_ratio": round((lfc[ref_genes] / ref_lfc[ref_genes]).median(), 3),
        "top_genes": len(top),
        "top_same_direction": int(same_sign[top].sum()),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    counts = pd.read_csv(config.PROCESSED / "counts_filtered.csv", index_col=0)
    meta = pd.read_csv(config.PROCESSED / "metadata.csv", index_col=0).loc[
        counts.columns
    ]

    names = None
    if config.GENE_NAMES.exists():
        names = pd.read_csv(config.GENE_NAMES, index_col="gene_id")["gene_name"]

    # PyDESeq2 wants samples in rows, genes in columns, integer counts.
    x = counts.T.round().astype(int)
    meta = meta.copy()
    meta["condition"] = meta["condition"].astype(str)
    for col, source in [
        ("batch", "char:batch"),
        ("sex", "char:infer_gender"),
        ("site", "char:site"),
    ]:
        meta[col] = meta[source].astype(str).where(meta[source].notna())

    # The paper does not test all FSHD biopsies against the controls. It splits them into
    # four groups by biomarker expression and tests each group against the controls.
    # TODO(student): the groups below are the authors' own (rnaseq.rank on GEO). They were
    # built from DUX4 target genes. Say in your report what that does to a test of DUX4
    # target genes between a group and the controls.
    rank = pd.to_numeric(meta[config.GROUP_COL], errors="coerce")
    meta["group"] = ("G" + rank.astype("Int64").astype(str)).where(rank.notna())
    meta.loc[meta["condition"] == config.REFERENCE_LEVEL, "group"] = (
        config.REFERENCE_LEVEL
    )

    blood = select(meta, "inliers").groupby(["site", "condition"])[BLOOD]
    blood = blood.agg(["count", "median", "min", "max"])
    blood.to_csv(OUT / "blood_by_condition.csv")
    print(blood.to_string())

    done: dict[str, dict] = {}  # run name -> {contrast: (res, n_control, n_case)}
    rows = []
    for name, (covariates, subset) in RUNS.items():
        out = OUT if name == "primary" else OUT / "sensitivity" / name
        out.mkdir(parents=True, exist_ok=True)
        selected = select(meta, subset)
        sub = selected.dropna(subset=FULL)
        dropped = selected.index.difference(sub.index)
        if len(dropped):
            print(
                f"{name}: dropped for missing covariates: {', '.join(map(str, dropped))}"
            )
        summary = [f"run: {name}  samples: {len(sub)}  genes: {x.shape[1]}", ""]
        done[name] = {}
        for factor in ("condition", "group"):
            lines, results = run(x.loc[sub.index], sub, factor, covariates, names, out)
            summary += lines + [""]
            done[name].update(results)
        (out / "summary.txt").write_text("\n".join(summary), encoding="utf-8")
        print("\n" + "\n".join(summary))

        if name == "primary":
            continue
        ref_name = REFERENCE_RUN.get(name, "primary")
        for tag in done["primary"]:
            row = {"run": name, "reference": ref_name, "contrast": tag}
            if tag not in done[name]:
                row["status"] = "not estimated in this run"
            elif tag not in done[ref_name]:
                row["status"] = "not estimated in reference"
            else:
                res, n_control, n_case = done[name][tag]
                ref, ref_n_control, ref_n_case = done[ref_name][tag]
                row.update(
                    status="ok",
                    n_control=n_control,
                    n_case=n_case,
                    ref_n_control=ref_n_control,
                    ref_n_case=ref_n_case,
                )
                row.update(compare(res, ref))
            rows.append(row)
        pd.DataFrame(rows).to_csv(OUT / "comparison.csv", index=False)

    table = pd.DataFrame(rows)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
