"""Step 3: differential expression with PyDESeq2 on raw counts.

Inputs:  data/processed/counts_filtered.csv and metadata.csv (from step 2)
Outputs (results/step3_de/):
  de_FSHD_vs_Control.csv        all FSHD biopsies against the controls
  de_G<k>_vs_Control.csv        each group of the paper (1 to 4) against the controls
  summary.txt                   design, samples per level, genes tested and significant
  ma_<tag>.png, volcano_<tag>.png, pvalue_histogram_<tag>.png

Decisions marked TODO(student) are yours. Record them in worksheets/step3_de.md.

Run:  uv run python scripts/step3_differential_expression.py
"""

from __future__ import annotations

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from pydeseq2.dds import DeseqDataSet  # noqa: E402
from pydeseq2.ds import DeseqStats  # noqa: E402

import config  # noqa: E402

OUT = config.RESULTS / "step3_de"
ALPHA = 0.05

#define the outliers samples
outliers = ["320016", "01-0037"]
#remove outliers from metadata.csv and from counts_filtered.csv
metadata = metadata.drop(index=outliers)
counts_filtered = counts_filtered.drop(columns=outliers)



# TODO(student): the design. False fits "~ <factor>" and ignores the batch; True fits
# "~ batch + <factor>". Look at results/step1_audit/crosstab_batch.csv and at your PCA
# before you decide, run both, and report what changes.
ADJUST_FOR_BATCH = False


def plots(res: pd.DataFrame, sig: pd.DataFrame, tag: str) -> None:
    label = res["gene_name"] if "gene_name" in res.columns else res.index.to_series()
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(np.log10(res["baseMean"] + 1), res["log2FoldChange"], s=4, c="lightgrey")
    ax.scatter(np.log10(sig["baseMean"] + 1), sig["log2FoldChange"], s=4, c="tab:red")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("log10(base mean + 1)")
    ax.set_ylabel(f"log2 fold change {tag}")
    ax.set_title(f"MA plot ({len(sig)} genes with padj < {ALPHA})")
    fig.tight_layout()
    fig.savefig(OUT / f"ma_{tag}.png", dpi=150)
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
    fig.savefig(OUT / f"volcano_{tag}.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(res["pvalue"].dropna(), bins=50)
    ax.set_xlabel("p-value")
    ax.set_ylabel("genes")
    ax.set_title(tag)
    fig.tight_layout()
    fig.savefig(OUT / f"pvalue_histogram_{tag}.png", dpi=150)
    plt.close(fig)


def run(x: pd.DataFrame, meta: pd.DataFrame, factor: str, names: pd.Series | None) -> list[str]:
    """Fit one model for `factor` and test each of its levels against the controls."""
    meta = meta[meta[factor].notna()]
    x = x.loc[meta.index]
    design = f"~ batch + {factor}" if ADJUST_FOR_BATCH else f"~ {factor}"
    dds = DeseqDataSet(counts=x, metadata=meta[[factor, "batch"]], design=design, quiet=True)
    dds.deseq2()
    lines = [f"design: {design}",
             "samples per level: " + ", ".join(f"{k}={v}" for k, v in meta[factor].value_counts().sort_index().items())]
    for level in sorted(set(meta[factor]) - {config.REFERENCE_LEVEL}):
        tag = f"{level}_vs_{config.REFERENCE_LEVEL}"
        stats = DeseqStats(dds, contrast=[factor, level, config.REFERENCE_LEVEL], alpha=ALPHA, quiet=True)
        stats.summary()
        res = stats.results_df
        if names is not None:
            res.insert(0, "gene_name", names.reindex(res.index))
        res = res.sort_values("padj")
        res.to_csv(OUT / f"de_{tag}.csv", float_format="%.5g")
        sig = res[res["padj"] < ALPHA]
        n_up = int((sig["log2FoldChange"] > 0).sum())
        lines.append(f"[{tag}] genes tested: {len(res)}  padj < {ALPHA}: {len(sig)}  "
                     f"up: {n_up}  down: {len(sig) - n_up}")
        plots(res, sig, tag)
    return lines


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    counts = pd.read_csv(config.PROCESSED / "counts_filtered.csv", index_col=0)
    meta = pd.read_csv(config.PROCESSED / "metadata.csv", index_col=0).loc[counts.columns]
    names = None
    if config.GENE_NAMES.exists():
        names = pd.read_csv(config.GENE_NAMES, index_col="gene_id")["gene_name"]

    # PyDESeq2 wants samples in rows, genes in columns, integer counts.
    x = counts.T.round().astype(int)
    meta = meta.copy()
    meta["condition"] = meta["condition"].astype(str)
    meta["batch"] = meta["char:batch"].astype(str)

    # The paper does not test all FSHD biopsies against the controls. It splits them into
    # four groups by biomarker expression and tests each group against the controls.
    # TODO(student): the groups below are the authors' own (rnaseq.rank on GEO). They were
    # built from DUX4 target genes. Say in your report what that does to a test of DUX4
    # target genes between a group and the controls.
    rank = pd.to_numeric(meta[config.GROUP_COL], errors="coerce")
    meta["group"] = ("G" + rank.astype("Int64").astype(str)).where(rank.notna())
    meta.loc[meta["condition"] == config.REFERENCE_LEVEL, "group"] = config.REFERENCE_LEVEL

    summary = [f"samples: {x.shape[0]}  genes: {x.shape[1]}", ""]
    summary += run(x, meta, "condition", names) + [""]
    summary += run(x, meta, "group", names)
    (OUT / "summary.txt").write_text("\n".join(summary), encoding="utf-8")
    print("\n".join(summary))
    print("\nNext: fill worksheets/step3_de.md. Run a second design before you conclude")


if __name__ == "__main__":
    main()
