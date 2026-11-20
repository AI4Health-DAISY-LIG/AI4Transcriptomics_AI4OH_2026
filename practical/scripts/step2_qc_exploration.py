"""Step 2: quality control and exploration of the raw counts.

Outputs (results/step2_qc/):
  library_sizes.png                 total counts per sample, by batch, colour and shape = condition
  library_sizes.csv                 total counts per sample with condition and batch
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

import config  # noqa: E402

OUT = config.RESULTS / "step2_qc"
AUDIT = config.RESULTS / "step1_audit"


def load_counts() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Raw counts (genes x samples) and metadata aligned on the sample columns."""
    counts = pd.read_csv(config.RAW / config.COUNTS_FILE, index_col=0)
    samples = pd.read_csv(AUDIT / "sample_table.csv")
        # Join on GEO `title`: it matches all 43 count columns, including 01-0022-1, whose
    # char:sample_id is 01-0022 (results/step1_audit/counts_header_check.txt,
    # worksheets/step1_audit.md). char:sample_id stays as a column, as the key to the
    # paper's Table S1.
    meta = samples.set_index("title").reindex(counts.columns)
    missing = meta["condition"].isna()
    if missing.any():
        raise ValueError(f"count columns without a GEO record: {list(counts.columns[missing])}")
    return counts, meta


def scatter_by(pcs: pd.DataFrame, values: pd.Series, var: str, ratio, title: str, path) -> None:
    """PC1 against PC2, coloured by one variable (colour scale if numeric with many values)."""
    fig, ax = plt.subplots(figsize=(6, 5))
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.notna().any() and numeric.nunique() > 8:
        na = numeric.isna()
        ax.scatter(pcs.loc[na, "PC1"], pcs.loc[na, "PC2"], c="lightgrey", s=35, label="no value")
        sc = ax.scatter(pcs.loc[~na, "PC1"], pcs.loc[~na, "PC2"], c=numeric[~na], cmap="viridis", s=35)
        fig.colorbar(sc, ax=ax, label=var.replace("char:", ""))
        if na.any():
            ax.legend(fontsize=7)
    else:
        levels = values.astype(object).where(values.notna(), "no value").astype(str)
        for level in sorted(levels.unique()):
            idx = levels == level
            ax.scatter(pcs.loc[idx, "PC1"], pcs.loc[idx, "PC2"], label=level, s=35, alpha=0.8)
        ax.legend(fontsize=7, title=var.replace("char:", ""))
    ax.set_xlabel(f"PC1 ({ratio[0]:.0%})")
    ax.set_ylabel(f"PC2 ({ratio[1]:.0%})")
    ax.set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    config.PROCESSED.mkdir(parents=True, exist_ok=True)
    counts, meta = load_counts()
    counts = counts.loc[:, ~meta["condition"].isna()]
    meta = meta.loc[counts.columns]
    print(f"counts: {counts.shape[0]} genes x {counts.shape[1]} samples")

    # Library sizes: total deposited counts per sample, i.e. reads assigned to genes after
    # the authors' gene filter, not sequencing depth. Shown per batch, because depth differs
    # by sequencing run; colour and shape show condition (worksheets/step2_qc.md).
    # Samples labelled: low mapping rate (32-0008, 01-0037, 32-0009) and the paper's
    # two PCA outliers (01-0037, 32-0016).
    watch = ["32-0008", "01-0037", "32-0009", "32-0016"]
    lib = counts.sum(axis=0)
    lib_table = meta[["condition", "char:batch"]].assign(lib_size=lib)
    lib_table.sort_values(["char:batch", "lib_size"]).to_csv(OUT / "library_sizes.csv")
    batches = sorted(lib_table["char:batch"].unique())
    x_pos = lib_table["char:batch"].map(batches.index)
    x_pos = x_pos + np.random.default_rng(0).uniform(-0.15, 0.15, len(x_pos))  # fixed jitter
    y_val = lib_table["lib_size"] / 1e6
    fig, ax = plt.subplots(figsize=(6, 4))
    styles = {"FSHD": ("o", "tab:red"), "Control": ("^", "tab:blue")}
    for condition, (marker, color) in styles.items():
        idx = lib_table["condition"] == condition
        ax.scatter(x_pos[idx], y_val[idx], marker=marker, s=45, color=color,
                   edgecolors="black", linewidths=0.5, alpha=0.85, label=condition)
    for sample in watch:
        ax.annotate(sample, (x_pos[sample], y_val[sample]), fontsize=6,
                    xytext=(4, 0), textcoords="offset points")
    ax.set_xticks(range(len(batches)))
    ax.set_xticklabels([f"batch {b}" for b in batches])
    ax.set_ylabel("counts (millions)")
    ax.set_title("Library size per sample, by batch")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "library_sizes.png", dpi=150)
    plt.close(fig)

    # Gene filtering.
    # TODO(student): choose and justify the rule. The default keeps a gene when at least
    # `min_samples` samples have `min_count` counts or more, with min_samples = smallest group.
    # Check what your rule does to the four biomarker genes of the paper (config.BIOMARKERS).
    min_count = 10
    min_samples = int(meta["condition"].value_counts().min())
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

    # log2 CPM and PCA on the most variable genes.
    cpm = filtered / filtered.sum(axis=0) * 1e6
    logcpm = np.log2(cpm + 1)
    n_top = 2000  # TODO(student): does the choice of n_top change the picture?
    top = logcpm.var(axis=1).sort_values(ascending=False).index[:n_top]
    x = logcpm.loc[top].T
    x = x - x.mean(axis=0)
    pca = PCA(n_components=10)
    names = [f"PC{i + 1}" for i in range(10)]
    pcs = pd.DataFrame(pca.fit_transform(x), index=x.index, columns=names)
    ratio = pca.explained_variance_ratio_
    pd.DataFrame({"component": names, "variance_explained": ratio}).to_csv(
        OUT / "pca_variance.csv", index=False)
    pcs.to_csv(OUT / "pca_scores.csv")

    variables = ["condition"] + config.TECHNICAL_VARS + config.BIOLOGICAL_VARS + config.PATHOLOGY_VARS
    for var in variables:
        if var not in meta.columns:
            continue
        safe = var.replace("char:", "").replace(" ", "_").replace(".", "_")
        scatter_by(pcs, meta[var], var, ratio,
                   f"PCA of log2 CPM, top {n_top} variable genes, coloured by {var}",
                   OUT / f"pca_by_{safe}.png")
    print("figures written to", OUT)
    print("\nNext: fill worksheets/step2_qc.md")


if __name__ == "__main__":
    main()
