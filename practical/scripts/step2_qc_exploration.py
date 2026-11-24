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

"""Step 2: quality control and PCA of the raw counts."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

import config

OUT = config.RESULTS / "step2_qc"
OUT.mkdir(parents=True, exist_ok=True)
config.PROCESSED.mkdir(parents=True, exist_ok=True)

# 1. counts and metadata
counts = pd.read_csv(config.RAW / config.COUNTS_FILE, index_col=0)
samples = pd.read_csv(config.RESULTS / "step1_audit" / "sample_table.csv")
# TODO(student): read results/step1_audit/counts_header_check.txt first.
# Only the count columns whose name matches a GEO sample_id are kept.
# If a column has no match, find out which biopsy it is and decide:
# rename it, or leave it out and say why.
key = "char:sample_id" if "char:sample_id" in samples.columns else "title"
meta = samples.set_index(key).reindex(counts.columns)
found = meta["condition"].notna()
print((~found).sum(), "count column(s) dropped")
counts = counts.loc[:, found]
meta = meta.loc[found]
print(counts.shape[0], "genes x", counts.shape[1], "samples")

# 2. library sizes
lib = counts.sum(axis=0)
colors = meta["condition"].map({"Control": "tab:blue", "FSHD": "tab:red"})
plt.figure(figsize=(12, 4))
plt.bar(range(len(lib)), lib.values / 1e6, color=colors.values)
plt.xticks(range(len(lib)), lib.index, rotation=90, fontsize=6)
plt.ylabel("counts (millions)")
plt.title("Library size per sample (blue: control, red: FSHD)")
plt.tight_layout()
plt.savefig(OUT / "library_sizes.png", dpi=150)
plt.close()

# 3. gene filter
# TODO(student): choose and justify the rule. A gene is kept if at least
# min_samples samples have min_count counts or more.
# Check what the rule does to the biomarker genes (config.BIOMARKERS).
min_count = 10
min_samples = int(meta["condition"].value_counts().min())
keep = (counts >= min_count).sum(axis=1) >= min_samples
filtered = counts.loc[keep]
with open(OUT / "gene_filter.txt", "w") as f:
    f.write(f"rule: >= {min_count} counts in >= {min_samples} samples\n")
    f.write(f"genes before: {len(counts)}\n")
    f.write(f"genes after: {len(filtered)}\n")
print("genes kept:", len(filtered), "/", len(counts))
filtered.to_csv(config.PROCESSED / "counts_filtered.csv")
meta.to_csv(config.PROCESSED / "metadata.csv")

# 4. log2 CPM and PCA on the most variable genes
cpm = filtered / filtered.sum(axis=0) * 1e6
logcpm = np.log2(cpm + 1)
n_top = 2000  # TODO(student): does n_top change the picture?
top = logcpm.var(axis=1).sort_values(ascending=False).index[:n_top]
x = logcpm.loc[top].T
x = x - x.mean(axis=0)
pca = PCA(n_components=10)
names = ["PC" + str(i + 1) for i in range(10)]
pcs = pd.DataFrame(pca.fit_transform(x), index=x.index, columns=names)
ratio = pca.explained_variance_ratio_
pd.DataFrame({"component": names, "variance_explained": ratio}).to_csv(OUT / "pca_variance.csv", index=False)
pcs.to_csv(OUT / "pca_scores.csv")

# 5. one PCA plot per variable
variables = ["condition"] + config.TECHNICAL_VARS + config.BIOLOGICAL_VARS + config.PATHOLOGY_VARS
for var in variables:
    if var not in meta.columns:
        continue
    values = meta[var]
    numbers = pd.to_numeric(values, errors="coerce")
    plt.figure(figsize=(6, 5))
    if numbers.notna().any() and numbers.nunique() > 8:
        plt.scatter(pcs["PC1"], pcs["PC2"], c=numbers, cmap="viridis")
        plt.colorbar(label=var)
    else:
        labels = values.fillna("no value").astype(str)
        for level in labels.unique():
            rows = labels == level
            plt.scatter(pcs.loc[rows, "PC1"], pcs.loc[rows, "PC2"], label=level)
        plt.legend(title=var)
    plt.xlabel(f"PC1 ({ratio[0]:.0%})")
    plt.ylabel(f"PC2 ({ratio[1]:.0%})")
    plt.title(f"PCA coloured by {var}", fontsize=9)
    plt.tight_layout()
    name = var.replace("char:", "").replace(" ", "_").replace(".", "_")
    plt.savefig(OUT / f"pca_by_{name}.png", dpi=150)
    plt.close()
print("figures written to", OUT)

