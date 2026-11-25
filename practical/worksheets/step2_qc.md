# Step 2: quality control and exploration

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

## Filtering

| Filter | Threshold I chose | Why | Samples and genes before | after |
|--------|-------------------|-----|--------------------------|-------|
| low-count genes | ≥10 counts in ≥9 samples | 9 is the control group size | 43 samples, 27506 genes | 43 samples, 19198 genes |

<!-- numbers from results/step2_qc/gene_filter.txt -->

Paper does not mention any filtering threshold.

All four biomarker genes (LEUTX, KHDC1L, PRAMEF2, TRIM43) pass the filter.

## Exploration

Figures (in `results/step2_qc/`):

- PCA coloured by condition: controls cluster on the left, FSHD spreads to the right along PC1.
- PCA coloured by each technical variable: batches I and II sit on the right with FSHD, while III and IV are mixed, so batch overlaps with condition.
- PCA coloured by the authors' group (`rnaseq.rank`): rank 1 to 4 move from left to right on PC1, with rank 4 at the top right and two rank 1 outliers at the bottom right.

<!-- variance from results/step2_qc/pca_variance.csv -->

PC1 explains 65% variance (condition / `rnaseq.rank`) and PC2 explains 11% (no single variable clearly explains it).

Compare with paper Figure 3: same overall pattern (controls on the left, groups 1 to 4 shifting right, and two group 1 outliers on the far right). Variance differs (our PC1 = 65%, PC2 = 11% vs paper PC1 = 38%, PC2 = 14%) because we used log2 CPM on top 2000 variable genes, and PC2 is flipped vertically.

## Samples I excluded, and why

-

## Decisions I made

-
