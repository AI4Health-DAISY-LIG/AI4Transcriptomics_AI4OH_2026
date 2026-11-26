# Step 2: quality control and exploration

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

## Filtering

| Filter | Threshold I chose | Why | Samples and genes before | after |
|--------|-------------------|-----|--------------------------|-------|
| low-count genes | >= 10 counts in >= 9 samples | 9 is the smallest group (controls) | 42 samples, 27506 genes | 42 samples, 19198 genes |

The paper does not state any filtering threshold.

The four biomarker genes (LEUTX, KHDC1L, PRAMEF2, TRIM43) pass the filter.

## Exploration

Figures (in `results/step2_qc/`):

- PCA coloured by condition: controls are on the left, FSHD samples spread to the right on PC1.
- PCA coloured by each technical variable: batches I and II are on the right with FSHD, batches III and IV are mixed. Batch overlaps with condition.
- PCA coloured by the authors' group (`rnaseq.rank`): rank 1 to 4 move from left to right on PC1. Rank 4 is at the top right and two rank 1 outliers are at the bottom right.

PC1 explains 65% of the variance (condition and `rnaseq.rank`). PC2 explains 11% (no clear variable).

Compare your PCA with Figure 3 of the paper. The same or not, and why:

- Similar pattern (controls left, groups 1 to 4 moving right, 2 group 1 outliers). Variance differs (65%/11% vs 38%/14%) because we used log2 CPM on top 2000 genes, and PC2 is flipped.

## Samples I excluded, and why

- `01-0022-1` was dropped by the script because its name in counts did not match `01-0022` in GEO metadata (`data/processed/metadata.csv` has 42 samples).
- None of the 42 remaining samples were excluded. All have normal library sizes (17M to 80M counts in `library_sizes.png`).

## Decisions I made

- Kept `min_count = 10` and `min_samples = 9` (controls count). Genes went from 27506 to 19198.
- Kept `n_top = 2000` for PCA.
- Left the default sample join in `scripts/step2_qc_exploration.py`, which kept 42 samples and dropped `01-0022-1`.
