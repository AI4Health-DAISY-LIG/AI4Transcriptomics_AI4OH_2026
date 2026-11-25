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

- PCA coloured by condition: 
- PCA coloured by each technical variable: 
- PCA coloured by the authors' group (`rnaseq.rank`): 

Name the variable that separates the samples on PC1, then on PC2: the condition,
a technical variable, the authors' group, or none you can name. Give the variance
PC1 and PC2 explain.

Compare your PCA with Figure 3 of the paper. The same or not, and why:

## Samples I excluded, and why

-

## Decisions I made

-
