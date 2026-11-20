# Step 2: quality control and exploration

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

## Filtering

| Filter | Threshold I chose | Why | Samples and genes before | after |
|--------|-------------------|-----|--------------------------|-------|
| Gene total count | total counts > 20 summed over the 43 samples | Same rule as the authors ([makeSE.R L164, commit 90cd343](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/makeSE.R#L164)), so the analysis starts from the same genes as theirs. The GEO count file was already filtered with it: the smallest gene total is 21 (`results/step2_qc/gene_filter.txt`) | 43 samples, 27,506 genes | 43 samples, 27,506 genes |

Does the paper state its thresholds? Quote them if it does, say so if it does not.

No. The Methods give no gene filter. The threshold appears only in the authors' code (total counts > 20, [makeSE.R L164](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/makeSE.R#L164)), and it was applied before deposit: every gene in the GEO count file has a total of at least 21 (`gene_filter.txt`). I do not add the starter's default filter (≥ 10 counts in ≥ 9 samples), which is not part of the paper's analysis. The starter's rule would remove 8,248 genes. It also cannot retain a gene whose counts are confined to a group with fewer than nine biopsies.

Do the four biomarker genes of the paper pass your filter?

Yes. Totals over the 43 samples: LEUTX 902, KHDC1L 2,212, PRAMEF2 1,723, TRIM43 2,550 (`gene_filter.txt`).

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

- None in the main analysis: all 43 deposited biopsies are kept. Three biopsies, all FSHD from batch II, have low mapping rates: 32-0008 (67.73%), 01-0037 (76.52%) and 32-0009 (85.60%), against 94.85% or more for every other biopsy, including the paper's other PCA outlier, 32-0016 (95.26%). Rates are calculated from the integer input and mapped read counts in GEO's R object and rounded to two decimals (`results/step2_qc/mapping_rates.csv`). A sensitivity run excluding these three biopsies is planned; it will keep the same 27,506 genes so that only the samples change.

## Decisions I made

- Apply the authors' gene filter (total counts > 20) instead of the starter's default, to reproduce their analysis; on the deposited file it removes no gene.
- Library sizes (total counts assigned to genes, not sequencing depth; `library_sizes.png`, `library_sizes.csv`) follow sequencing batch: 59.8–80.2 million per biopsy in batch III, at most 47.7 million in the other batches. Within batches III and IV, controls and FSHD biopsies have similar sizes. Depth differences are handled by normalisation (CPM in the PCA, size factors in DESeq2), so no biopsy is excluded for depth; the smallest is 32-0008 (17.6 million).