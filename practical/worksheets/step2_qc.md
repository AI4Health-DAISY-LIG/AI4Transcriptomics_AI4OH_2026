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

Figures (in `results/step2_qc/`). PCA with the authors' method: rlog of the GEO counts (`r/rlog_from_geo.R`), 500 most variable genes, centred (`scripts/step2_qc_exploration.py`); shape = condition, watched samples labelled.

- PCA coloured by condition: `pca_by_condition.png`
- PCA coloured by each technical variable: `pca_by_batch.png`, `pca_by_instrument_model.png` and `pca_by_tissue.png` (one value for all 43 samples, so they cannot separate anything), plus `pca_by_library_size_millions.png` and `pca_by_mapping_rate_pct.png`
- PCA coloured by the authors' group (`rnaseq.rank`): `pca_by_paper_group.png`, with the five categories of Figure 3 (controls are in group 1 by construction of the authors' clustering)

Name the variable that separates the samples on PC1, then on PC2: the condition,
a technical variable, the authors' group, or none you can name. Give the variance
PC1 and PC2 explain.

PC1 (37.3%, `pca_variance.csv`) is separated by the authors' group: it runs from the controls through groups 1 to 4, and the group explains 61% of PC1's variance, against 30% for condition alone and 36% for batch (`pca_associations.csv`). PC1 also follows the disease measures: inflammation (Spearman ρ = 0.70), DUX4 score (0.67), pathology (0.63), T1 (0.49), fat fraction (0.47) and STIR (0.44). Its link with batch mostly reflects which groups each batch holds, since 4 of the 6 group 4 biopsies are in batch I (`results/step1_audit/crosstabs.md`, `sample_table.csv`). Library size and mapping rate correlate weakly and negatively with PC1 (ρ = −0.39 and −0.41), also in line with batch.

PC2 (13.6%) is driven by four biopsies at the top: 32-0016, 01-0037, 32-0008 and 32-0009 (`pca_scores.csv`). Condition explains nothing on PC2 (R² = 0.004); batch explains 33%, through batch II. None I can name for the whole axis: three of the four are the low-mapping batch II libraries, but 32-0016 (batch IV) maps normally (95.26%, `mapping_rates.csv`) and has the highest fat fraction among them (0.43), which fits the paper's explanation of it as a fat-rich biopsy.

Compare your PCA with Figure 3 of the paper. The same or not, and why:

The same. Replaying the authors' PCA on their saved rlog matrix (`r/pca_replay.R`) gives 37.94% / 13.77%, the published 38% / 14%. Recomputing rlog from the GEO counts with DESeq2 1.50.2 (`r/rlog_from_geo.R`) gives 37.32% / 13.61%, with the same layout (PC1 and PC2 score correlations above 0.9997 with the replay) and the same two outliers named in the legend, 01-0037 and 32-0016 (`rlog_geo_comparison.txt`). The small difference comes from 9 sparse genes whose rlog values differ between DESeq2 versions; without them, both PCAs give 38.38% / 13.91%. 32-0008 and 32-0009 also sit high on PC2 in the authors' own matrix (`pca_replay_scores.csv`), so they are in Figure 3 too, but the paper does not name them.

The course default, log2(CPM + 1) on the 2,000 most variable genes, gives 63.8% / 11.2% and is not comparable, for three reasons (`pca_logcpm_comparison.csv`): the transformation alone, on Figure 3's 500 genes, raises PC1 to 46.1%; letting log-CPM choose its own 500 genes raises it to 58.1%, because only 197 of them are Figure 3's genes; and taking 2,000 genes gives 63.8%.

Compare your PCA with Figure 3 of the paper. The same or not, and why:

## Samples I excluded, and why

- None in the main analysis: all 43 deposited biopsies are kept. Three biopsies, all FSHD from batch II, have low mapping rates: 32-0008 (67.73%), 01-0037 (76.52%) and 32-0009 (85.60%), against 94.85% or more for every other biopsy, including the paper's other PCA outlier, 32-0016 (95.26%). Rates are calculated from the integer input and mapped read counts in GEO's R object and rounded to two decimals (`results/step2_qc/mapping_rates.csv`). A sensitivity run excluding these three biopsies is planned; it will keep the same 27,506 genes so that only the samples change.

## Decisions I made

- Apply the authors' gene filter (total counts > 20) instead of the starter's default, to reproduce their analysis; on the deposited file it removes no gene.
- Library sizes (total counts assigned to genes, not sequencing depth; `library_sizes.png`, `library_sizes.csv`) follow sequencing batch: 59.8–80.2 million per biopsy in batch III, at most 47.7 million in the other batches. Within batches III and IV, the control and FSHD ranges overlap. Normalisation accounts for these differences in overall count scale (DESeq2 size factors in the rlog PCA and the differential-expression analysis; CPM only in the exploratory log-CPM comparison), so no biopsy is excluded for depth; the smallest is 32-0008 (17.6 million).
- PCA: reproduce the authors' method (DESeq2 rlog, 500 most variable genes) rather than the starter's log2(CPM + 1) on 2,000 genes, using a locked R environment for the rlog (`r/pixi.toml`, R 4.5.3, DESeq2 1.50.2); Python recomputes the PCA from the R matrix and stops if it disagrees with R. The log-CPM PCA is kept only as a comparison (`pca_logcpm_comparison.csv`).
