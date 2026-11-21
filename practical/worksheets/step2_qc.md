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

PC1 (37.3%, `pca_variance.csv`) is separated by the authors' group: it runs from the controls through groups 1 to 4, and the group explains 61% of PC1's variance, against 30% for condition alone and 36% for batch (`pca_associations.csv`). PC1 also follows the disease measures: inflammation (Spearman ρ = 0.70), DUX4 score (0.67), pathology (0.63), T1 (0.49), fat fraction (0.47) and STIR (0.44). Its link with batch mostly reflects which groups each batch holds, since 4 of the 6 group 4 biopsies are in batch I (`paper_group_by_batch.csv`). Library size and mapping rate correlate weakly and negatively with PC1 (ρ = −0.39 and −0.41), also in line with batch.

PC2 (13.6%) is driven by four biopsies at the top: 32-0016, 01-0037, 32-0008 and 32-0009 (`pca_scores.csv`). Condition explains nothing on PC2 (R² = 0.004); batch explains 33%, through batch II. None I can name for the whole axis: three of the four are the low-mapping batch II libraries, but 32-0016 (batch IV) maps normally (95.26%, `mapping_rates.csv`) and has the highest fat fraction among them (0.43), which fits the paper's explanation of it as a fat-rich biopsy. Without the three low-mapping libraries, PC2's link with batch drops from R² 0.33 to 0.08, and 32-0016 stays at the top (see "Samples I excluded"). PC2 is not a single factor (`results/step2_qc/pca_top_genes.csv`). In the main run, 4 of the 8 genes with the largest PC2 weight are sex-chromosome genes: XIST, which is expressed in women, has the opposite sign to the Y-chromosome genes EIF1AY, RPS4Y1 and USP9Y. The other four are MYH8 (a developmental myosin that regenerating muscle fibres make), EGFL6, DMBT1 and ADIPOQ (a fat-cell gene). Without the three libraries, 7 of the top 8 PC2 genes are sex-chromosome genes, and the eighth is HBG2, a haemoglobin gene (`results/step2_qc/sensitivity_40/pca_top_genes.csv`). So PC2 mixes sex with tissue composition, and which tissue signal shows up depends on which samples are included.

Compare your PCA with Figure 3 of the paper. The same or not, and why:

The same. Replaying the authors' PCA on their saved rlog matrix (`r/pca_replay.R`) gives 37.94% / 13.77%, the published 38% / 14%. Recomputing rlog from the GEO counts with DESeq2 1.50.2 (`r/rlog_from_geo.R`) gives 37.32% / 13.61%, with the same layout (PC1 and PC2 score correlations above 0.9997 with the replay) and the same two outliers named in the legend, 01-0037 and 32-0016 (`rlog_geo_comparison.txt`). The small difference comes from 9 sparse genes whose rlog values differ between DESeq2 versions; without them, both PCAs give 38.38% / 13.91%. 32-0008 and 32-0009 also sit high on PC2 in the authors' own matrix (`pca_replay_scores.csv`), so they are in Figure 3 too, but the paper does not name them.

The course default, log2(CPM + 1) on the 2,000 most variable genes, gives 63.8% / 11.2% and is not comparable, for three reasons (`pca_logcpm_comparison.csv`): the transformation alone, on Figure 3's 500 genes, raises PC1 to 46.1%; letting log-CPM choose its own 500 genes raises it to 58.1%, because only 197 of them are Figure 3's genes; and taking 2,000 genes gives 63.8%.

Compare your PCA with Figure 3 of the paper. The same or not, and why:

## Samples I excluded, and why

- **Samples I excluded:** none in the main analysis; all 43 biopsies are kept. As a sensitivity check, I removed the three libraries with a mapping rate below 90% (32-0008, 01-0037 and 32-0009, all FSHD biopsies from batch II; `results/step2_qc/mapping_rates.csv`). I then refitted the PCA on the same rlog values and the same 500 genes (`results/step2_qc/sensitivity_40/pca_comparison.txt`). PC1 holds. It explains 38.65% of the variance (37.32% with all 43), its sample scores correlate with the main run at r = 0.996, and it still tracks condition (R² 0.30) and the paper's groups (R² 0.75, up from 0.61). PC2 changes more. It falls from 13.61% to 10.76%, its scores correlate with the main run at r = 0.90, and its link with batch drops from R² 0.33 to 0.08. 32-0016 is still the sample with the highest PC2 score. Re-selecting the 500 most variable genes on the 40 samples gives nearly the same result (37.23% and 10.65%; 432 of the 500 genes are shared). So the low-mapping libraries do not drive the main FSHD axis, but they explain most of PC2's apparent link with batch.
- No biopsy is removed because of sex, but one sex label contradicts the expression data. I inferred each biopsy's sex from its expression (`results/step2_qc/sex_check.csv`, `sex_check.png`). The method compares XIST, which only women express, with the mean of six Y-chromosome genes, both in log2(CPM + 1). A biopsy is called male when the Y mean is above XIST, a rule I fixed before looking at the data. The biopsies split into two clear groups: the closest one to the boundary is still 2.65 log2 units from it. The nine controls come out as 6 women and 3 men (01-0044, 01-0046 and 01-0049 are the men), as the paper reports. GEO's label agrees with the expression for 33 of the 34 FSHD biopsies. The exception is 32-0004, which GEO and the paper's Table S1 both list as female (`results/step1_audit/paper_table_s1_id_check.txt`). Its library has almost no XIST (0.03, against 3.01 to 6.73 for the women) and Y-gene expression at the male level (5.64, against 4.97 to 5.88 for the other men). Either its sex is wrong in both sources, or the library comes from another person. The counts cannot tell which. A swap between two of the 43 biopsies would show a second mismatch, and there is none. I keep 32-0004, as the paper does. It is one of the six group 4 biopsies, so if the library belongs to someone else, its MRI and pathology values do not describe the tissue that was sequenced, which matters for Step 4. It goes on the list of questions for the authors.

## Decisions I made

- Apply the authors' gene filter (total counts > 20) instead of the starter's default, to reproduce their analysis; on the deposited file it removes no gene.
- Library sizes (total counts assigned to genes, not sequencing depth; `library_sizes.png`, `library_sizes.csv`) follow sequencing batch: 59.8–80.2 million per biopsy in batch III, at most 47.7 million in the other batches. Within batches III and IV, the control and FSHD ranges overlap. Normalisation accounts for these differences in overall count scale (DESeq2 size factors in the rlog PCA and the differential-expression analysis; CPM only in the exploratory log-CPM comparison), so no biopsy is excluded for depth; the smallest is 32-0008 (17.6 million).
- PCA: reproduce the authors' method (DESeq2 rlog, 500 most variable genes) rather than the starter's log2(CPM + 1) on 2,000 genes, using a locked R environment for the rlog (`r/pixi.toml`, R 4.5.3, DESeq2 1.50.2); Python recomputes the PCA from the R matrix and stops if it disagrees with R. The log-CPM PCA is kept only as a comparison (`pca_logcpm_comparison.csv`).
- Sex: where a model needs sex (Step 3), use the sex inferred from expression for all 43 biopsies (`results/step2_qc/sex_check.csv`, column `inferred_sex`), not GEO's label. GEO gives no sex for the controls, and one FSHD label contradicts the expression data.
