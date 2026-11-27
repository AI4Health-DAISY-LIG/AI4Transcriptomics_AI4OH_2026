# Blood, inflammation and the PCA outliers in FSHD muscle

GSE115650 · 27 November 2026 · Exploratory Track B results

## What the blood and neutrophil adjustments showed

The question in the primary Track B comparison is how whole FSHD muscle biopsies differ from controls. Inflammation, repair and changes in tissue composition can contribute to that difference. The primary model therefore stays as specified: batch, inferred sex and condition. Site remains a separate sensitivity. Track A, the reproduction, is unaffected.

- **Haemoglobin.** The percentage of reads assigned to globin markers is a proxy for erythroid RNA in a biopsy. Blood introduced during collection is one possible source. In the 41-biopsy analysis linked to course issue 24, median haemoglobin percentages were 0.087% in Rochester controls, 0.089% in Rochester FSHD and 0.367% in Seattle FSHD. There were no Seattle controls. Rochester used a modified Bergström needle and Seattle a regular one. Needle, site and other site-specific procedures cannot be separated here. Adjusting for site accounts for mean site differences. A haemoglobin covariate can still account for variation between biopsies at the same site. [Site-stratified results](https://github.com/Tan9uy/AI4Transcriptomics_AI4OH_2026/blob/ce23bcf0320fcfe6f06f7bd838bc546458ef0c73/practical/results/step3_de/blood_by_condition.csv) · [Wang et al., Methods](https://doi.org/10.1093/hmg/ddy364)
- **Neutrophils.** A neutrophil-marker RNA percentage can reflect inflammatory cells in diseased muscle as well as blood in the specimen. Histology records inflammatory infiltrates of any cell type. In the earlier comparison of the 32 sequenced FSHD biopsies with histology, the median neutrophil-marker percentage was 0.00379% in the 11 with infiltrates and 0.00222% in the 21 without (unadjusted Mann–Whitney p = 0.18). The highest RNA marker value occurred in 32-0014, which had histology grade 0. These measurements cannot be substituted for each other. The RNA and histology also come from different pieces of a biopsy. [Course issue 24](https://github.com/AI4Health-DAISY-LIG/AI4Transcriptomics_AI4OH_2026/issues/24) · [Paper](https://doi.org/10.1093/hmg/ddy364)
- **What changed after adjustment.** The batch-, sex- and site-adjusted model linked to course issue 24 used **32 FSHD and 9 controls**. Adding both markers changed the number of significant genes from 824 to 901. Of these, 765 were shared, all in the same direction. Genome-wide log2 fold-change ranks correlated at 0.996. This describes agreement between the two disease-effect estimates. It does not measure whether the proxies are still associated with diagnosis. The earlier partial-rank checks gave ρ = 0.197 for haemoglobin and 0.220 for neutrophils after accounting for batch, sex and site. [Model comparison](https://github.com/Tan9uy/AI4Transcriptomics_AI4OH_2026/blob/ce23bcf0320fcfe6f06f7bd838bc546458ef0c73/practical/results/step3_de/comparison.csv)
- **Track B decision.** Account for systematic differences between collection sites in the sensitivity models. Keep haemoglobin and neutrophil RNA as descriptive measures. The prespecified histology contrasts retain immune and repair signals that may form part of the disease response. Per-biopsy haemoglobin addresses a separate question about blood variation within sites and remains an optional follow-up.

## Two genes that illustrate the distinction

| Gene | Earlier exploratory checks using the published issue-24 design | Results in this Track B analysis | Interpretation |
|---|---|---|---|
| AREG (amphiregulin) | Overall rank was about 77 without either marker, 25–26 with haemoglobin only, 50–52 with neutrophils only and 10 with both. The largest single-marker rank gain came from haemoglobin. The joint model reported log2FC +6.27. | FSHD–control log2FC +4.35 in the primary model and +5.33 in Rochester alone. The within-FSHD DUX4-score association was weak. | AREG is a candidate repair-associated signal. Mouse muscle experiments support a role for amphiregulin in repair. Its rank in bulk human biopsies does not identify the cell producing it or put it upstream of DUX4. |
| MTRNR2L1 | Overall FSHD–control log2FC changed from −1.67 with neither marker to −3.02 with haemoglobin alone and −2.90 with both. The overall effect did not pass FDR 0.05. | FSHD–control log2FC −1.06 (q = 0.57), −1.56 with site (q = 0.55), −1.99 in Rochester alone (q ≈ 0.38). Within the 30 retained FSHD biopsies, the DUX4-target score slope was −0.535 per unit (q = 0.023). With all 34 FSHD biopsies the slope was −0.267, just below the prespecified 50% effect-retention threshold. | These are distinct comparisons. The within-FSHD association is sensitive to sample selection. The whole-biopsy measurement cannot identify the cell of origin, mitochondrial RNA release, or production of a protective peptide. |

The single-marker fits reported above were earlier exploratory checks using the published issue-24 design. The published comparison adds both markers together. The existing Track B primary, site, Rochester and score fits answer different questions. The four additional Track B fits below use the project’s own counts, cohorts and covariates. The marker-adjusted findings from issue 24 come from a different cohort and design. [Issue-24 worksheet (Tan9uy)](https://github.com/Tan9uy/AI4Transcriptomics_AI4OH_2026/blob/ce23bcf0320fcfe6f06f7bd838bc546458ef0c73/practical/worksheets/step3_de.md) · [Track B results](https://github.com/pierrebaumierresearch/AI4Transcriptomics_AI4OH_2026/tree/c5f71d287da9882ee77830301062886ca430dd0c/practical/results/step3_de/track_b) · [Amphiregulin and muscle repair](https://pubmed.ncbi.nlm.nih.gov/24315098/)

MTRNR2L1 is a nuclear relative of the mitochondrial 16S rRNA gene MT-RNR2. NCBI annotates MTRNR2L1 as a pseudogene on chromosome 17. Sequence similarity also makes locus-specific read assignment relevant. There is independent evidence that mitochondrial function changes during development: pyruvate dehydrogenase transiently localizes to the nucleus around human embryonic genome activation. In human muscle-cell models, DUX4 decreased complex I-linked respiration in myoblasts and increased it in myotubes. These findings support a stage-dependent mitochondrial response to DUX4. They leave the involvement of MTRNR2L1 or a neural developmental programme in FSHD unresolved. [Locus annotation](https://www.ncbi.nlm.nih.gov/gene/100462977) · [Nagaraj et al.](https://pubmed.ncbi.nlm.nih.gov/28086092/) · [Heher et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8899416/)

## PCA label correction and interpretation

| Biopsy | Role in the analysis | Evidence at biopsy |
|---|---|---|
| 32-0009 | PCA label and low-mapping exclusion | STIR 2, pathology 10, inflammation 3, DUX4-target score 10.63. Strong evidence of active injury despite low mapping. |
| 32-0008 | PCA label and low-mapping exclusion | STIR 2, pathology 2, inflammation 0, DUX4-target score 7.54. Elevated DUX4 targets without recorded inflammatory infiltrates. |
| 01-0037 | PCA label and low-mapping exclusion | STIR 0, pathology 0, inflammation 0, DUX4-target score 0.024. One of the two group 1 outliers named in the paper. |
| 32-0016 | Fourth PCA label, retained in primary Track B | STIR 2, MRI fat fraction 42.8%, DUX4-target score 0.51, histology unavailable. The paper’s other group 1 outlier. |
| 32-0004 | Fourth Track B exclusion, restored in the 30-to-34 FSHD sensitivity | RNA expression disagrees with recorded sex. The linked record gives STIR 2, pathology 6, inflammation 1 and DUX4-target score 16.25. The discrepancy limits confidence in linking the RNA to that participant’s phenotype. |

**Caption correction.** “Labels identify the three low-mapping libraries (32-0008, 01-0037 and 32-0009) and the fat-rich biopsy 32-0016.” Only the first three were removed for low mapping. The fourth exclusion in the Track B primary analysis was 32-0004. The four labels are not a biological group. An actively affected muscle is plausible for 32-0009. One baseline biopsy cannot date a flare or establish a common acute episode. [GEO-derived sample table](https://github.com/pierrebaumierresearch/AI4Transcriptomics_AI4OH_2026/blob/c5f71d287da9882ee77830301062886ca430dd0c/practical/results/step1_audit/sample_table.csv) · [Wang et al., Figure 3](https://doi.org/10.1093/hmg/ddy364)

## Expression by histological inflammation

**Objective.** Distinguish expression differences already present in FSHD muscle without observed infiltrates from differences associated with infiltrates within FSHD. Groups were defined from histology (`char:inflam`) independently of the expression matrix. Controls form their own reference category. All nine lack a GEO inflammation grade.

| Population | Grade 0 | Grade 1 | Grade 2 | Grade 3 | Missing |
|---|---:|---:|---:|---:|---:|
| All 34 sequenced FSHD biopsies | 21 | 7 | 2 | 2 | 2 |
| 30 FSHD biopsies retained in primary Track B | 19 | 6 | 2 | 1 | 2 |

The paper’s 11/33 refers to the histology cohort. The sequenced cohort has 11/32 with recorded infiltrates. The authors also found active myopathic changes in 18/33 when their definition included necrosis and regeneration. An inflammation grade of zero therefore does not imply inactive disease. [Histology results](https://pmc.ncbi.nlm.nih.gov/articles/PMC6337697/) · [GEO-derived sample table](https://github.com/pierrebaumierresearch/AI4Transcriptomics_AI4OH_2026/blob/c5f71d287da9882ee77830301062886ca430dd0c/practical/results/step1_audit/sample_table.csv)

**Models fitted.** All three contrasts were estimated on the same 37 biopsies using a model with batch and inferred sex (`~ batch + sex + infiltrates`) and its site-adjusted sensitivity (`~ batch + sex + site + infiltrates`). `sex` denotes expression-inferred sex. The normalization and DE settings follow the existing Track B fits. The group factor is control (9), FSHD grade 0 (19), or FSHD grade 1–3 (9). Six of the nine biopsies with infiltrates were grade 1. Two FSHD biopsies lack a grade. The site-adjusted design has full rank (8/8 columns). Both sites contribute to the within-FSHD contrast. Rochester provides the only direct case–control comparison within a site. The analysis assumes a common additive site effect.

| Contrast, positive log2FC indicates higher expression in the first group | n | Question |
|---|---:|---|
| FSHD grade 0 vs control | 19 vs 9 | Which differences are already detectable in FSHD biopsies with no histological infiltrates? |
| FSHD grade 1–3 vs FSHD grade 0 | 9 vs 19 | Which changes accompany histological inflammation among people with FSHD? |
| FSHD grade 1–3 vs control | 9 vs 9 | How does the overall expression profile of FSHD biopsies with infiltrates differ from control muscle? |

**Observed contrasts.** The first two numeric columns count genes at Benjamini–Hochberg q < 0.05 within each contrast. The final columns give the site-adjusted log2 fold change and q for the two genes discussed above. The contrasts share samples and were explored in the same dataset.

| Contrast | Significant without site | Significant with site | MTRNR2L1 with site | AREG with site |
|---|---:|---:|---:|---:|
| FSHD grade 0 vs controls | 938 | 1,243 | −1.93 (q = 0.389) | +5.33 (q = 0.00112) |
| FSHD grades 1–3 vs FSHD grade 0 | 4,040 | 4,635 | −4.71 (q = 0.00279) | −0.91 (q = 0.579) |
| FSHD grades 1–3 vs controls | 7,514 | 7,227 | −6.64 (q = 0.000854) | +4.42 (q = 0.0164) |

The direct comparison within FSHD is the informative one for infiltrates. MTRNR2L1 expression was lower in biopsies with infiltrates both before site adjustment (−3.42, SE 1.25, q = 0.0366) and after it (−4.71, SE 1.25, q = 0.00279). AREG was elevated against controls even in grade-0 biopsies. Its within-FSHD contrast had SE 1.10 after adjustment, with little evidence of a difference by infiltrate status. These are bulk RNA associations. Multiple cell types, disease severity and sampling may contribute. Both focal genes passed all three recorded optimizer checks in these fits.

**Why keep all 19 grade-0 biopsies?** All 19 contribute information about expression and variation in this group. The two case–control contrasts share the same nine controls and have different precision. Their numbers of significant genes therefore cannot establish which FSHD group differs more from controls. The direct 9-vs-19 contrast estimates the expression difference associated with infiltrates within FSHD.

**What the earlier analysis did.** The primary Track B disease comparison used 30 FSHD biopsies and 9 controls. The four published RNA groups contributed 9, 5, 11 and 5 retained FSHD biopsies, each compared with the same nine controls. The deposited RNA labels determined membership. G1 happened to contain nine retained biopsies. The separate Step 4 bootstrap resampled FSHD biopsies to estimate uncertainty in MRI and pathology correlations. [Step 3 cohort and contrasts](https://github.com/pierrebaumierresearch/AI4Transcriptomics_AI4OH_2026/blob/c5f71d287da9882ee77830301062886ca430dd0c/practical/worksheets/step3_de.md#L30) · [Step 4 bootstrap](https://github.com/pierrebaumierresearch/AI4Transcriptomics_AI4OH_2026/blob/c5f71d287da9882ee77830301062886ca430dd0c/practical/scripts/step4_track_b.py#L245-L253)

**Optional equal-size sensitivity.** To examine how sample size and the choice of biopsies affect the results, repeatedly draw 9 of the 19 grade-0 biopsies without replacement within each draw. Keep the same nine controls and nine FSHD biopsies with infiltrates. Preserve batch and sex representation where feasible, use the same model settings and verify that each design is estimable. Summarise how much the log2 fold changes vary across draws, particularly for AREG and MTRNR2L1. The draws overlap and share controls, so they measure subset sensitivity and cannot count as independent replications. This check is proposed.

A literal **low versus high severity** comparison would group grades 0–1 (25) and 2–3 (3), plus the same 9 controls. All three grade 2–3 biopsies retained in Track B are from male participants. Keep their gene-level estimates or count plots descriptive. Restoring excluded libraries would give only four grade 2–3 cases among all sequenced biopsies with histology.

Site enters the histology sensitivity, the DUX4-score sensitivity and the RNA-group sensitivity. Neutrophil-marker RNA is descriptive because immune recruitment may be part of the tissue response. Individual haemoglobin adjustment remains a separate, optional analysis of within-site variation.

## What changed with site adjustment

| Site in the 37-biopsy analysis | Controls | FSHD grade 0 | FSHD grades 1–3 | Median haemoglobin in FSHD | Median neutrophil-marker RNA in FSHD |
|---|---:|---:|---:|---:|---:|
| Rochester | 9 | 12 | 3 | 0.0889% | 0.00185% |
| Seattle | 0 | 7 | 6 | 0.2488% | 0.00336% |

These percentages are shares of assigned gene counts. Seattle FSHD biopsies had 2.80 times the Rochester median haemoglobin and 1.82 times the median neutrophil signal. Haemoglobin and neutrophil RNA also rose together within each site (Spearman ρ = 0.77 in Rochester and 0.78 in Seattle). Circulating blood could contribute to both measurements. A site term captures average differences in collection, tissue and other site-associated factors. It cannot isolate the effect of biopsy needle or blood specifically. A site-average haemoglobin variable duplicates the information in the two-level site factor. [Per-biopsy marker measurements](https://github.com/Tan9uy/AI4Transcriptomics_AI4OH_2026/blob/ce23bcf0320fcfe6f06f7bd838bc546458ef0c73/practical/results/step2_qc/tissue_composition_ratios.csv) · [Histology and identifiers](https://github.com/pierrebaumierresearch/AI4Transcriptomics_AI4OH_2026/blob/c5f71d287da9882ee77830301062886ca430dd0c/practical/results/step1_audit/sample_table.csv)

- **Overall FSHD–control expression, already tested.** The original 39-biopsy comparison has a saved site-adjusted fit and a Rochester-only fit. Adding site gave 583 significant genes versus 2,577 in the primary model. Of those 2,577, 2,416 kept the same direction and at least half the primary effect magnitude. Reuse these results for the site question. [Fit definitions](https://github.com/pierrebaumierresearch/AI4Transcriptomics_AI4OH_2026/blob/c5f71d287da9882ee77830301062886ca430dd0c/practical/scripts/step3_differential_expression.py#L96-L127) · [Paired results](https://github.com/pierrebaumierresearch/AI4Transcriptomics_AI4OH_2026/blob/c5f71d287da9882ee77830301062886ca430dd0c/practical/results/step3_de/track_b/summary.txt#L62-L71)
- **Within-FSHD DUX4-target score, retrospective sensitivity.** The model was `~ batch + sex + site + score` on the original 30 biopsies. There are 15 cases per site, with overlapping score ranges (Rochester 0.14–15.98, Seattle 0.04–21.23). The design is full rank (7/7). MTRNR2L1 had a slope of −0.535 without site (SE 0.138, q = 0.0229) and −0.534 with site (SE 0.140, q = 0.0258) per score unit. AREG had slopes of +0.185 (q = 0.566) and +0.180 (q = 0.547), respectively. The model gave 144 discoveries without site and 163 with site. The score measures expression of four DUX4 targets.
- **RNA groups, retrospective sensitivity.** The model was `~ batch + sex + site + group` on the original 39 biopsies. All four FSHD groups occur at both sites, and the design is full rank (10/10). At q < 0.05, G1, G2, G3 and G4 against controls gave 1,002, 518, 639 and 5,409 discoveries without site, and 28, 230, 283 and 4,446 with site. G4 vs G1 increased from 458 to 749. For MTRNR2L1 the G4–G1 estimate weakened from −6.80 (q = 0.0396) to −4.86 (q = 0.143). AREG had little G4–G1 separation in either fit (−0.37 and −0.25, both q > 0.92). G4 contains only two Rochester and three Seattle biopsies, all male. Full rank permits estimation but leaves this contrast sensitive to individual biopsies. These two retrospective fits are exploratory additions to Track B. The original reproduction and prespecified gene-selection criteria retain their original definitions.
- **Neutrophil-marker RNA versus histological infiltrates, retrospective check completed.** Among all 32 sequenced FSHD biopsies with histology, the earlier pooled Mann–Whitney p-value was 0.18. Comparing ranks within each site gave a two-sided p-value of 0.327. Among the 28 retained FSHD biopsies with histology, the corresponding site-stratified p-value was 0.293. In the latter set, median marker percentages for grade 0 versus grades 1–3 were 0.00187% versus 0.00185% in Rochester (12 vs 3) and 0.00220% versus 0.00546% in Seattle (7 vs 6). The estimates are too imprecise to establish different effects by site. The exact test sums the within-site rank sums for the infiltrate-positive group, permutes histology labels within site at fixed group sizes and counts both tails by distance from the null mean. This check accounts for site only and uses the same per-biopsy marker and histology sources linked above.

## Reproducing the follow-up

From the practical directory, the four additional fits can be selected in one run. The full Track B workflow also includes them. The original robustness criteria and Translator selection still use their published fit lists.

```bash
uv run --locked python scripts/step3_differential_expression.py --fit score_site --fit groups_site --fit inflammation --fit inflammation_site
```

The descriptive site check uses the public per-biopsy marker percentages from the analysis linked in course issue 24. Downloading its pinned CSV into the ignored data directory is required only to regenerate the descriptive report. The four expression fits use this project’s Step 1 and Step 2 inputs directly.

```bash
mkdir -p data/raw
curl --fail --location --output data/raw/tissue_composition_ratios.csv https://raw.githubusercontent.com/Tan9uy/AI4Transcriptomics_AI4OH_2026/ce23bcf0320fcfe6f06f7bd838bc546458ef0c73/practical/results/step2_qc/tissue_composition_ratios.csv
uv run --locked python scripts/check_blood_markers_by_site.py data/raw/tissue_composition_ratios.csv
```

The report is saved as `results/step3_de/track_b/blood_marker_site_check.txt`. It records the marker input checksum, site-specific medians and correlations, and both exact stratified tests.

**Numerical checks and provenance.** The per-gene tables cover all 27,506 deposited count features. The four design matrices were full rank. The unscaled site-plus-score fit encountered a singular matrix during dispersion estimation. Centering the score and dividing by its within-cohort SD (5.5262) resolved the numerical failure. The Wald contrast divides the fitted coefficient and its SE by the same SD to express both per original GEO score unit. All 163 discoveries in that fit passed the three reported PyDESeq2 optimizer checks. Some categorical contrasts include discoveries with a failed optimizer flag (17 of 4,635 in the site-adjusted within-FSHD inflammation contrast). The reported totals include these genes. The saved per-gene flags allow them to be identified. For Cook’s replacement, PyDESeq2 reports the flags from the initial fit.

The four full-matrix fits were executed independently with Python 3.12.14 and the course-locked PyDESeq2 0.5.4 on the public counts. All 43 library totals were cross-checked against the project’s saved QC. Gene names, types and chromosomes were recovered by versioned gene ID from the existing primary result table for all 27,506 features. These columns label the output and do not enter the fits. The local reproduction command above uses the project’s original gene annotation files. Source baseline: [project commit c5f71d2](https://github.com/pierrebaumierresearch/AI4Transcriptomics_AI4OH_2026/tree/c5f71d287da9882ee77830301062886ca430dd0c).
