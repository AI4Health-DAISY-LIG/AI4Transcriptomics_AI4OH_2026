# Step 2: quality control and exploration

_Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence._

## Filtering

| Filter      | Threshold I chose | Why                                                                                                                                                                | Samples and genes before | after                    |
| ----------- | ----------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------ | ------------------------ |
| min_samples | 9                 | Size of the smallest group (9 controls).                                                                                                                           | 43 samples, 27,506 genes | 43 samples, 19,258 genes |
| min_count   | 10                | Standard DESeq2 default. Per `filtering_genes_effect.csv`, any threshold between **9 and 14** keeps all four biomarker genes; above ~14 `LEUTX` starts to be lost. | (same rule, see above)   | -                        |

**Does the paper state its thresholds?** No. It only says counts come from GENCODE v24 (`summarizeOverlaps`, IntersectionStrict), the same annotation we use, so our counts are comparable.

**Do the four biomarker genes pass?** Yes, by construction (threshold chosen after looking at them). Samples with at least 10 counts: `KHDC1L` 20, `PRAMEF2` 17, `TRIM43` 14, `LEUTX` 10. `LEUTX` is the tight one: it would be lost at `min_count = 15` (8 samples). These genes are expressed only in a minority of FSHD samples (groups 3 and 4 hold 19), so a filter requiring expression in many samples would remove them.

**What the filter removes:** 8,248 genes (30%). It runs on all 43 samples, before the outlier exclusion, so the outliers do not decide which genes are kept.

## Sequencing depth

| Group     | n   | Mean (M reads) | Median (M) | Range (M) |
| :-------- | :-- | :------------- | :--------- | :-------- |
| Batch I   | 5   | 43.5           | 42.3       | 40.6–47.7 |
| Batch II  | 5   | 30.6           | 30.8       | 17.6–41.5 |
| Batch III | 22  | **71.8**       | 71.8       | 59.8–80.2 |
| Batch IV  | 11  | **27.4**       | 28.0       | 20.0–30.2 |
| Site 01   | 25  | 58.3           | 69.3       | 25.9–80.2 |
| Site 32   | 18  | 44.2           | 41.3       | 17.6–75.2 |
| Control   | 9   | 59.6           | 74.3       | 26.5–80.2 |
| FSHD      | 34  | 50.5           | 46.3       | 17.6–79.2 |

- **Depth follows batch, not condition.** Batch III is 2.6× deeper than batch IV, and their ranges do not overlap (59.8M vs ≤30.2M). Batch therefore carries a depth effect on top of any prep effect, and the two cannot be separated.
- **Controls are split across both extremes:** 6 in batch III (deep), 3 in batch IV (shallow). The mean (59.6M) sits well below the median (74.3M) for this reason.
- **Impact:** DESeq2 size factors correct library size, not detection power. Low-expressed genes (the DUX4 biomarkers) are easier to miss in the 3 shallow controls, so their zeros are less informative than the zeros of the deep ones. All samples have at least 17.6M reads, which is enough for gene-level DE.
- The paper does not report depth, so no comparison is possible.

## Tissue composition analysis

![](../results/step2_qc/tissue_composition_ratios.png)

Some samples are contaminated by blood `32-0014`, `32-0005` and `32-0008`. All are FSHD samples, so I believe that the patient muscle is more fragile and can be more impacted by the needle during the biopsy.

## PCA exploration

Setup: 1,500 most variable genes, three runs.

| Run               | Transformation              | Samples | Notes                                                         |
| :---------------- | :-------------------------- | :------ | :------------------------------------------------------------ |
| `RLOG2_PCA`       | rlog (Python approximation) | 43      | Closest to the paper's DESeq2 rlog. PC1 = 37% of variance.    |
| `VST_PCA`         | VST (PyDESeq2)              | 43      | PC1 = 51% of variance, but driven by two samples (see below). |
| `INLIERS_VST_PCA` | VST                         | 41      | Reference run for the rest of the analysis.                   |

### Outliers

| Sample                          | rlog PC1 / PC2 | VST PC1 / PC2 | Paper (Wang et al. 2019)                                                                                       |
| :------------------------------ | :------------- | :------------ | :------------------------------------------------------------------------------------------------------------- |
| `01-0037`                       | 28.1 / 28.2    | 128.0 / −65.0 | PCA outlier in FSHD group 1. Biomarkers very low or undetectable; ECM, inflammatory and immune genes enriched. |
| `32-0016`                       | 33.8 / 24.8    | 133.2 / −50.5 | Same as above, plus high adipocyte-gene content and high T1 signal (STIR+ but no DUX4 targets).                |
| _Next most extreme (`32-0009`)_ | 37.1 / 10.9    | 71.9 / 15.6   | -                                                                                                              |

- **Both transformations flag the same two samples**, as in the paper.
- **VST:** both sit at PC1 ≈ 130, nearly twice the next sample (71.9). The 51% variance on PC1 is therefore mostly these two samples, not a biological gradient.
- **rlog:** they are not isolated on PC1 (`32-0009` is higher) but are the top two on PC2 (28.2 and 24.8, next is 18.2). The outlier signal is weaker in rlog than in VST.
- **Watch list, not excluded:** `32-0014` is extreme on one PC in every run (rlog PC3 32.9; VST PC3 −42.7; inliers PC2 −46.0, about 2× the next sample). Not flagged by the paper. It's because of blood contamination foubd before.

### What drives the PCs (inliers, 41 samples)

Effect size is R² (numeric) or η² (categorical). FSHD-only variables use 31–32 samples, the others 41.

![](../results/step2_qc/INLIERS_VST_PCA/PC1_PC2_top_variables.png)

**Across methods:** `inflam` and `batch` are top-2 on PC1 in rlog and VST; `rnaseq.rank` and `rnaseq.score` are top-2 on PC2 in both. Beyond that the order differs, and the two outliers distort PC1 in both. Removing them gives the cleaner picture above.

### Reading of each figure

| Variable                                                | Reading                                                                                                               |
| :------------------------------------------------------ | :-------------------------------------------------------------------------------------------------------------------- |
| active                                                  | Same trend in all three runs, with some overlap. Moves with PC1.                                                      |
| age                                                     | No visible structure within FSHD. Controls missing.                                                                   |
| batch                                                   | III and IV overlap; I is offset. No controls in I–II, so batch effect and disease cannot be told apart.               |
| condition                                               | Controls form a compact group with little overlap with FSHD, except the low-DUX4 FSHD samples.                        |
| fat fraction                                            | Higher fat shifts expression, but it explains ~17–24% of PC1 (vs 63% for the DUX4 score). Not available for controls. |
| gender / infer_gender                                   | Two clusters on PC2, clearest with inferred sex and without outliers. Keep `infer_gender` as covariate.               |
| inflam, path.score, rnaseq.rank, rnaseq.score, stir, t1 | All follow PC1 (table above). One shared disease-severity axis.                                                       |
| instrument_model, tissue                                | Constant. Nothing to see.                                                                                             |

<details>
<summary><strong>All PCA figures</strong> (rlog · VST · VST inliers)</summary>

#### **active**

![](../results/step2_qc/RLOG2_PCA/pca_by_active.png) ![](../results/step2_qc/VST_PCA/pca_by_active.png) ![](../results/step2_qc/INLIERS_VST_PCA/pca_by_active.png)

#### **age**

![](../results/step2_qc/RLOG2_PCA/pca_by_age.png) ![](../results/step2_qc/VST_PCA/pca_by_age.png) ![](../results/step2_qc/INLIERS_VST_PCA/pca_by_age.png)

#### **batch**

![](../results/step2_qc/RLOG2_PCA/pca_by_batch.png) ![](../results/step2_qc/VST_PCA/pca_by_batch.png) ![](../results/step2_qc/INLIERS_VST_PCA/pca_by_batch.png)

#### **condition**

![](../results/step2_qc/RLOG2_PCA/pca_by_condition.png) ![](../results/step2_qc/VST_PCA/pca_by_condition.png) ![](../results/step2_qc/INLIERS_VST_PCA/pca_by_condition.png)

#### **fat fraction**

![](../results/step2_qc/RLOG2_PCA/pca_by_fat_fraction.png) ![](../results/step2_qc/VST_PCA/pca_by_fat_fraction.png) ![](../results/step2_qc/INLIERS_VST_PCA/pca_by_fat_fraction.png)

#### **gender**

![](../results/step2_qc/RLOG2_PCA/pca_by_gender.png) ![](../results/step2_qc/VST_PCA/pca_by_gender.png) ![](../results/step2_qc/INLIERS_VST_PCA/pca_by_gender.png)

#### **infer gender**

![](../results/step2_qc/RLOG2_PCA/pca_by_infer_gender.png) ![](../results/step2_qc/VST_PCA/pca_by_infer_gender.png) ![](../results/step2_qc/INLIERS_VST_PCA/pca_by_infer_gender.png)

#### **inflam**

![](../results/step2_qc/RLOG2_PCA/pca_by_inflam.png) ![](../results/step2_qc/VST_PCA/pca_by_inflam.png) ![](../results/step2_qc/INLIERS_VST_PCA/pca_by_inflam.png)

#### **instrument model**

![](../results/step2_qc/RLOG2_PCA/pca_by_instrument_model.png) ![](../results/step2_qc/VST_PCA/pca_by_instrument_model.png) ![](../results/step2_qc/INLIERS_VST_PCA/pca_by_instrument_model.png)

#### **path score**

![](../results/step2_qc/RLOG2_PCA/pca_by_path_score.png) ![](../results/step2_qc/VST_PCA/pca_by_path_score.png) ![](../results/step2_qc/INLIERS_VST_PCA/pca_by_path_score.png)

#### **rnaseq rank**

![](../results/step2_qc/RLOG2_PCA/pca_by_rnaseq_rank.png) ![](../results/step2_qc/VST_PCA/pca_by_rnaseq_rank.png) ![](../results/step2_qc/INLIERS_VST_PCA/pca_by_rnaseq_rank.png)

#### **rnaseq score**

![](../results/step2_qc/RLOG2_PCA/pca_by_rnaseq_score.png) ![](../results/step2_qc/VST_PCA/pca_by_rnaseq_score.png) ![](../results/step2_qc/INLIERS_VST_PCA/pca_by_rnaseq_score.png)

#### **stir**

![](../results/step2_qc/RLOG2_PCA/pca_by_stir.png) ![](../results/step2_qc/VST_PCA/pca_by_stir.png) ![](../results/step2_qc/INLIERS_VST_PCA/pca_by_stir.png)

#### **t1**

![](../results/step2_qc/RLOG2_PCA/pca_by_t1.png) ![](../results/step2_qc/VST_PCA/pca_by_t1.png) ![](../results/step2_qc/INLIERS_VST_PCA/pca_by_t1.png)

#### **tissue**

![](../results/step2_qc/RLOG2_PCA/pca_by_tissue.png) ![](../results/step2_qc/VST_PCA/pca_by_tissue.png) ![](../results/step2_qc/INLIERS_VST_PCA/pca_by_tissue.png)

</details>
Raw PC coordinates and per-variable effect-size tables for the three runs are in `results/step2_qc/`.
 
## Comparison with the paper (Wang et al. 2019)
 
**Where my results differ:** the Python rlog is only an approximation of DESeq2's, and 1,500 genes is a different gene set from the paper's. The exact PC axes therefore differ, but the conclusions about outliers and group structure agree.
 
## Samples I excluded, and why
 
| Sample    | Why                                                                                                                                                             |
| :-------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `01-0037` | Outlier in rlog and VST PCA; paper reports undetectable DUX4 biomarkers and an ECM/immune signature instead.                                                     |
| `32-0016` | Same, plus adipocyte-rich profile (paper: high T1 signal, STIR+ without DUX4 targets). Fat content would also distort bulk expression.                           |
 
Result: 43 → 41 samples (9 controls, 32 FSHD). **Cost:** the paper reads these two as early ECM/immune disease, so excluding them removes real biology. This is justified here because the analysis is built on the four DUX4 biomarkers, but it is tested in a sensitivity run with both samples included.
 
## Decisions I made
 
- **Exclude** `32-0016` and `01-0037` (see above). Reference PCA = VST on the 41 inliers.
- **Gene set for PCA:** 1,500 most variable genes. The rlog plot is then close to the paper's and the choice reduces noise without dropping informative genes.
- **Covariates:** keep `infer_gender` (PC2 η² 0.51, clean clusters). I connot adjust for `site`, `age` or `fat.fraction`: collinear with condition or missing for controls.
- **Batch:** it carries depth and disease effects that cannot be separated without controls in batches I–II. I will test with and without it.
- **Sensitivity designs for Step 3:** with/without the two outliers, with/without batch, batches III–IV only, and with/without `32-0014`.
- **Limits to report:** control age and kinship are unknown, depth differs 2.6× between batches III and IV, and site, batch and condition overlap.
