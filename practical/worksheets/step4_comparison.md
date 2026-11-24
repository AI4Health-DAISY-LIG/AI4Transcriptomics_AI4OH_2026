# Step 4: comparison with the paper

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Step 4 has two tracks, as Step 3:

- **Track A, the reproduction.** The authors' score, groups, Figure 2, correlations and DUX4-target counts are recomputed in R from their DESeq2 object on GEO, following their code (`r/score_reproduction.R`, `results/step4_comparison/track_a/`).
- **Course parts and Track B, independent.** The course parts are the Python score (`biomarker_score()`) and the gene claims (`scripts/step4_compare_with_paper.py`, `paper/paper_claims.csv`, `results/step4_comparison/`). Track B is an alternative DUX4 score and a tissue-state analysis of the Translator panel genes (`r/vst_track_b.R`, `scripts/step4_track_b.py`, `results/step4_comparison/track_b/`).

## The biomarker score

How the paper computes the score, in your words, with the sentence you rely on:

The Methods say the score is obtained by "taking the sum of the four biomarkers' log2 difference from the control average". The four biomarkers are LEUTX, KHDC1L, PRAMEF2 and TRIM43, on DESeq2's regularised log (rlog) scale. For each biopsy, take each marker's rlog value minus its mean in the nine controls, then add the four differences. The authors' code adds three details the paper does not give (`misc_tools.R` L84–L129 and L185–L198):
- the rlog is fitted on the four-gene subset only, with a local dispersion fit;
- the size factors come from all genes;
- the score deposited on GEO is that sum minus its minimum over the 43 biopsies.

The four groups come from hierarchical clustering of the four rlog values (Euclidean distance, complete linkage, k = 4), followed by a relabelling by hand (`addColDataToDDS.R` L14–L31): cluster 2 → group 1; rows 20–24 of the score-sorted table → group 2; cluster 1 → group 3; clusters 3 and 4 → group 4.

What you could not apply as stated, and the substitute you chose:

PyDESeq2 has no rlog. Substitute: a blind PyDESeq2 VST (design `~1`), fitted on the full count matrix of the 43 biopsies; the paper's formula is then applied to the four markers. A full-matrix fit gives a stable dispersion trend; fitting four genes, as the authors did for rlog, rests the trend on four points. PyDESeq2's size factors equal the deposited ones (max relative difference 4.7e-15). Sensitivity: `log2(normalised count + 0.5)`, the authors' own commented-out alternative (`misc_tools.R` L106–L108). The exact rlog score is reproduced in R (Track A).

| Item | Paper or GEO | Me | Agree? |
|------|--------------|-----|--------|
| Biopsies in groups 1, 2, 3 and 4 | Paper and GEO: 10, 5, 13, 6 | Track A: 10, 5, 13, 6, the same biopsies as GEO, and the 9 controls in group 1 (`track_a/groups_crosstab.txt`). The raw tree cut is 24/13/5/1 (controls included); the authors' relabelling, applied as written, gives the published groups. Clustering the Python VST values instead puts 3 group 3 biopsies with the controls (`score_clusters.txt`) | Yes (Track A); the clusters depend on the transformation |
| My score against the deposited score (Spearman) | GEO `rnaseq.score`, 34 FSHD | Python VST: 0.996; log2(normalised + 0.5): 0.994; values differ by up to 5.1 (VST and rlog scales differ). Track A: ρ = 1, max difference 4.7e-05 | Ranks yes; values only in Track A |
| Biomarker expression in the controls | Results: zero or < 0.1 RPKM | Max FPKM in the 9 controls: PRAMEF2 0.000, TRIM43 0.020, KHDC1L 0.019, LEUTX 0.016; 0 of 36 values ≥ 0.1 (`track_a/controls_fpkm.csv`) | Yes |
| STIR-positive biopsies that express DUX4 targets | 20/23 (87%); 91% without 32-0016 | 20/23 and 20/22 with the paper's definition (groups 2–4); 15/23 if "expressing" means at least one marker above relative rlog 2 (Figure 1 dot size) | Yes, with the paper's definition |
| Mean score in MRI-normal and in STIR-positive muscles | 1.4 (SD 2.4) and 7.7 (SD 5.7), P < 0.001, test not named | 1.383 and 7.666 (n = 10 and 23; one biopsy, 01-0022-1, is STIR 0 but T1 3, in neither group). SD: population 2.45 and 5.73, sample 2.58 and 5.85 (`track_a/figure2.csv`). Two-sided candidate tests: Welch p = 0.00017, exact Wilcoxon rank-sum p = 0.00064 (no ties), normal approximation p = 0.0012 | Means yes; SDs partly (5.7 only); the test cannot be identified |
| Score against STIR rating (Spearman) | Results 0.61 (P < 0.001); Fig. S7 0.54 | 0.555 (n = 34, p = 0.00065); group: 0.544 against 0.63 | No |
| Score against pathology score (Spearman) | Results 0.63 (P < 0.001) | 0.584 (n = 32, p = 0.00045) | No |

All correlations: `track_a/correlations.csv`. The same pattern holds for T1 (0.392 against 0.49) and fat fraction (0.226 against 0.34). Two labelled diagnostics:
- **Without 32-0016**, which the paper itself sets aside for the 91% figure: the values move up (STIR 0.590, T1 0.456).
- **On the 32 biopsies with pathology**: fat fraction is recovered to rounding (0.342 and 0.321 against 0.34 and 0.32), while T1 (0.503 and 0.485) and STIR (0.635 and 0.653) move past the paper's values.

No subset recovers all the reported values, so the subset the authors used cannot be identified from the deposited data. GEO STIR and T1 match the whole-muscle ratings of Table S1 (Step 1). Which level the GEO fat fraction belongs to is unclear.

The paper gives several of its correlations a p-value. Count the correlations it reports, and say whether the authors corrected for that number.

The Results text gives 12 correlations with a p-value:
- 7 between the RNA measures and pathology or MRI (pathology 1; T1, fat fraction and STIR 2 each, with score and with group);
- 5 between MRI and pathology.

The Figure 2 test adds a thirteenth p-value. Supplementary Figures S6 and S7 add a full correlation matrix. The Methods mention Spearman correlations only, and the paper names no correction for multiple testing anywhere.

## Differential expression

Summarise what the paper reports: which comparisons, which method and software versions, which thresholds, which genes. Then the same items from your own `results/`.

The paper compares each RNA group with the nine controls: the two group 1 outliers (Table S4), group 1 (S5), group 2 (S6), group 3 (S8) and group 4 (S10).
- **Method:** DESeq2 and GOSeq (Bioconductor); the paper names no versions, and the authors' scripts load R 3.4.
- **Threshold:** padj < 0.05 and |log2FC| ≥ 2.
- **Genes highlighted:** 47 of the 54 DUX4 candidates differ in group 3 and all 54 in group 4. Figures 4 and 5 show 18 immune, extracellular-matrix and other genes that separate the two outliers from the controls and, according to the captions, rise progressively in groups 2–4.

| Item | Paper | Me | Agree? |
|------|-------|-----|--------|
| Samples used | Each group with the 9 controls; all 43 biopsies | Track A: the same subsets. Track B: 39 biopsies (3 low-mapping libraries and 32-0004 excluded) | Yes (Track A) |
| Filtering | None in the Methods; code: total counts > 20, applied before deposit | Same deposited matrix, no extra filter | Yes |
| DE method and version | DESeq2, version not stated (R 3.4 in the scripts) | Track A: DESeq2 1.50.2, R 4.5.3; Track B: PyDESeq2 0.5.4 | Track A: same tables (Step 3) |
| Design and covariates | `~ pheno_type`, fitted on each subset, no covariates | Track A: the same. Track B: `~ batch + sex + group`, all groups in one model | Yes (Track A) |
| DUX4 target genes that differ in group 3 | 47 of 54 | 47 of 54 (Track A, `track_a/dux4_targets_de.csv`). The authors' 55 candidates include TRIM43CP, absent from the deposited counts, hence 54 | Yes |
| DUX4 target genes that differ in group 4 | 54 of 54 | 54 of 54 | Yes |
| Other genes the paper highlights, and their direction | 18 genes of Figures 4–5, up in the outliers and progressively in groups 2–4 | All 18 are up in every contrast. At padj < 0.05: Track A 17/18 in S4 (CTLA4 is not in Table S4), 17/18 in S6, 18/18 in S8 and S10; Track B 10, 11 and 15 of 18 in groups 2, 3 and 4 (`claims_vs_results_*.csv`). Group medians on the figures' scale, log10(FPKM + 0.5): all 18 above the controls in groups 2–4, but only 10/18 ordered G2 ≤ G3 ≤ G4; the rise is mostly at group 4 (`figure_genes_group_medians.csv`) | Direction yes; "progressive" only for 10/18 |

## What is missing to reproduce the paper

List each missing element and classify it:

- stated in the paper but not matching the deposited data;
- not stated anywhere, decided by me (say what you chose);
- stated and matching, but I could not apply it (say why).

**Stated in the paper but not matching the deposited data**
- **Spearman correlations.** The reported values are all higher than the deposited data give: STIR 0.61 against 0.555, T1 0.49 against 0.392, fat fraction 0.34 against 0.226, pathology 0.63 against 0.584. The STIR value also differs between the Results (0.61) and Fig. S7 (0.54).
- **Figure 2 SD for MRI-normal muscles:** 2.4, against 2.45 (population) or 2.58 (sample).
- **The Methods attribute the four groups to hierarchical clustering.** The clustering alone gives other groups (raw cut 24/13/5/1). The published groups need a relabelling by hand that the paper does not describe.

**Not stated anywhere, decided by me**
- **The test behind Figure 2's P < 0.001.** I report Welch's t-test and the exact Wilcoxon rank-sum test, both two-sided and both below 0.001, and the normal approximation (0.0012); none is confirmed.
- **The GEO score convention** (the sum minus its cohort minimum) and the rlog details (four-gene subset, local fit). Both are taken from the authors' code.
- **"Expressing DUX4 targets" = groups 2–4.** This is consistent with the 20/23; I report the per-marker definition next to it.
- **MRI level of the GEO variables.** STIR and T1 match the whole-muscle ratings of Table S1; the level of the fat fraction is unclear.
- **Software versions.** The paper names none. Track A runs DESeq2 1.50.2, and it reproduces the score to 5e-05 and every DE list exactly.

**Stated and matching, but I could not apply it**
- **rlog in Python.** PyDESeq2 has none; the VST substitute preserves the ranks (ρ = 0.996) but not the values, and changes the clusters. Track A applies rlog in R.

## What the design allows the paper to claim

The paper compares FSHD biopsies with control biopsies. List the differences between the two sets of biopsies other than the disease, and say for each one whether you could adjust for it with the deposited data.

- **Collection and banking.** The controls are banked Rochester biopsies, collected under the same conditions according to the Methods, but storage times are not reported. Not adjustable.
- **Site.** All controls are from Rochester; the FSHD biopsies come from Rochester and Seattle. Partly adjustable: a site term (learned from the FSHD biopsies) and a Rochester-only fit in Step 3 Track B.
- **Muscle.** All controls are quadriceps, but only 6 of the 34 FSHD biopsies; the muscle is not in GEO (Table S1 only). Not adjusted. With Table S1, a muscle term would rest on those 6 biopsies.
- **Age.** Controls mean 35 years (Results) against about 54 for FSHD; GEO gives no age for the controls. Not adjustable.
- **Sex.** 6 women and 3 men among controls (Results); group 4 is all men by inferred sex; GEO gives no sex for the controls. Adjustable with sex inferred from expression (Step 3 Track B).
- **Sequencing batch.** Batches I and II contain only FSHD biopsies. Adjustable in part (Step 3 Track B), but partly confounded with condition.
- **Library quality.** Three FSHD libraries map below 90%. Adjustable by exclusion (Track B and the sensitivity runs).
- **Relatedness.** The controls are spouses, parents or siblings of patients. Not recorded, not adjustable.
- **Selection of the biopsied muscle.** STIR-positive muscles were preferentially biopsied. Not adjustable; it limits the claims to this kind of muscle.
- **Tissue composition** (fat, fibrosis, inflammation). It is part of the disease; not adjusted, and not adjustable without separating it from the disease.

Then take one conclusion of the paper and say how far your analysis supports it.

**Conclusion: STIR-positive muscles are more likely to express DUX4 targets, and the STIR rating is associated with the RNA score.**

**Supported:**
- 20/23 STIR-positive biopsies are in groups 2–4.
- The STIR–score association is clearly present: ρ = 0.555, 95% bootstrap interval [0.27, 0.77], n = 34. It is weaker than the reported 0.61.
- A broad DUX4-target score without the four markers gives a similar positive association: Yao 2014 Table S1, 137 features after removing the five features carrying the four marker names, ρ = 0.569. The change against the deposited score is +0.013 [−0.115, +0.154]: the interval spans zero and allows a modest change in either direction, so the two scores are not shown to be equivalent.
- The same holds without the three low-mapping libraries, and without 32-0004 (`track_b/delta_rho.csv`).

**Limits:**
- 34 biopsies, in muscles chosen for their MRI signal, at one time point.
- No correction for the number of correlations tested.
- An association between MRI and RNA cannot say which comes first.
- The sharper MRI-normal against STIR-positive contrast (Figure 2) reproduces in its means, but its test cannot be identified.

## What I would ask the authors

Three questions at most.

1. **Correlations.** Which biopsies and which MRI level (whole muscle or biopsy region) gave the Spearman values in the Results (STIR 0.61, T1 0.49, fat fraction 0.34, pathology 0.63)? The deposited data give 0.555, 0.392, 0.226 and 0.584, and Fig. S7 gives 0.54 for STIR.
2. **Figure 2.** Which test gave P < 0.001, and is the SD a population or a sample SD?
3. **32-0004.** Is it recorded as a woman while its expression profile is male because of a labelling error or a sample exchange?

## Choices I made and would defend

- **Two tracks.** Track A reproduces the score, groups and statistics exactly with the authors' object and code, so any remaining gap lies in the paper's reported values, not in the method. Track B uses independent choices.
- **Python substitute for rlog:** a blind, full-matrix VST, with the authors' own `log2(n + 0.5)` alternative as a check. Both preserve the ranks (0.996 and 0.994).
- **Candidate tests reported, none chosen.** The paper does not name its test, and the threshold P < 0.001 does not discriminate.
- **The paper's definition first, mine next to it**: groups 2–4 for "expressing DUX4 targets", with the per-marker count beside it.
- **Alternative DUX4 score (Track B).**
  - Construction: the mean of control-median-centred VST values over the Yao 2014 Table S1 genes, without the four markers. The authors' own Table S3 mapping is a second panel, and floored z-scores a sensitivity.
  - Comparison: against the deposited score as baseline, with paired bootstrap intervals (5,000 resamples, seed 20261124) and fixed wording (weaker, stronger, difference unresolved) instead of a binary rule.
  - Result: every change is "difference unresolved" and small (`track_b/delta_rho.csv`, `panels.csv`).
- **Tissue state of the Translator panel (Track B, exploratory).** The 30 panel genes, plus POSTN chosen from the literature, are correlated with MRI and pathology, with and without adjustment for the DUX4 score. The panel's two groups were fixed from Step 3 before this analysis (`track_b/panel_genes.csv`, `tissue_state.csv`).
  - **SOX11** follows the DUX4 score most closely (ρ = 0.76). Beyond the score, it still follows the pathology score: partial ρ 0.48 [0.12, 0.74] on 32 biopsies, and 0.33 [−0.05, 0.66] on the 28 of Track B. It does not follow STIR: −0.12 [−0.40, 0.22].
  - **RUNX1** is most tied to inflammation, even beyond the score: partial ρ 0.55 [0.21, 0.74].
  - **RUNX1 and POSTN** are strongly co-expressed: ρ 0.79 [0.62, 0.88].
  - In bulk tissue this fits a shared remodelling compartment. It shows neither the cells involved nor any regulation.
  - Every table here is descriptive and unadjusted for its 31 genes × 7 variables.
