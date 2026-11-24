# Step 3: differential expression with PyDESeq2

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Step 3 has two tracks:

- **Track A, the reproduction.** It reruns the authors' DE in R with their own DESeq2 object from GEO (`r/de_reproduction.R`, `results/step3_de/track_a/`).
- **Track B, the independent analysis.** This is the PyDESeq2 analysis that this worksheet asks for (`scripts/step3_differential_expression.py`, `results/step3_de/track_b/`).

The design and sample rules of Track B were fixed before the fits were run. Two fits were added after the first run, before any gene-level result was examined: Rochester-only and the G3 check.

## Design

Input matrix (which file, which filter): `data/processed/counts_filtered.csv`, the GEO raw counts (27,506 genes × 43 biopsies). There is no extra filter: the authors had already removed genes with a total count ≤ 20 before deposit (`results/step2_qc/gene_filter.txt`). Size factors are estimated by PyDESeq2 for each fit.

Design formula: `~ batch + sex + condition` (primary, FSHD vs control). The groups are fitted in one joint model, `~ batch + sex + group`, with controls as the reference.

Covariates included and why, or why not:
- **Batch (4 runs).** Library size and the instrument follow batch (`results/step2_qc/`). Batch explains 36% of PC1 and 33% of PC2 (`pca_associations.csv`).
- **Sex, inferred from expression** (`results/step2_qc/sex_check.csv`). GEO gives no sex for the controls. Sex is imbalanced: controls are 6 women and 3 men, while G4 is 0 women and 5 men. Sex genes sit on PC2. The recorded GEO labels stay unchanged in the sample table, and X/Y genes are not interpreted.
- **Site, sensitivity only** (`~ batch + sex + site + condition`). All 9 controls are from Rochester (banked biopsies), so the site effect is learned from FSHD biopsies alone. `32` = Seattle (authors' follow-up workbook); `01` = Rochester, inferred from that roster and the paper's 18 participants per site (`scripts/config.py`, `SITE_BY_PREFIX`).
- **Library size and mapping rate are not covariates.** Size factors and batch absorb them. The three libraries below 90% mapping are excluded instead.

Covariates you wanted and could not include, and why:
- **Age:** missing for every control in GEO.
- **Biopsied muscle:** all controls are quadriceps, so it cannot be separated from condition.
- **Relatedness of controls to patients:** not recorded in GEO.
- **Tissue composition (fat, inflammation, blood):** deliberately not adjusted, because it is part of the disease.

Number of samples per group: 39 biopsies in the primary set, after excluding 32-0008, 01-0037 and 32-0009 (mapping < 90%) and 32-0004 (recorded sex disagrees with expression), as set in `scripts/config.py`, `TRACK_B_EXCLUDE`. Control 9; FSHD 30, made of G1 9, G2 5, G3 11, G4 5 (`results/step3_de/track_b/preflight.txt`).

## Results

Track B, padj < 0.05 (`results/step3_de/track_b/summary.txt`). "Genes tested" means genes with a finite p-value; the genes with an adjusted p-value after independent filtering are given in brackets.

| Comparison | Genes tested | padj < 0.05 | of which up | of which down |
|------------|--------------|-------------|-------------|---------------|
| FSHD vs Control | 27,471 (20,021) | 2,577 | 2,078 | 499 |
| Group 1 vs Control | 27,407 (19,968) | 1,002 | 866 | 136 |
| Group 2 vs Control | 27,407 (19,968) | 518 | 463 | 55 |
| Group 3 vs Control | 27,407 (20,498) | 639 | 564 | 75 |
| Group 4 vs Control | 27,407 (22,090) | 5,409 | 3,212 | 2,197 |

Figures: MA plot, volcano plot and p-value histogram, in `results/step3_de/track_b/primary/` and in the two ladder fits.

Track A, the reproduction of the paper's tables, is exact (`results/step3_de/track_a/published_comparison.txt`):
- S4 (the two G1 outliers) 5,680, S5 (G1) 393, S6 (G2) 303, S8 (G3) 816, S10 (G4) 2,826 genes at padj < 0.05 and |log2FC| ≥ 2;
- the same genes as the published tables (0 missing, 0 extra), with shared log2FC within 0.0004;
- the notebook-only contrasts also match: G1 without its two outliers 40, G2 active 540.

Top 10 genes by adjusted p-value, primary FSHD vs Control (`primary/de_FSHD_vs_Control.csv.gz`):

| Gene | Chromosome | log2FC | padj |
|---|---|---|---|
| COMP | chr19 | 7.52 | 1.5e-09 |
| TNMD | chrX | 5.55 | 1.1e-07 |
| PRG4 | chr1 | 4.38 | 1.7e-05 |
| CDKN1A | chr6 | 2.78 | 1.8e-05 |
| MS4A14 | chr11 | 1.92 | 2.1e-05 |
| EPB41L3 | chr18 | 1.48 | 2.1e-05 |
| CTGF | chr6 | 2.14 | 2.6e-05 |
| FOSB | chr19 | 2.86 | 3.7e-05 |
| TNFAIP6 | chr2 | 1.66 | 3.7e-05 |
| FOS | chr14 | 3.16 | 3.7e-05 |

These are extracellular matrix and fibrosis genes (COMP, TNMD, PRG4, CTGF), stress and inflammation genes (CDKN1A, TNFAIP6) and immediate early genes (FOS, FOSB). They match the top of the paper's S5 (G1): PRG4 1st, COMP 2nd, TNMD 4th, CDKN1A 7th. A test pooling all FSHD biopsies favours changes shared by most of them, and those are not the DUX4 targets, which are switched on in only some biopsies. FOS and FOSB respond within minutes to tissue handling, so part of their signal may come from how the banked controls and the study biopsies were collected. TNMD is on chrX and is not interpreted.

## Checks

- **Describe the p-value histogram: a flat part with a spike near 0, or another shape?**
  A large spike in the first bin (about 4,200 genes with p < 0.02), then a flat part from about 0.3 to 1 (about 350–450 genes per bin). Between them, a gentle decline from 0.02 to about 0.3 shows many small effects (`primary/pvalue_histogram_FSHD_vs_Control.png`). There is no hump near 1 and no U shape, so the model does not look misspecified. The unadjusted model has the same shape with a higher spike, about 6,500 genes (`ladder_none/`).

- **Run the analysis with and without batch in the design. Report how many genes are significant in both, and in one only.**
  On the same 39 biopsies, `~ condition` gives 5,509 genes and `~ batch + condition` gives 4,315. That is 4,166 genes in both, 1,343 only without batch and 149 only with batch. Adding batch removes about a quarter of the calls, but the effects of the genes kept are stable. Of the 2,577 primary genes, 2,560 keep their direction and at least half their log2FC without any adjustment, and 2,565 without sex (`sensitivity_comparison.csv`, `robustness.csv.gz`). Adding sex brings the count down to 2,577. The full ladder of models and the other sensitivity runs (site, all 43, 32-0004 restored, 32-0016 omitted, Rochester only) are in `sensitivity_comparison.csv`. Effects hold in 2,416–2,577 of the 2,577 genes; site mostly widens the standard errors (median SE ratio 1.15).

- **Group 1 biopsies show almost no expression of the biomarkers. Count the genes that differ between group 1 and the controls, and list the explanations you can think of other than the disease.**
  Track B, `~ batch + sex + group`: 1,002 genes (866 up, 136 down); 1,945 without batch and sex (`groups_none`). Track A, the paper's design: 393 at |log2FC| ≥ 2, but only 40 once the two G1 outliers, 01-0037 and 32-0016, are removed (`track_a/summary.txt`). Explanations other than the disease:
  - **Two outlier libraries.** 01-0037 maps at 76.5%; 32-0016 is a fat-rich biopsy. They carry most of the paper's G1 signal.
  - **How the biopsies were collected and stored.** The controls are banked Rochester biopsies, so needle, handling and ischaemia time may differ; the immediate early genes are a sign of this. Site and banking are confounded with condition.
  - **Muscle and age.** All controls are quadriceps and younger (mean 35 vs 53.6 years), while G1 biopsies come from several muscles.
  - **Sex and batch imbalance**, adjusted only in Track B.
  - **Relatedness:** controls are spouses, parents or siblings of patients.
  - **Tissue composition** (fat, fibrous tissue, immune cells): it changes bulk expression whatever its cause.
  - **A different DE design.** The paper fits each group separately with controls, so for G1 the dispersions come from 19 biopsies; Track B fits all groups together.

  G1 can also reflect the disease without detectable DUX4, because DUX4 is sporadic. The shared fibrosis signal (PRG4, COMP, TNMD) points that way, but these data cannot separate it from the explanations above.

## Decisions I made

- **Two tracks.** Track A reproduces the paper exactly in R, because the paper used DESeq2 in R; that also separates design differences from software differences. Track B is the PyDESeq2 analysis with the best adjustment the data allow.
- **Track B sample rules, fixed after QC and before any DE result** (`scripts/config.py`):
  - exclude the three libraries mapped below 90%;
  - exclude 32-0004, a provenance flag, not a proven swap;
  - keep the biological extreme 32-0016.

  Each exclusion is undone in a sensitivity run: `all43`, `restore_32_0004` and `omit_32_0016`.
- **Primary design `~ batch + sex + condition` with inferred sex.** Site is a sensitivity only, because all controls come from one site.
- **Robustness rule:** same sign and at least 50% of the primary |log2FC| in each adjusted sensitivity. Significance is required in the primary only, because smaller fits lose power. Result: 2,408 robust genes over the five planned sensitivities (2,319 on chr1–22), and 2,385 with Rochester (2,297 on chr1–22).
- **Rochester-only fit added**, to remove the site difference between cases and controls. Its effects hold for 2,436 of 2,577 genes. They do not depend on Cook's replacement: 2,440 without it (`rochester_norefit`).
- **G3 check added.** With 32-0008 and 32-0009 put back, G3 rises from 639 to 1,669 genes (from 235 to 415 at |log2FC| ≥ 2; published 816), while G1 and G2 drop. These two low-mapping libraries account for part of the gap with the published G3 count.
- **DUX4-score models** (`score*`, FSHD only), for the bonus. They are described in `worksheets/bonus_translator.md`.
