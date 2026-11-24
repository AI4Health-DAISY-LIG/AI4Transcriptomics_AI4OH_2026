# Bonus: Biolink and Translator

## Gene selection rule

Fixed on 24 Nov 2026, before the selection script was run. Agreed between the student, the lead assistant, Fable and Assistant 2. Applied by `scripts/bonus_select_genes.py`, which writes to `results/bonus/translator/`. The reference lists are in `bonus_lists/`, with their sources, versions and checksums in `bonus_lists/sources.tsv`.

This is exploratory, for the bonus. No comparison of biopsies taken at a single time can show what causes DUX4 expression.

### Two views, both from `results/step3_de/track_b/`

**View A: FSHD differential expression that tracks the DUX4-target score. This view supplies the required gene query.**

- Primary FSHD vs control, `~ batch + sex + condition`, n = 39 (`primary/de_FSHD_vs_Control.csv.gz`).
- Eligible genes meet all of the following:
  - padj < 0.05 and |log2FC| ≥ 1;
  - robust: same sign and ≥ 50% of |log2FC| in the five planned sensitivities and in Rochester-only (`robustness.csv.gz`, `robust_with_rochester`);
  - the same sign and padj < 0.05 in the score model, so the gene follows the DUX4-target score within FSHD.
- Rochester is in the gate because its effects do not depend on Cook's replacement: 2,440 vs 2,436 primary effects hold with and without it (`sensitivity_comparison.csv`).

**View B: regulatory genes associated with the DUX4-target score within FSHD. This view supplies the optional second query.**

- Fitted on FSHD biopsies only, `~ batch + sex + score`, n = 30 (`score/de_score_slope.csv.gz`).
- Eligible genes meet all of the following:
  - padj < 0.05;
  - |log2FC across the G1→G4 score span| ≥ 1, where the span is frozen at 14.98423 (from the G1 median 0.35207 to the G4 median 15.33630);
  - robust (`score_robustness.csv.gz`, `robust_b`): same sign and ≥ 50% of the slope without sex, with all 34 FSHD biopsies, without 32-0006 and without 01-0025; and the same sign as G4 vs G1.
- B is circular by construction: the score comes from DUX4 targets, and DUX4 induces MBD3L proteins that relieve NuRD silencing.

### Exclusions, both views

- Not on chr1–22.
- Not protein-coding (GENCODE v24 `gene_type`).
- One of the 4 group biomarkers: LEUTX, KHDC1L, PRAMEF2, TRIM43.
- **Immediate early genes:** Tullai et al. 2007, NHLBI ESBL copy, 46 symbols, a screen derived from growth-factor-stimulated T98G cells. Matched in upper case, with CCN1 = CYR61 and CCN2 = CTGF (the GENCODE v24 names).
- **Composition genes:**
  - haemoglobins: HBA1, HBA2, HBB, HBD, HBE1, HBG1, HBG2, HBM, HBQ1, HBZ;
  - every `HLA-` symbol;
  - adipocyte genes: ADIPOQ, LEP, PLIN1, PLIN4, CIDEA, CIDEC, FABP4;
  - neutrophil genes: S100A8, S100A9, S100A12, FCGR3B, CXCR1, CXCR2, CSF3R, ELANE, MPO, CAMP, LTF;
  - immunoglobulins are already excluded as not protein-coding.
- **View B pick only:** the 213 robust DUX4-induced genes of Yao et al. 2014, Supplementary Table S1, matched on the Ensembl ID without version. For the IDs absent from GENCODE v24 (retired or renamed since 2014), the match falls back on the symbol, so that DUXA and five PRAMEF genes are not missed. The match status of each gene is in `results/bonus/translator/yao2014_matching.csv`. These genes are flagged in the tables, not hidden.

### Functional filters

1. **Transcription factors:** Lambert et al. 2018, TF names v1.01, 1,639 symbols.
2. **Chromatin regulators:** EpiFactors v2.1, all entries except those whose only functions are RNA modification (with or without scaffold; 84 entries) or protein modification of non-histone proteins (4 entries). That leaves 713 of 801.
3. **Muscle regeneration:** GO:0043403, human annotations, descendants through `is_a` and `part_of` (QuickGO, 45 annotations).
4. **Known DUX4 regulators, reported in both views and never selected:** SMCHD1, DNMT3B, LRIF1, CHD4, HDAC2, MTA2, MBD2, RBBP4, CHAF1A, CHAF1B, MBD3L2, MAPK14, MAPK11.

Symbols that match no gene in our annotation are counted and reported.

### Ranking, influence check, choice

- **Ranking** within each pool: padj, then |Wald statistic| (higher first), then gene ID.
- **Influence, view A:** refit the primary model without the biopsy that has the largest Cook's distance for the gene. If the gene loses its sign or more than 50% of its |log2FC|, skip it, log the skip, and take the next gene.
- **Influence, view B:** covered by its robustness checks.
- **View A gene:** the first gene that passes, taken from the regeneration pool (3), then the transcription-factor and chromatin pool (1 ∪ 2), then all eligible genes.
- **View B gene:** the first gene from the chromatin pool (2), then the transcription-factor pool (1). If none passes, report "none" and do not relax the rule.
- **Reported:** the top 10 of every pool, with the reason for each skip, and the known DUX4 regulators in both views.

### Exploration: does the top move between runs?

This section is for looking, not choosing; it does not change the two designated genes. For each run, the script takes that run's own top 10: the gene filters above, plus padj < 0.05 and |log2FC| ≥ 1 measured in that run (for view B, across the score span), ranked by that run's padj.
- View A runs: primary, ladder_batch, site, all43, restore_32_0004, omit_32_0016, rochester and ladder_none.
- View B runs: score and its four checks.

The union of these top-10 lists, with each gene's rank in every run and whether it passes the rule, is in `top10_by_run_A.csv` and `top10_by_run_B.csv`. Genes that rise into the top after a correction are inspected, and may be queried in Translator as exploration. The report lists every gene queried and its result count, and states how the two highlighted genes were chosen.

### Translator queries (README, bonus)

- **Disease query:** "What drugs may impact conditions related to ...", with facioscapulohumeral muscular dystrophy.
- **Gene query:** the view A gene. If it is higher in FSHD, use "What chemicals may decrease the activity of a gene?"; if it is lower, use "increase".
- **Optional view B query:** the direction comes from the sign of the score slope, stated explicitly. It is worded as "what modulates a gene associated with DUX4-high muscle".
- Neither gene is changed after seeing Translator output.

## Queries and results

To fill in: date, link ("Share Result Set"), number of results after "Sync New Results", first results, the evidence behind one relationship, and what it tells us.
