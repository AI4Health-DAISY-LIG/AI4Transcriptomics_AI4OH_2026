# Step 1: fetch and audit

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Keep it short and factual. Take each number from `results/step1_audit/`.

## What I downloaded

| File | Size | What it holds | Why I need it (or not) |
|------|------|---------------|------------------------|
| `GSE115650_all.soft.txt` (GEO series record) | 13.5 MB | Series and sample metadata, 43 GSM records | Builds `sample_table.csv`: condition, batch, sex, age, pathology, MRI, authors' score and group |
| `GSE115650_biopsy_genecounts.csv.gz` | 1.7 MB | Raw gene counts, GENCODE v24 gene ids x 43 samples | Input of steps 2 and 3: DESeq2 needs raw integer counts |
| `gencode.v24.annotation.gtf.gz` | 38.7 MB | GENCODE v24 annotation, the release named in the paper | The counts file has gene ids only; gives the symbols (60,554 genes in `gene_names.csv`) |
| FPKM (8.2 MB), DESeq2-normalised counts (6.7 MB), `biopsy.dds.rda` (14.9 MB) | not downloaded | Normalised values and an R object | Not valid input for DESeq2, and the `.rda` needs R |

## Sample table

Number of samples: 43
Condition groups and counts: FSHD 34, Control 9
Variables available per sample (list the columns of `sample_table.csv`): GEO fields (gsm, title, source_name_ch1, organism, molecule, extraction and library protocols, data_processing, instrument_model, library_selection/source/strategy, relation to BioSample/SRA) and the characteristics `subject status`, `age`, `gender`, `path.score`, `inflam`, `active`, `stir`, `t1`, `fat.fraction`, `tissue`, `rnaseq.score`, `rnaseq.rank`, `sample_id`, `batch`; plus `condition` added by the script.
Variables that have a value for one condition only (see `missing_values.csv`): `age`, `gender`, `path.score`, `inflam`, `active`, `stir`, `t1`, `fat.fraction`, `rnaseq.score`, `rnaseq.rank` are filled for FSHD only, empty for all 9 controls. `path.score`, `inflam` and `active` are also missing for 2 FSHD samples (32-0011, 32-0016), consistent with the paper ("in three subjects, the histopathologic samples were inadequate").

## Do the counts match the sample table?

What `counts_header_check.txt` shows, and what you decided: 43 count columns and 43 GSM records. Matched by `title`: 43/43. Matched by `sample_id`: 42/43. The unmatched column `01-0022-1` is the GEO **title** of the record whose `sample_id` is `01-0022` (batch III, authors' group 3): same GSM record. Decision: keep this sample and rename the column `01-0022` in step 2. The meaning of the suffix "-1" (re-run, second library?) is not documented; asked in issue #22.

## Cross-tabulations: condition against the other variables

| Variable | Type | Balanced across conditions? |
|----------|------|-----------------------------|
| `batch` (`crosstab_batch.csv`) | technical | **No.** Batches I (5 FSHD) and II (5 FSHD) contain no control; controls are only in III (6) and IV (3). Batch is partly confounded with the condition: **confounder**. Batch also tracks the authors' groups (batch I: 4 of 5 biopsies are group 4). |
| `instrument_model` (`crosstab_instrument_model.csv`) | technical | Constant (HiSeq 2500 for all 43): not a confounder, but carries no information. |
| `tissue` (`crosstab_tissue.csv`) | technical | Constant ("muscle biopsy"): hides the muscle type, which differs between groups (see next section). |
| `gender` (`crosstab_gender.csv`) | biological | Cannot be assessed from GEO: missing for all controls. The paper gives 6 F / 3 M for controls (67 % F) against 12 F / 22 M for sequenced FSHD (35 % F): unbalanced. |
| `age` (`crosstab_age.csv`) | biological | Cannot be assessed from GEO: missing for all controls. Paper: controls mean 35 years (19-56); sequenced FSHD mean 53.6 (20-75) from GEO: unbalanced. |
| recruitment site (prefix of the sample id, not a GEO variable) | technical | **No.** `01-` = Rochester, `32-` = Seattle (18 subjects per site in the paper). All 9 controls are `01-` (Rochester biobank); FSHD are 16 `01-` and 18 `32-`. Site is a **confounder** that GEO does not record. |

## What the paper says about the samples and GEO does not

Read the paper's description of the cohort and of the controls. List the
variables the paper describes that you cannot find in the GEO record.

- Muscle biopsied: all controls are quadriceps; FSHD are quadriceps (6), gastrocnemius (14), tibialis anterior (13), hamstring (1) (Results; Table S1). GEO `tissue` is "muscle biopsy" for all. The authors write that some differences "might reflect differences in the type of muscle".
- Sex and age of the controls (6 F / 3 M, mean 35 years, 19-56).
- Recruitment site (University of Washington, Seattle, or University of Rochester) and biopsy needle (modified Bergstrom in Rochester, regular Bergstrom in Seattle).
- Controls are unaffected spouses, parents or siblings of FSHD individuals: some controls may be related to some patients, and the relation is not given.
- Genotype (FSHD1 or FSHD2), short 4q allele size and clinical severity score (CSS) (Table S1).
- The paper does not give the accession: it writes "GSExxxxx" in the Methods.

## FAIR score

Total: 21/30 (see `fair_score.md`). The two weakest criteria and why: R1 (0, no licence or usage terms anywhere) and R6 (0, nothing tells which deposited file feeds which figure). F2 and R2 are also weak for the reproduction: no sex or age for the controls, and no DESeq2 version, filter, design or thresholds.

## Decisions I made

- Keep sample `01-0022-1` and rename it `01-0022` (same GEO record, matched by title).
- Treat batch as a technical confounder: run step 3 with and without batch in the design.
- Site and muscle type cannot be adjusted for with the deposited variables (site is only in the id prefix and is fully confounded with being a control; muscle type is absent). I report them as limits.
- Sex and age cannot be used as covariates: they are missing for all controls.