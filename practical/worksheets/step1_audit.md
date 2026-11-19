# Step 1: fetch and audit

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Keep it short and factual. Take each number from `results/step1_audit/`.

## What I downloaded

| File | Size | What it holds | Why I need it (or not) |
|------|------|---------------|------------------------|
|GSE115650_all.soft.txt |13.5 MB |GEO series record and sample metadata |To identify the biopsies and audit their metadata |
|GSE115650_biopsy_genecounts.csv.gz |1.7 MB |Raw gene counts for each biopsy |To check that the count columns match the sample records and provide inputs for analysis |
|gencode.v24.annotation.gtf.gz |38.7 MB |GENCODE v24 gene annotation |To map gene IDs to gene symbols and types |


## Sample table

Number of samples: 43
Condition groups and counts: 34 FSHD and 9 controls
Variables available per sample (list the columns of `sample_table.csv`): gsm, title, geo_accession, status, submission_date, last_update_date, type,
channel_count, source_name_ch1, organism_ch1, taxid_ch1, char:subject status,
char:age, char:gender, char:path.score, char:inflam, char:active, char:stir,
char:t1, char:fat.fraction, char:tissue, char:rnaseq.score, char:rnaseq.rank,
char:sample_id, char:batch, growth_protocol_ch1, molecule_ch1,
extract_protocol_ch1, data_processing, platform_id, contact_name,
contact_email, contact_laboratory, contact_department, contact_institute,
contact_address, contact_city, contact_state, contact_zip/postal_code,
contact_country, instrument_model, library_selection, library_source,
library_strategy, relation, supplementary_file_1, series_id, data_row_count,
condition
Variables that have a value for one condition only (see `missing_values.csv`): char:age, char:gender, char:path.score, char:inflam, char:active, char:stir,
char:t1, char:fat.fraction, char:rnaseq.score, char:rnaseq.rank


## Do the counts match the sample table?

What `counts_header_check.txt` shows, and what you decided:

All 43 sample columns in the count table match the GEO title field, while only 42 match char:sample_id. The single mismatch is record GSM3186680, which is labeled 01-0022-1 in title (matching the count column) but 01-0022 in char:sample_id. As no other sample uses 01-0022 and both title and char:sample_id are unique across the 43 records, we'll consider it as an inconsistent label for one biopsy, not a duplicate donor.
In Step 2, the merge key will be changed to title so the concerned biopsy will not be left unmatched and silently dropped. The choice was made so that the original data are not modified and the sample_id column will be kept in the future merged table for traceability and verification purposes.
Moreover, paper's Table S1 confirms that both keys are needed: it lists 01-0022 without the "-1" suffix, so title links GEO to the count matrix (43/43) while char:sample_id links GEO to the paper's tables (34/34 FSHD) (in results:`paper_table_s1_id_check.txt`)


## Cross-tabulations: condition against the other variables

Paste or link the tables from `results/step1_audit/crosstab_*.csv`. For each
variable say in one sentence whether it is balanced across conditions and
whether it is technical (batch, instrument, tissue of sampling) or biological
(sex, age). A technical variable that is unbalanced across conditions is a
confounder: name it as such.

Source: `results/step1_audit/crosstabs.md`

**Batch** (`char:batch`)

| condition | I | II | III | IV | All |
|:----------|--:|---:|----:|---:|----:|
| Control   | 0 |  0 |   6 |  3 |   9 |
| FSHD      | 5 |  5 |  16 |  8 |  34 |
| All       | 5 |  5 |  22 | 11 |  43 |

Batch is a partial confounder because Batches I and II contain only FSHD samples, meaning technical batch effects in those batches could look like disease effects (course 4).

**Instrument** (`instrument_model`)

| condition | Illumina HiSeq 2500 | All |
|:----------|--------------------:|----:|
| Control   |                   9 |   9 |
| FSHD      |                  34 |  34 |
| All       |                  43 |  43 |


**Tissue** (`char:tissue`)

| condition | muscle biopsy | All |
|:----------|--------------:|----:|
| Control   |             9 |   9 |
| FSHD      |            34 |  34 |
| All       |            43 |  43 |


Instrument (Illumina HiSeq 2500) and tissue (muscle biopsy) are technical and present respectively a single value in all samples, and are thus balanced across conditions.

**Sex** (`char:gender`)

| condition | female | male | missing | All |
|:----------|-------:|-----:|--------:|----:|
| Control   |      0 |    0 |       9 |   9 |
| FSHD      |     12 |   22 |       0 |  34 |
| All       |     12 |   22 |       9 |  43 |


**Age** (`char:age`, binned in years)

| condition | [0, 30) | [30, 45) | [45, 60) | [60, 120) | missing | All |
|:----------|--------:|---------:|---------:|----------:|--------:|----:|
| Control   |       0 |        0 |        0 |         0 |       9 |   9 |
| FSHD      |       2 |        5 |       14 |        13 |       0 |  34 |
| All       |       2 |        5 |       14 |        13 |       9 |  43 |


We cannot tell whether sex or age (biological) are balanced across conditions as both are missing for all 9 control samples.


## What the paper says about the samples and GEO does not

Read the paper's description of the cohort and of the controls. List the
variables the paper describes that you cannot find in the GEO record.

Checked with `scripts/check_paper_sample_ids.py` (output: `results/step1_audit/paper_table_s1_id_check.txt`).

| Detail described in the paper | Paper source | What GEO provides |
|---|---|---|
| Clinical Severity Score (CSS) and D4Z4 repeat-array size | Supplementary Table S1 | Neither value is provided per participant. |
| FSHD subtype and mosaic status | Supplementary Table S1 | The series description gives overall subtype counts. Sample records lack subtype and mosaic status. |
| Specific muscle biopsied and side | Supplementary Table S1 | `char:tissue` says `muscle biopsy` for all samples. The series description specifies quadriceps for controls. |
| Recruitment site | Methods | The paper reports 18 participants recruited at each site, Seattle and Rochester. GEO lacks a participant-level recruitment site. |
| Strength of the biopsied muscle | Methods, Clinical assessments | The paper describes measurement with a hand-held dynamometer. GEO lacks the measurements. |
| Control age and sex | Results, Study cohort | The paper reports six women and three men, with a mean age of 35 years (range 19–56). It gives these figures for the group as a whole. Age and sex are blank for each control sample in GEO. |
| Control relationship to people with FSHD and genetic confirmation | Results, Study cohort | GEO lacks these details. The paper describes unaffected spouses, parents, or siblings and reports genetic confirmation for at-risk controls. |

GEO records a sequencing batch (I to IV) for every sample. The paper does not
report these batches.
Additional discrepancies.** The article and GEO both report 34 FSHD1 and 2 FSHD2 participants among the 36 recruited, and report that 34% were women. Supplementary Table S1
lists three participants labelled FSHD2, one separately labelled mosaic,
and 13 women among 36 participants (36%). The GEO sample table contains
12 women among the 34 sequenced FSHD participants (35%).

## FAIR score

Total: /30 (see `fair_score.md`). The two weakest criteria and why:

## Decisions I made

-
