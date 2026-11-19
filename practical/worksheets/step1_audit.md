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

All 43 sample columns in the count table match the GEO title field, while only 42 match char:sample_id. The single mismatch is record GSM3186680, which is labeled 01-0022-1 in title (matching the count column) but 01-0022 in char:sample_id. As no other sample uses 01-0022, we'll consider it as an inconsistent label for one biopsy, not a duplicate donor.
In Step 2, the merge key will be changed to title so the concerned biopsy will no be left unmatched and silently dropped.


## Cross-tabulations: condition against the other variables

Paste or link the tables from `results/step1_audit/crosstab_*.csv`. For each
variable say in one sentence whether it is balanced across conditions and
whether it is technical (batch, instrument, tissue of sampling) or biological
(sex, age). A technical variable that is unbalanced across conditions is a
confounder: name it as such.

FSHD and control biopsies are unevenly distributed across batches, so technical batch effects could look like disease effects in PCA (Course 4)

## What the paper says about the samples and GEO does not

Read the paper's description of the cohort and of the controls. List the
variables the paper describes that you cannot find in the GEO record.

-

## FAIR score

Total: /30 (see `fair_score.md`). The two weakest criteria and why:

## Decisions I made

-
