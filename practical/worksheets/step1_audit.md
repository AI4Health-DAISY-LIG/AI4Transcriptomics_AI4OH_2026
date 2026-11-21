# Step 1: fetch and audit

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Keep it short and factual. Take each number from `results/step1_audit/`.

## What I downloaded

| File                             | Size   | What it holds                  | Why I need it (or not)                   |
|----------------------------------|--------|--------------------------------|------------------------------------------|
| GSE115650_all.soft.txt           |13.5 MB |GEO metadata and sample records | To obtain sample information and metadata|
|----------------------------------|--------|--------------------------------|------------------------------------------|
|GSE115650_biopsy_genecounts.csv.gz|1.7 MB  |Raw gene counts                 |Main input for RNA-seq analysis           |
|----------------------------------|--------|--------------------------------|------------------------------------------|
|gencode.v24.annotation.gtf.gz     |38.7 MB |Gene annotation                 |To map gene identifiers to gene names     |
|----------------------------------|--------|--------------------------------|------------------------------------------|
|gene_names.csv                    |3.7 MB  |Gene names extracted from annotation|To identify genes in the count matrix |
|-----------------------------------------------------------------------------------------------------------------------|


## Sample table

Number of samples: 43
Condition groups and counts: 34 FSHD and 9 Control.
Variables available per sample (list the columns of `sample_table.csv`): 
gsm, title, geo_accession, status, submission_date, last_update_date, type, channel_count, source_name_ch1, organism_ch1, taxid_ch1, char:subject status, char:age, char:gender, char:path.score, char:inflam, char:active, char:stir, char:t1, char:fat.fraction, char:tissue, char:rnaseq.score, char:rnaseq.rank, char:sample_id, char:batch, growth_protocol_ch1, molecule_ch1, extract_protocol_ch1, data_processing, platform_id, contact_name, contact_email, contact_laboratory, contact_department, contact_institute, contact_address, contact_city, contact_state, contact_zip/postal_code, contact_country, instrument_model, library_selection, library_source, library_strategy, relation, supplementary_file_1, series_id, data_row_count, condition.

Variables that have a value for one condition only (see `missing_values.csv`): 

Age and gender are missing for all 9 Control samples. The variables path.score, inflam, and active are missing for 2 FSHD samples. The other listed variables are available for all samples.
## Do the counts match the sample table?

What `counts_header_check.txt` shows, and what you decided: 
The count matrix contains 43 sample columns, matching the 43 GEO samples. All 43 columns matched GEO records using the sample title. However, matching by char:sample_id matched only 42 samples: the count file contains 01-0022-1, while GEO contains 01-0022.

I kept the samples because the title-based matching was complete and documented the sample ID discrepancy for further investigation.

## Cross-tabulations: condition against the other variables

Paste or link the tables from `results/step1_audit/crosstab_*.csv`. For each
variable say in one sentence whether it is balanced across conditions and
whether it is technical (batch, instrument, tissue of sampling) or biological
(sex, age). A technical variable that is unbalanced across conditions is a
confounder: name it as such.

Age: Age is available only for FSHD samples and is missing for all 9 controls, so age cannot be compared between conditions.

Batch: Batch is unbalanced across conditions. Controls are present only in batches III and IV, while FSHD samples are present in all four batches. Batch is a technical variable and a potential confounder.

Gender: Gender information is available only for FSHD samples and is missing for all controls; therefore, comparison between conditions is not possible.

Instrument: All 43 samples were sequenced using Illumina HiSeq 2500. This technical variable is balanced and has no variation.

Tissue: All 43 samples are muscle biopsies. This variable is balanced and has no variation.

## What the paper says about the samples and GEO does not

Read the paper's description of the cohort and of the controls. List the
variables the paper describes that you cannot find in the GEO record.

The GEO metadata states that muscle biopsy samples were obtained from 34 of 36 FSHD subjects and from 9 control quadriceps muscle biopsies from unaffected individuals. Age and gender information are missing for all Control samples in the GEO metadata.

Further comparison with the paper and supplementary tables is needed to identify additional cohort or biopsy characteristics described in the publication but absent from GEO.

-

## FAIR score

Total: /30 (see `fair_score.md`). The two weakest criteria and why:
Not completed yet. The FAIR score will be assessed at the end of the project after reviewing the dataset, GEO record, paper and reproducibility issues.

## Decisions I made
Kept all 43 samples because sample titles matched all GEO records.

Documented the discrepancy between sample IDs 01-0022-1 and 01-0022 rather than changing it without evidence.

Identified batch as a potential confounder because its distribution differs between FSHD and Control samples.

Documented missing demographic information for Control samples.

Deferred the FAIR score until the end of the project, when all reproducibility issues can be assessed.

-
