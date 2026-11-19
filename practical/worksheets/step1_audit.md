# Step 1: fetch and audit

_Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence._

Keep it short and factual. Take each number from `results/step1_audit/`.

## What I downloaded

| File                          | Size   | What it holds                                  | Why I need it (or not)                                                                                |
| ----------------------------- | ------ | ---------------------------------------------- | ----------------------------------------------------------------------------------------------------- |
| counts_header_check.txt       | 392 B  | Header validation results for the count matrix | Needed to confirm file header format is correct and is associated with the correct sample             |
| crosstab_age.csv              | 125 B  | Age group cross-tabulation statistics          | Useful for demographic analysis as a covariate, but not usable here since we lack control group data. |
| crosstab_batch.csv            | 78 B   | Sample counts per batch by condition           | Useful for QC and potentially for correcting batch effects if applicable during execution.            |
| crosstab_gender.csv           | 81 B   | Gender distribution cross-tabulation           | Like age, this data is not useful here since we lack control group information.                       |
| crosstab_instrument_model.csv | 67 B   | Machine used for mRNA sequencing analysis      | Every sample was processed on the same machine, so no need to use this as a covariate.                |
| crosstab_tissue.csv           | 61 B   | Tissue type cross-tabulation                   | Data is too coarse to be useful as a covariate; only one general tissue type is listed.               |
| crosstabs.md                  | 1.6 KB | Summary of all cross-tabulations               | Quick reference guide for all crosstab files.                                                         |
| missing_values.csv            | 332 B  | Missing value counts per sample metadata field | Needed to assess data completeness and identify gaps in the dataset.                                  |
| sample_table.csv              | 89 KB  | Sample metadata dataset                        | Used for mapping samples to geo records and as a source of covariate information.                     |
| series_record.txt             | 4.3 KB | Geo-series metadata                            | Helps locate supplementary files and provides general study info (PubMed ID, abstract, authors).      |

## Sample table

Number of samples: 43\
Condition groups and counts: There 2 condition (control : 9 and FSHD : 34)\
Variables available per sample (list the columns of `sample_table.csv`):

| column_name             |
| ----------------------- |
| gsm                     |
| title                   |
| geo_accession           |
| status                  |
| submission_date         |
| last_update_date        |
| type                    |
| channel_count           |
| source_name_ch1         |
| organism_ch1            |
| taxid_ch1               |
| char:subject status     |
| char:age                |
| char:gender             |
| char:path.score         |
| char:inflam             |
| char:active             |
| char:stir               |
| char:t1                 |
| char:fat.fraction       |
| char:tissue             |
| char:rnaseq.score       |
| char:rnaseq.rank        |
| char:sample_id          |
| char:batch              |
| growth_protocol_ch1     |
| molecule_ch1            |
| extract_protocol_ch1    |
| data_processing         |
| platform_id             |
| contact_name            |
| contact_email           |
| contact_laboratory      |
| contact_department      |
| contact_institute       |
| contact_address         |
| contact_city            |
| contact_state           |
| contact_zip/postal_code |
| contact_country         |
| instrument_model        |
| library_selection       |
| library_source          |
| library_strategy        |
| relation                |
| supplementary_file_1    |
| series_id               |
| data_row_count          |
| condition               |

Variables that have a value for one condition only (see `missing_values.csv`):

- char:age (Control: 0)
- char:gender (Control: 0)
- char:path.score (Control: 0)
- char:inflam (Control: 0)
- char:active (Control: 0)
- char:stir (Control: 0)
- char:t1 (Control: 0)
- char:fat.fraction (Control: 0)
- char:rnaseq.score (Control: 0)
- char:rnaseq.rank (Control: 0)

## Do the counts match the sample table?

What `counts_header_check.txt` shows, and what you decided:

The header does not match the sample ID for one sample: `01-0022` should be `01-0022-1`. Since the title is correctly set to `01-0022-1` for this sample, I chose to keep it and correct the sample ID to `01-0022-1` for the rest of the analysis.

## Cross-tabulations: condition against the other variables

Paste or link the tables from `results/step1_audit/crosstab_*.csv`. For each
variable say in one sentence whether it is balanced across conditions and
whether it is technical (batch, instrument, tissue of sampling) or biological
(sex, age). A technical variable that is unbalanced across conditions is a
confounder: name it as such.

| Variable                                                                 | Type       | Balance Status | Assessment                                                                                                                                                           |
| :----------------------------------------------------------------------- | :--------- | :------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [Age](../results/step1_audit/crosstab_age.csv)                           | Biological | Unknown        | The balance status is unknown because there is no data available for the control group.                                                                              |
| [Batch](../results/step1_audit/crosstab_batch.csv)                       | Technical  | Unbalanced     | This variable acts as a confounder because controls are absent in the first two batches and the last two batches show a 2.6x higher sample count for the FSHD group. |
| [Gender](../results/step1_audit/crosstab_gender.csv)                     | Biological | Unknown        | The balance status is unknown due to missing data for the control group.                                                                                             |
| [Instrument Model](../results/step1_audit/crosstab_instrument_model.csv) | Technical  | Balanced       | This variable is balanced because every sample was sequenced on the same machine.                                                                                    |
| [Tissue](../results/step1_audit/crosstab_tissue.csv)                     | Technical  | Balanced       | This variable is considered balanced as all samples are muscle biopsies, although they originate from different specific muscles.                                    |

## What the paper says about the samples and GEO does not

Read the paper's description of the cohort and of the controls. List the
variables the paper describes that you cannot find in the GEO record.

- There are two types of FSHD (34 FSHD1, 2 FSHD2)
- The Gender of the control is composed of 6 females, 3 males but no mapping to sample is given. Also the male/female ratio differs a lot between the groups : 34% females in the FSHD group vs 66% females in the control.
- The detailed biopsy tissue type for the samples of the FSHD (additional files)
- Regarding the age the mean age between the FSHD and the control respectively 54 vs 35 years old differ.
- The control are related (spouses, parents or siblings of FSHD individuals).

## FAIR score

Total: 16/30 (see `fair_score.md`). The two weakest criteria and why:

The weakest criteria for me is F2:

> **F2:** Sample metadata are rich; condition, donor, sex, age, tissue, batch are present for every sample.

Because a lot of those variables are important for the analysis to add as covariates or confounders.

The second weakest criterion is R4.

> **R4:** The sample table lets you rebuild the groups compared in the paper without guessing.

All the analyses were based on these four groups, but we cannot recompute them because they did not provide filtering rules and theasholds used.

## Decisions I made

-
