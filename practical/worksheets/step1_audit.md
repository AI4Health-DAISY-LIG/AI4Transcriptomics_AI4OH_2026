# Step 1: fetch and audit

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Keep it short and factual. Take each number from `results/step1_audit/`.

## What I downloaded

| File | Size | What it holds | Why I need it (or not) |
|------|------|---------------|------------------------|
| gencode.v24.annotation.gtf.gz     | 38.7 MB | GENCODE file with the annotations of the genes                      | Conversion of the genes data to genes names to interpret the results |
| GSE115650_all.soft.txt            | 13.5 MB | Information about the 43 samples: condition, batch, sex, age        | Interpretation of the results and comparison of the groups           |
| GSE115650_biopsy_genecounts.csv.gz| 1.7 MB  | Compressed CSV file with one gene per row and one sample per column | Expression data we will analyze in next steps                        |

## Sample table

Number of samples: 43
Condition groups and counts: FSHD 34 | Control 9
Variables available per sample (list the columns of `sample_table.csv`): 49 | gsm,title,geo_accession,status,submission_date,last_update_date,type,channel_count,source_name_ch1,organism_ch1,taxid_ch1,char:subject status,char:age,char:gender,char:path.score,char:inflam,char:active,char:stir,char:t1,char:fat.fraction,char:tissue,char:rnaseq.score,char:rnaseq.rank,char:sample_id,char:batch,growth_protocol_ch1,molecule_ch1,extract_protocol_ch1,data_processing,platform_id,contact_name,contact_email,contact_laboratory,contact_department,contact_institute,contact_address,contact_city,contact_state,contact_zip/postal_code,contact_country,instrument_model,library_selection,library_source,library_strategy,relation,supplementary_file_1,series_id,data_row_count,condition
Variables that have a value for one condition only (see `missing_values.csv`): 10 | char:age,char:gender,char:path.score,char:inflam,char:active,char:stir,char:t1,char:fat.fraction,char:rnaseq.score,char:rnaseq.rank 

## Do the counts match the sample table?

What `counts_header_check.txt` shows, and what you decided: "counts_header_check.txt" shows that the counts file matches GEO (43 columns, 43 records), but one sample id is different: 01-0022-1 in our file and 01-0022 in GEO.
If we look at the first columns, sample 22 comes after 21 and before 23, 24, 25 and so on, so it looks like a naming error.
For the next steps we will use the title, which matches every column (43/43).

## Cross-tabulations: condition against the other variables

Paste or link the tables from `results/step1_audit/crosstab_*.csv`. For each
variable say in one sentence whether it is balanced across conditions and
whether it is technical (batch, instrument, tissue of sampling) or biological
(sex, age). A technical variable that is unbalanced across conditions is a
confounder: name it as such.

## What the paper says about the samples and GEO does not

Read the paper's description of the cohort and of the controls. List the
variables the paper describes that you cannot find in the GEO record.

-

## FAIR score

Total: /30 (see `fair_score.md`). The two weakest criteria and why:

## Decisions I made

-
