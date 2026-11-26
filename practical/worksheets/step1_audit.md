# Step 1: fetch and audit

## What I downloaded

| File | Size | What it holds | Why I need it (or not) |
|------|------|---------------|------------------------|
| gencode.v24.annotation.gtf.gz     | 38.7 MB | GENCODE file with the annotations of the genes                      | Conversion of the genes data to genes names to interpret the results |
| GSE115650_all.soft.txt            | 13.5 MB | Information about the 43 samples: condition, batch, sex, age        | Interpretation of the results and comparison of the groups           |
| GSE115650_biopsy_genecounts.csv.gz| 1.7 MB  | Compressed CSV file with one gene per row and one sample per column | Expression data we will analyze in next steps                        |

## Sample table

Number of samples: 43

Condition groups and counts: FSHD 34 | Control 9

Variables available per sample : 49 | gsm,title,geo_accession,status,submission_date,last_update_date,type,channel_count,source_name_ch1,organism_ch1,taxid_ch1,char:subject status,char:age,char:gender,char:path.score,char:inflam,char:active,char:stir,char:t1,char:fat.fraction,char:tissue,char:rnaseq.score,char:rnaseq.rank,char:sample_id,char:batch,growth_protocol_ch1,molecule_ch1,extract_protocol_ch1,data_processing,platform_id,contact_name,contact_email,contact_laboratory,contact_department,contact_institute,contact_address,contact_city,contact_state,contact_zip/postal_code,contact_country,instrument_model,library_selection,library_source,library_strategy,relation,supplementary_file_1,series_id,data_row_count,condition

Variables that have a value for one condition only : 10 | char:age,char:gender,char:path.score,char:inflam,char:active,char:stir,char:t1,char:fat.fraction,char:rnaseq.score,char:rnaseq.rank 

## Do the counts match the sample table?

File "counts_header_check.txt" shows that the counts file matches GEO (43 columns, 43 records), but one sample id is different: 01-0022-1 in our file and 01-0022 in GEO.
If we look at the first columns, sample 22 comes after 21 and before 23, 24, 25 and so on, so it looks like a naming error.

For the next steps we widecide to use the title, which matches every column (43/43).

## Cross-tabulations: condition against the other variables

[Age](../results/step1_audit/crosstab_age.csv) - Biological variable - Unbalanced because the age is missing for the Control group. 

[Batch](../results/step1_audit/crosstab_batch.csv) - Technical variable - Unbalanced because batch 1 and 2 contain only FSHD (5 each) when Controls are all in batch 3 and 4 and the size of each batch is not the same for all. It is a confounder

[Gender](../results/step1_audit/crosstab_gender.csv) - Biological variable - Unbalanced because the gender is missing for the Control group.

[Instrument_model](../results/step1_audit/crosstab_instrument_model.csv) - Technical variable - Balanced because the instrument was the same for each condition (Illumina HiSeq 2500).

[Tissue](../results/step1_audit/crosstab_tissue.csv) - Technical variable - Balanced because the tissue sampling was the same for each condition (muscle biopsy).

## What the paper says about the samples and GEO does not

Difference between :
- Cohort size and exclusion of the 2 person
- Details about the Controls (sex, age, familial link)
- Details about the muscle used for the Biopsy
- Distinction between FSHD1 and FSHD2
- Score CSS + Muscle strenght (dynamometer) + Fonctionnal test (30 feet walk)
- MRI - 10 FSHD with a mean of pathology score of 3.5/7 with 50% low to moderate fibrosis --> Detailled as 0/4 in STIR and 0/5 in T1 MRI score
- (No MRI for the control group but "blank cells" on GEO )

## FAIR score

Total: /30 (see [fair_score.md](../worksheets/fair_score.md)). The two weakest criteria and why:

## Decisions I made


