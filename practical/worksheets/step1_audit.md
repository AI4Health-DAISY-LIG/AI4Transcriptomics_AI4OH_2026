# Step 1: fetch and audit

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Keep it short and factual. Take each number from `results/step1_audit/`.

## What I downloaded

| File | Size | What it holds | Why I need it (or not) |
|------|------|---------------|------------------------|
| | | | |
File: gencode.v24.annotation.gtf.gz;  GSE115650_all.soft.txt;    GSE115650_biopsy_genecounts.csv.gz
size:          36.9MB                      12.8 MB                   1.61MB
content:  GENCODE v24 gene annotation |  GEO series and sample metadata  |  gene-level counts for 43 biopsies
reason: Used to map Ensembl gene IDs to gene symbols |  Used to build the sample table | Used as the input for QC, PCA and differential expression.                                              

## Sample table

Number of samples: 43
Condition groups and counts: Control-9, FSHD-34, All-43
Variables available per sample (list the columns of `sample_table.csv`): 

gsm,title,geo_accession,status,submission_date,last_update_date,type,channel_count,source_name_ch1,organism_ch1,taxid_ch1,char:subject status,char:age,char:gender,char:path.score,char:inflam,char:active,char:stir,char:t1,char:fat.fraction,char:tissue,char:rnaseq.score,char:rnaseq.rank,char:sample_id,char:batch,growth_protocol_ch1,molecule_ch1,extract_protocol_ch1,data_processing,platform_id,contact_name,contact_email,contact_laboratory,contact_department,contact_institute,contact_address,contact_city,contact_state,contact_zip/postal_code,contact_country,instrument_model,library_selection,library_source,library_strategy,relation,supplementary_file_1,series_id,data_row_count,condition.

Variables that have a value for one condition only (see `missing_values.csv`): 

age, gender, path.score, inflam, active, stir, t1, fat.fraction, rnaseq.score, rnaseq.rank

## Do the counts match the sample table?
What `counts_header_check.txt` shows, and what you decided: 

Yes. The counts_header_check dislplayed the number of columns and GEO records in file identified by title and sample-id.
Using title to identify the count matrix, it contains 43 samples that match the GEO sccession numbers. Using sample-id, the count column 01-0022-1 matches the title of GSM3186680, but its char:sample_id is 01-0022.
Thus, I decided to use 'title' to verify 43 samples in next steps.

## Cross-tabulations: condition against the other variables

Paste or link the tables from `results/step1_audit/crosstab_*.csv`. For each
variable say in one sentence whether it is balanced across conditions and
whether it is technical (batch, instrument, tissue of sampling) or biological
(sex, age). A technical variable that is unbalanced across conditions is a
confounder: name it as such.

# crosstable_age

condition,"[0, 30)","[30, 45)","[45, 60)","[60, 120)",missing,All
Control,0,0,0,0,9,9
FSHD,2,5,14,13,0,34
All,2,5,14,13,9,43

age is a biological variable, not balanced in FSHD group. its balance condition is not known in control group.

# crosstable_batch

condition,I,II,III,IV,All
Control,0,0,6,3,9
FSHD,5,5,16,8,34
All,5,5,22,11,43

Batch is a technical variable and not balanced across samples. All samples in control group are from batch 3,4, most Fshd samples are from batch 3, 4 too.

# crosstable_gender

condition,female,male,missing,All
Control,0,0,9,9
FSHD,12,22,0,34
All,12,22,9,43

Gender is a biological variable. It's not balanced in FSHD samples, unknown in control.

# crosstable_instrument model

condition,Illumina HiSeq 2500,All
Control,9,9
FSHD,34,34
All,43,43

The instrument model is a technical variable, balanced. The instrument model has only one condition, Illumina HiSeq 2500.

# crosstabletable_tissue

condition,muscle biopsy,All
Control,9,9
FSHD,34,34
All,43,43

Tissue is an biological variable, balanced. All tissue is muscle biopsy.

Confounder: batch

## What the paper says about the samples and GEO does not

Read the paper's description of the cohort and of the controls. List the
variables the paper describes that you cannot find in the GEO record.

The paper reports the clinical severity score (CSS), genotype information and the specific muscle biopsied that are not included in samples and GEO.
-

## FAIR score

Total: 27/30 (see `fair_score.md`). The two weakest criteria and why:

F2: Sample metadata are rich: condition, donor, sex, age, tissue, batch are present for every sample. Some age and gender values in sample are missing.
R6: A reader can tell which deposited file feeds which figure of the paper. 
The files are named by data type, but the deposit does not map each figure to its input files.

## Decisions I made
Using 'title' to match the count columns.
?
-
