# Step 1: fetch and audit

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Keep it short and factual. Take each number from `results/step1_audit/`.

## What I downloaded

| File | Size | What it holds | Why I need it (or not) |
|Gene_names.cvs|3.38MB|Gene identity, names, and types|identify genes studied in the RNA-seq analysis with role in the disease|
|Gencode.v24.annotation.gtf.gz|36.9MB|gene ID names, genomic position, transcript isoform, exon position, gene biotype |to audit the biopsy gene counts|
|GSE115650_all.soft.txt|12.8MB |every FSHD sample characteristics, site, age , gender, clinical variables, batch, etc |to audit characteristics of patients and detect potential confounders |
|GSE11650_biopsy_genecounts.csv.gz|1.61MB|gene counts of the 43 samples| I need it for all tasks (PCA, comparisons, DUX4 score)|
|Counts_header_check.txt|1KB|holds columns in counts file, GSM records in GEO, columns matched, etc| I need it as quality check (samples match the metadata,file imported correctly, no samples omitted or duplicated)|
|Crosstabs.md|2KB|distribution of samples across categories of techincal and biological variables|I need it to check for differences between groups and potential confounders/covariates|
|missing_values.csv|1KB|missing values for each variable|I need it before analysis as quality check 
|sample_table.csv|89KB|holds experimental design for each patient and control sample| I need it to detect batch effect and other experimental design biases|
|series_record.txt|5KB|holds list of contributors, samples ID, emails, links to supplementary files, etc| most information are redundant|


## Sample table

Number of samples: 43 (34 FSHD and 9 Controls)

Condition groups and counts: 
condition,I,II,III,IV,All
Control,0,0,6,3,9
FSHD,5,5,16,8,34
All,5,5,22,11,43

Variables available per sample (list the columns of `sample_table.csv`): 
Variables that have a value for one condition only (see `missing_values.csv`): 
char:path.score,0,32,43
char:inflam,0,32,43
char:active,0,32,43
2 FSHD patients missing the above 3 variables.

## Do the counts match the sample table?

What `counts_header_check.txt` shows, and what you decided: 

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
