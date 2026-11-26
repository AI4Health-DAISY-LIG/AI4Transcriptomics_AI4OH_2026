# Step 1: fetch and audit

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Keep it short and factual. Take each number from `results/step1_audit/`.

## What I downloaded

| File | Size | What it holds | Why I need it (or not) |
|Gene_names.cvs|3.38MB|Gene identity, names, and types|identify genes studied in the RNA-seq analysis with role in the disease; to audit the biopsy gene counts; quality check before running all tasks|
|Gencode.v24.annotation.gtf.gz|36.9MB|gene ID names, genomic position, transcript isoform, exon position, gene biotype |to audit the annotation and the biopsy gene counts; quality check before running all tasks|
|GSE115650_all.soft.txt|12.8MB |every FSHD sample characteristics, site, age , gender, clinical variables, batch, etc |to audit characteristics of patients and detect potential confounders |
|GSE11650_biopsy_genecounts.csv.gz|1.61MB|gene counts of the 43 samples| I need it for all tasks (PCA, comparisons, DUX4 score)|
|Counts_header_check.txt|1KB|holds columns in counts file, GSM records in GEO, columns matched, etc| I need it as quality check (samples match the metadata,file imported correctly, no samples omitted or duplicated)|
|Crosstabs.md|2KB|distribution of samples across categories of technical and biological variables|I need it to check for differences between groups and potential confounders/covariates|
|missing_values.csv|1KB|missing values for each variable|I need it before analysis as quality check 
|sample_table.csv|89KB|holds experimental design for each patient and control sample with variables such as sociodemographics and clinical| I need it to make comparisons between groups, detect batch effect and other experimental design biases, and answer the tasks|
|series_record.txt|5KB|holds list of contributors, samples ID, emails, links to supplementary files, etc| I don't think it's really helpful because most information are redundant|


## Sample table

Number of samples: 43 (34 FSHD and 9 Controls)

Condition groups and counts: 
condition,I,II,III,IV,All
Control,0,0,6,3,9
FSHD,5,5,16,8,34
All,5,5,22,11,43

Variables available per sample (list the columns of `sample_table.csv`): (I'm having difficulties writing a code to extract the list so I opened the file in excel and copied the first title row => should write code to extract list)
gsm	title;	geo_accession;	status	submission_date;	last_update_date;	type;	channel_count;	source_name_ch1;	organism_ch1;	taxid_ch1;	char:subject status;	char:age;	char:gender;	char:path.score;	char:inflam;	char:active;	char:stir;	char:t1;	char:fat.fraction;	char:tissue;	char:rnaseq.score;	char:rnaseq.rank;	char:sample_id;	char:batch;	growth_protocol_ch1;	molecule_ch1;	extract_protocol_ch1;	data_processing;	platform_id;	contact_name;	contact_email;	contact_laboratory;	contact_department;	contact_institute;	contact_address;	contact_city;	contact_state;	contact_zip/postal_code;	contact_country;	instrument_model;	library_selection;	library_source;	library_strategy;	relation;	supplementary_file_1;	series_id;	data_row_count;	condition
Controls miss the following variable: char:age;	char:gender;	char:path.score;	char:inflam;	char:active;	char:stir;	char:t1;	char:fat.fraction;	char:rnaseq.score;	char:rnaseq.rank;


Variables that have a value for one condition only (see `missing_values.csv`): 
Variable,Control,FSHD,n_samples_total
char:subject status,9,34,43
char:age,0,34,43
char:gender,0,34,43
char:path.score,0,32,43
char:inflam,0,32,43
char:active,0,32,43
char:stir,0,34,43
char:t1,0,34,43
char:fat.fraction,0,34,43
char:tissue,9,34,43
char:rnaseq.score,0,34,43
char:rnaseq.rank,0,34,43
char:sample_id,9,34,43
char:batch,9,34,43

ALL Controls miss the following variables: age, gender, path.score, inflam, active, stir, t1, fat.fraction, rna.seq, rnaseq,rank
2 FSHD patients missing 3 variables (path.score, inflam, active).

## Do the counts match the sample table? 
Yes

What `counts_header_check.txt` shows, and what you decided: 
It shows if the columns in counts file match the columns in the GEO records:
columns in counts file: 43
GSM records in GEO: 43
[title] columns matched: 43
[title] columns without a GEO record: 
[title] GEO records without a column: 
[char:sample_id] columns matched: 42
[char:sample_id] columns without a GEO record: 01-0022-1
[char:sample_id] GEO records without a column: 01-0022
first columns: 01-0021, 01-0022-1, 01-0023, 01-0024, 01-0025, 01-0026, 01-0027, 01-0029
Sample_id 01-0022-1 doesn't have a GEO record ; and sample_id 01-0022 in GEO records doesn't have a column. 

char:sample_id 01-0022 in sample_table.csv should be replaced with 01-0022-1


## Cross-tabulations: condition against the other variables

Paste or link the tables from `results/step1_audit/crosstab_*.csv`. For each
variable say in one sentence whether it is balanced across conditions and
whether it is technical (batch, instrument, tissue of sampling) or biological
(sex, age). A technical variable that is unbalanced across conditions is a
confounder: name it as such.

## char:batch  >> Technical; Unbalanced between controls and FSHD so considered as confounder

| condition   |   I |   II |   III |   IV |   All |
|:------------|----:|-----:|------:|-----:|------:|
| Control     |   0 |    0 |     6 |    3 |     9 |
| FSHD        |   5 |    5 |    16 |    8 |    34 |
| All         |   5 |    5 |    22 |   11 |    43 |

## instrument_model >> Technical; Balanced >> Not a confounder in this case

| condition   |   Illumina HiSeq 2500 |   All |
|:------------|----------------------:|------:|
| Control     |                     9 |     9 |
| FSHD        |                    34 |    34 |
| All         |                    43 |    43 |

## char:tissue >> Biological; All samples from muscle biopsy without specific details about muscle>> this Var Is not a confounder as such. 

| condition   |   muscle biopsy |   All |
|:------------|----------------:|------:|
| Control     |               9 |     9 |
| FSHD        |              34 |    34 |
| All         |              43 |    43 |

## char:gender >> Biological; I don't have records of gender among controls although in paper results, authors mentioned 6 females among controls. If I take results of paper, this variable is not balanced! 

| condition   |   female |   male |   missing |   All |
|:------------|---------:|-------:|----------:|------:|
| Control     |        0 |      0 |         9 |     9 |
| FSHD        |       12 |     22 |         0 |    34 |
| All         |       12 |     22 |         9 |    43 |


## char:age >> Biological; no records about age for controls; cannot assess balance from table; In results section of the paper, authors mention mean age of controls 35 (range 19-56); controls are younger and no control taken from group 60 and above whereas 13 FSHD patients are taken 60 and above. Not balanced but cannot assess from crosstabs.

| condition   |   [0, 30) |   [30, 45) |   [45, 60) |   [60, 120) |   missing |   All |
|:------------|----------:|-----------:|-----------:|------------:|----------:|------:|
| Control     |         0 |          0 |          0 |           0 |         9 |     9 |
| FSHD        |         2 |          5 |         14 |          13 |         0 |    34 |
| All         |         2 |          5 |         14 |          13 |         9 |    43 |


## What the paper says about the samples and GEO does not

Read the paper's description of the cohort and of the controls. List the
variables the paper describes that you cannot find in the GEO record.

Variables not found in the GEO record (Age and Gender are recorded for patients and not recorded for controls):
-FSHD type (1 or 2) knowing that 2 patients were FSHD2.
-Parental link between controls and patients (authors mention controls were chosen from spouses, parents, or siblings)
-Muscle biopsied (authors reported that biopsies of patients are taken from different muscles of lower extremities such as quadriceps, tibialis anterior, gastrocnemius, hamstring; whereas all control biospies are taken from quadriceps)
-Biospy side (authors reported that Biopsies of patients are taken from left or right side of the body; not reported for controls)
-10-point Clinical severity score CSS reporded for patients in the paper 
-Size in KB of chromosome Short4q reported for patients
-Other MRI measures: Quantitative DIXON fat/water maps; Fat saturated T2
-Needle technique to extract biopsies (Biopsies were obtained under sterile condition using either UCH modified Bergstrom needles (Rochester) or regular Bergstrom needles (Seattle))
-The muscle biopsy pathology score is found in the GEO records however the 4 variables to compute the score are not found (variability in fiber size, percent of centrally located nuclei, interstitial fibrosis and muscle fiber)


## FAIR score

Total: 25/30 (see `fair_score.md`). The two weakest criteria and why:
F2 (important variabes are missing and cannot be assessed to reproduce all paper results, compare groups for balance, and confounders )
R1 (we don't have a legal permission to reuse the data)

## Decisions I made

-I should inspect the raw gene counts to decide if any filtering is needed
-What about missing values in this case, do I leave them empty? Age, gender, ...
-How should I normalize the data for PCA and plots?
-DE takes raw conts
-Log1p?
-Should I perform feature selection top 2000 genes ?
-How many principal components ? screeplot
-What covariates should I consider: batch
- I should plot PCA 4 times colored by condition, batch, age , sex?




