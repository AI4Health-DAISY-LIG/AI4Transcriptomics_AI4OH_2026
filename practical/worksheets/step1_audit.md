# Step 1: fetch and audit

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Keep it short and factual. Take each number from `results/step1_audit/`.

## What I downloaded

| File | Size | What it holds | Why I need it (or not) |
|------|------|---------------|------------------------|
| `GSE115650_all.soft.txt` | 13.5 MB | Complete series and sample metadata from GEO for GSE115650 | Needed to extract sample attributes, batch, conditions (Control vs FSHD), and clinical scores |
| `GSE115650_biopsy_genecounts.csv.gz` | 1.7 MB | Raw RNA-seq read counts for 43 muscle biopsy samples across 60,554 genes | Needed as primary numerical input for QC, filtering, and PyDESeq2 differential expression |
| `gencode.v24.annotation.gtf.gz` | 38.7 MB | Human genome annotation (GENCODE release 24, hg38) | Needed to map Ensembl gene IDs to official gene symbols (e.g. LEUTX, TRIM43) because the count file lacks gene symbols |

## Sample table

Number of samples: 43
Condition groups and counts: FSHD = 34, Control = 9
Variables available per sample (list the columns of `sample_table.csv`): ['gsm', 'title', 'geo_accession', 'status', 'submission_date', 'last_update_date', 'type', 'channel_count', 'source_name_ch1', 'organism_ch1', 'taxid_ch1', 'char:subject status', 'char:age', 'char:gender', 'char:path.score', 'char:inflam', 'char:active', 'char:stir', 'char:t1', 'char:fat.fraction', 'char:tissue', 'char:rnaseq.score', 'char:rnaseq.rank', 'char:sample_id', 'char:batch', 'growth_protocol_ch1', 'molecule_ch1', 'extract_protocol_ch1', 'data_processing', 'platform_id', 'contact_name', 'contact_email', 'contact_laboratory', 'contact_department', 'contact_institute', 'contact_address', 'contact_city', 'contact_state', 'contact_zip/postal_code', 'contact_country', 'instrument_model', 'library_selection', 'library_source', 'library_strategy', 'relation', 'supplementary_file_1', 'series_id', 'data_row_count', 'condition']
Variables that have a value for one condition only (see `missing_values.csv`): 
char:age, char:gender, char:path.score, char:inflam, 
char:active, char:stir, char:t1, char:fat.fraction, 
char:rnaseq.score, char:rnaseq.rank - all FSHD only, missing for Controls.

## Do the counts match the sample table?

What `counts_header_check.txt` shows, and what you decided: 
The file shows 43 columns in the counts file and 43 GSM records in GEO. 
Matching by sample title gives 43/43 matches. Matching by char:sample_id 
gives only 42/43 - one sample is named "01-0022-1" in the counts file but 
"01-0022" in GEO. These are the same sample; the "-1" suffix is a typo in 
the counts file. Decision: use the title column for matching so all 43 samples are kept.

## Cross-tabulations: condition against the other variables

Paste or link the tables from `results/step1_audit/crosstab_*.csv`. For each
variable say in one sentence whether it is balanced across conditions and
whether it is technical (batch, instrument, tissue of sampling) or biological
(sex, age). A technical variable that is unbalanced across conditions is a
confounder: name it as such.

### char:batch (Technical)

| condition   |   I |   II |   III |   IV |   All |
|:------------|----:|-----:|------:|-----:|------:|
| Control     |   0 |    0 |     6 |    3 |     9 |
| FSHD        |   5 |    5 |    16 |    8 |    34 |
| All         |   5 |    5 |    22 |   11 |    43 |

Technical variable. Severely unbalanced - batches I and II have zero controls. This makes batch a **confounder**: it is mixed up with condition.

### instrument_model (Technical)

| condition   |   Illumina HiSeq 2500 |   All |
|:------------|----------------------:|------:|
| Control     |                     9 |     9 |
| FSHD        |                    34 |    34 |
| All         |                    43 |    43 |

Technical variable. Perfectly balanced - all 43 samples used the same instrument. Not a confounder.

### char:tissue (Technical)

| condition   |   muscle biopsy |   All |
|:------------|----------------:|------:|
| Control     |               9 |     9 |
| FSHD        |              34 |    34 |
| All         |              43 |    43 |

Technical variable. Constant across all samples — all are muscle biopsies. Not a confounder.

### char:gender (Biological)

| condition   |   female |   male |   missing |   All |
|:------------|---------:|-------:|----------:|------:|
| Control     |        0 |      0 |         9 |     9 |
| FSHD        |       12 |     22 |         0 |    34 |
| All         |       12 |     22 |         9 |    43 |

Biological variable. Gender is completely missing for all 9 controls so balance cannot be assessed; sex cannot be used as a covariate in the FSHD vs Control model.

### char:age (Biological)

| condition   |   [0, 30) |   [30, 45) |   [45, 60) |   [60, 120) |   missing |   All |
|:------------|----------:|-----------:|-----------:|------------:|----------:|------:|
| Control     |         0 |          0 |          0 |           0 |         9 |     9 |
| FSHD        |         2 |          5 |         14 |          13 |         0 |    34 |
| All         |         2 |          5 |         14 |          13 |         9 |    43 |

Biological variable. Age is completely missing for all 9 controls so balance cannot be assessed; age cannot be used as a covariate in the FSHD vs Control model.

## What the paper says about the samples and GEO does not

Read the paper's description of the cohort and of the controls. List the
variables the paper describes that you cannot find in the GEO record.

- Age and sex of the 9 control subjects - the paper mentions unaffected individuals were biopsied but their demographics are not in GEO.
- The exact muscle that was biopsied for each sample (e.g. gastrocnemius, tibialis anterior, quadriceps, hamstring). GEO just says "muscle biopsy" for everyone.
- D4Z4 repeat count for each FSHD patient - the paper describes patients genetically confirmed with D4Z4 repeat contraction on chromosome 4q35, but this is not in GEO.
- Clinical Severity Score (CSS) - a 10-grade validated score measured per subject, not deposited in GEO.
- Muscle strength measurements taken with a hand-held dynamometer - mentioned in Methods but absent from GEO metadata.

## FAIR score

Total: /30 (see `fair_score.md`). The two weakest criteria and why:

## Decisions I made

-
