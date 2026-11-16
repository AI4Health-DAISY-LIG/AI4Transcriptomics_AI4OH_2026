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

Number of samples: 
Condition groups and counts: 
Variables available per sample (list the columns of `sample_table.csv`): 
Variables that have a value for one condition only (see `missing_values.csv`): 

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
