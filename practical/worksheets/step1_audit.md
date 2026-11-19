# Step 1: fetch and audit

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Keep it short and factual. Take each number from `results/step1_audit/`.

## What I downloaded

| File | Size | What it holds | Why I need it (or not) |
|------|------|---------------|------------------------|
| GSE115650_all.soft.txt | 13.5 MB | GEO series record: metadata of the 43 samples | Tells which sample is FSHD or control. It also gives batch, sex, age and the authors' scores per sample |
| GSE115650_biopsy_genecounts.csv.gz | 1.7 MB | Raw gene counts, genes x 43 samples (GENCODE v24 IDs) | Main input. Raw counts, not normalised values, because I normalise myself and PyDESeq2 requires raw counts |
| gencode.v24.annotation.gtf.gz | 38.7 MB | GENCODE v24 annotation: gene ID to gene symbol | The counts use gene IDs while the paper names genes by symbol (e.g. ZSCAN4, LEUTX); needed to compare with the paper |

## Sample table

Number of samples: 43
Condition groups and counts: FSHD 34, Control 9
Variables available per sample (list the columns of `sample_table.csv`): subject status, age, gender, path.score, inflam, active, stir, t1, fat.fraction, tissue, rnaseq.score, rnaseq.rank, sample_id, batch (the char: columns; the other columns are GEO administrative fields)
Variables that have a value for one condition only (see `missing_values.csv`): age, gender, path.score, inflam, active, stir, t1, fat.fraction, rnaseq.score and rnaseq.rank have values for FSHD only (0 of 9 controls). path.score, inflam and active are also missing for 2 FSHD samples (32 of 34)

## Do the counts match the sample table?

What `counts_header_check.txt` shows, and what you decided: The counts file has 43 columns and GEO has 43 samples. Matching on `title` pairs all 43 columns; matching on `char:sample_id` pairs only 42: column `01-0022-1` has no sample_id match, and sample_id `01-0022` has no column. Both names belong to the same GEO record (GSM3186680, title `01-0022-1`, sample_id `01-0022`, FSHD, batch III), so this is a naming difference between two GEO fields, not a missing sample.
Decision: match counts to samples on `title` (43/43) to keep all 43 samples.

## Cross-tabulations: condition against the other variables

Paste or link the tables from `results/step1_audit/crosstab_*.csv`. For each
variable say in one sentence whether it is balanced across conditions and
whether it is technical (batch, instrument, tissue of sampling) or biological
(sex, age). A technical variable that is unbalanced across conditions is a
confounder: name it as such.

See results/step1_audit/crosstabs.md for the tables.

- instrument_model (technical): constant, all 43 samples on Illumina HiSeq 2500, so it cannot confound.
- tissue (technical, sampling): constant, all 43 are muscle biopsies, so it cannot confound.
- batch (technical): unbalanced. Controls are only in batches III (6) and IV (3); batches I and II hold 10 FSHD samples and no control. Batch is a confounder of condition.
- gender (biological): recorded for FSHD only (12 female, 22 male), missing for all 9 controls. Balance cannot be checked, and sex cannot be used as a covariate.
- age (biological): recorded for FSHD only (mostly 45 and over), missing for all 9 controls. Balance cannot be checked, and age cannot be used as a covariate.

## What the paper says about the samples and GEO does not

Read the paper's description of the cohort and of the controls. List the
variables the paper describes that you cannot find in the GEO record.

-

## FAIR score

Total: /30 (see `fair_score.md`). The two weakest criteria and why:

## Decisions I made

- Join counts to the sample table on `title`, not `sample_id`: 43/43 vs 42/43 matches; the mismatch (01-0022-1 vs 01-0022) is the same sample, GSM3186680.
- Treat batch as a confounder to examine in step 2 (PCA coloured by batch) and step 3 (design with and without batch).
-
