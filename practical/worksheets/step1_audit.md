# Step 1: fetch and audit

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Keep it short and factual. Take each number from `results/step1_audit/`.

## What I downloaded

| File | Size | What it holds | Why I need it (or not) |
|GSE115650_all.soft.txt|13.5 MB|GEO study and sample metadata|To retrieve sample information|
|GSE115650_biopsy_genecounts.csv.gz |GSE115650_biopsy_genecounts.csv.gz | 1.7 MB |Raw gene counts |Needed for gene expression analysis |
|GENCODE v24 annotation|38.7 MB|Gene IDs and gene names|Needed to annotate genes|

## Sample table

Number of samples: 43
Condition groups and counts: 34 FSHD, 9 Control 
Variables available per sample (list the columns of `sample_table.csv`): GEO sample identifiers, title, status, source, organism and sample characteristics (subject status, age, gender, path.score, inflam, active, stir, t1, fat.fraction, tissue, rnaseq.score, rnaseq.rank, sample_id and batch).
Variables that have a value for one condition only (see `missing_values.csv`): Age, gender and several MRI/pathology variables are available only for FSHD samples.

## Do the counts match the sample table?
The count file contains 43 samples, matching the 43 samples in the GEO sample table. However, one sample identifier differs between the two files: 01-0022-1 in the count file and 01-0022 in GEO.
What `counts_header_check.txt` shows, and what you decided: 
All 43 count columns match the GEO sample titles. When matching by sample_id, 42 samples match and one does not. I kept the original identifiers and documented the discrepancy rather than changing the sample ID without confirmation.

## Cross-tabulations: condition against the other variables

Paste or link the tables from `results/step1_audit/crosstab_*.csv`. For each
variable say in one sentence whether it is balanced across conditions and
whether it is technical (batch, instrument, tissue of sampling) or biological
(sex, age). A technical variable that is unbalanced across conditions is a
confounder: name it as such.

Available in results/step1_audit/.

Batch: The samples not equally distributed between the 2 groups. All controls belong to batches III and IV, while FSHD samples are found in all four batches. Since batch is a technical factor, it could affect the results and act as a confounder.

Instrument: no difference between the groups because all samples were analysed using the same instrument model.

Tissue: In GEO, all samples are described as muscle biopsies. However, the paper mentions that control samples came from the quadriceps, while FSHD samples came from different muscles.

Age: Available for FSHD samples but missing for all controls, so the two groups cannot be compared using the GEO data. Age is a biological variable.

Gender: same problem: information is available for FSHD samples but missing for controls. Gender is also a biological variable.

## What the paper says about the samples and GEO does not

Read the paper's description of the cohort and of the controls. List the
variables the paper describes that you cannot find in the GEO record.

-The paper gives some additional info that is not available in the GEO:
- The FSHD group included 34 patients with FSHD1 and 2 with FSHD2.
- The average age of the FSHD group was 54 years, and the average age of the controls was 35 years.
- The paper gives the clinical severity score (CSS), with an average of 5 for the FSHD group.
- It explains that the nine controls were unaffected spouses, parents or siblings of FSHD patients.
- It gives the gender distribution of the controls (6 females and 3 males).
- It specifies that all control biopsies were taken from the quadriceps, while FSHD biopsies came from different muscles.
- It also describes the recruitment centres and the participants' eligibility criteria, including age, muscle weakness and ability to walk independently.

## FAIR score

Total: /30 (see `fair_score.md`). The two weakest criteria and why:
R1 – Licence: I could not find clear information about the licence or the conditions for reusing the dataset.

F2 – Sample metadata: Some information is missing from the GEO record. In particular, age and gender are not available for the control samples, which makes the sample description incomplete.

## Decisions I made

I used the GEO metadata to identify the samples and their conditions.
kept the original sample IDs and noted the difference between 01-0022-1 and 01-0022 instead of changing it myself.
I did not try to fill in the missing age and gender information for the controls.
I noted the batch imbalance because it could affect the analysis.
I kept the difference between the 36 recruited patients and the 34 FSHD samples with RNA-seq data and I recorded the downloaded files and their purpose in my audit.
