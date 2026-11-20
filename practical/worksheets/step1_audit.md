# Step 1: fetch and audit

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Keep it short and factual. Take each number from `results/step1_audit/`.

## What I downloaded

| File | Size | What it holds | Why I need it (or not) |
|------|------|---------------|------------------------|
|GSE115650_all.soft.txt |13.5 MB |GEO series record and sample metadata |To identify the biopsies and audit their metadata |
|GSE115650_biopsy_genecounts.csv.gz |1.7 MB |Raw gene counts for each biopsy |To check that the count columns match the sample records and provide inputs for analysis |
|gencode.v24.annotation.gtf.gz |38.7 MB |GENCODE v24 gene annotation |To map gene IDs to gene symbols and types |
| Wang et al. (2019), article PDF | 5.3 MB | Published methods, results, and figures | To compare the study design, methods, and results with GEO and my analysis. Downloaded by hand from the publisher; kept in `doc/`, not in git. |
| `ddy364_supp.zip` | 13.7 MB | Supplementary Tables S1–S10 and Figures S1–S7 | Table S1 supports the participant ID and metadata checks against GEO. Table S2 and the differential expression tables provide reference results for Step 4. Downloaded by hand from the publisher; kept in `doc/`, not in git. |
| `GSE115650_biopsy.dds.rda.gz` | 15.6 MB | DESeq2 object containing gene counts, FPKM, sample size factors, gene annotations (symbols, types) and per-sample TopHat alignment summaries | Inspected separately for the FAIR audit and to calculate mapping rates for all 43 samples. Downloaded from GEO over HTTPS with `curl`, not by the Step 1 script |
| `GSE115650_biosy_normalized_genecounts.csv.gz` | 7.1 MB | DESeq2-normalized gene counts (filename misspelled on GEO) | To check Supplementary Table S2 against GEO. Downloaded from GEO over HTTPS with `curl`, not by the Step 1 script |
| `GSE115650_biopsy_FPKM.csv.gz` | 8.6 MB | FPKM per gene and sample | To check Supplementary Table S2 and the inputs of Figures 4–5. Downloaded from GEO over HTTPS with `curl`, not by the Step 1 script |


## Sample table

Number of samples: 43
Condition groups and counts: 34 FSHD and 9 controls
Variables available per sample (list the columns of `sample_table.csv`): gsm, title, geo_accession, status, submission_date, last_update_date, type,
channel_count, source_name_ch1, organism_ch1, taxid_ch1, char:subject status,
char:age, char:gender, char:path.score, char:inflam, char:active, char:stir,
char:t1, char:fat.fraction, char:tissue, char:rnaseq.score, char:rnaseq.rank,
char:sample_id, char:batch, growth_protocol_ch1, molecule_ch1,
extract_protocol_ch1, data_processing, platform_id, contact_name,
contact_email, contact_laboratory, contact_department, contact_institute,
contact_address, contact_city, contact_state, contact_zip/postal_code,
contact_country, instrument_model, library_selection, library_source,
library_strategy, relation, supplementary_file_1, series_id, data_row_count,
condition
Variables that have a value for one condition only (see `missing_values.csv`): char:age, char:gender, char:path.score, char:inflam, char:active, char:stir,
char:t1, char:fat.fraction, char:rnaseq.score, char:rnaseq.rank


## Do the counts match the sample table?

What `counts_header_check.txt` shows, and what you decided:

All 43 sample columns in the count table match the GEO title field, while only 42 match char:sample_id. The single mismatch is record GSM3186680, which is labeled 01-0022-1 in title (matching the count column) but 01-0022 in char:sample_id. No other GEO sample uses 01-0022, and both title and char:sample_id are unique across the 43 records, so GEO holds one sample for this participant. The authors' code explains the suffix: it renames sequencing library 0952-1 to 01-0022-1 ([renameSample.R](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/renameSample.R#L17-L19)) and removes a library named 0952 before analysis ([makeSE.R](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/makeSE.R#L162)). Whether library 0952 comes from the same participant or the same biopsy is not stated in the code or the paper.
In Step 2, the merge key will be changed to title so the concerned biopsy will not be left unmatched and silently dropped. The choice was made so that the original data are not modified and the sample_id column will be kept in the future merged table for traceability and verification purposes.
Moreover, paper's Table S1 confirms that both keys are needed: it lists 01-0022 without the "-1" suffix, so title links GEO to the count matrix (43/43) while char:sample_id links GEO to the paper's tables (34/34 FSHD) (in results:`paper_table_s1_id_check.txt`)


## Cross-tabulations: condition against the other variables

Paste or link the tables from `results/step1_audit/crosstab_*.csv`. For each
variable say in one sentence whether it is balanced across conditions and
whether it is technical (batch, instrument, tissue of sampling) or biological
(sex, age). A technical variable that is unbalanced across conditions is a
confounder: name it as such.

Source: `results/step1_audit/crosstabs.md`

**Batch** (`char:batch`)

| condition | I | II | III | IV | All |
|:----------|--:|---:|----:|---:|----:|
| Control   | 0 |  0 |   6 |  3 |   9 |
| FSHD      | 5 |  5 |  16 |  8 |  34 |
| All       | 5 |  5 |  22 | 11 |  43 |

Batch is a partial confounder because Batches I and II contain only FSHD samples, meaning technical batch effects in those batches could look like disease effects (course 4).
The authors' code shows that each batch is one sequencing run, and the run names associate batches I–II with one instrument identifier and III–IV with another ([ReadMe.txt](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/ReadMe.txt#L3-L9)). Batches are therefore nested within instrument: instrument is determined by batch, so adjusting for batch also absorbs instrument effects, and the two cannot be estimated separately in the same model.

**Instrument** (`instrument_model`)

| condition | Illumina HiSeq 2500 | All |
|:----------|--------------------:|----:|
| Control   |                   9 |   9 |
| FSHD      |                  34 |  34 |
| All       |                  43 |  43 |


**Tissue** (`char:tissue`)

| condition | muscle biopsy | All |
|:----------|--------------:|----:|
| Control   |             9 |   9 |
| FSHD      |            34 |  34 |
| All       |            43 |  43 |


Instrument and tissue are technical. GEO gives every sample the same sequencer model (Illumina HiSeq 2500) and the same broad label (muscle biopsy), so these fields cannot confound condition as recorded. They hide real differences, however: the run names show two instrument identifiers, tied to batch (see above), and the paper reports four different muscles among FSHD biopsies but quadriceps for all controls (Table S1, series description).

**Sex** (`char:gender`)

| condition | female | male | missing | All |
|:----------|-------:|-----:|--------:|----:|
| Control   |      0 |    0 |       9 |   9 |
| FSHD      |     12 |   22 |       0 |  34 |
| All       |     12 |   22 |       9 |  43 |


**Age** (`char:age`, binned in years)

| condition | [0, 30) | [30, 45) | [45, 60) | [60, 120) | missing | All |
|:----------|--------:|---------:|---------:|----------:|--------:|----:|
| Control   |       0 |        0 |        0 |         0 |       9 |   9 |
| FSHD      |       2 |        5 |       14 |        13 |       0 |  34 |
| All       |       2 |        5 |       14 |        13 |       9 |  43 |


Sex and age (biological) are missing in GEO for all 9 controls, so their balance cannot be checked sample by sample. The paper's group summaries show that controls were younger on average (35 vs 53.6 years for the sequenced FSHD participants) and more often female (6/9 vs 12/34). GEO lacks the individual control values needed to adjust for either variable.


## What the paper says about the samples and GEO does not

Read the paper's description of the cohort and of the controls. List the
variables the paper describes that you cannot find in the GEO record.

Evidence comes from `results/step1_audit/paper_table_s1_id_check.txt`
(generated by `scripts/check_paper_sample_ids.py`),
`results/step1_audit/sample_table.csv`, `results/step1_audit/series_record.txt`,
and the publisher's article and supplement (Supplementary Table S1,
Results “Study cohort”, p. 477, and Methods, pp. 483–484).

| Detail described in the paper | Paper source | What GEO provides |
|---|---|---|
| Clinical Severity Score (CSS) and D4Z4 repeat-array size | Supplementary Table S1 | Neither value is provided per participant. |
| FSHD subtype and mosaic status | Supplementary Table S1 | The series description gives overall subtype counts. Sample records lack subtype and mosaic status. |
| Specific muscle biopsied and side | Supplementary Table S1 | `char:tissue` says `muscle biopsy` for all samples. The series description specifies quadriceps for controls. |
| Recruitment site | Methods | The paper reports 18 participants recruited at each site, Seattle and Rochester, with different biopsy needles (modified Bergstrom at Rochester, regular Bergstrom at Seattle); all controls were banked at Rochester. GEO lacks a participant-level site or needle field. |
| Strength of the biopsied muscle | Methods, Clinical assessments | The paper describes measurement with a hand-held dynamometer. GEO lacks the measurements. |
| MRI measurements by anatomical region | Supplementary Table S1 and Figure S7 | GEO STIR and T1 match the whole-muscle ratings in Table S1 for all 34 sequenced FSHD biopsies, and GEO fat fraction matches Table S1, whose fat-fraction column leaves the region unspecified (`paper_table_s1_id_check.txt`, section 3). GEO provides no separate biopsy-region STIR, quantitative T2, or fat-fraction fields. |
| Control age and sex | Results, Study cohort | The paper reports six women and three men, with a mean age of 35 years (range 19–56). It gives these figures for the group as a whole. Age and sex are blank for each control sample in GEO. |
| Control relationship to people with FSHD and genetic confirmation | Results, Study cohort | GEO lacks these details. The paper describes unaffected spouses, parents, or siblings and reports genetic confirmation for at-risk controls. |

GEO records a sequencing batch (I to IV) for every sample. The paper does not report these batches. The authors' code links each batch to one Illumina sequencing run, dated 2 Oct 2015 (I), 8 Feb 2016 (II), 8 Sep 2016 (III) and 19 Jan 2017 (IV); the run folder names carry one instrument identifier for batches I–II and another for III–IV ([ReadMe.txt](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/ReadMe.txt#L3-L9)).

Additional discrepancies: The article and GEO both report 34 FSHD1 and 2 FSHD2 participants among the 36 recruited, and report that 34% were women.
The Methods state that all procedures, including muscle biopsies, were identical at both sites, but also report different biopsy needles and MRI skin markers at each site.
Supplementary Table S1 lists three participants labelled FSHD2, one separately labelled mosaic, and 13 women among 36 participants (36%). The GEO sample table contains 12 women among the 34 sequenced FSHD participants (35%).

Library selection: The authors' alignment summary lists 49 library IDs. Their pre-filter count object `se.rda` contains 48 libraries (37 FSHD and 11 controls), with `2374` already absent. Applying the five removals listed in `makeSE.R` (`0952`, `32-0008-1`, `32-0002b`, `2374-1` and `2374-2`) yields exactly the 43 library IDs in the saved `sanitized.dds.rda` and the GEO count object. The paper does not explain these selection decisions. The retained `32-0008` mapped at 67.7%, and the excluded `32-0008-1` mapped at 96.2%. The files checked do not establish whether the suffixes denote resequencing or distinct biopsies. [Library selection code](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/makeSE.R#L151-L164) · [Alignment summaries](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/data/tophat_summary.rda)

## FAIR score

Total: 22/30 (see `fair_score.md`). The two weakest criteria and why: R1 (0/2), no dataset-specific licence or reuse terms were identified; R2 (1/2), because the paper and GEO sample metadata misprint the TopHat version and omit read layout and strandedness. These gaps were recovered from SRA RunInfo and the authors' code, but the software versions remain inconsistent between the paper and the authors' README.

## Decisions I made

- Use GEO `title` to join the count matrix to the sample table in Step 2.
  It matches all 43 count columns, including `01-0022-1`. Keep
  `char:sample_id` to link the 34 FSHD participants to Supplementary
  Table S1, where that participant is `01-0022`. Resolve the join in code
  and preserve the deposited identifiers. See `counts_header_check.txt`
  and `paper_table_s1_id_check.txt`.
- Start the main analysis with all 43 deposited biopsies. In Step 2,
  examine library size and mapping rate by batch and flag the three
  low-mapping biopsies (`32-0008`, `01-0037`, `32-0009`). All three are
  FSHD biopsies from Batch II. I will rerun key results in Steps 2–4
  without these three to see how sensitive the conclusions are.
- Keep versioned Ensembl IDs as gene identifiers and use GENCODE v24 to
  annotate gene symbols and types. All 27,506 deposited gene IDs matched
  this annotation release.
- Treat batch as partially confounded with condition. Examine it in
  the PCA in Step 2 and compare differential-expression designs with
  and without batch in Step 3.
- Consult the authors' code to fill gaps in the paper and identify which
  findings come from it. [Authors' code at commit 90cd343](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/tree/90cd343b3d60a696a35041a33015131ab9908164)
