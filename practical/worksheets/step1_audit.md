# Step 1: Fetch and Audit

## Downloaded Files

| File                                           | Size | Contents                                                                       | Required?                                                                                                                        |
| :--------------------------------------------- | :--- | :----------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------------------------- |
| `GSE115650_all.soft.txt`                       | 13M  | GEO Soft format metadata (titles, characteristics, platform)                   | **Yes** — maps samples to conditions (Control vs. FSHD).                                                                         |
| `GSE115650_biopsy_genecounts.csv.gz`           | 1.7M | Raw gene-level expression count matrix from muscle biopsies                    | **Yes** — primary input matrix for differential expression analysis.                                                             |
| `gencode.v24.annotation.gtf.gz`                | 37M  | GENCODE v24 gene annotation                                                    | **Yes** — maps gene symbol IDs to gene names and biotypes.                                                                       |
| `GSE115650_biopsy.dds.rda.gz`                  | 15M  | Pre-built R `DESeqDataSet` object from the original study                      | **Benchmark** — useful for inspecting the authors' original design formula, sample filtering, and intermediate `DESeq2` objects. |
| `GSE115650_biosy_normalized_genecounts.csv.gz` | 6.8M | Author-normalized gene count matrix _(note original `biosy` typo in filename)_ | **Benchmark** — useful for verifying how closely our normalized counts match the paper's processed matrix.                       |
| `GSE115650_biopsy_FPKM.csv.gz`                 | 8.3M | FPKM-normalized expression matrix                                              | **Benchmark** — not used for DE modeling, but useful for cross-checking expression levels against figures/tables in the paper.   |

## Sample Table Overview

- **Total samples in metadata:** 43 (Control: 9, FSHD: 34) _(Note: 42 samples in expression matrix after audit)_
- **Variables present only in FSHD (`Control = 0` in `missing_values.csv`):** `char:age`, `char:gender`, `char:path.score`, `char:inflam`, `char:active`, `char:stir`, `char:t1`, `char:fat.fraction`, `char:rnaseq.score`, `char:rnaseq.rank`.

<details>
<summary><strong>All 49 columns in (without infer variables) <code>sample_table.csv</code></strong> (click to expand)</summary>

`gsm`, `title`, `geo_accession`, `status`, `submission_date`, `last_update_date`, `type`, `channel_count`, `source_name_ch1`, `organism_ch1`, `taxid_ch1`, `char:subject status`, `char:age`, `char:gender`, `char:path.score`, `char:inflam`, `char:active`, `char:stir`, `char:t1`, `char:fat.fraction`, `char:tissue`, `char:rnaseq.score`, `char:rnaseq.rank`, `char:sample_id`, `char:batch`, `growth_protocol_ch1`, `molecule_ch1`, `extract_protocol_ch1`, `data_processing`, `platform_id`, `contact_name`, `contact_email`, `contact_laboratory`, `contact_department`, `contact_institute`, `contact_address`, `contact_city`, `contact_state`, `contact_zip/postal_code`, `contact_country`, `instrument_model`, `library_selection`, `library_source`, `library_strategy`, `relation`, `supplementary_file_1`, `series_id`, `data_row_count`, `condition`,

</details>

## Counts vs. Sample Table Check (`counts_header_check.txt`)

- **Discrepancy:** One sample header is labeled `01-0022` in the count matrix instead of `01-0022-1`.
- **Resolution:** Retained the sample and standardized its ID to `01-0022-1` to match the metadata `title`.

## Cross-Tabulations (`results/step1_audit/crosstab_*.csv`)

Values are shown as **Count (Row %)**.

### Technical Variables

**`char:batch`** — Technical variable; unbalanced across conditions (Controls are absent from Batches I & II), making it a potential confounder.

| Condition   | I         | II        | III        | IV        | Total  |
| :---------- | :-------- | :-------- | :--------- | :-------- | :----- |
| **Control** | 0 (0.0%)  | 0 (0.0%)  | 6 (66.7%)  | 3 (33.3%) | **9**  |
| **FSHD**    | 5 (15.2%) | 5 (15.2%) | 15 (45.5%) | 8 (24.2%) | **33** |
| **All**     | **5**     | **5**     | **21**     | **11**    | **42** |

**`char:site`** — Technical variable; significantly unbalanced across conditions (all Controls come from Site 01), making it a primary confounder.

| Condition   | 01         | 32         | Total  |
| :---------- | :--------- | :--------- | :----- |
| **Control** | 9 (100.0%) | 0 (0.0%)   | **9**  |
| **FSHD**    | 15 (45.5%) | 18 (54.5%) | **33** |
| **All**     | **24**     | **18**     | **42** |

**`instrument_model` & `char:tissue`** — Technical variables; balanced across conditions (100% `Illumina HiSeq 2500` and 100% `muscle biopsy` for all 9 Control and 33 FSHD samples).

### Biological Variables

**`char:gender` & `char:infer_gender`** — Biological variable; 100% missing for Controls in GEO, and unbalanced using inferred sex (66.7% female in Control vs. 33.3% in FSHD).

| Condition   | GEO: Female | GEO: Male  | GEO: Missing | Inferred: Female | Inferred: Male | Total  |
| :---------- | :---------- | :--------- | :----------- | :--------------- | :------------- | :----- |
| **Control** | 0 (0.0%)    | 0 (0.0%)   | 9 (100.0%)   | 6 (66.7%)        | 3 (33.3%)      | **9**  |
| **FSHD**    | 12 (36.4%)  | 21 (63.6%) | 0 (0.0%)     | 11 (33.3%)       | 22 (66.7%)     | **33** |
| **All**     | **12**      | **21**     | **9**        | **17**           | **25**         | **42** |

**`char:age`** — Biological variable; 100% missing for Controls in GEO (missingness strongly associated with condition, adjusted $p = 5.61 \times 10^{-9}$), preventing direct balance assessment from metadata alone.

| Condition   | <30      | 30–44     | 45–59      | 60+        | Missing    | Total  |
| :---------- | :------- | :-------- | :--------- | :--------- | :--------- | :----- |
| **Control** | 0 (0.0%) | 0 (0.0%)  | 0 (0.0%)   | 0 (0.0%)   | 9 (100.0%) | **9**  |
| **FSHD**    | 2 (6.1%) | 5 (15.2%) | 13 (39.4%) | 13 (39.4%) | 0 (0.0%)   | **33** |
| **All**     | **2**    | **5**     | **13**     | **13**     | **9**      | **42** |

## Statistical Summary of Potential Confounders

All the informations can be found in the following files :
[crosstabs.md](../results/step1_audit/crosstabs.md)
[statistical_tests.csv](../results/step1_audit/statistical_tests.csv)

| Variable       | Type       | Balance Status              | Statistical Evidence                                                  | Assessment                                                                                  |
| :------------- | :--------- | :-------------------------- | :-------------------------------------------------------------------- | :------------------------------------------------------------------------------------------ |
| **Site**       | Technical  | **Unbalanced (Confounder)** | Fisher's exact $p = 0.0054$ (adj. $p = 0.0090$), Cramér's $V = 0.430$ | Strongly associated with condition. Controls originate exclusively from Site 01.            |
| **Batch**      | Technical  | **Potentially unbalanced**  | $\chi^2 = 3.587$, adj. $p = 0.3097$, Cramér's $V = 0.112$             | Controls are absent from Batches I–II. Lacks statistical significance but warrants caution. |
| **Age**        | Biological | **Potentially unbalanced**  | Paper mean: 54 yrs (FSHD) vs. 35 yrs (Control)                        | Likely biological confounder, but 100% missing in GEO for Controls prevents adjustment.     |
| **Gender**     | Biological | **Potentially unbalanced**  | Fisher's exact $p = 0.1242$, Cramér's $V = 0.234$ (inferred)          | 66.7% female in Control vs. 33.3% in FSHD. Non-significant but notable imbalance.           |
| **Tissue**     | Technical  | **Broadly balanced**        | 100% labeled `muscle biopsy`                                          | Broad label matches, though specific sampled muscle types vary (per supplemental files).    |
| **Instrument** | Technical  | **Balanced**                | 100% `Illumina HiSeq 2500`                                            | No platform-related confounding.                                                            |

## Information Reported in Paper but Missing from GEO

- **FSHD subtypes:** Cohort comprises 34 FSHD1 and 2 FSHD2 patients (unannotated in GEO).
- **Control demographics & kinship:** Controls include 6 females and 3 males (unmapped to sample IDs), have a younger mean age (35 vs. 54 years), and are relatives (spouses, parents, or siblings) of FSHD patients.
- **Specific muscle sampled:** Exact anatomical biopsy locations for FSHD samples are only available in supplemental files.

## Pre-Sequencing & Wet-Lab Sources of Variation

Unmeasured anatomical and wet-lab differences between Control and FSHD biopsies directly enter the estimated condition effect.

| Stage & Factor                        | Protocol / Reporting Status                                                                                                                                                                                            | Impact on Expression Counts                                                                                                                                                              | Adjustable?                                                    |
| :------------------------------------ | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :------------------------------------------------------------- |
| **Anatomical muscle & MRI targeting** | Controls: _quadriceps_. FSHD: _tibialis anterior_, _gastrocnemius_, _vastus lateralis_, or _biceps femoris_ (MRI-guided to T2-STIR hyperintense regions).                                                              | Fiber-type mix (slow vs. fast myosin), metabolism, and vascularization are collinear with condition; STIR targeting enriches FSHD samples for active disease (71% express DUX4 targets). | **No**                                                         |
| **Biopsy needle & tissue splitting**  | Modified Bergstrom needle (Rochester, Site 01) vs. regular Bergstrom (Seattle, Site 32). 75–150 mg split into 3 pieces (2 frozen in $\text{LN}_2$ within 30s for RNA, 1 for histology; stored at $-80^\circ\text{C}$). | Needle differences alter tissue yield and local trauma across sites; spatial patchiness between RNA and histology pieces dilutes pathology correlations for focal _DUX4_ targets.        | **Partly** (restrict to Site 01 in sensitivity checks)         |
| **Cell-type composition**             | Bulk tissue mix (myofibers, fat, fibroblasts, immune/endothelial cells); high adipocyte signal noted in `32-0016` and `32-0006`.                                                                                       | Fat replacement and fibrosis shift cell proportions—upregulated genes may reflect tissue remodeling rather than intracellular regulation.                                                | **Diagnostic only** (marker-gene scores & FSHD `fat.fraction`) |
| **RNA extraction, RIN & library kit** | **Not reported** (no extraction protocol, RIN scores, or poly-A vs. rRNA-depletion details).                                                                                                                           | Degradation introduces $3'$ bias and lowers gene detection; library chemistry dictates coverage of non-coding, mitochondrial, and histone RNAs.                                          | **No** (library size & detected genes as QC proxies)           |
| **Sequencing depth & `batch`**        | Depth and physical meaning of `char:batch` (prep date, library prep, or flow cell) are **not reported**.                                                                                                               | Depth sets the detection floor (6/9 Controls have zero reads across the 4 DUX4 biomarkers).                                                                                              | **Partly** (library size in QC; `batch` tested cautiously)     |
| **Kinship & physiological state**     | Relatives/spouses share genetics/environment; fasting, circadian time, exercise, and medications are **not reported**.                                                                                                 | Sample non-independence makes permutation tests optimistic; age tracks PC4; acute state alters immediate-early genes.                                                                    | **No**                                                         |

### Downstream Interpretation Rules

- **Composite condition effect:** Estimated condition fold-changes reflect **disease + muscle type + site/needle + age + relatedness**, alongside modeled covariates (`infer_gender` and partially `batch`).
- **Signal vs. background:** Large, disease-specific signals (DUX4 targets in Groups 3–4 vs. near-zero Control baselines) remain robust, whereas subtle shifts in genes sensitive to muscle type, age, or site cannot be isolated.
- **Internal baseline:** Group 1 (FSHD samples lacking DUX4 biomarker expression) helps flag genes driven by muscle type, site, or non-DUX4 disease background rather than DUX4 activity.
- **Bounding bias:** Pre-sequencing variation cannot be corrected post hoc; it will be bounded via sensitivity designs (with/without `batch`, Site 01 only, Batches III–IV only) and reported as a study limitation.

## FAIR Score (`fair_score.md`)

**Total Score:** 17 / 30

1. **F2 (Weakest):** Essential sample metadata (`sex`, `age`, `batch`) are incomplete across conditions, limiting covariate adjustment.
2. **R4 (Second weakest):** Filtering rules and thresholds used to define the paper's four comparison groups are omitted, preventing exact group reconstruction.

## Key Analytical Decisions

- **Sample ID alignment:** Renamed count matrix column `01-0022` to `01-0022-1` to match metadata.
- **Sex covariate:** Included `char:infer_gender` as a model covariate to compensate for missing Control sex metadata.
- **Batch handling:** Will evaluate `char:batch` impact cautiously, as the absence of Controls in Batches I and II limits reliable batch correction.
