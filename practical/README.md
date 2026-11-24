# Practical: reproduce a published analysis (GSE115650)

*From Nucleus to Neural Networks: Foundations of AI in Transcriptomics*,
AI4Health@DAISY group, LIG, Université Grenoble Alpes.
Deadlines, submission and grading are in the [README at the root](../README.md).

## The task

Redo parts of one published analysis from the counts its authors deposited on
GEO, and find out how far you can reproduce the published results.

You cannot reproduce this paper from its Methods section alone. In places the
paper is silent or unclear, and in others it does not match the deposited data.
Finding those places is part of the work. Each time you meet one, make a
decision, write it down and give your reason. Do not tune your analysis until
your numbers match the paper. Analyse the data soundly, then show with evidence
where your results differ from the paper and why.

You hand in two things:

| What | Where | Public? |
|------|-------|---------|
| Your code, your results and your filled worksheets | your fork, through your pull request | yes |
| Your report | UNESS | no |

## Dataset

| | |
|---|---|
| GEO series | [GSE115650](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115650) |
| Assay | bulk RNA-seq of skeletal muscle needle biopsies, Illumina HiSeq 2500, TruSeq RNA v2 |
| Samples | 43 biopsies: 34 from people with FSHD, 9 from unaffected controls |
| Processing by the authors | Trimmomatic 0.32, TopHat / Bowtie2 on hg38, GENCODE v24, counts by `summarizeOverlaps` |
| Variables per sample | subject status, sample_id, batch, tissue; for some samples: sex, age, MRI scores, pathology scores, an RNA-seq score and rank |
| Files on GEO | raw gene counts (1.6 MB), FPKM, normalised counts, an R object |

The raw count file is your input. The FPKM and normalised files let you check
your own normalisation against the authors'. Step 1 also downloads the GENCODE
v24 annotation (38.7 MB), because no file on GEO gives the gene symbols.

## The paper

Wang LH et al., *MRI-informed muscle biopsies correlate MRI with pathology and
DUX4 target gene expression in FSHD*, Hum Mol Genet 2019, 28(3):476-486.
DOI [10.1093/hmg/ddy364](https://doi.org/10.1093/hmg/ddy364), open access in
[PMC6337697](https://pmc.ncbi.nlm.nih.gov/articles/PMC6337697/). More links in
[paper/README.md](paper/README.md).

The authors score each biopsy from four DUX4 target genes, split the FSHD
biopsies into four groups with that score, relate the score to MRI signal and to
histology, and compare each group with the controls.

## The four tasks

The report is built on four tasks. Each one redoes a part of the paper that
depends only on files you can download. You do not redo the MRI scoring, the
histology or the immunostaining: use the values the authors deposited in GEO as
given.

| Task | Part of the paper | Script |
|---|---|---|
| 1. Recompute the four-gene DUX4 score and compare it with the score deposited in GEO | Figure 1, Methods | step 4, part A |
| 2. Redo the principal component analysis | Figure 3 | step 2 |
| 3. Redo the comparison of the score with MRI | Figure 2, Results | step 4, part A |
| 4. Run a differential expression analysis with PyDESeq2, compare at least two design formulas and justify the one you keep | Results, last section | step 3, then step 4, part B |

## The four scripts

| Script | Writes into | Feeds |
|--------|-------------|-------|
| `scripts/step1_fetch_and_audit.py` | `data/`, `results/step1_audit/` | all four tasks |
| `scripts/step2_qc_exploration.py` | `data/processed/`, `results/step2_qc/` | task 2, and the filtered counts of task 4 |
| `scripts/step3_differential_expression.py` | `results/step3_de/` | task 4 |
| `scripts/step4_compare_with_paper.py` | `results/step4_comparison/` | tasks 1, 3 and 4 |

**Step 1, fetch and audit.** The script downloads the GEO record, the count file
and the annotation, one request per file. It writes the sample table,
cross-tabulates the condition against batch, instrument, tissue, sex and age,
counts the samples that have a value for each variable, and checks the column
names of the count file against the GEO samples.

**Step 2, QC and exploration.** Filter low-count genes, plot library sizes,
compute a PCA on log2 CPM and colour it by condition, batch, sex, age and the
authors' variables. In this fork, the PCA uses the authors' method (DESeq2
rlog, 500 most variable genes; see "Running this fork's analysis"), and log2
CPM is kept only as a comparison.

**Step 3, differential expression.** PyDESeq2 on the filtered raw counts: all
FSHD biopsies against the controls, then each of the authors' four groups
against the controls. You choose the design, with or without batch.

**Step 4, comparison with the paper.** Part A: you write the biomarker score
from the paper's Methods, and the script compares it with the score the authors
deposited and with the MRI and pathology variables. Part B: you fill
`paper/paper_claims.csv` with the genes the paper reports, and the script joins
them with your results.

Each script has decision points marked `TODO(student)`: a threshold, a design,
a sample to keep or drop. Treat the default value as a starting point. Part 1 of
the report asks for the decisions that most affected your results, with the
evidence behind them.

## The worksheets

`worksheets/` holds one form per step and a FAIR grid. They are your lab
notebook: fill each one as you work, with the numbers your scripts produce and
the decisions you take, and commit it with your code.

| Worksheet | Fill it after | It prepares |
|-----------|---------------|-------------|
| `worksheets/step1_audit.md` | step 1 | the issues of Part 1 |
| `worksheets/fair_score.md` | step 1 | a metric for Part 2 |
| `worksheets/step2_qc.md` | step 2 | task 2 |
| `worksheets/step3_de.md` | step 3 | task 4 |
| `worksheets/step4_comparison.md` | step 4 | tasks 1, 3 and 4, and Part 2 |

The worksheets are not graded as such. They help you write the report, and the
report cites them: when you state a number or a decision in the report, link to
the line of the worksheet that records it, next to the code that produces it.
The grade goes to the reasoning in your report.

## The report

You write one report and submit it as a PDF on UNESS. It does not go into git.
It is short: one page per part, three pages in total.

| Part | Question |
|------|----------|
| 1. Your reasoning and the issues you met | What you did in the four tasks, how far each one reproduces the paper, and the two or three issues that mattered most |
| 2. Measuring and improving reproducibility | At least three metrics you define and compute, one reproducibility measure you code in your fork, and suggestions for the authors |
| 3. Your approach to reproducible and open science | What this exercise changed for you, and the tools and strategies you will keep |
| Bonus (optional) | Your Translator queries and what they tell you |

1. Download the report template from UNESS and write your report in `doc/`.
   The template is on UNESS only.
2. Support each factual claim with evidence: a permanent link to the lines of
   your fork that produce what you describe, or to the table or figure in
   `results/` that holds the number, or to the line of a worksheet that records
   your decision. The template explains how to get such a link. A claim without
   evidence does not count.
3. Export the report as a PDF and upload it on UNESS before the deadline.

## Bonus: Biolink and Translator

This part is optional and earns up to 2 bonus points. You take genes from your
own results and ask the NCATS Biomedical Data Translator what the knowledge
graphs already say about them. Translator describes each relationship with the
[Biolink Model](https://biolink.github.io/biolink-model/): a subject, a
predicate and an object, each with a category.

Use the production instance, and no other: **https://ui.transltr.io/**

1. **Log in with your GitHub account.** The site asks you to log in before you
   can type a search term. Create a Project for the practical and add each query
   to it, so that you find your queries again in your Query History.
2. **Run a disease query.** Smart Query, "What drugs may impact conditions
   related to ...", with facioscapulohumeral muscular dystrophy.
3. **Run a gene query from your own results.** Pick one gene from a table of
   `results/step3_de/` and say which table and why this gene. For a gene that
   is higher in FSHD, use "What chemicals may decrease the activity of a gene?".
   For a gene that is lower, use the question with "increase".
4. **Optional.** Try a Pathfinder Query between one of your genes and the
   disease. This query type is in beta.
5. **Save each query.** Click "Share Result Set" and keep the link. Translator
   does not promise to keep results, so write down the date, the number of
   results and the first results as well. Results load in several waves: click
   "Sync New Results" until no new result arrives before you count.
6. **Open the evidence.** In one result, click the term that connects two
   entities and read what supports the relationship.

Write it up in the bonus part of the report: your queries, the links, the results,
and what they tell you about your analysis and about the paper. A query that
returns nothing is a result too: report it and say what it means.

Translator is a research tool. Do not read its results as medical advice.

## Two folders that stay on your machine

`.gitignore` lists two folders. Git never commits what you put there.

| Folder | Who creates it | What goes there | Why git ignores it |
|--------|----------------|-----------------|--------------------|
| `data/` | step 1 | `data/raw/`: the GEO record, the count file, the GENCODE annotation. `data/processed/`: the filtered counts, the sample metadata and the gene names that the scripts derive | Anyone re-creates it by running step 1, and the files are too large for a repository |
| `doc/` | you (`mkdir doc`) | the PDF of the paper, its supplementary files, your reading notes, your report | The publisher distributes the paper, not you. Your notes and your report are yours and stay private |

Everything else in this folder goes into your fork: the scripts, `results/`,
`worksheets/`, `paper/paper_claims.csv`.

## How to run

You need [uv](https://docs.astral.sh/uv/) and git. uv installs Python 3.12 for
you (see `.python-version`). Work in your fork, inside this folder:

```bash
cd practical
```

```bash
uv sync
```

```bash
uv run python scripts/step1_fetch_and_audit.py
```

```bash
uv run python scripts/step2_qc_exploration.py
```

```bash
uv run python scripts/step3_differential_expression.py
```

```bash
uv run python scripts/step4_compare_with_paper.py
```

The four scripts run as they are, in a few minutes on a laptop. Commit your
scripts, `results/` and `worksheets/` each time you stop working.

### Running this fork's analysis

Step 2 here reproduces the paper's PCA with the authors' method (DESeq2 rlog), so it
needs R through a locked [Pixi](https://pixi.sh) environment in `r/`, and two inputs
beyond the default Step 1 downloads. From `practical/`:

```bash
uv sync --locked
uv run --locked python scripts/step1_fetch_and_audit.py --optional   # adds FPKM, normalised counts and the GEO R object
```

```bash
mkdir -p data/raw/authors
curl -fL -o data/raw/authors/sanitized.rlg.rda \
  https://raw.githubusercontent.com/TapscottLab/FSHD-Biopsy-RNA-seq/90cd343b3d60a696a35041a33015131ab9908164/data/sanitized.rlg.rda
sha256sum data/raw/authors/sanitized.rlg.rda   # 70b1758029229dbb89c33f3f0abc2dea5165d43f61e0da102566994c2dc18239
```

```bash
pixi install --manifest-path r/pixi.toml
pixi run --manifest-path r/pixi.toml Rscript r/pca_replay.R      # Figure 3 replay, sample link, mapping rates
pixi run --manifest-path r/pixi.toml Rscript r/rlog_from_geo.R   # rlog recomputed from the GEO counts
uv run --locked python scripts/step2_qc_exploration.py
uv run --locked python scripts/step2_qc_exploration.py --sensitivity   # PCA without the 3 low-mapping libraries
pixi run --manifest-path r/pixi.toml Rscript r/de_reproduction.R             # Step 3 Track A: the authors' DE, per contrast
uv run --locked python scripts/step3_compare_published.py doc/supp          # Track A vs supplementary tables S4-S10
uv run --locked python scripts/step3_differential_expression.py               # Step 3 Track B: independent DE, all planned fits
# Bonus: first extract the reference-list snapshot into data/raw/bonus_lists/ (bonus_lists/README.md)
uv run --locked python scripts/bonus_select_genes.py                          # Bonus: Translator genes, rule in worksheets/bonus_translator.md
```

Optional, with the paper's supplementary archive unzipped in `doc/supp/`:

```bash
uv run --locked python scripts/check_paper_sample_ids.py "doc/supp/Supplemental Table S1.docx"
```

## A pace that fits the milestones

| Date | Milestone in the root README | Where you should be |
|------|------------------------------|---------------------|
| 25 November | fork, first commit, pull request opened | step 1 runs on your machine and its results are committed |
| 4 December | first complete draft of the code | the four scripts run with your own decisions, the worksheets are filled |
| 16 December | final deadline | code final in your pull request, report submitted as a PDF on UNESS |

## Layout

```
practical/
|-- README.md
|-- pyproject.toml, uv.lock   dependencies (uv)
|-- .python-version           3.12
|-- .gitignore                excludes data/, doc/, .venv/
|-- scripts/
|   |-- geo.py                fetch helpers: one request per file, SOFT parser
|   |-- config.py             accession, files, condition column, variables
|   |-- step1_fetch_and_audit.py
|   |-- step2_qc_exploration.py
|   |-- step3_differential_expression.py
|   `-- step4_compare_with_paper.py
|-- results/                  tables and figures the scripts write, one folder per step
|-- worksheets/               one form per step and the FAIR grid: fill and commit them
|-- paper/                    links to the paper, paper_claims.csv
|-- data/                     NOT in git: downloads and derived tables
`-- doc/                      NOT in git: the paper, your notes, your report
```

## Before you trust a result

- Read `results/step1_audit/counts_header_check.txt` before any join between the
  counts and the sample table. In this fork, Step 2 joins on GEO title and stops with an error if a count column has no GEO record.
- Look at `missing_values.csv`. Some variables exist for one condition only.
- Compare the paper's description of the biopsies with the fields GEO records
  for them.
- Look at the batch cross-tabulation before you choose a design.
- The counts use GENCODE v24 identifiers and the paper names genes by symbol.
- The paper's groups come from the expression of DUX4 target genes. Keep that in
  mind when a test between a group and the controls returns DUX4 target genes.
