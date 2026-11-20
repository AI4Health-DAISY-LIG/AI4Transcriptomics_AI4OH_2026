# FAIR score of the GEO series

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Score each criterion 0 (absent), 1 (partly), 2 (fully met). Write one line of
evidence per criterion: the GEO field, the file name or the sentence of the
paper you relied on. The total is out of 30. The number matters less than the
evidence, so do not round up.

| Principle | Criterion | Score (0 to 2) | Evidence |
|-----------|-----------|----------------|----------|
| F1 | The dataset has a persistent identifier (GEO accession, DOI of the paper) |2|GEO series GSE115650 with one GSM per sample; paper DOI 10.1093/hmg/ddy364; BioProject PRJNA475678 and SRA SRP150311 (`series_record.txt`) |
| F2 | Sample metadata are rich: condition, donor, sex, age, tissue, batch are present for every sample |1|Condition, sample ID, tissue and batch are present for all 43 samples; no separate donor field; age and sex are missing for the 9 controls (`missing_values.csv`) |
| F3 | The metadata point to the data files and the paper, and the paper points back to the accession |1|GEO links the files and PMID 30312408 (`series_record.txt`); the paper's RNA-Seq analysis section prints "GSExxxxx" in place of the accession |
| F4 | The series is indexed and searchable (GEO search, PubMed link) |2|Indexed in GEO DataSets as 200115650; NCBI links it to PMID 30312408 in both directions (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=gds&id=30312408) |
| A1 | The files can be retrieved by a standard protocol (HTTPS or FTP) without a login |2|Step 1 retrieves the SOFT record and the count file over HTTPS without login (`scripts/geo.py`); GEO also lists FTP links |
| A2 | The metadata stay available even if the files were withdrawn (GEO keeps the record) |1|For published accessions, GEO retains the record and adds a withdrawal notice. Record contents may be modified or deleted. [GEO FAQ, "How can I delete my records?"](https://www.ncbi.nlm.nih.gov/geo/info/faq.html) |
| I1 | The data use an open, documented format (CSV or TSV, not a binary from one tool) |1|Raw counts, normalized counts and FPKM are available as CSV, and sample metadata as SOFT. Per-sample TopHat alignment summaries (input and mapped reads) are stored only in `GSE115650_biopsy.dds.rda.gz`, an R binary, and are absent from the open-format files ([makeSE.R](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/makeSE.R#L18-L37)) |
| I2 | Gene identifiers follow a named reference (which genome build, which annotation release?) |2|`data_processing` names hg38 and GENCODE v24; all 27,506 count IDs are unique, versioned and found in GENCODE v24 (`data/processed/gene_names.csv`) |
| I3 | The metadata use controlled vocabularies, or at least consistent spelling, for condition and tissue |2|`char:subject status` uses one consistent label per condition; all 43 `char:tissue` entries are "muscle biopsy"; no ontology terms (`sample_table.csv`) |
| R1 | A licence or usage terms are stated |0|No dataset-specific licence or reuse terms in the GEO record; NCBI's general terms are not dataset-specific; the article is "All rights reserved"; the authors' code repository has no licence file |
| R2 | The processing provenance is stated: aligner, annotation, counting method, software versions |1|GEO and the paper name Trimmomatic 0.32, TopHat/Bowtie2, hg38, GENCODE v24 and summarizeOverlaps IntersectionStrict, but misprint the TopHat version ("21.0") and omit read layout and strandedness. These gaps were recovered elsewhere: SRA RunInfo lists all 43 runs as single-end, and the authors' code gives TopHat 2.1.0, Bowtie2 2.2.6 and unstranded settings. Versions remain uncertain: the paper and the `makeSE.R` header give R 3.4, the README R 3.5 / Bioconductor 3.7 ([SRA RunInfo](https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=SRP150311); [authors' code](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/tree/90cd343b3d60a696a35041a33015131ab9908164)) |
| R3 | Raw counts are deposited, not only normalised values |2|`GSE115650_biopsy_genecounts.csv.gz` is the count file, separate from the normalised file; all 1,182,758 entries (27,506 genes × 43 samples) are non-negative integers; raw reads in SRA (SRP150311) |
| R4 | The sample table lets you rebuild the groups compared in the paper without guessing |2|`condition` and `char:rnaseq.rank` give the 9 controls and FSHD groups 1–4 with the paper's sizes (10/5/13/6), and `source_name_ch1` agrees. The published DE tables use the same 9 controls: S4 compares the two group 1 PCA outliers (01-0037, 32-0016) with them, and S5, S6, S8 and S10 compare groups 1–4, with sample columns matching GEO's group memberships. Among the 34 sequenced FSHD biopsies, "normal-appearing" MRI needs STIR = 0 and T1 = 0 (10), and STIR+ is STIR > 0 (23) (`sample_table.csv`). The labels can be reused, but the stated clustering method alone does not reproduce them (see "What blocks the reproduction") |
| R5 | The deposited files and the paper agree on the number of samples |2|43 count columns and 43 GEO records (`counts_header_check.txt`). The paper's 36 recruited FSHD participants include 2 without an RNA-seq sample, leaving 34 sequenced FSHD participants plus 9 controls (`paper_table_s1_id_check.txt`) |
| R6 | A reader can tell which deposited file feeds which figure of the paper |1|Neither the paper nor GEO maps files to figures. All 54 genes × 43 biopsies in each of the three Supplementary Table S2 sheets match the corresponding GEO CSV (raw counts identical; normalized counts and the S2 "RPKM" sheet, which matches the GEO "FPKM" file, to rounding precision; check run outside this repository). The authors' code shows that Figures 4–5 plot log10(FPKM + 0.5), an offset missing from the legends, and that Figure 3 uses their saved rlog object, not a GEO file ([figure code](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/makeManuscripts_supFigs.R#L75-L141)) |
| **Total** | | **22/30** | Findable 6/8, Accessible 3/4, Interoperable 5/6, Reusable 8/12 |

## What blocks the reproduction

List here what you could not find in GEO or in the paper and needed to decide
for yourself (thresholds, covariates, versions, sample exclusions). This list
feeds Part 2 of your report.

- **Control covariates.** The paper gives age and sex only as group
  summaries for the controls, and GEO leaves both fields blank for all
  nine. The deposited metadata cannot support individual age or sex
  adjustment in the FSHD–control comparison.
- **Batches.** GEO lists batches I–IV and the authors' code identifies
  them as sequencing runs. The paper gives no batch adjustment strategy.
  Controls occur only in batches III and IV, so within-batch
  FSHD–control comparisons are possible in those two batches. I will
  examine batch in Step 2 and compare DE designs in Step 3.
- **Original library selection.** The alignment summary lists 49 libraries.
  `2374` is already absent from the 48-library pre-filter count object.
  The five removals listed in the authors' code then yield exactly the
  43 deposited libraries. The paper does not explain these decisions.
  I will start from the 43 libraries in the deposited count matrix.
  [Library selection code](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/makeSE.R#L151-L164)
- **Read layout and strandedness.** SRA RunInfo identifies all 43 runs
  as single-end. The authors' code supplies the unstranded settings.
  [SRA RunInfo](https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=SRP150311) · [Alignment settings](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/do_tophat.sh#L15-L22)
- **Mapping QC.** GEO stores the per-sample TopHat alignment summaries
  in `GSE115650_biopsy.dds.rda.gz`. The open-format GEO files omit
  these fields. I will document how I extract and use them in Step 2.
- **Gene filter.** The paper gives no threshold. The smallest total count
  per gene in the deposited matrix is 21, consistent with the authors'
  rule of more than 20 counts across samples. I will document the
  additional filtering rule I choose in Step 2.
  [Authors' filter](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/makeSE.R#L161-L166)
- **Published PCA.** The paper does not state the input or settings for
  Figure 3. The authors' code runs DESeq2's `plotPCA` (500 most variable
  genes by default) on an rlog transform saved as an R object, not as a
  GEO file. A check run outside this repository (NumPy on that saved rlog
  matrix, 500 most variable genes) gave 37.94% for PC1 and 13.77% for PC2,
  matching the published 38% and 14%. Reproducing this in my code, and
  comparing it with a fresh PCA from the GEO counts, is part of Step 2.
  [PCA code](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/makeManuscripts_supFigs.R#L75-L93)
- **DE design and thresholds.** The Methods name DESeq2 and omit the
  design and cutoffs. The authors' code fits `~ pheno_type` separately
  for each comparison and reports adjusted p < 0.05 and absolute log2
  fold change ≥ 2. The published tables use the nine GEO controls: S4
  compares the two group 1 PCA outliers with them, and S5, S6, S8 and S10
  compare groups 1–4, with sample columns matching GEO's group
  memberships. The code's historical-control comparisons are separate
  analyses. I will document my design and thresholds in Step 3.
  [Authors' DE helper](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/misc_tools.R#L320-L359)
- **Biomarker score.** The authors fit DESeq2's rlog (`fitType="local"`)
  on the four biomarker genes only and sum each gene's difference from
  the control mean
  ([getRankScore](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/misc_tools.R#L84-L129)).
  The locked PyDESeq2 release lacks rlog, so I will document the
  transformation I use in Step 4.
- **RNA groups.** The Methods describe hierarchical clustering into four
  groups, but the authors' code relabels the clusters by hand: it assigns
  five samples to group 2 by their row position in a score-sorted table
  and merges two clusters into group 4. The deposited labels can be
  reused, but the stated method alone does not reproduce them.
  [Group relabelling](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/makeManuscripts_supFig_scatterPlot.R#L17-L22)
- **MRI measure level.** GEO's STIR and T1 are the whole-muscle ratings
  (they match Table S1's whole-muscle columns), which supports the
  Figure 2 comparison. The biopsy-region correlations in Supplementary
  Figure S7 require measurements absent from GEO.
- **Figure 2 test.** The paper reports P < 0.001 for MRI-normal vs STIR+
  scores but does not name the test; the Biostatistical analysis section
  mentions only Spearman correlations. I will choose and name the test in
  Step 4.
- **Software versions.** The paper reports R 3.4 and Bioconductor 3.5, and
  the header of `makeSE.R` names R 3.4.0, but the authors' README lists
  R 3.5 and Bioconductor 3.7 as requirements. The versions may have varied
  between scripts, so the exact environment that generated the published
  results is uncertain. I will report the versions used for my runs.
  [makeSE.R header](https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/makeSE.R#L1-L5)
