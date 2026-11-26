# FAIR score of the GEO series

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Score each criterion 0 (absent), 1 (partly), 2 (fully met). Write one line of
evidence per criterion: the GEO field, the file name or the sentence of the
paper you relied on. The total is out of 30. The number matters less than the
evidence, so do not round up.

| Principle | Criterion | Score (0 to 2) | Evidence |
|-----------|-----------|----------------|----------|
| F1 | The dataset has a persistent identifier (GEO accession, DOI of the paper) | 2 | GEO accession GSE115650; paper DOI 10.1093/hmg/ddy364 |
| F2 | Sample metadata are rich: condition, donor, sex, age, tissue, batch are present for every sample | 1 | Condition, donor id (`sample_id`), tissue and batch present for all 43; sex and age missing for the 9 controls (`missing_values.csv`); recruitment site only implicit in the id prefix; muscle type absent |
| F3 | The metadata point to the data files and the paper, and the paper points back to the accession | 1 | GEO record gives `pubmed_id = 30312408` and the supplementary files; the paper's Methods write "deposited to GEO series GSExxxxx": no accession in the paper |
| F4 | The series is indexed and searchable (GEO search, PubMed link) | 2 | Indexed in GEO, linked to PubMed 30312408 and to BioProject PRJNA475678 / SRA SRP150311 (`series_record.txt`) |
| A1 | The files can be retrieved by a standard protocol (HTTPS or FTP) without a login | 2 | Downloaded by script over HTTPS from ftp.ncbi.nlm.nih.gov, no login |
| A2 | The metadata stay available even if the files were withdrawn (GEO keeps the record) | 2 | GEO keeps series and sample records; each sample also links to BioSample and SRA (`relation` field) |
| I1 | The data use an open, documented format (CSV or TSV, not a binary from one tool) | 2 | Counts, FPKM and normalised counts are gzipped CSV; the extra `biopsy.dds.rda` is an R binary but is not needed |
| I2 | Gene identifiers follow a named reference (which genome build, which annotation release?) | 2 | `data_processing`: genome build hg38, GENCODE version 24; ids are versioned GENCODE ids, no symbol column (symbols rebuilt from the GTF) |
| I3 | The metadata use controlled vocabularies, or at least consistent spelling, for condition and tissue | 1 | Condition has 2 consistent values; tissue is the uninformative "muscle biopsy"; `source_name_ch1` labels controls "Group1_Control" although controls have no group; file name typo "biosy_normalized_genecounts" |
| R1 | A licence or usage terms are stated | 0 | No licence or usage terms in the series record or the sample records |
| R2 | The processing provenance is stated: aligner, annotation, counting method, software versions | 1 | Trimmomatic 0.32, TopHat 2.1.0 / Bowtie2 2.2.6, GENCODE v24, `summarizeOverlaps` IntersectionStrict are stated; DESeq2 version, gene filter, design and thresholds are not; R 3.4 / Bioconductor 3.5 only in the paper |
| R3 | Raw counts are deposited, not only normalised values | 2 | `GSE115650_biopsy_genecounts.csv.gz` holds raw integer counts |
| R4 | The sample table lets you rebuild the groups compared in the paper without guessing | 1 | Authors' groups deposited (`rnaseq.rank`, 10/5/13/6 as in the paper), but built by a clustering the paper does not specify; the "group 1 outliers" (01-0037, 32-0016) and the muscle type are not flagged |
| R5 | The deposited files and the paper agree on the number of samples | 2 | 43 count columns = 34 FSHD + 9 controls, as in the paper |
| R6 | A reader can tell which deposited file feeds which figure of the paper | 0 | Neither GEO nor the paper map files to figures; the authors' code repository is cited but I did not use it |
| **Total** | | **21/30** | |

## What blocks the reproduction

List here what you could not find in GEO or in the paper and needed to decide
for yourself (thresholds, covariates, versions, sample exclusions). This list
feeds Part 2 of your report.

- Gene filtering rule before DESeq2: not stated.
- DESeq2 design and covariates (batch? site?): not stated; DESeq2 version not stated.
- Thresholds to call a gene differentially expressed (padj, log2 fold change): not stated.
- How the four groups were cut from the hierarchical clustering (distance, linkage, cut): not stated; I use the deposited `rnaseq.rank`.
- Biomarker score: "sum of the four biomarkers' log2 difference from the control average" after rlog; the exact computation (per gene mean of controls, rlog settings) must be decided, and PyDESeq2 has no rlog.
- Muscle type and recruitment site per sample: not in GEO (muscle type only in Table S1, site only in the id prefix).
- Age and sex of the controls: not in GEO (only summary figures in the paper).
- Meaning of the "-1" suffix of sample `01-0022-1` (issue #22).