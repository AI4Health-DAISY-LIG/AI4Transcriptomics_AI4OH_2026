# FAIR score of the GEO series

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Score each criterion 0 (absent), 1 (partly), 2 (fully met). Write one line of
evidence per criterion: the GEO field, the file name or the sentence of the
paper you relied on. The total is out of 30. The number matters less than the
evidence, so do not round up.

| Principle | Criterion | Score (0 to 2) | Evidence |
|-----------|-----------|----------------|----------|
| F1 | The dataset has a persistent identifier (GEO accession, DOI of the paper) |2 |persistent DOI https://doi.org/10.1093/hmg/ddy364 and persistent GEO accession https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115650 |
| F2 | Sample metadata are rich: condition, donor, sex, age, tissue, batch are present for every sample |1 |sample metadata are relatively more rich for patients than controls with many variables missing for both (sample_table.csv)|
| F3 | The metadata point to the data files and the paper, and the paper points back to the accession |2 | Taken from the paper (The RNA seq data and gene counts have been deposited to GEO series GSExxxxx. Codes for preprocessing pipeline and reproducible analysis can be found in TapscottLab’s GitHub repository,https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq); https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115650|
| F4 | The series is indexed and searchable (GEO search, PubMed link) | 2|https://www.ncbi.nlm.nih.gov/geo/ |
| A1 | The files can be retrieved by a standard protocol (HTTPS or FTP) without a login |2 |Not logged in. https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE115650 |
| A2 | The metadata stay available even if the files were withdrawn (GEO keeps the record) | 2|https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html |
| I1 | The data use an open, documented format (CSV or TSV, not a binary from one tool) |2 |https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE115650 |
| I2 | Gene identifiers follow a named reference (which genome build, which annotation release?) |2 |Taken from the paper: "To profile the gene expression, the gene features and annotation were first obtained from GENCODE,version 24." |
| I3 | The metadata use controlled vocabularies, or at least consistent spelling, for condition and tissue |2 |"sample_table.csv" |
| R1 | A licence or usage terms are stated |0 | CC BY or CC0 for data Not stated anywhere|
| R2 | The processing provenance is stated: aligner, annotation, counting method, software versions |2 |Taken from paper methods: "The preprocessing of RNA-Seq data started with filtering unqualified reads, followed by trimming Illumina universal adapter using Trimmomatic-0.32. Subsequently, the preprocessed reads were aligned to genome build hg38 using Tophat-21.0/Bowtie2-2.2.6. To profile the gene expression, the gene features and annotation were first obtained from GENCODE, version 24. Gene counts were then estimated by using Bioconductor function GenomicAlignment::summarizeOverlap with IntersectionStrict mode, which restrict read counts to those that fall completely within the range of exons.The bioinformatics analysis and visualization were mostly performed using R/3.4 and Bioconductor 3.5 packages." |
| R3 | Raw counts are deposited, not only normalised values |2 |Series RNA-seq raw counts matrix Raw count matrices may be suitable for input into differential expression analysis tools. GSE115650_raw_counts_GRCh38.p13_NCBI.tsv.gz |
| R4 | The sample table lets you rebuild the groups compared in the paper without guessing |2 |in "sample_table.csv" Column source_name_ch1|
| R5 | The deposited files and the paper agree on the number of samples |2 |34 patients and 9 controls in deposited files, similar to paper |
| R6 | A reader can tell which deposited file feeds which figure of the paper |1 | Not explicitly stated by the authors in the GEO series or the paper; in GEO https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE115650 we have statements such as Raw count matrices may be suitable for input into differential expression analysis tools. GSE115650_raw_counts_GRCh38.p13_NCBI.tsv.gz|
| **Total** || **26/30** | |

## What blocks the reproduction

List here what you could not find in GEO or in the paper and needed to decide
for yourself (thresholds, covariates, versions, sample exclusions). This list
feeds Part 2 of your report.

-Although data are available but I didn't find any statement from the authors' to the licence to legally reuse the data.
-I didn't find some technical and biological variables that could be used as covariates such as FSHD type, age, gender of controls, biopsy needle type, 

