# FAIR score of the GEO series

*Worksheet. Fill it as you work and commit it with your code. Your report on UNESS cites it as evidence.*

Score each criterion 0 (absent), 1 (partly), 2 (fully met). Write one line of
evidence per criterion: the GEO field, the file name or the sentence of the
paper you relied on. The total is out of 30. The number matters less than the
evidence, so do not round up.

| Principle | Criterion | Score (0 to 2) | Evidence |
|-----------|-----------|----------------|----------|
| F1 | The dataset has a persistent identifier (GEO accession, DOI of the paper) | 2 | GSE115650 on GEO, DOI in paper. |
| F2 | Sample metadata are rich: condition, donor, sex, age, tissue, batch are present for every sample | 1 | Age and sex missing for all 9 controls (missing_values.csv). |
| F3 | The metadata point to the data files and the paper, and the paper points back to the accession | 2 | GEO links to PubMed; paper cites GSE115650. |
| F4 | The series is indexed and searchable (GEO search, PubMed link) | 2 | GSE115650 search on GEO works. |
| A1 | The files can be retrieved by a standard protocol (HTTPS or FTP) without a login | 2 | Downloaded in step 1 without login. |
| A2 | The metadata stay available even if the files were withdrawn (GEO keeps the record) | 2 | GEO keeps records permanently. |
| I1 | The data use an open, documented format (CSV or TSV, not a binary from one tool) | 2 | Counts file is a plain CSV. |
| I2 | Gene identifiers follow a named reference (which genome build, which annotation release?) | 1 | Build not on GEO; found hg38 and GENCODE v24 in paper Methods. |
| I3 | The metadata use controlled vocabularies, or at least consistent spelling, for condition and tissue | 1 | Condition labels consistent but sample 01-0022 has a typo in counts file. |
| R1 | A licence or usage terms are stated | 1 | No explicit licence on GEO. |
| R2 | The processing provenance is stated: aligner, annotation, counting method, software versions | 1 | Tools named in Methods but no scripts or exact versions provided. |
| R3 | Raw counts are deposited, not only normalised values | 2 | Raw counts file is on GEO. |
| R4 | The sample table lets you rebuild the groups compared in the paper without guessing | 1 | Groups in rnaseq.rank but score formula only in paper. |
| R5 | The deposited files and the paper agree on the number of samples | 2 | Both say 43 samples. |
| R6 | A reader can tell which deposited file feeds which figure of the paper | 1 | Not stated on GEO, had to read Methods. |
| **Total** | | **23/30** | |

## What blocks the reproduction

- Controls have no age or sex in GEO so I could not adjust for them.
- The paper does not say whether batch was included in the model; I will test both ways.
- No filtering threshold is given in the paper; I picked one myself.
- The score formula the paper uses is not directly available in PyDESeq2; I will use log2 CPM as a substitute.
- One sample ID is written differently in the counts file and GEO; I matched it by hand.
