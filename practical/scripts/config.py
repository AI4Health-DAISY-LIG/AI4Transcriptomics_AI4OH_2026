"""Settings for GSE115650 (bulk RNA-seq, FSHD muscle biopsies, Wang et al. 2019)."""

from pathlib import Path

ACC = "GSE115650"

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

# Series-level supplementary files on GEO. Only the raw gene counts feed steps 2 and 3.
COUNTS_FILE = "GSE115650_biopsy_genecounts.csv.gz"  # raw counts, 1.6 MB
OPTIONAL_FILES = [
    "GSE115650_biopsy_FPKM.csv.gz",  # FPKM, 8.2 MB (not for DESeq2)
    "GSE115650_biosy_normalized_genecounts.csv.gz",  # DESeq2-normalised, 6.7 MB (the typo is GEO's)
    "GSE115650_biopsy.dds.rda.gz",  # R object (DESeqDataSet), 14.9 MB, needs R
]

# The counts use GENCODE v24 gene identifiers and no file on GEO gives the gene symbols.
# Step 1 downloads the annotation the authors name (38.7 MB, one request) and extracts
# one row per gene into data/processed/gene_names.csv.
ANNOTATION_FILE = "gencode.v24.annotation.gtf.gz"
ANNOTATION_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/"
    + ANNOTATION_FILE
)
GENE_NAMES = PROCESSED / "gene_names.csv"

# Column of the sample table that defines the condition, and short labels for it.
CONDITION_COL = "char:subject status"
CONDITION_LABELS = {
    "Control": "Control",
    "Facioscapulohumeral muscular dystrophy (FSHD)": "FSHD",
}
REFERENCE_LEVEL = "Control"

# Variables to cross-tabulate against the condition in step 1 and to colour the PCA in step 2.
# Decide in your report which ones are technical and which ones are biological.
TECHNICAL_VARS = ["char:batch", "instrument_model", "char:tissue"]
BIOLOGICAL_VARS = ["char:gender", "char:age"]
# Pathology, MRI and RNA-seq score variables deposited by the authors.
PATHOLOGY_VARS = [
    "char:stir",
    "char:active",
    "char:inflam",
    "char:path.score",
    "char:rnaseq.rank",
    "char:rnaseq.score",
    "char:fat.fraction",
    "char:t1",
]

# Add two infer variables from th meta data
INFER_VARS = ["char:infer_gender", "char:site"]

# The four DUX4 target genes the paper uses as biomarkers, by symbol as the paper names them.
BIOMARKERS = ["LEUTX", "KHDC1L", "PRAMEF2", "TRIM43"]
# The authors' group (1 to 4) and score for each FSHD biopsy, as deposited on GEO.
GROUP_COL = "char:rnaseq.rank"
SCORE_COL = "char:rnaseq.score"


# Outliers found by the PCA methods in step2 QC
OUTLIERS = ["01-0037", "32-0016"]
