"""Settings for GSE115650 (bulk RNA-seq, FSHD muscle biopsies, Wang et al. 2019)."""

from pathlib import Path

ACC = "GSE115650"

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

# Series-level supplementary files on GEO. The raw gene counts are the input of the Python
# steps 2 and 3; the R object feeds the replay and reproduction scripts in r/.
COUNTS_FILE = "GSE115650_biopsy_genecounts.csv.gz"          # raw counts, 1.6 MB
OPTIONAL_FILES = [
    "GSE115650_biopsy_FPKM.csv.gz",                          # FPKM, 8.2 MB (not for DESeq2)
    "GSE115650_biosy_normalized_genecounts.csv.gz",          # DESeq2-normalised, 6.7 MB (the typo is GEO's)
    "GSE115650_biopsy.dds.rda.gz",                           # R object (DESeqDataSet), 14.9 MB, needs R
]

# The counts use GENCODE v24 gene identifiers and no CSV file on GEO gives the gene symbols
# (the R object carries them).
# Step 1 downloads the annotation the authors name (38.7 MB, one request) and extracts
# one row per gene into data/processed/gene_names.csv.
ANNOTATION_FILE = "gencode.v24.annotation.gtf.gz"
ANNOTATION_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/" + ANNOTATION_FILE
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
PATHOLOGY_VARS = ["char:stir", "char:active", "char:inflam", "char:path.score",
                  "char:rnaseq.rank", "char:rnaseq.score", "char:fat.fraction", "char:t1"]

# The four DUX4 target genes the paper uses as biomarkers, by symbol as the paper names them.
BIOMARKERS = ["LEUTX", "KHDC1L", "PRAMEF2", "TRIM43"]
# The authors' group (1 to 4) and score for each FSHD biopsy, as deposited on GEO.
GROUP_COL = "char:rnaseq.rank"
# Biopsies excluded in the sensitivity runs: mapping rate below 90%, all FSHD from batch II
# (32-0008 67.73%, 01-0037 76.52%, 32-0009 85.60%; every other biopsy 94.85% or more,
# results/step2_qc/mapping_rates.csv). The main analysis keeps all 43 biopsies.
SENSITIVITY_EXCLUDE = ["32-0008", "01-0037", "32-0009"]
# Step 3, Track B (own analysis, not the reproduction): sample rules fixed after QC and before
# any DE result (23 Nov). (a) 32-0004: unresolved disagreement between recorded sex
# and the expression profile (results/step2_qc/sex_check.csv); a provenance flag, not a proven
# swap. (b) mapping rate below 90% (SENSITIVITY_EXCLUDE). Biological extremes (32-0016) are kept.
TRACK_B_EXCLUDE = SENSITIVITY_EXCLUDE + ["32-0004"]
# Recruitment site from the participant ID prefix. 32 = Seattle in the authors' follow-up
# workbook (FredHutch/RWellstone_FSHD_muscle_biopsy @ 27e705c, Sheet1, Visit 1); 01 = Rochester,
# inferred from that roster and the paper's 18 participants per site. The nine controls were
# banked at the University of Rochester (paper, Materials and Methods).
SITE_BY_PREFIX = {"01": "Rochester", "32": "Seattle"}
SCORE_COL = "char:rnaseq.score"
