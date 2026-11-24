"""Step 1 (FAIR R6): compare the paper's Supplementary Table S2 with the three GEO matrices.

Table S2 gives, for the DUX4 candidate biomarkers and all 43 biopsies, the raw read counts,
the normalised counts and "RPKM". Each sheet is aligned with the matching GEO file on the
versioned gene ID and the biopsy name, and the largest absolute difference is reported.

Inputs
  doc/supp/Supplemental Table S2.xlsx            from the paper's supplementary archive
  data/raw/GSE115650_biopsy_genecounts.csv.gz    raw counts
  data/raw/GSE115650_biosy_normalized_genecounts.csv.gz   normalised counts (GEO's spelling)
  data/raw/GSE115650_biopsy_FPKM.csv.gz          FPKM
Output
  results/step1_audit/table_s2_check.csv         one row per sheet

Run, from practical/:
  uv run --locked python scripts/check_table_s2.py "doc/supp/Supplemental Table S2.xlsx"
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pandas as pd

import config

SHEETS = {"raw_counts": config.COUNTS_FILE,
          "normalized_counts": "GSE115650_biosy_normalized_genecounts.csv.gz",
          "RPKM": "GSE115650_biopsy_FPKM.csv.gz"}


def read_sheet(path: Path, sheet: str) -> pd.DataFrame:
    """Table S2 sheets have a few title rows above the header row that starts with gene_id."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")      # openpyxl: unsupported worksheet extension
        raw = pd.read_excel(path, sheet_name=sheet, header=None)
    header = raw.index[raw.iloc[:, 0].astype(str) == "gene_id"][0]
    table = raw.iloc[header + 1:].copy()
    table.columns = raw.iloc[header]
    table = table[table["gene_id"].notna()].set_index("gene_id").drop(columns="gene_name")
    return table.apply(pd.to_numeric, errors="coerce")


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    workbook = Path(sys.argv[1])
    rows = []
    for sheet, geo_file in SHEETS.items():
        s2 = read_sheet(workbook, sheet)
        geo = pd.read_csv(config.RAW / geo_file, index_col=0)
        genes = s2.index.intersection(geo.index)
        samples = s2.columns.intersection(geo.columns)
        diff = (s2.loc[genes, samples] - geo.loc[genes, samples]).abs()
        rows.append({"sheet": sheet, "geo_file": geo_file, "genes_in_s2": len(s2),
                     "genes_matched": len(genes), "samples_matched": len(samples),
                     "max_abs_difference": float(diff.max().max())})
    out = config.RESULTS / "step1_audit" / "table_s2_check.csv"
    table = pd.DataFrame(rows)
    table.to_csv(out, index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
