"""Step 3, Track A: compare the reproduced DE tables with the published supplementary tables.

For each published contrast (S4, S5, S6, S8, S10), it compares
  1. gene membership: the recomputed genes at padj < 0.05 and |log2FC| >= 2 vs the genes in the table;
  2. the shared genes' statistics: baseMean, log2FoldChange, lfcSE, padj;
  3. the per-sample normalised counts in the table vs raw count / full-cohort size factor,
     which shows which size factors the published fits used.

Input:
  results/step3_de/track_a/de_<contrast>.csv.gz  written by r/de_reproduction.R
  results/step3_de/track_a/size_factors.csv      idem
  data/raw/GSE115650_biopsy_genecounts.csv.gz    GEO raw counts
  Supplemental Table S4/S5/S6/S8/S10.xlsx        the paper's supplementary archive, in doc/ (git-ignored)

Output:
  results/step3_de/track_a/published_comparison.csv
  results/step3_de/track_a/published_comparison.txt

Run:  uv run --locked python scripts/step3_compare_published.py doc/supp
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import config

TRACK_A = config.RESULTS / "step3_de" / "track_a"
TABLES = {"S4": "S4_G1_outliers", "S5": "S5_G1", "S6": "S6_G2", "S8": "S8_G3", "S10": "S10_G4"}
STATS = ["baseMean", "log2FoldChange", "lfcSE", "padj"]


def read_published(path: Path) -> pd.DataFrame:
    """First sheet; row 0 is a title, row 1 the header, the first column the gene ID."""
    table = pd.read_excel(path, sheet_name=0, header=1)
    table = table.rename(columns={table.columns[0]: "gene_id"}).set_index("gene_id")
    return table[table.index.astype(str).str.startswith("ENSG")]


def sample_columns(table: pd.DataFrame) -> dict[str, str]:
    """Normalised-count columns such as Control_01.0041 -> GEO title 01-0041."""
    return {c: c.split("_", 1)[1].replace(".", "-") for c in table.columns
            if c.startswith(("Control_", "FSHD_"))}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("supp_dir", type=Path, help="folder with the Supplemental Table files")
    supp = parser.parse_args().supp_dir

    counts = pd.read_csv(config.RAW / "GSE115650_biopsy_genecounts.csv.gz", index_col=0)
    size_factors = pd.read_csv(TRACK_A / "size_factors.csv", index_col="sample")[
        "size_factor_deposited"]
    normalised = counts / size_factors[counts.columns]

    rows = []
    for table_id, contrast in TABLES.items():
        published = read_published(supp / f"Supplemental Table {table_id}.xlsx")
        recomputed = pd.read_csv(TRACK_A / f"de_{contrast}.csv.gz", index_col=0)
        called = recomputed[(recomputed["padj"] < 0.05) & (recomputed["log2FoldChange"].abs() >= 2)]
        shared = published.index.intersection(called.index)
        row = {"table": table_id, "contrast": contrast,
               "published": len(published), "reproduced": len(called), "shared": len(shared),
               "published_only": len(published.index.difference(called.index)),
               "reproduced_only": len(called.index.difference(published.index))}
        for stat in STATS:
            a = published.loc[shared, stat].astype(float)
            b = called.loc[shared, stat].astype(float)
            row[f"{stat}_max_abs_diff"] = float((a - b).abs().max())
        cols = sample_columns(published)
        pub_norm = published.loc[shared, list(cols)].astype(float).rename(columns=cols)
        recomputed_norm = normalised.loc[shared, pub_norm.columns]
        rel = ((pub_norm - recomputed_norm).abs() / recomputed_norm.where(recomputed_norm > 0)).stack()
        row["norm_samples"] = len(cols)
        row["norm_max_rel_diff"] = float(rel.max())
        rows.append(row)

    result = pd.DataFrame(rows)
    result.to_csv(TRACK_A / "published_comparison.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", None):
        text = result.to_string(index=False, float_format=lambda v: f"{v:.3g}")
    (TRACK_A / "published_comparison.txt").write_text(
        "Reproduced (padj < 0.05, |log2FC| >= 2) vs published supplementary tables\n"
        "norm_max_rel_diff: published normalised counts vs raw / full-cohort size factor\n\n"
        + text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
