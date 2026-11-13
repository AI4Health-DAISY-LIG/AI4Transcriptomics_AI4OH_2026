"""Step 1: fetch GSE115650 by code and audit it.

Outputs (results/step1_audit/):
  sample_table.csv            one row per GSM with every GEO field and a `condition` column
  series_record.txt           the !Series_* fields as fetched
  crosstab_<variable>.csv     condition against each variable in config.TECHNICAL_VARS + BIOLOGICAL_VARS
  crosstabs.md                the same tables in one Markdown file for the report
  missing_values.csv          for each GEO variable, how many samples of each condition have a value
  counts_header_check.txt     do the column names of the counts file match the sample table?
  data/processed/gene_names.csv   gene identifier, symbol and type from GENCODE v24

Run:  uv run python scripts/step1_fetch_and_audit.py [--optional]
"""

from __future__ import annotations

import argparse
import gzip
import re

import pandas as pd

import config
import geo

OUT = config.RESULTS / "step1_audit"


def gene_names_from_gtf(gtf_gz, dest) -> pd.DataFrame:
    """One row per gene of a GENCODE GTF: versioned id, id without version, symbol, type."""
    rows = []
    with gzip.open(gtf_gz, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 9)
            if fields[2] != "gene":
                continue
            attr = dict(re.findall(r'(\S+) "([^"]*)"', fields[8]))
            rows.append({"gene_id": attr["gene_id"], "gene_id_noversion": attr["gene_id"].split(".")[0],
                         "gene_name": attr.get("gene_name", ""), "gene_type": attr.get("gene_type", "")})
    table = pd.DataFrame(rows)
    table.to_csv(dest, index=False)
    return table


def main(optional: bool) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    config.PROCESSED.mkdir(parents=True, exist_ok=True)

    # 1. Series and sample records (one request).
    soft = geo.fetch_soft(config.ACC, config.RAW, targ="all")
    series, samples = geo.parse_soft(soft)
    (OUT / "series_record.txt").write_text(
        "\n".join(f"{k} = {v}" for k, vals in series.items() for v in vals), encoding="utf-8"
    )

    # 2. Sample table with a clean condition column.
    samples["condition"] = samples[config.CONDITION_COL].map(config.CONDITION_LABELS)
    if samples["condition"].isna().any():
        raise SystemExit(f"unmapped condition values: {samples[config.CONDITION_COL].unique()}")
    samples.to_csv(OUT / "sample_table.csv", index=False)
    print(samples["condition"].value_counts().to_string(), "\n")

    # 3. Supplementary files and annotation (one request each).
    geo.download(geo.suppl_url(config.ACC, config.COUNTS_FILE), config.RAW / config.COUNTS_FILE)
    if optional:
        for name in config.OPTIONAL_FILES:
            geo.download(geo.suppl_url(config.ACC, name), config.RAW / name)
    print("files listed by GEO:", geo.series_supplementary_files(series))
    gtf = geo.download(config.ANNOTATION_URL, config.RAW / config.ANNOTATION_FILE)
    if not config.GENE_NAMES.exists():
        names = gene_names_from_gtf(gtf, config.GENE_NAMES)
        print(f"gene names: {len(names)} genes written to {config.GENE_NAMES.name}")

    # 4. Condition against each variable.
    lines = ["# Cross-tabulations: condition against sample variables", ""]
    for var in config.TECHNICAL_VARS + config.BIOLOGICAL_VARS:
        if var not in samples.columns:
            lines += [f"## {var}", "", "_not in the sample table_", ""]
            continue
        col = samples[var]
        if var == "char:age":
            col = pd.to_numeric(col, errors="coerce")
            col = pd.cut(col, bins=[0, 30, 45, 60, 120], right=False).astype(object)
        table = pd.crosstab(samples["condition"], col.fillna("missing"), margins=True)
        safe = var.replace("char:", "").replace(" ", "_").replace(".", "_")
        table.to_csv(OUT / f"crosstab_{safe}.csv")
        lines += [f"## {var}", "", table.to_markdown(), ""]
    (OUT / "crosstabs.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))

    # 5. Which variables have a value, for how many samples of each condition?
    chars = [c for c in samples.columns if c.startswith("char:")]
    present = samples[chars].notna().groupby(samples["condition"]).sum().T
    present["n_samples_total"] = len(samples)
    present.to_csv(OUT / "missing_values.csv")
    print("samples with a value, per variable and condition:\n", present.to_string(), "\n")

    # 6. Do the counts columns match the samples? GEO lists 43 GSMs; the CSV has its own names.
    with gzip.open(config.RAW / config.COUNTS_FILE, "rt") as fh:
        header = fh.readline().rstrip("\n").split(",")
    cols = [c.strip('"') for c in header[1:]]
    report = [f"columns in counts file: {len(cols)}", f"GSM records in GEO: {len(samples)}"]
    for key in ["title", "char:sample_id"]:
        ids = set(samples[key].astype(str)) if key in samples.columns else set()
        report += [
            f"[{key}] columns matched: {sum(c in ids for c in cols)}",
            f"[{key}] columns without a GEO record: " + ", ".join(c for c in cols if c not in ids),
            f"[{key}] GEO records without a column: " + ", ".join(sorted(ids - set(cols))),
        ]
    report.append("first columns: " + ", ".join(cols[:8]))
    (OUT / "counts_header_check.txt").write_text("\n".join(report), encoding="utf-8")
    print("\n".join(report))

    print("\nNext: fill worksheets/step1_audit.md and worksheets/fair_score.md")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--optional", action="store_true",
                        help="also download the FPKM, normalised and R object files")
    main(parser.parse_args().optional)
