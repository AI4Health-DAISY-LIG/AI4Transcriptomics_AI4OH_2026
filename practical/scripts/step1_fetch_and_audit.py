
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

import gzip
import re

import pandas as pd

import config
import geo

DOWNLOAD_OPTIONAL = False  

OUT = config.RESULTS / "step1_audit"
OUT.mkdir(parents=True, exist_ok=True)
config.PROCESSED.mkdir(parents=True, exist_ok=True)

# 1. series and sample records
soft = geo.fetch_soft(config.ACC, config.RAW, targ="all")
series, samples = geo.parse_soft(soft)
with open(OUT / "series_record.txt", "w", encoding="utf-8") as f:
    for key, values in series.items():
        for value in values:
            f.write(f"{key} = {value}\n")

# 2. table
samples["condition"] = samples[config.CONDITION_COL].map(config.CONDITION_LABELS)
if samples["condition"].isna().any():
    raise SystemExit("unmapped condition values: " + str(samples[config.CONDITION_COL].unique()))
samples.to_csv(OUT / "sample_table.csv", index=False)
print(samples["condition"].value_counts(), "\n")

# 3. files and annotation
geo.download(geo.suppl_url(config.ACC, config.COUNTS_FILE), config.RAW / config.COUNTS_FILE)
if DOWNLOAD_OPTIONAL:
    for name in config.OPTIONAL_FILES:
        geo.download(geo.suppl_url(config.ACC, name), config.RAW / name)
print("files listed by GEO:", geo.series_supplementary_files(series))

gtf = geo.download(config.ANNOTATION_URL, config.RAW / config.ANNOTATION_FILE)
if not config.GENE_NAMES.exists():
    rows = []
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 9)
            if parts[2] != "gene":
                continue
            attrs = dict(re.findall(r'(\S+) "([^"]*)"', parts[8]))
            rows.append({
                "gene_id": attrs["gene_id"],
                "gene_id_noversion": attrs["gene_id"].split(".")[0],
                "gene_name": attrs.get("gene_name", ""),
                "gene_type": attrs.get("gene_type", ""),
            })
    names = pd.DataFrame(rows)
    names.to_csv(config.GENE_NAMES, index=False)
    print(len(names), "genes written to", config.GENE_NAMES.name)

# 4. condition/variable
md = ["# Condition against sample variables", ""]
for var in config.TECHNICAL_VARS + config.BIOLOGICAL_VARS:
    if var not in samples.columns:
        md += [f"## {var}", "", "not in the sample table", ""]
        continue
    values = samples[var]
    if var == "char:age":
        values = pd.to_numeric(values, errors="coerce")
        values = pd.cut(values, bins=[0, 30, 45, 60, 120], right=False).astype(object)
    table = pd.crosstab(samples["condition"], values.fillna("missing"), margins=True)
    name = var.replace("char:", "").replace(" ", "_").replace(".", "_")
    table.to_csv(OUT / f"crosstab_{name}.csv")
    md += [f"## {var}", "", table.to_markdown(), ""]
with open(OUT / "crosstabs.md", "w", encoding="utf-8") as f:
    f.write("\n".join(md))
print("\n".join(md))

# 5. how many samples of each condition have a value for each variable
chars = [c for c in samples.columns if c.startswith("char:")]
present = samples[chars].notna().groupby(samples["condition"]).sum().T
present["n_samples_total"] = len(samples)
present.to_csv(OUT / "missing_values.csv")
print(present, "\n")

# 6. counts columns match samples?
with gzip.open(config.RAW / config.COUNTS_FILE, "rt") as f:
    header = f.readline().rstrip("\n").split(",")
cols = [c.strip('"') for c in header[1:]]
lines = [f"columns in counts file: {len(cols)}", f"GSM records in GEO: {len(samples)}"]
for key in ["title", "char:sample_id"]:
    ids = set(samples[key].astype(str)) if key in samples.columns else set()
    lines.append(f"[{key}] columns matched: {sum(c in ids for c in cols)}")
    lines.append(f"[{key}] columns without a GEO record: " + ", ".join(c for c in cols if c not in ids))
    lines.append(f"[{key}] GEO records without a column: " + ", ".join(sorted(ids - set(cols))))
lines.append("first columns: " + ", ".join(cols[:8]))
with open(OUT / "counts_header_check.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
print("\n".join(lines))
