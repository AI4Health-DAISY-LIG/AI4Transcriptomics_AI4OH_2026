
"""Download the GSE115650 data and check the sample information.

Run:
    uv run python scripts/step1_fetch_and_audit.py
    uv run python scripts/step1_fetch_and_audit.py --optional
"""

import argparse
import gzip
import re

import pandas as pd

import config
import geo


# Download the extra GEO files only when --optional is used.
parser = argparse.ArgumentParser()
parser.add_argument("--optional", action="store_true")
args = parser.parse_args()

OUT = config.RESULTS / "step1_audit"
OUT.mkdir(parents=True, exist_ok=True)
config.PROCESSED.mkdir(parents=True, exist_ok=True)


# 1. Download and read the GEO records
soft = geo.fetch_soft(config.ACC, config.RAW, targ="all")
series, samples = geo.parse_soft(soft)

with open(OUT / "series_record.txt", "w", encoding="utf-8") as f:
    for key, values in series.items():
        for value in values:
            f.write(f"{key} = {value}\n")


# 2. Add the condition of each sample
samples["condition"] = samples[config.CONDITION_COL].map(
    config.CONDITION_LABELS
)

if samples["condition"].isna().any():
    raise SystemExit(
        "Some condition values could not be identified: "
        + str(samples[config.CONDITION_COL].unique())
    )

samples.to_csv(OUT / "sample_table.csv", index=False)

print("Samples by condition:")
print(samples["condition"].value_counts())


# 3. Download the count file and the gene annotation
counts_path = config.RAW / config.COUNTS_FILE

geo.download(
    geo.suppl_url(config.ACC, config.COUNTS_FILE),
    counts_path
)

if args.optional:
    for name in config.OPTIONAL_FILES:
        geo.download(geo.suppl_url(config.ACC, name), config.RAW / name)

print("\nFiles available from GEO:")
print(geo.series_supplementary_files(series))

gtf = geo.download(
    config.ANNOTATION_URL,
    config.RAW / config.ANNOTATION_FILE
)

# Extract gene IDs, names and types from the GTF file.
if not config.GENE_NAMES.exists():
    genes = []

    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.split("\t", 9)

            if parts[2] != "gene":
                continue

            attributes = dict(
                re.findall(r'(\S+) "([^"]*)"', parts[8])
            )

            gene_id = attributes["gene_id"]

            genes.append({
                "gene_id": gene_id,
                "gene_id_noversion": gene_id.split(".")[0],
                "gene_name": attributes.get("gene_name", ""),
                "gene_type": attributes.get("gene_type", "")
            })

    gene_names = pd.DataFrame(genes)
    gene_names.to_csv(config.GENE_NAMES, index=False)

    print(
        "\nGenes saved:",
        len(gene_names),
        "->",
        config.GENE_NAMES.name
    )


# 4. Compare sample variables between conditions
markdown = ["# Condition against sample variables", ""]

for var in config.TECHNICAL_VARS + config.BIOLOGICAL_VARS:

    if var not in samples.columns:
        markdown.extend([
            f"## {var}",
            "",
            "Variable not found in the sample table.",
            ""
        ])
        continue

    values = samples[var]

    # Group ages into ranges to make the table easier to read.
    if var == "char:age":
        values = pd.to_numeric(values, errors="coerce")
        values = pd.cut(
            values,
            bins=[0, 30, 45, 60, 120],
            right=False
        ).astype(object)

    table = pd.crosstab(
        samples["condition"],
        values.fillna("missing"),
        margins=True
    )

    name = (
        var.replace("char:", "")
        .replace(" ", "_")
        .replace(".", "_")
    )

    table.to_csv(OUT / f"crosstab_{name}.csv")

    markdown.extend([
        f"## {var}",
        "",
        table.to_markdown(),
        ""
    ])

with open(OUT / "crosstabs.md", "w", encoding="utf-8") as f:
    f.write("\n".join(markdown))

print("\n".join(markdown))


# 5. Check missing values in the sample information
char_columns = [
    col for col in samples.columns
    if col.startswith("char:")
]

present = samples[char_columns].notna().groupby(
    samples["condition"]
).sum().T

present["n_samples_total"] = len(samples)

present.to_csv(OUT / "missing_values.csv")

print("\nAvailable values for each variable:")
print(present)


# 6. Check whether the count file columns match the GEO samples
with gzip.open(counts_path, "rt") as f:
    header = f.readline().rstrip("\n").split(",")

columns = [col.strip('"') for col in header[1:]]

lines = [
    f"Columns in counts file: {len(columns)}",
    f"GEO samples: {len(samples)}"
]

for key in ["title", "char:sample_id"]:

    if key not in samples.columns:
        lines.append(f"\n[{key}] not found in GEO metadata.")
        continue

    sample_ids = set(samples[key].astype(str))
    count_ids = set(columns)

    matched = sum(col in sample_ids for col in columns)

    lines.append(f"\n[{key}] columns matched: {matched}")

    missing_geo = [col for col in columns if col not in sample_ids]
    missing_counts = sorted(sample_ids - count_ids)

    lines.append(
        f"[{key}] columns without a GEO record: "
        + ", ".join(missing_geo)
    )

    lines.append(
        f"[{key}] GEO records without a column: "
        + ", ".join(missing_counts)
    )

lines.append("First columns: " + ", ".join(columns[:8]))

with open(OUT / "counts_header_check.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(lines))

print("\nCounts file check:")
print("\n".join(lines))

print("\nStep 1 finished.")
