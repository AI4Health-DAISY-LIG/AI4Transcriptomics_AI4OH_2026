import polars as pl

import config  # noqa: E402

OUT = config.RESULTS / "step2_qc"
AUDIT = config.RESULTS / "step1_audit"

# Map gene IDs to gene names for biomarkers of interest
bio = (
    pl.read_csv(config.PROCESSED / "gene_names.csv")
    .filter(pl.col("gene_name").is_in(config.BIOMARKERS))
    .select("gene_id", "gene_name")
)


def display_genes_and_sample_size_for_biomarkers():
    """Read counts_filtered.csv and check if the 4 biomarkers remain after filtering, plus report how many samples are kept."""

    # Load raw count matrix and rename the first column to 'gene_id'
    filter_counts = pl.read_csv(config.PROCESSED / "counts_filtered.csv").rename(
        {"": "gene_id"}
    )

    long = filter_counts.join(bio, on="gene_id")
    real_shape = long.shape
    print(f"Biomarker genes kept: {long["gene_name"].to_list()}")

    # -2 => remove the gene_id and gene_name column
    print(f"sample count kept: {real_shape[1]-2}")


if __name__ == "__main__":
    # Load raw count matrix and rename the first column to 'gene_id'
    raw = pl.read_csv(config.RAW / config.COUNTS_FILE).rename({"": "gene_id"})

    # Parse sample groups from the audit table, extracting group labels from source names
    groups = pl.read_csv(AUDIT / "sample_table.csv").select(
        "title", group=pl.col("source_name_ch1").str.replace(r"_(fe)?male$", "")
    )

    # Join raw counts with gene names based on gene_id, then melt into long format
    long = (
        raw.join(bio, on="gene_id")
        .unpivot(
            index=["gene_id", "gene_name"], variable_name="title", value_name="count"
        )
        .join(groups, on="title")
    )

    thresholds = [9, 10, 11, 12, 13, 14, 15, 20]

    # For each threshold, count how many samples have counts >= threshold, aggregated by gene and group
    out = (
        long.group_by("gene_name", "group")
        .agg(
            [(pl.col("count") >= t).sum().alias(f"n_ge_{t}") for t in thresholds]
            + [pl.len().alias("n_samples")]
        )
        .sort("gene_name", "group")
    )
    print(out)

    # For each threshold, count how many samples have counts >= threshold, aggregated by gene only
    total = long.group_by("gene_name").agg(
        [(pl.col("count") >= t).sum().alias(f"n_ge_{t}") for t in thresholds]
    )
    print(total)
    total.write_csv(OUT / "filtering_genes_effect.csv")
