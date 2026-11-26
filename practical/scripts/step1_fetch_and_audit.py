"""Step 1: fetch GSE115650 by code and audit it.

Outputs (results/step1_audit/):
  sample_table.csv            one row per GSM with GEO fields, `condition`, `char:infer_gender`, and `char:site`
  series_record.txt           the !Series_* fields as fetched
  counts_header_check.txt     do the column names of the counts file match the sample table?
  crosstab_<variable>.csv     condition against each variable in config.TECHNICAL_VARS + BIOLOGICAL_VARS + INFER_VARS
  crosstabs.md                descriptive cross-tabulations, row percentages, and statistical test report
  statistical_tests.csv       FDR-adjusted association and missingness test results across audited variables
  missing_values.csv          for each `char:*` variable, how many samples of each condition have a value
  (plus sex-verification artifacts written to results/step1_audit/ by `verify_sample_sex`)

Run:  uv run python scripts/step1_fetch_and_audit.py [--optional]
"""

from __future__ import annotations

import argparse
import gzip
import re

import pandas as pd
import numpy as np
import config
import geo
from infer_sex import verify_sample_sex
from step2_qc_exploration import load_counts
from helpers.statistics import (
    categorical_test,
    continuous_test,
    benjamini_hochberg,
)

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
            rows.append(
                {
                    "gene_id": attr["gene_id"],
                    "gene_id_noversion": attr["gene_id"].split(".")[0],
                    "gene_name": attr.get("gene_name", ""),
                    "gene_type": attr.get("gene_type", ""),
                }
            )
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
        "\n".join(f"{k} = {v}" for k, vals in series.items() for v in vals),
        encoding="utf-8",
    )

    # 2. Sample table with a clean condition column.
    samples["condition"] = samples[config.CONDITION_COL].map(config.CONDITION_LABELS)
    if samples["condition"].isna().any():
        raise SystemExit(
            f"unmapped condition values: {samples[config.CONDITION_COL].unique()}"
        )

    samples.to_csv(OUT / "sample_table.csv", index=False)

    # 3.  Supplementary files and annotation (one request each).
    geo.download(
        geo.suppl_url(config.ACC, config.COUNTS_FILE), config.RAW / config.COUNTS_FILE
    )
    if optional:
        for name in config.OPTIONAL_FILES:
            geo.download(geo.suppl_url(config.ACC, name), config.RAW / name)
    print("files listed by GEO:", geo.series_supplementary_files(series))
    gtf = geo.download(config.ANNOTATION_URL, config.RAW / config.ANNOTATION_FILE)
    if not config.GENE_NAMES.exists():
        names = gene_names_from_gtf(gtf, config.GENE_NAMES)
        print(f"gene names: {len(names)} genes written to {config.GENE_NAMES.name}")

    # 4. Do the counts columns match the samples? GEO lists 43 GSMs; the CSV has its own names.
    with gzip.open(config.RAW / config.COUNTS_FILE, "rt") as fh:
        header = fh.readline().rstrip("\n").split(",")
    cols = [c.strip('"') for c in header[1:]]
    report = [
        f"columns in counts file: {len(cols)}",
        f"GSM records in GEO: {len(samples)}",
    ]
    for key in ["title", "char:sample_id"]:
        ids = set(samples[key].astype(str)) if key in samples.columns else set()
        report += [
            f"[{key}] columns matched: {sum(c in ids for c in cols)}",
            f"[{key}] columns without a GEO record: "
            + ", ".join(c for c in cols if c not in ids),
            f"[{key}] GEO records without a column: "
            + ", ".join(sorted(ids - set(cols))),
        ]
    report.append("first columns: " + ", ".join(cols[:8]))
    (OUT / "counts_header_check.txt").write_text("\n".join(report), encoding="utf-8")
    print("\n".join(report))

    # Need to align the count matrix column with the samples id
    counts, samples = load_counts()

    gene_names = pd.read_csv(config.GENE_NAMES, index_col=0)

    # 5. Add infer col for gender and site

    # Add a col "char:infer_gender" with the right gender based on Xist and Y-chromosome
    samples, sex_audit = verify_sample_sex(
        counts=counts,
        meta=samples,
        gene_names=gene_names,
        gender_col="char:gender",
        output_dir=OUT,
    )

    print(samples["condition"].value_counts().to_string(), "\n")

    # Infer a site id based on the first 2 number of the title
    samples["char:site"] = [str(title)[:2] for title in samples["title"]]

    samples.to_csv(OUT / "sample_table.csv", index=True)

    # 6. Cross-tabulations and statistical tests.

    lines = [
        "# Cross-tabulations and statistical analysis",
        "",
        "## Statistical methodology",
        "",
        "- Categorical variables: Fisher's exact test for sparse 2x2 "
        "tables; otherwise Pearson's chi-square test.",
        "- Continuous variables: Mann-Whitney U for two conditions; "
        "Kruskal-Wallis for more than two.",
        "- Categorical effect size: bias-corrected Cramér's V.",
        "- Missingness: separately tested against condition.",
        "- Multiple testing: Benjamini-Hochberg FDR correction.",
        "- Statistical significance: adjusted p-value < 0.05.",
        "",
        "Missing observations are excluded from association tests "
        "but retained as a category in descriptive cross-tabs.",
        "",
    ]

    test_results = []

    variables = config.TECHNICAL_VARS + config.BIOLOGICAL_VARS + config.INFER_VARS

    for var in variables:

        lines += [f"## {var}", ""]

        if var not in samples.columns:
            lines += ["_Not in the sample table._", ""]
            continue

        raw_col = samples[var]
        condition = samples["condition"]

        # A. Descriptive cross-tabulation

        if var == "char:age":
            numeric_col = pd.to_numeric(raw_col, errors="coerce")

            display_col = pd.cut(
                numeric_col,
                bins=[0, 30, 45, 60, 120],
                right=False,
                labels=["<30", "30-44", "45-59", "60+"],
            ).astype(object)

            display_col = display_col.where(numeric_col.notna(), "missing")

        else:
            numeric_col = None
            display_col = raw_col.astype("object").fillna("missing")

        table = pd.crosstab(
            condition,
            display_col,
            margins=True,
        )

        safe = var.replace("char:", "").replace(" ", "_").replace(".", "_")

        table.to_csv(OUT / f"crosstab_{safe}.csv")

        lines += [
            "### Counts",
            "",
            table.to_markdown(),
            "",
        ]

        # Percentages within each condition
        percentage_table = (
            pd.crosstab(
                condition,
                display_col,
                normalize="index",
            )
            * 100
        )

        lines += [
            "### Row percentages (%)",
            "",
            percentage_table.round(1).to_markdown(),
            "",
        ]

        # B. Statistical association test

        if var == "char:age":

            result = continuous_test(raw_col, condition)

            result.update(
                {
                    "variable": var,
                    "analysis": "Condition vs continuous variable",
                    "effect_size": np.nan,
                    "effect_measure": "Not calculated",
                }
            )

            test_results.append(result)

        else:

            valid = pd.DataFrame(
                {
                    "condition": condition,
                    "value": raw_col,
                }
            ).dropna()

            observed_table = pd.crosstab(
                valid["condition"],
                valid["value"],
            )

            result = categorical_test(observed_table)

            result.update(
                {
                    "variable": var,
                    "analysis": "Condition vs observed categories",
                    "effect_size": result["cramers_v"],
                    "effect_measure": "Cramér's V",
                }
            )

            test_results.append(result)

        # C. Test whether missingness differs by condition

        missing_table = pd.crosstab(
            condition,
            raw_col.isna().map(
                {
                    True: "missing",
                    False: "present",
                }
            ),
        )

        if missing_table.shape == (2, 2):
            missing_result = categorical_test(missing_table)

            missing_result.update(
                {
                    "variable": var,
                    "analysis": "Condition vs missingness",
                    "effect_size": missing_result["cramers_v"],
                    "effect_measure": "Cramér's V",
                }
            )

            test_results.append(missing_result)

            lines += [
                "### Missingness by condition",
                "",
                missing_table.to_markdown(),
                "",
            ]

    # 7. Multiple testing

    results_df = pd.DataFrame(test_results)

    if not results_df.empty:

        results_df["p_adjusted"] = benjamini_hochberg(results_df["p_value"])

        results_df["significant"] = results_df["p_adjusted"] < 0.05

        results_df["p_value"] = results_df["p_value"].map(
            lambda x: f"{x:.4g}" if pd.notna(x) else "NA"
        )

        results_df["p_adjusted"] = results_df["p_adjusted"].map(
            lambda x: f"{x:.4g}" if pd.notna(x) else "NA"
        )

        results_df["statistic"] = results_df["statistic"].map(
            lambda x: f"{x:.4g}" if pd.notna(x) else "NA"
        )

        results_df["effect_size"] = results_df["effect_size"].map(
            lambda x: f"{x:.4f}" if pd.notna(x) else "NA"
        )

        results_df["min_expected"] = results_df["min_expected"].map(
            lambda x: f"{x:.2f}" if pd.notna(x) else "NA"
        )

        results_df.to_csv(
            OUT / "statistical_tests.csv",
            index=False,
        )

        lines += [
            "# Statistical test results",
            "",
            results_df.to_markdown(index=False),
            "",
            "## Interpretation",
            "",
            "Adjusted p-values below 0.05 indicate statistically "
            "significant associations after FDR correction.",
            "",
            "Cramér's V measures the strength of categorical "
            "associations. Statistical significance should be "
            "interpreted alongside effect size and sample size.",
            "",
        ]

    # Save the complete Markdown report.
    (OUT / "crosstabs.md").write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print("Cross-tabulations and statistical tests completed.")

    # 9. Which variables have a value, for how many samples of each condition?
    chars = [c for c in samples.columns if c.startswith("char:")]
    present = samples[chars].notna().groupby(samples["condition"]).sum().T
    present["n_samples_total"] = len(samples)
    present.to_csv(OUT / "missing_values.csv")
    print(
        "samples with a value, per variable and condition:\n", present.to_string(), "\n"
    )

    print("\nNext: fill worksheets/step1_audit.md and worksheets/fair_score.md")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--optional",
        action="store_true",
        help="also download the FPKM, normalised and R object files",
    )
    main(parser.parse_args().optional)
