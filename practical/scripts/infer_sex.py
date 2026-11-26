from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans


def verify_sample_sex(
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    gene_names: pd.DataFrame,
    gender_col: str = "char:gender",
    condition_col: str = "condition",
    control_label: str = "control",
    output_dir: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Infer sample sex from XIST and Y-chromosome gene expression.
    Compare inferred sex with metadata and highlight mismatches in the plot.
    """

    # Keep samples present in both count data and metadata
    samples = counts.columns.intersection(meta.index)

    if len(samples) == 0:
        raise ValueError("No sample IDs are shared between counts and metadata.")

    print(samples)

    count_data = counts.loc[:, samples].copy()
    sample_meta = meta.loc[samples].copy()

    # Map gene IDs to gene symbols
    gene_map = gene_names["gene_name"].dropna().to_dict()
    symbols = count_data.index.map(gene_map)

    # Fallback for Ensembl IDs with version numbers
    if symbols.isna().any() and "gene_id_noversion" in gene_names.columns:
        fallback_map = (
            gene_names.dropna(subset=["gene_id_noversion", "gene_name"])
            .drop_duplicates("gene_id_noversion")
            .set_index("gene_id_noversion")["gene_name"]
            .to_dict()
        )
        fallback = count_data.index.astype(str).str.split(".").str[0].map(fallback_map)
        symbols = symbols.fillna(fallback)

    count_data["gene_name"] = symbols
    count_data = count_data.dropna(subset=["gene_name"])
    count_data["gene_name"] = count_data["gene_name"].astype(str).str.strip()
    count_data = count_data[count_data["gene_name"] != ""]

    # Combine rows that map to the same gene symbol
    count_data = count_data.groupby("gene_name").sum(numeric_only=True)

    # Normalize counts to log2(CPM + 1)
    library_sizes = count_data.sum(axis=0).replace(0, np.nan)
    cpm = count_data.div(library_sizes, axis=1) * 1_000_000
    log_cpm = np.log2(cpm + 1).fillna(0)

    # Select sex-associated genes
    x_marker = "XIST"
    y_markers = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY"]
    available_y = [gene for gene in y_markers if gene in log_cpm.index]

    if x_marker not in log_cpm.index or not available_y:
        raise ValueError(
            f"Missing sex markers. XIST found: {x_marker in log_cpm.index}; "
            f"Y markers found: {available_y}"
        )

    results = pd.DataFrame(index=samples)
    results["XIST_logCPM"] = log_cpm.loc[x_marker]
    results["Y_mean_logCPM"] = log_cpm.loc[available_y].mean(axis=0)

    # Infer Female/Male using two clusters
    features = results[["XIST_logCPM", "Y_mean_logCPM"]]

    if len(results) < 2:
        raise ValueError("At least two samples are required for clustering.")

    model = KMeans(n_clusters=2, random_state=42, n_init=10)
    results["cluster"] = model.fit_predict(features)

    # Higher XIST and lower Y expression indicates the female-like cluster
    cluster_score = (
        results.groupby("cluster")["XIST_logCPM"].mean()
        - results.groupby("cluster")["Y_mean_logCPM"].mean()
    )
    female_cluster = cluster_score.idxmax()

    results["char:infer_gender"] = np.where(
        results["cluster"] == female_cluster, "Female", "Male"
    )
    results.drop(columns="cluster", inplace=True)

    # Read and standardize reported gender
    if gender_col in sample_meta.columns:
        reported = (
            sample_meta[gender_col]
            .astype("string")
            .str.strip()
            .str.lower()
            .map({"female": "Female", "male": "Male"})
        )
        results["reported_gender"] = reported
    else:
        print(f"Warning: '{gender_col}' was not found in metadata.")
        results["reported_gender"] = np.nan

    results["has_reported_gender"] = results["reported_gender"].notna()
    results["sex_mismatch"] = results["has_reported_gender"] & (
        results["char:infer_gender"] != results["reported_gender"]
    )

    # Keep condition information
    if condition_col in sample_meta.columns:
        results["condition"] = sample_meta[condition_col]
        results["is_control"] = (
            results["condition"]
            .fillna("")
            .astype(str)
            .str.strip()
            .str.lower()
            .str.contains(control_label.lower(), regex=False)
        )
    else:
        results["condition"] = np.nan
        results["is_control"] = False

    # Add inferred gender to metadata
    meta_updated = meta.copy()
    meta_updated["char:infer_gender"] = results["char:infer_gender"].reindex(meta.index)

    # Print summary
    mismatches = results[results["sex_mismatch"]]
    unknown = results[~results["has_reported_gender"]]

    print(f"Checked {len(results)} samples.")
    print("Inferred gender:", results["char:infer_gender"].value_counts().to_dict())
    print(f"Samples without reported gender: {len(unknown)}")
    print(f"Gender mismatches: {len(mismatches)}")

    if not mismatches.empty:
        print("\nGender mismatches:")
        print(
            mismatches[
                ["reported_gender", "char:infer_gender", "XIST_logCPM", "Y_mean_logCPM"]
            ]
        )
    else:
        print("No gender mismatches found.")

    # Plot and save the image
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        results.to_csv(output_dir / "sex_verification_results.csv")

        fig, ax = plt.subplots(figsize=(9, 6))

        # Samples with known gender and matching metadata
        for gender, color in [("Female", "orange"), ("Male", "green")]:
            subset = results[
                (results["char:infer_gender"] == gender)
                & (~results["sex_mismatch"])
                & results["has_reported_gender"]
            ]

            ax.scatter(
                subset["XIST_logCPM"],
                subset["Y_mean_logCPM"],
                color=color,
                label=f"{gender} (match)",
                alpha=0.75,
                s=55,
            )

        # Samples with missing reported gender
        unknown = results[~results["has_reported_gender"]]
        ax.scatter(
            unknown["XIST_logCPM"],
            unknown["Y_mean_logCPM"],
            facecolors="none",
            edgecolors="gray",
            s=75,
            linewidths=1.5,
            label="Gender not reported",
        )

        # Gender mismatches: red X markers
        if not mismatches.empty:
            ax.scatter(
                mismatches["XIST_logCPM"],
                mismatches["Y_mean_logCPM"],
                color="red",
                edgecolors="black",
                marker="X",
                s=150,
                label="GENDER MISMATCH",
                zorder=5,
            )

            # Label mismatched samples
            for sample_id, row in mismatches.iterrows():
                ax.annotate(
                    str(sample_id),
                    (row["XIST_logCPM"], row["Y_mean_logCPM"]),
                    xytext=(6, 6),
                    textcoords="offset points",
                    color="red",
                    fontsize=8,
                    fontweight="bold",
                )

        ax.set_xlabel("XIST expression [log2(CPM + 1)]")
        ax.set_ylabel(
            f"Mean Y-gene expression ({', '.join(available_y)}) " "[log2(CPM + 1)]"
        )
        ax.set_title("Inferred Gender vs Reported Gender")
        ax.legend()
        ax.grid(True, linestyle="--", alpha=0.4)
        fig.tight_layout()

        fig.savefig(output_dir / "sex_verification_qc.png", dpi=150)
        plt.close(fig)

    return meta_updated, results
