"""Describe published per-biopsy blood markers by site and histological inflammation.

Input: Tan9uy/AI4Transcriptomics_AI4OH_2026 at ce23bcf0320fcfe6f06f7bd838bc546458ef0c73,
practical/results/step2_qc/tissue_composition_ratios.csv. Percentages are fractions of
assigned gene counts multiplied by 100. The exact test conditions on site and the
number of FSHD biopsies with infiltrates at each site.
"""

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

import config


HB = "hemoglobin_pct"
NEUT = "neutrophil_pct"


def load_samples(marker_file: Path) -> pd.DataFrame:
    meta = pd.read_csv(config.RESULTS / "step1_audit" / "sample_table.csv").set_index(
        "title", verify_integrity=True
    )
    markers = pd.read_csv(marker_file).set_index("sample_id", verify_integrity=True)
    if set(meta.index) != set(markers.index):
        raise ValueError("The GEO titles and marker sample IDs must match.")
    markers = markers.loc[meta.index]
    if not meta["condition"].eq(markers["condition"]).all():
        raise ValueError("Condition labels differ between the GEO and marker tables.")
    if not np.isfinite(markers[[HB, NEUT]]).all().all():
        raise ValueError("The two marker percentages must be finite for every biopsy.")
    data = markers[["condition", HB, NEUT]].copy()
    data["site"] = data.index.str[:2].map(config.SITE_BY_PREFIX)
    if data["site"].isna().any():
        raise ValueError("An unsupported site prefix was found.")
    grades = pd.to_numeric(meta["char:inflam"], errors="raise")
    if not grades.dropna().isin([0, 1, 2, 3]).all():
        raise ValueError("Inflammation grade must be 0–3 or missing.")
    data["grade"] = grades
    return data


def exact_site_test(cases: pd.DataFrame) -> float:
    """Two-sided rank-sum test, permuting grades within each site at fixed sizes."""
    # Doubling midranks gives integer sums, including when marker values tie.
    total_observed = total_center = 0
    combined = Counter({0: 1})
    for _, site in cases.groupby("site", sort=True):
        ranks = (2 * site[NEUT].rank(method="average")).astype(int)
        positive = site["grade"].gt(0)
        k = int(positive.sum())
        total_observed += int(ranks[positive].sum())
        total_center += k * (len(site) + 1)
        within = Counter(sum(selection) for selection in combinations(ranks, k))
        next_combined = Counter()
        for left, left_count in combined.items():
            for right, right_count in within.items():
                next_combined[left + right] += left_count * right_count
        combined = next_combined
    extreme = sum(
        count for value, count in combined.items()
        if abs(value - total_center) >= abs(total_observed - total_center)
    )
    return extreme / sum(combined.values())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("marker_csv", type=Path, help="Published tissue_composition_ratios.csv")
    args = parser.parse_args()
    data = load_samples(args.marker_csv)
    track_b = data.drop(index=config.TRACK_B_EXCLUDE)
    # Healthy controls have no recorded histology grade and remain a separate group.
    inflammation = track_b.loc[track_b["condition"].eq("Control") | track_b["grade"].notna()]
    lines = [
        "GSE115650: per-biopsy blood markers",
        f"Marker CSV SHA-256: {sha256(args.marker_csv.read_bytes()).hexdigest()}",
        "Marker percentages use assigned gene counts as the denominator.",
        "",
        "cohort\tsite\tcondition\tn\tmedian_Hb_pct\tmedian_neutrophil_pct\tSpearman_rho_Hb_neut",
    ]
    for cohort, samples in [("Track_B", track_b), ("inflammation", inflammation)]:
        for (site, cond), group in samples.groupby(["site", "condition"], sort=True):
            rho = group[HB].corr(group[NEUT], method="spearman")
            lines.append(
                f"{cohort}\t{site}\t{cond}\t{len(group)}\t{group[HB].median():.8f}"
                f"\t{group[NEUT].median():.8f}\t{rho:.4f}"
            )

    lines.extend(["", "cohort\tsite\tinfiltrates\tn\tmedian_neutrophil_pct"])
    tests = []
    for cohort, samples in [("all_scored_FSHD", data), ("Track_B_scored_FSHD", track_b)]:
        cases = samples.loc[samples["condition"].eq("FSHD") & samples["grade"].notna()]
        for (site, positive), group in cases.groupby(["site", cases["grade"].gt(0)], sort=True):
            label = "present" if positive else "absent"
            lines.append(f"{cohort}\t{site}\t{label}\t{len(group)}\t{group[NEUT].median():.8f}")
        tests.append(f"{cohort}\t{len(cases)}\t{exact_site_test(cases):.6f}")
    lines.extend(["", "cohort\tn\texact_two_sided_site_stratified_p", *tests])

    out = config.RESULTS / "step3_de" / "track_b" / "blood_marker_site_check.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(out)


if __name__ == "__main__":
    main()
