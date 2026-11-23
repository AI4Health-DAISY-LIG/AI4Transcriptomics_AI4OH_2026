import pandas as pd
import numpy as np
from scipy import stats


def cramers_v(table: pd.DataFrame) -> float:
    """Calculate bias-corrected Cramér's V."""
    observed = table.to_numpy(dtype=float)

    if observed.size == 0 or observed.sum() == 0:
        return np.nan

    n = observed.sum()
    r, k = observed.shape

    if min(r, k) < 2:
        return np.nan

    chi2 = stats.chi2_contingency(observed, correction=False).statistic

    phi2 = chi2 / n
    phi2_corrected = max(0, phi2 - ((k - 1) * (r - 1)) / (n - 1))

    r_corrected = r - ((r - 1) ** 2) / (n - 1)
    k_corrected = k - ((k - 1) ** 2) / (n - 1)

    denominator = min(r_corrected - 1, k_corrected - 1)

    if denominator <= 0:
        return np.nan

    return float(np.sqrt(phi2_corrected / denominator))


def categorical_test(table: pd.DataFrame) -> dict:
    """
    Test association between condition and a categorical variable.

    Fisher exact for 2x2 tables with small expected frequencies.
    Chi-square otherwise.
    """
    observed = table.to_numpy(dtype=int)

    if min(observed.shape) < 2 or observed.sum() == 0:
        return {
            "test": "Not applicable",
            "statistic": np.nan,
            "p_value": np.nan,
            "cramers_v": np.nan,
            "min_expected": np.nan,
            "warning": "Insufficient categories",
        }

    chi2, chi_p, dof, expected = stats.chi2_contingency(observed, correction=False)

    min_expected = expected.min()

    if observed.shape == (2, 2) and min_expected < 5:
        odds_ratio, p_value = stats.fisher_exact(observed)

        return {
            "test": "Fisher exact",
            "statistic": odds_ratio,
            "p_value": p_value,
            "cramers_v": cramers_v(table),
            "min_expected": min_expected,
            "warning": "",
        }

    warning = ""

    if min_expected < 5:
        warning = (
            "Sparse expected frequencies: chi-square approximation "
            "may be unreliable."
        )

    return {
        "test": "Chi-square",
        "statistic": chi2,
        "p_value": chi_p,
        "cramers_v": cramers_v(table),
        "min_expected": min_expected,
        "warning": warning,
    }


def continuous_test(
    values: pd.Series,
    conditions: pd.Series,
) -> dict:
    """Compare a continuous variable between conditions."""

    data = pd.DataFrame(
        {
            "value": pd.to_numeric(values, errors="coerce"),
            "condition": conditions,
        }
    ).dropna()

    groups = [group["value"].to_numpy() for _, group in data.groupby("condition")]

    if len(groups) < 2 or any(len(g) == 0 for g in groups):
        return {
            "test": "Not applicable",
            "statistic": np.nan,
            "p_value": np.nan,
            "warning": "Insufficient observations",
        }

    if len(groups) == 2:
        result = stats.mannwhitneyu(
            groups[0],
            groups[1],
            alternative="two-sided",
        )

        return {
            "test": "Mann-Whitney U",
            "statistic": result.statistic,
            "p_value": result.pvalue,
            "warning": "",
        }

    result = stats.kruskal(*groups)

    return {
        "test": "Kruskal-Wallis",
        "statistic": result.statistic,
        "p_value": result.pvalue,
        "warning": "",
    }


def benjamini_hochberg(p_values):
    """Benjamini-Hochberg FDR correction, preserving NaNs."""

    p_values = np.asarray(p_values, dtype=float)
    adjusted = np.full(len(p_values), np.nan)

    valid = np.isfinite(p_values)

    if not valid.any():
        return adjusted

    p = p_values[valid]
    n = len(p)

    order = np.argsort(p)
    ranked = p[order]

    corrected = ranked * n / np.arange(1, n + 1)
    corrected = np.minimum.accumulate(corrected[::-1])[::-1]
    corrected = np.clip(corrected, 0, 1)

    restored = np.empty(n)
    restored[order] = corrected

    adjusted[valid] = restored

    return adjusted
