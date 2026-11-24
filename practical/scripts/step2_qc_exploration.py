"""Step 2: quality control and exploration of the raw counts.

Outputs (results/step2_qc/):
  library_sizes.png                 total counts per sample, by batch, colour and shape = condition
  library_sizes.csv                 total counts per sample with condition and batch
  gene_filter.txt                   rule used, smallest gene total, genes before and after,
                                    and whether the four biomarkers pass
  pca_variance.csv                  variance explained by PC1 to PC10 (%), authors' PCA method
                                    on the rlog of the GEO counts
  pca_scores.csv                    PC1 to PC10 per sample, signs aligned to the Figure 3 replay
  pca_by_<variable>.png             PC1 vs PC2, colour = variable, shape = condition
  pca_logcpm_comparison.csv         PC1 and PC2 variance of log2(CPM + 1) PCAs, for comparison
  pca_associations.csv              association of PC1 to PC3 with each plotted variable
  paper_group_by_batch.csv          number of samples per paper group and batch
  sex_check.csv                     sex inferred from XIST and 6 Y genes vs GEO label
  sex_check.png                     XIST vs Y genes, colour = GEO label, shape = condition
  sample_similarity.csv             each sample's most similar other sample (duplicate check)
  pca_genes.txt                     the 500 genes of the PCA, reused by the sensitivity run
  pca_top_genes.csv                 the 20 genes with the largest weight on PC1 to PC3, with symbols
  data/processed/counts_filtered.csv   filtered raw counts (genes x samples), input of step 3
  data/processed/metadata.csv          one row per sample column of the counts, with condition

Inputs written by the R scripts, which must run first (worksheets/step2_qc.md):
  data/processed/rlog_geo.csv.gz             rlog of the GEO counts (r/rlog_from_geo.R)
  results/step2_qc/rlog_geo_pca_variance.csv R's PCA variance, used as a cross-check
  results/step2_qc/pca_replay_scores.csv     Figure 3 replay scores, for axis orientation
  results/step2_qc/pca_replay_genes.txt      the 500 genes of Figure 3 (r/pca_replay.R)
  results/step2_qc/mapping_rates.csv         mapping rate per sample (r/pca_replay.R)

Run, from practical/:
  pixi run --manifest-path r/pixi.toml Rscript r/pca_replay.R
  pixi run --manifest-path r/pixi.toml Rscript r/rlog_from_geo.R
  uv run python scripts/step2_qc_exploration.py

Sensitivity run (after the main run), without the low-mapping biopsies of
config.SENSITIVITY_EXCLUDE: same rlog values, same 500 genes, PCA refitted on the remaining
samples. Writes to results/step2_qc/sensitivity_40/ and leaves data/processed/ untouched:
  uv run python scripts/step2_qc_exploration.py --sensitivity

Figures only: redraws every pca_by_<variable>.png of both runs from the saved scores and
variances (pca_scores.csv, pca_variance.csv), without the counts; nothing else changes:
  uv run python scripts/step2_qc_exploration.py --figures

Figure conventions: shape = condition (FSHD orange circle, control blue triangle); a variable
with a single value (instrument, GEO tissue label) is stated in the subtitle and the points
are coloured by condition; continuous variables share the colour limits of the 43-biopsy
cohort in both runs; ordered ratings (STIR, T1, active, inflammation) use an ordered palette
with integer levels; the mapping rate is shown in five ordered bins (BINNED); missing values
are grey.
"""

from __future__ import annotations

import argparse

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import config  # noqa: E402

OUT = config.RESULTS / "step2_qc"
AUDIT = config.RESULTS / "step1_audit"

# Samples labelled on the figures: low mapping rate (32-0008, 01-0037, 32-0009) and the
# paper's two PCA outliers (01-0037, 32-0016). Only those present in the data are labelled.
WATCH = ["32-0008", "01-0037", "32-0009", "32-0016"]
LABEL_LEFT = {"32-0009"}   # labels placed left of the point (the right frame crosses them)
MARKERS = {"FSHD": "o", "Control": "^"}
CONDITION_COLOURS = {"FSHD": "tab:orange", "Control": "tab:blue"}
ORDINAL = {"stir", "t1", "active", "inflam"}       # ordered clinical ratings
# Mapping rate: most libraries lie between 94.85% and 97.20%, so a continuous scale over the
# full range (67.73-97.20%) shows them in nearly one colour. Binned instead: below 90% (the
# three libraries of config.SENSITIVITY_EXCLUDE), then the quartiles of the libraries at 90% or
# more, with edges fixed on the full cohort.
BINNED = {"mapping_rate_pct": 90.0}
CONSTANT_NOTE = {"instrument_model": "All {n} samples were sequenced on an {value}",
                 "tissue": "GEO records every sample as '{value}'"}
N_TOP = 500  # genes in the PCA, as DESeq2's plotPCA default used for Figure 3


def load_counts() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Raw counts (genes x samples) and metadata aligned on the sample columns."""
    counts = pd.read_csv(config.RAW / config.COUNTS_FILE, index_col=0)
    samples = pd.read_csv(AUDIT / "sample_table.csv")
    # Join on GEO `title`: it matches all 43 count columns, including 01-0022-1, whose
    # char:sample_id is 01-0022 (results/step1_audit/counts_header_check.txt,
    # worksheets/step1_audit.md). char:sample_id stays as a column, as the key to the
    # paper's Table S1.
    meta = samples.set_index("title").reindex(counts.columns)
    missing = meta["condition"].isna()
    if missing.any():
        raise ValueError(f"count columns without a GEO record: {list(counts.columns[missing])}")
    return counts, meta


def pca_top(mat: pd.DataFrame, n_top: int, genes=None):
    """PCA as DESeq2's plotPCA: the n_top genes with the highest variance (or the given genes),
    centred, not scaled. Returns scores (samples x PC1..PC10), variance explained (%) and genes."""
    if genes is None:
        genes = mat.var(axis=1).sort_values(ascending=False).index[:n_top]
    x = mat.loc[list(genes)].T
    x = x - x.mean(axis=0)
    pca = PCA(n_components=10).fit(x)
    scores = pd.DataFrame(pca.transform(x), index=x.index,
                          columns=[f"PC{i + 1}" for i in range(10)])
    return scores, 100 * pca.explained_variance_ratio_, list(genes)


def scatter_by(pcs: pd.DataFrame, values: pd.Series, conditions: pd.Series, name: str,
               pct, title: str, path, limits: dict | None = None) -> None:
    """PC1 against PC2: shape = condition; colour = the variable (see the module docstring for
    constant, continuous and ordered variables); watched samples labelled. `limits` gives the
    colour range of a continuous variable, or the levels of an ordered one, taken from the full
    cohort so that the main and sensitivity figures share them."""
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(7.5, 5))
    numeric = pd.to_numeric(values, errors="coerce")
    present = values.dropna().astype(str).unique()
    constant = len(present) == 1
    binned = name in BINNED and name in (limits or {})
    ordinal = (name in ORDINAL or binned) and bool(numeric.notna().any())
    continuous = not ordinal and bool(numeric.notna().any()) and numeric.nunique() > 8
    style = dict(s=45, edgecolors="black", linewidths=0.5)
    handles = []

    def condition_handles():
        return [Line2D([0], [0], marker=MARKERS[c], color="w", markerfacecolor=col,
                       markeredgecolor="black", markersize=7, label=c)
                for c, col in CONDITION_COLOURS.items()]

    if name == "condition" or constant:
        for condition, marker in MARKERS.items():
            idx = conditions == condition
            ax.scatter(pcs.loc[idx, "PC1"], pcs.loc[idx, "PC2"],
                       color=CONDITION_COLOURS[condition], marker=marker, **style)
        handles = condition_handles()
        legend_title = "condition"
        if constant and name in CONSTANT_NOTE:
            title = title.replace(f"colour = {name}", "colour = condition") + "\n" + \
                CONSTANT_NOTE[name].format(n=len(pcs), value=present[0])
    elif continuous:
        vmin, vmax = (limits or {}).get(name, (numeric.min(), numeric.max()))
        sc = None
        for condition, marker in MARKERS.items():
            idx = conditions == condition
            na, ok = idx & numeric.isna(), idx & numeric.notna()
            ax.scatter(pcs.loc[na, "PC1"], pcs.loc[na, "PC2"], c="lightgrey", marker=marker, **style)
            sc = ax.scatter(pcs.loc[ok, "PC1"], pcs.loc[ok, "PC2"], c=numeric[ok], cmap="viridis",
                            vmin=vmin, vmax=vmax, marker=marker, **style)
        fig.colorbar(sc, ax=ax, label=name)
        if numeric.isna().any():
            handles.append(Line2D([0], [0], marker="s", color="w", markerfacecolor="lightgrey",
                                  markeredgecolor="black", label="no value"))
        handles += [Line2D([0], [0], marker=m, color="w", markerfacecolor="white",
                           markeredgecolor="black", label=c) for c, m in MARKERS.items()]
        legend_title = name
    else:
        if binned:
            edges, labels = limits[name]
            cmap = plt.get_cmap("viridis", len(labels))
            colours = {lab: cmap(i) for i, lab in enumerate(labels)}
            levels = pd.cut(numeric, edges, labels=labels, include_lowest=True) \
                .astype(object).where(numeric.notna(), "no value").astype(str)
            order = list(colours)
        elif ordinal:
            levels_all = (limits or {}).get(name, sorted(numeric.dropna().unique()))
            cmap = plt.get_cmap("viridis", len(levels_all))
            colours = {f"{int(v)}": cmap(i) for i, v in enumerate(levels_all)}
            levels = numeric.map(lambda v: f"{int(v)}" if pd.notna(v) else "no value")
            order = list(colours)
        else:
            levels = values.astype(object).where(values.notna(), "no value").astype(str)
            order = sorted(set(levels) - {"no value"})
            colours = dict(zip(order, plt.cm.tab10.colors))
        colours["no value"] = "lightgrey"
        # Colours stay fixed by the full cohort; the legend lists only the levels present.
        order = [lv for lv in order if (levels == lv).any()]
        if (levels == "no value").any():
            order = order + ["no value"]
        for condition, marker in MARKERS.items():
            for level in order:
                sel = (conditions == condition) & (levels == level)
                ax.scatter(pcs.loc[sel, "PC1"], pcs.loc[sel, "PC2"], color=colours[level],
                           marker=marker, **style)
        handles += [Line2D([0], [0], marker="s", color="w", markerfacecolor=colours[lv],
                           markeredgecolor="black", label=lv) for lv in order]
        handles += [Line2D([0], [0], marker=m, color="w", markerfacecolor="white",
                           markeredgecolor="black", label=c) for c, m in MARKERS.items()]
        legend_title = name
    if continuous:   # the colour bar sits on the right: shape legend inside the axes
        ax.legend(handles=handles, fontsize=7, title=legend_title, title_fontsize=7, loc="best")
    else:            # outside the axes, so that it never hides a point
        ax.legend(handles=handles, fontsize=7, title=legend_title, title_fontsize=7,
                  loc="upper left", bbox_to_anchor=(1.01, 1), borderaxespad=0)
    for sample in WATCH:
        if sample in pcs.index:
            left = sample in LABEL_LEFT
            ax.annotate(sample, (pcs.loc[sample, "PC1"], pcs.loc[sample, "PC2"]), fontsize=6,
                        xytext=(-4 if left else 4, 2), textcoords="offset points",
                        ha="right" if left else "left")
    ax.set_xlabel(f"PC1 ({pct[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({pct[1]:.1f}%)")
    ax.set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def pca_variables(meta: pd.DataFrame, lib: pd.Series) -> dict[str, pd.Series]:
    """The variables the PCA figures are coloured by."""
    # Paper groups as in Figure 3 (controls are in group 1 by construction of the authors'
    # clustering).
    rank = pd.to_numeric(meta[config.GROUP_COL], errors="coerce")
    paper_group = rank.map(lambda r: f"Group{int(r)}_FSHD" if pd.notna(r) else None)
    paper_group[meta["condition"] == "Control"] = "Group1_Control"
    variables = {
        "condition": meta["condition"],
        "paper_group": paper_group,
        "batch": meta["char:batch"],
        "instrument_model": meta["instrument_model"],
        "tissue": meta["char:tissue"],
        "library_size_millions": lib.reindex(meta.index) / 1e6,
    }
    rates_file = OUT / "mapping_rates.csv"
    if rates_file.exists():
        rates = pd.read_csv(rates_file, index_col="geo_title")["mapping_rate_pct"]
        variables["mapping_rate_pct"] = rates.reindex(meta.index)
    for var in config.BIOLOGICAL_VARS + config.PATHOLOGY_VARS:
        if var in meta.columns and var != config.GROUP_COL:
            variables[var.replace("char:", "").replace(".", "_")] = meta[var]
    return variables


def colour_limits(variables: dict[str, pd.Series]) -> dict:
    """Full-cohort colour range of each numeric variable and levels of each ordered rating."""
    limits = {}
    for name, values in variables.items():
        numeric = pd.to_numeric(values, errors="coerce").dropna()
        if numeric.empty:
            continue
        if name in BINNED:
            high = numeric[numeric >= BINNED[name]]
            # The low bin is strictly below the QC threshold (a library at exactly 90% is kept
            # by the QC rule, so it goes in the next bin): its upper edge sits just below 90.
            edges = [float(numeric.min()), float(np.nextafter(BINNED[name], -np.inf))] + \
                [round(float(q), 2) for q in high.quantile([0.25, 0.5, 0.75])] + \
                [float(numeric.max())]
            bins = pd.cut(numeric, edges, include_lowest=True)
            # Each label gives the observed range of its bin in the full cohort.
            labels = [f"{numeric[bins == b].min():.2f}-{numeric[bins == b].max():.2f}%"
                      for b in bins.cat.categories]
            labels[0] = f"below {BINNED[name]:.0f}% ({labels[0]})"
            limits[name] = (edges, labels)
        elif name in ORDINAL:
            limits[name] = sorted(numeric.unique())
        else:
            limits[name] = (float(numeric.min()), float(numeric.max()))
    return limits


def full_cohort_limits() -> dict:
    meta = load_sample_table()
    lib = pd.read_csv(OUT / "library_sizes.csv", index_col=0)["lib_size"]
    return colour_limits(pca_variables(meta, lib))


def load_sample_table() -> pd.DataFrame:
    """GEO sample table indexed by title (the count-matrix column names), 43 biopsies."""
    meta = pd.read_csv(AUDIT / "sample_table.csv").set_index("title")
    return meta[meta["condition"].notna()]


def pca_figures(pcs: pd.DataFrame, pct, meta: pd.DataFrame, lib: pd.Series, out, label: str,
                limits: dict) -> dict[str, pd.Series]:
    variables = pca_variables(meta, lib)
    for name, values in variables.items():
        scatter_by(pcs, values, meta["condition"], name, pct,
                   f"PCA, {label}, top {N_TOP} genes: colour = {name}",
                   out / f"pca_by_{name}.png", limits)
    return variables


def describe_pca(pcs: pd.DataFrame, pct, meta: pd.DataFrame, lib: pd.Series, out, label: str,
                 limits: dict) -> None:
    """Figures (colour = variable, shape = condition), paper group by batch, and the association
    of PC1 to PC3 with each plotted variable, written to `out`."""
    variables = pca_figures(pcs, pct, meta, lib, out, label, limits)
    # Paper groups against batch: which groups each sequencing run holds (worksheets/step2_qc.md).
    pd.crosstab(variables["paper_group"], meta["char:batch"]).to_csv(out / "paper_group_by_batch.csv")

    # Association of PC1 to PC3 with each plotted variable, to name what separates the
    # samples: Spearman rho for numeric variables, and for categorical ones the share of the
    # PC's variance explained by the levels (R^2). Descriptive only: batch and group overlap.
    rows = []
    for name, values in variables.items():
        ok = values.notna()
        numeric = pd.to_numeric(values[ok], errors="coerce")
        for pc in ["PC1", "PC2", "PC3"]:
            y = pcs.loc[ok[ok].index, pc]
            if numeric.notna().all():
                rho, p = spearmanr(y, numeric)
                rows.append({"pc": pc, "variable": name, "measure": "spearman_rho",
                             "value": round(rho, 3), "p_value": float(f"{p:.3g}"),
                             "n": int(ok.sum())})
            else:
                within = ((y - y.groupby(values[ok]).transform("mean")) ** 2).sum()
                r2 = 1 - within / ((y - y.mean()) ** 2).sum()
                rows.append({"pc": pc, "variable": name, "measure": "r2_between_levels",
                             "value": round(r2, 3), "p_value": None, "n": int(ok.sum())})
    pd.DataFrame(rows).to_csv(out / "pca_associations.csv", index=False)


def write_pca(pcs: pd.DataFrame, pct, out) -> None:
    names = [f"PC{i + 1}" for i in range(10)]
    pd.DataFrame({"component": names, "variance_explained_pct": pct.round(4)}).to_csv(
        out / "pca_variance.csv", index=False)
    pcs.round(6).to_csv(out / "pca_scores.csv")


def orient(pcs: pd.DataFrame, reference: pd.DataFrame) -> pd.DataFrame:
    """Flip PC1 and PC2 when they point opposite to the reference scores (shared samples)."""
    shared = pcs.index.intersection(reference.index)
    for pc in ["PC1", "PC2"]:
        if np.corrcoef(pcs.loc[shared, pc], reference.loc[shared, pc])[0, 1] < 0:
            pcs[pc] = -pcs[pc]
    return pcs

def write_top_genes(mat: pd.DataFrame, genes, pcs: pd.DataFrame, out, n: int = 20) -> None:
    """The n genes with the largest weight on PC1 to PC3, with symbols. A gene's weight is
    proportional to its covariance with the PC scores; scaled to unit length, it equals the
    PCA loading, with the sign of the oriented axis."""
    x = mat.loc[list(genes), pcs.index].T
    x = x - x.mean(axis=0)
    symbols = (pd.read_csv(config.GENE_NAMES).drop_duplicates("gene_id")
               .set_index("gene_id")["gene_name"])
    rows = []
    for pc in ["PC1", "PC2", "PC3"]:
        w = x.T @ pcs[pc]
        w = w / np.sqrt((w ** 2).sum())
        for rank, gene in enumerate(w.abs().nlargest(n).index, start=1):
            rows.append({"pc": pc, "rank": rank, "gene_id": gene,
                         "gene_name": symbols.get(gene, ""), "loading": round(float(w[gene]), 4)})
    pd.DataFrame(rows).to_csv(out / "pca_top_genes.csv", index=False)

Y_GENES = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y"]


def sex_check(counts: pd.DataFrame, meta: pd.DataFrame, out) -> pd.DataFrame:
    """Sex inferred from expression, compared with the GEO label (FSHD only; controls have none).

    XIST is expressed from the inactive X, so only in women; the six Y-chromosome genes only
    in men. Rule, fixed before looking at the data: male when the mean log2(CPM + 1) of the
    Y genes is above XIST's log2(CPM + 1). `margin` is that difference; values near 0 would
    be ambiguous.
    """
    symbol_to_id = (pd.read_csv(config.GENE_NAMES).drop_duplicates("gene_name")
                    .set_index("gene_name")["gene_id"])
    logcpm = np.log2(counts / counts.sum(axis=0) * 1e6 + 1)
    table = pd.DataFrame({"condition": meta["condition"],
                          "geo_gender": meta["char:gender"],
                          "XIST": logcpm.loc[symbol_to_id["XIST"]]})
    for gene in Y_GENES:
        table[gene] = logcpm.loc[symbol_to_id[gene]]
    table["Y_mean"] = table[Y_GENES].mean(axis=1)
    table["margin"] = table["Y_mean"] - table["XIST"]
    table["inferred_sex"] = np.where(table["margin"] > 0, "male", "female")
    labelled = table["geo_gender"].notna()
    table["agrees_with_geo"] = (table["inferred_sex"] == table["geo_gender"]).where(labelled)
    table.round(3).sort_values(["condition", "inferred_sex", "margin"]).to_csv(
        out / "sex_check.csv")

    fig, ax = plt.subplots(figsize=(5, 4))
    colours = {"female": "tab:orange", "male": "tab:green"}
    for condition, marker in MARKERS.items():
        for label, colour in [*colours.items(), ("not given", "lightgrey")]:
            gender = table["geo_gender"].fillna("not given")
            idx = (table["condition"] == condition) & (gender == label)
            if idx.any():
                ax.scatter(table.loc[idx, "XIST"], table.loc[idx, "Y_mean"], marker=marker,
                           s=45, color=colour, edgecolors="black", linewidths=0.5,
                           label=f"{condition}, GEO: {label}")
    # Biopsies whose GEO label disagrees: redrawn on top in their GEO colour, ringed in red
    # and labelled with an arrow, so they stay visible inside the other sex's cluster.
    wrong = table["agrees_with_geo"] == False  # noqa: E712
    for sample in table.index[wrong]:
        x, y = table.loc[sample, "XIST"], table.loc[sample, "Y_mean"]
        ax.scatter(x, y, marker=MARKERS[table.loc[sample, "condition"]], s=60, zorder=4,
                   color=colours.get(table.loc[sample, "geo_gender"], "lightgrey"),
                   edgecolors="black", linewidths=0.5)
        ax.annotate(f"{sample} (GEO: {table.loc[sample, 'geo_gender']})", (x, y),
                    xytext=(60, -60), textcoords="offset points", fontsize=8, color="red",
                    arrowprops={"arrowstyle": "->", "color": "red"})
    ax.scatter(table.loc[wrong, "XIST"], table.loc[wrong, "Y_mean"], s=160, facecolors="none",
               edgecolors="red", linewidths=1.2, zorder=5, label="label disagrees")
    ax.set_xlabel("XIST, log2(CPM + 1)")
    ax.set_ylabel("mean of 6 Y genes, log2(CPM + 1)")
    ax.set_title("Sex from expression vs GEO label")
    ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(out / "sex_check.png", dpi=150)
    plt.close(fig)

    fshd = table[labelled]
    disagree = list(fshd.index[fshd["agrees_with_geo"] == False])  # noqa: E712
    controls = table.loc[table["condition"] == "Control", "inferred_sex"].value_counts()
    print(f"sex check: {int(fshd['agrees_with_geo'].sum())}/{len(fshd)} GEO labels agree "
          f"with expression; disagree: {', '.join(disagree) or 'none'}; smallest |margin| "
          f"{table['margin'].abs().min():.2f}; controls inferred: "
          f"{controls.get('female', 0)} female, {controls.get('male', 0)} male")
    return table


def similarity_check(counts: pd.DataFrame, out) -> None:
    """Looks for duplicated libraries: identical count columns, and each sample's most similar
    other sample (Pearson r of log2(CPM + 1) over genes with >= 10 counts in at least half the
    samples). A library sequenced twice or put under two names would stand out above the
    pairs of different people."""
    identical = int(counts.T.duplicated().sum())
    expressed = (counts >= 10).sum(axis=1) >= counts.shape[1] / 2
    logcpm = np.log2(counts / counts.sum(axis=0) * 1e6 + 1).loc[expressed]
    r = logcpm.corr()
    r = r.where(~np.eye(len(r), dtype=bool))  # ignore each sample's correlation with itself
    table = pd.DataFrame({"most_similar": r.idxmax(), "r": r.max().round(4)})
    table.sort_values("r", ascending=False).to_csv(out / "sample_similarity.csv")
    print(f"similarity check: {identical} identical count columns; {int(expressed.sum())} genes; "
          f"highest pair r {table['r'].max():.4f} ({table['r'].idxmax()} - "
          f"{table.loc[table['r'].idxmax(), 'most_similar']}); median nearest r "
          f"{table['r'].median():.4f}")


def sensitivity_run() -> None:
    """PCA without config.SENSITIVITY_EXCLUDE: same rlog values and same 500 genes as the main
    run, refitted on the remaining samples, so that only the samples change. A variant that
    re-selects the 500 genes on the remaining samples is reported alongside."""
    out = OUT / "sensitivity_40"
    out.mkdir(parents=True, exist_ok=True)
    genes_file = OUT / "pca_genes.txt"
    if not genes_file.exists():
        raise SystemExit("missing results/step2_qc/pca_genes.txt: run the main Step 2 first")
    counts, meta = load_counts()
    keep = [s for s in meta.index if s not in config.SENSITIVITY_EXCLUDE]
    meta = meta.loc[keep]
    lib = counts[keep].sum(axis=0)
    rlog = pd.read_csv(config.PROCESSED / "rlog_geo.csv.gz", index_col=0)[keep]
    main_genes = genes_file.read_text(encoding="utf-8").split()
    main_pcs = pd.read_csv(OUT / "pca_scores.csv", index_col=0)
    main_pct = pd.read_csv(OUT / "pca_variance.csv")["variance_explained_pct"].to_numpy()

    pcs, pct, _ = pca_top(rlog, N_TOP, main_genes)
    pcs = orient(pcs, main_pcs)
    write_pca(pcs, pct, out)
    write_top_genes(rlog, main_genes, pcs, out)
    label = f"without {', '.join(config.SENSITIVITY_EXCLUDE)}"
    describe_pca(pcs, pct, meta, lib, out, label, full_cohort_limits())

    _, pct_new, genes_new = pca_top(rlog, N_TOP)
    corr = [np.corrcoef(pcs[pc], main_pcs.loc[keep, pc])[0, 1] for pc in ["PC1", "PC2"]]
    assoc = {name: pd.read_csv(d / "pca_associations.csv").set_index(["pc", "variable"])["value"]
             for name, d in [("main", OUT), ("sensitivity", out)]}
    lines = [
        f"excluded: {', '.join(config.SENSITIVITY_EXCLUDE)} ({len(keep)} samples left)",
        f"same {N_TOP} genes as the main PCA (results/step2_qc/pca_genes.txt), same rlog values",
        f"PC1 {pct[0]:.2f}% (main {main_pct[0]:.2f}%), PC2 {pct[1]:.2f}% (main {main_pct[1]:.2f}%)",
        f"score correlation with the main PCA on the {len(keep)} shared samples: "
        f"PC1 {corr[0]:.4f}, PC2 {corr[1]:.4f}",
        "PC2 highest: " + ", ".join(f"{s} ({v:.1f})" for s, v in pcs["PC2"].nlargest(3).items()),
    ]
    for pc, var in [("PC1", "condition"), ("PC1", "paper_group"), ("PC1", "batch"),
                    ("PC2", "condition"), ("PC2", "batch")]:
        lines.append(f"{pc} R^2 {var}: {assoc['sensitivity'][(pc, var)]:.3f} "
                     f"(main {assoc['main'][(pc, var)]:.3f})")
    # Batch R^2 of the main run's PC2 scores on the remaining samples, before refitting:
    # separates the effect of removing samples from the effect of refitting the PCA.
    y = main_pcs.loc[keep, "PC2"]
    within = ((y - y.groupby(meta["char:batch"]).transform("mean")) ** 2).sum()
    r2_kept = 1 - within / ((y - y.mean()) ** 2).sum()
    lines.append(f"PC2 R^2 batch, main run's PC2 scores of the {len(keep)} samples, "
                 f"not refitted: {r2_kept:.3f}")
    lines.append(f"variant, {N_TOP} genes re-selected on the {len(keep)} samples: PC1 {pct_new[0]:.2f}%, "
                 f"PC2 {pct_new[1]:.2f}%, {len(set(genes_new) & set(main_genes))} genes shared "
                 "with the main PCA")
    (out / "pca_comparison.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print("written to", out)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    config.PROCESSED.mkdir(parents=True, exist_ok=True)
    counts, meta = load_counts()
    counts = counts.loc[:, ~meta["condition"].isna()]
    meta = meta.loc[counts.columns]
    print(f"counts: {counts.shape[0]} genes x {counts.shape[1]} samples")

    # Library sizes: total deposited counts per sample, i.e. reads assigned to genes after
    # the authors' gene filter, not sequencing depth. Shown per batch, because depth differs
    # by sequencing run; colour and shape show condition (worksheets/step2_qc.md).
    # Watched samples (WATCH) are labelled when present.
    lib = counts.sum(axis=0)
    lib_table = meta[["condition", "char:batch"]].assign(lib_size=lib)
    lib_table.sort_values(["char:batch", "lib_size"]).to_csv(OUT / "library_sizes.csv")
    batches = sorted(lib_table["char:batch"].unique())
    x_pos = lib_table["char:batch"].map(batches.index)
    x_pos = x_pos + np.random.default_rng(0).uniform(-0.15, 0.15, len(x_pos))  # fixed jitter
    y_val = lib_table["lib_size"] / 1e6
    fig, ax = plt.subplots(figsize=(6, 4))
    styles = {"FSHD": ("o", "tab:red"), "Control": ("^", "tab:blue")}
    for condition, (marker, color) in styles.items():
        idx = lib_table["condition"] == condition
        ax.scatter(x_pos[idx], y_val[idx], marker=marker, s=45, color=color,
                   edgecolors="black", linewidths=0.5, alpha=0.85, label=condition)
    for sample in [w for w in WATCH if w in lib_table.index]:
        ax.annotate(sample, (x_pos[sample], y_val[sample]), fontsize=6,
                    xytext=(4, 0), textcoords="offset points")
    ax.set_xticks(range(len(batches)))
    ax.set_xticklabels([f"batch {b}" for b in batches])
    ax.set_ylabel("counts (millions)")
    ax.set_title("Library size per sample, by batch")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "library_sizes.png", dpi=150)
    plt.close(fig)

    # Gene filtering: the authors' rule, total counts above 20 across the samples
    # (makeSE.R line 164 at commit 90cd343:
    # https://github.com/TapscottLab/FSHD-Biopsy-RNA-seq/blob/90cd343b3d60a696a35041a33015131ab9908164/scripts/makeSE.R#L164).
    # The paper states no threshold. The deposited count file was already filtered with
    # this rule, so every gene should pass; gene_filter.txt records the evidence (smallest
    # gene total, genes before and after, and the four biomarkers of config.BIOMARKERS).
    min_total = 20
    totals = counts.sum(axis=1)
    keep = totals > min_total
    filtered = counts.loc[keep]
    # For comparison only: the starter's default rule, not applied (worksheets/step2_qc.md).
    starter_min_samples = int(meta["condition"].value_counts().min())
    starter_keep = (counts >= 10).sum(axis=1) >= starter_min_samples
    names = pd.read_csv(config.GENE_NAMES)
    symbol_to_id = dict(zip(names["gene_name"], names["gene_id"]))
    marker_lines = [
        f"  {symbol} ({symbol_to_id[symbol]}): total {int(totals[symbol_to_id[symbol]])}, "
        f"kept {bool(keep[symbol_to_id[symbol]])}"
        for symbol in config.BIOMARKERS
    ]
    (OUT / "gene_filter.txt").write_text(
        f"rule: total counts > {min_total} across {counts.shape[1]} samples "
        f"(authors' makeSE.R line 164, commit 90cd343)\n"
        f"smallest gene total in the deposited file: {int(totals.min())}\n"
        f"genes before: {len(counts)}\ngenes after: {len(filtered)}\n"
        f"starter default, not applied (>= 10 counts in >= {starter_min_samples} samples, "
        f"the smallest group): would keep {int(starter_keep.sum())}, "
        f"remove {int((~starter_keep).sum())}\n"
        "biomarkers:\n" + "\n".join(marker_lines) + "\n",
        encoding="utf-8",
    )
    print(f"genes kept: {len(filtered)} / {len(counts)} (smallest gene total {int(totals.min())})")
    filtered.to_csv(config.PROCESSED / "counts_filtered.csv")
    meta.to_csv(config.PROCESSED / "metadata.csv")

    # PCA with the authors' method (DESeq2 plotPCA on rlog, top 500 genes; makeManuscripts_supFigs.R
    # L82 at commit 90cd343), applied to the rlog recomputed from the GEO counts by
    # r/rlog_from_geo.R. Python recomputes the PCA from that matrix and must agree with R.
    rlog_file = config.PROCESSED / "rlog_geo.csv.gz"
    if not rlog_file.exists():
        raise SystemExit(f"missing {rlog_file}: run r/pca_replay.R, then r/rlog_from_geo.R")
    rlog = pd.read_csv(rlog_file, index_col=0)[meta.index]
    pcs, pct, rlog_genes = pca_top(rlog, N_TOP)
    r_pct = pd.read_csv(OUT / "rlog_geo_pca_variance.csv")["variance_explained_pct"].to_numpy()
    if not np.allclose(pct[:2], r_pct[:2], atol=0.01):
        raise SystemExit(f"Python PCA {pct[:2].round(2)} differs from R {r_pct[:2].round(2)}")
    # Orient PC1 and PC2 as in the Figure 3 replay, so the figures read like the paper's.
    pcs = orient(pcs, pd.read_csv(OUT / "pca_replay_scores.csv", index_col="geo_title"))
    write_pca(pcs, pct, OUT)
    (OUT / "pca_genes.txt").write_text("\n".join(rlog_genes) + "\n", encoding="utf-8")
    write_top_genes(rlog, rlog_genes, pcs, OUT)
    print(f"PCA (rlog of GEO counts, top {N_TOP} genes): PC1 {pct[0]:.2f}%, PC2 {pct[1]:.2f}% "
          f"(R: {r_pct[0]:.2f}%, {r_pct[1]:.2f}%)")

    # Comparison with the course default, log2(CPM + 1): on the 500 genes of Figure 3, on its
    # own top 500 and on its own top 2,000 genes.
    replay_genes = (OUT / "pca_replay_genes.txt").read_text(encoding="utf-8").split()
    logcpm = np.log2(filtered / filtered.sum(axis=0) * 1e6 + 1)
    rows = [{"input": "rlog of GEO counts", "genes": f"own top {N_TOP}",
             "pc1_pct": round(pct[0], 2), "pc2_pct": round(pct[1], 2),
             "figure3_genes_included": len(set(rlog_genes) & set(replay_genes))}]
    for label, n_top, genes in [("Figure 3's 500 genes", 500, replay_genes),
                                ("own top 500", 500, None),
                                ("own top 2000 (course default)", 2000, None)]:
        _, p, used = pca_top(logcpm, n_top, genes)
        rows.append({"input": "log2(CPM + 1)", "genes": label,
                     "pc1_pct": round(p[0], 2), "pc2_pct": round(p[1], 2),
                     "figure3_genes_included": len(set(used) & set(replay_genes))})
    pd.DataFrame(rows).to_csv(OUT / "pca_logcpm_comparison.csv", index=False)

    describe_pca(pcs, pct, meta, lib, OUT, "rlog of GEO counts",
                 colour_limits(pca_variables(meta, lib)))
    sex_check(counts, meta, OUT)
    similarity_check(counts, OUT)
    print("figures written to", OUT)
    print("\nNext: fill worksheets/step2_qc.md")


def figures_only() -> None:
    """Redraw the PCA figures of both runs from the saved scores and variances."""
    limits = full_cohort_limits()
    meta_all = load_sample_table()
    lib = pd.read_csv(OUT / "library_sizes.csv", index_col=0)["lib_size"]
    for out, label in [(OUT, "rlog of GEO counts"),
                       (OUT / "sensitivity_40",
                        f"without {', '.join(config.SENSITIVITY_EXCLUDE)}")]:
        pcs = pd.read_csv(out / "pca_scores.csv", index_col=0)
        pct = pd.read_csv(out / "pca_variance.csv")["variance_explained_pct"].to_numpy()
        pca_figures(pcs, pct, meta_all.loc[pcs.index], lib, out, label, limits)
        print(f"redrawn: {out} ({len(pcs)} samples)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step 2: QC and exploration")
    parser.add_argument("--sensitivity", action="store_true",
                        help="PCA without config.SENSITIVITY_EXCLUDE (after the main run)")
    parser.add_argument("--figures", action="store_true",
                        help="redraw the PCA figures of both runs from the saved scores")
    args = parser.parse_args()
    if args.sensitivity:
        sensitivity_run()
    elif args.figures:
        figures_only()
    else:
        main()
