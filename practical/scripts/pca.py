import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.dds import DeseqDataSet
import config  # noqa: E402

OUT = config.RESULTS / "step2_qc"
AUDIT = config.RESULTS / "step1_audit"

filtered = pd.read_csv(config.PROCESSED / "counts_filtered.csv", index_col=0)


def fit_dispersion_trend(base_mean, raw_disp):
    """
    Replicates DESeq2's parametric dispersion-mean trend: alpha_tr(mu) = a1 / mu + a0
    using Gamma-likelihood optimization on genes with sufficient mean expression.
    """
    valid = (base_mean > 1.0) & (raw_disp > 1e-6)
    mu_v = base_mean[valid]
    disp_v = raw_disp[valid]

    def gamma_deviance(params):
        a0, a1 = np.exp(params[0]), np.exp(params[1])
        pred = a1 / mu_v + a0
        return np.sum((disp_v - pred) / pred - np.log(disp_v / pred))

    res = minimize(gamma_deviance, x0=[np.log(0.05), np.log(1.0)], method="Nelder-Mead")
    a0, a1 = np.exp(res.x[0]), np.exp(res.x[1])
    return a1 / np.maximum(base_mean, 1e-8) + a0


def deseq2_rlog_exact_python(counts_df):
    """
    Pure-Python replication of DESeq2::rlog(dds, blind=True).
    For Transparency : this function was genereated by gemini to be a reference for the top genes selection.
    I do not fully understand each steps but my goal was not to use it as the primary method but just as a reference.
    """
    counts = counts_df.astype(float).copy()
    k_mat = counts.values
    ln2 = np.log(2.0)

    # 1. DESeq2 size factors (median-of-ratios)
    nonzero = k_mat[(k_mat > 0).all(axis=1)]
    log_geo_means = np.log(nonzero).mean(axis=1)
    size_factors = np.median(nonzero / np.exp(log_geo_means)[:, None], axis=0)

    norm_counts = k_mat / size_factors[None, :]
    base_mean = norm_counts.mean(axis=1)

    # 2. Gene-wise dispersion + parametric trend alpha(mu) = a1 / mu + a0
    s_inv_mean = np.mean(1.0 / size_factors)
    gene_var = norm_counts.var(axis=1, ddof=1)
    raw_disp = np.maximum(
        (gene_var - base_mean * s_inv_mean) / (base_mean**2 + 1e-8), 1e-6
    )
    disp_trend = fit_dispersion_trend(base_mean, raw_disp)

    # 3. Geometric intercept (beta_0)
    log_norm_05 = np.log2(norm_counts + 0.5)
    beta_0 = log_norm_05.mean(axis=1)

    # 4. Empirical Bayes prior variance (sigma_b^2)
    raw_lfc = log_norm_05 - beta_0[:, None]
    noise_var = (1.0 / np.maximum(norm_counts, 0.5) + disp_trend[:, None]) / (ln2**2)

    high = base_mean > 10.0
    excess_var = np.maximum(
        np.var(raw_lfc[high], axis=1, ddof=1) - np.mean(noise_var[high], axis=1), 0.0
    )
    prior_var = float(np.quantile(excess_var, 0.85))
    prior_var = max(prior_var, 0.05)

    # 5. Vectorized Newton-Raphson MAP estimation
    beta_ij = np.zeros_like(k_mat)
    alpha_mat = disp_trend[:, None]
    sf_mat = size_factors[None, :]
    b0_mat = beta_0[:, None]

    for _ in range(25):
        mu = sf_mat * (2.0 ** (b0_mat + beta_ij))
        denom = 1.0 + alpha_mat * mu
        w = mu / denom
        grad = ln2 * (k_mat - mu) / denom - (beta_ij / prior_var)
        hess = -(ln2**2) * w - (1.0 / prior_var)
        step = grad / hess
        beta_ij -= step
        if np.max(np.abs(step)) < 1e-6:
            break

    rlog_mat = b0_mat + beta_ij
    return pd.DataFrame(rlog_mat, index=counts.index, columns=counts.columns), prior_var


def python_pydeseq2_vst(counts_df):
    """Simpler version for selecting the best genes based on the method Variance Stabilizing Transformation"""
    dds = DeseqDataSet(
        counts=counts_df.T,
        metadata=pd.DataFrame(index=counts_df.columns),
        design="~1",
        quiet=True,
    )

    # Run Variance Stabilizing Transformation
    dds.vst(fit_type="parametric")

    # get the result for each genes
    vst_df = pd.DataFrame(
        dds.layers["vst_counts"].T,
        index=counts_df.index,
        columns=counts_df.columns,
    )
    return vst_df


if __name__ == "__main__":

    rlog_py, p_var = deseq2_rlog_exact_python(filtered)
    print(f"Genes retained: {rlog_py.shape[0]} | estimated prior_var: {p_var:.4f}")

    ntop = 2500
    top_genes = rlog_py.var(axis=1).sort_values(ascending=False).index[:ntop]
    x = rlog_py.loc[top_genes].T
    x = x - x.mean(axis=0)

    pca = PCA(n_components=10)
    names = [f"PC{i + 1}" for i in range(10)]
    pcs = pd.DataFrame(pca.fit_transform(x), index=x.index, columns=names)
    ratio = pca.explained_variance_ratio_ * 100

    print(
        f"Pure-Python rlog PCA (ntop={ntop}) -> PC1: {ratio[0]:.2f}%, PC2: {ratio[1]:.2f}%"
    )

    plt.figure(figsize=(8, 6))
    plt.scatter(pcs["PC1"], pcs["PC2"], s=65, alpha=0.85, edgecolors="black")
    for sample_id, row in pcs.iterrows():
        plt.annotate(sample_id, (row["PC1"], row["PC2"]), fontsize=7, alpha=0.8)
    plt.xlabel(f"PC1: {ratio[0]:.0f}% variance")
    plt.ylabel(f"PC2: {ratio[1]:.0f}% variance")
    plt.title(f"Pure-Python rlog PCA (ntop={ntop})")
    plt.axhline(0, color="gray", linestyle="--", linewidth=0.5)
    plt.axvline(0, color="gray", linestyle="--", linewidth=0.5)
    plt.tight_layout()
    plt.show()
