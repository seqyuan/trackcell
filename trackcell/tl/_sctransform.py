"""
Core SCTransform / VST implementation.

Ported from Seurat 4.4.0 (sctransform::vst) with AnnData-oriented helpers.
Reference: Hafemeister & Satija, Genome Biology 2019.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Optional, Sequence, Union

import numba
import numpy as np
import pandas as pd
import statsmodels.api as sm
from joblib import Parallel, cpu_count, delayed
from patsy import dmatrix
from scipy import optimize
from scipy.interpolate import interp1d
from scipy.sparse import csr_matrix, issparse
from scipy.special import digamma, polygamma
from ._r_density import r_density_interp
from ._r_sample import make_vst_rng, r_sample_available
from ._sctransform_v2 import (
    LOG10_SLOPE,
    apply_vst_flavor_settings,
    build_poisson_gene_tables,
    fit_nb_offset,
    fit_offset_model,
    mark_suspicious_theta,
    resolve_min_variance,
)

ArrayLike = Union[np.ndarray, csr_matrix]
KNOWN_CELL_ATTR = {
    "umi",
    "gene",
    "log_umi",
    "log_gene",
    "umi_per_gene",
    "log_umi_per_gene",
}
RESIDUAL_TYPES = {"pearson", "deviance", "none"}
SCALE_MODELS = {"linear", "poisson", "negbinom"}


def _as_genes_by_cells(
    x: ArrayLike,
    genes: Sequence[str],
    cells: Sequence[str],
) -> csr_matrix:
    """Return UMI matrix as genes x cells CSR float64."""
    if issparse(x):
        mat = x.T.tocsr()
    else:
        mat = csr_matrix(np.asarray(x, dtype=np.float64).T)
    if mat.shape != (len(genes), len(cells)):
        raise ValueError(
            f"Expression matrix shape {mat.shape} does not match "
            f"{len(genes)} genes x {len(cells)} cells."
        )
    return mat.astype(np.float64, copy=False)


def _expand_genes_gbc(
    mat: csr_matrix,
    genes_subset: pd.Index,
    genes_full: pd.Index,
) -> csr_matrix:
    """Expand a genes-by-cells matrix from ``genes_subset`` order to ``genes_full``."""
    if len(genes_subset) == len(genes_full) and genes_subset.equals(genes_full):
        return mat
    pos = genes_full.get_indexer(genes_subset)
    valid = pos >= 0
    out = csr_matrix((len(genes_full), mat.shape[1]), dtype=mat.dtype)
    if valid.any():
        out[pos[valid], :] = mat[valid, :]
    return out


def _as_list(value: Optional[Sequence[str]], default: Sequence[str]) -> list[str]:
    if value is None:
        return list(default)
    return list(value)


def _build_model_str(
    latent_var: Sequence[str],
    batch_var: Optional[str] = None,
) -> str:
    latent_var = _as_list(latent_var, ("log_umi",))
    lv = " + ".join(latent_var)
    if batch_var:
        return f"y ~ ({lv}):{batch_var} + {batch_var} - 1"
    return f"y ~ {lv}"


def _get_model_formula(model_str: str) -> str:
    return model_str.replace("y ~", "~", 1)


def _build_design_matrix(cell_attr: pd.DataFrame, model_str: str) -> pd.DataFrame:
    return dmatrix(_get_model_formula(model_str), cell_attr, return_type="dataframe")


def _coef_columns(model_pars: pd.DataFrame, regressor_data: pd.DataFrame) -> list[str]:
    return [col for col in model_pars.columns if col != "theta" and col in regressor_data.columns]


def _batch_coef_columns(
    columns: Sequence[str],
    batch_var: str,
    batch_level: Any,
) -> list[str]:
    token = f"{batch_var}[{batch_level}]"
    return [col for col in columns if token in col]


def _predict_mu(
    model_pars: pd.DataFrame,
    genes: Sequence[str],
    regressor_data: pd.DataFrame,
) -> np.ndarray:
    coef_cols = _coef_columns(model_pars, regressor_data)
    coefs = model_pars.loc[genes, coef_cols].to_numpy(dtype=np.float64)
    design = regressor_data[coef_cols].to_numpy(dtype=np.float64)
    return np.exp(coefs @ design.T)


def _make_cell_attr(
    umi: csr_matrix,
    cells: pd.Index,
    cell_attr: Optional[pd.DataFrame] = None,
    latent_var: Optional[Sequence[str]] = None,
    batch_var: Optional[str] = None,
    latent_var_nonreg: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    latent_var = _as_list(latent_var, ("log_umi",))
    needed = set(latent_var)
    if batch_var:
        needed.add(batch_var)
    if latent_var_nonreg:
        needed.update(latent_var_nonreg)

    if cell_attr is None:
        cell_attr = pd.DataFrame(index=cells)
    else:
        cell_attr = cell_attr.reindex(cells).copy()

    missing = [name for name in needed if name not in cell_attr.columns]
    unknown = [name for name in missing if name not in KNOWN_CELL_ATTR]
    if unknown:
        raise ValueError(
            "Unknown cell attributes requested via latent_var/batch_var/latent_var_nonreg: "
            f"{unknown}. Provide them in cell_attr."
        )

    if any(name in missing for name in ("umi", "log_umi", "umi_per_gene", "log_umi_per_gene")):
        umi_per_cell = np.asarray(umi.sum(axis=0)).ravel()
        cell_attr["umi"] = umi_per_cell
        cell_attr["log_umi"] = np.log10(umi_per_cell + 1e-9)
    if any(name in missing for name in ("gene", "log_gene", "umi_per_gene", "log_umi_per_gene")):
        gene_per_cell = np.asarray((umi > 0).sum(axis=0)).ravel()
        cell_attr["gene"] = gene_per_cell
        cell_attr["log_gene"] = np.log10(gene_per_cell + 1e-9)
    if any(name in missing for name in ("umi_per_gene", "log_umi_per_gene")):
        cell_attr["umi_per_gene"] = cell_attr["umi"] / np.maximum(cell_attr["gene"], 1)
        cell_attr["log_umi_per_gene"] = np.log10(cell_attr["umi_per_gene"] + 1e-9)

    for name in needed:
        values = cell_attr[name]
        if pd.api.types.is_numeric_dtype(values):
            if pd.isna(values).any() or np.isinf(values.to_numpy(dtype=np.float64, na_value=np.nan)).any():
                raise ValueError(f'cell attribute "{name}" contains NA or infinite values.')
        elif pd.isna(values).any():
            raise ValueError(f'cell attribute "{name}" contains NA values.')

    if batch_var:
        cell_attr[batch_var] = cell_attr[batch_var].astype("category")
    return cell_attr


def _design_log_umi(log_umi: np.ndarray, index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {"Intercept": np.ones(len(log_umi)), "log_umi": log_umi},
        index=index,
    )


def _make_cell_attr_simple(umi: csr_matrix, cells: pd.Index) -> pd.DataFrame:
    return _make_cell_attr(umi, cells)


@numba.jit(cache=True, forceobj=True, nogil=True)
def _row_gmean(values: np.ndarray, gmean_eps: float = 1.0) -> np.ndarray:
    return np.exp(np.log(values + gmean_eps).mean(1)) - gmean_eps


def _row_gmean_sparse(umi: csr_matrix, gmean_eps: float = 1.0) -> np.ndarray:
    """Per-gene geometric mean (sctransform::row_gmean), genes × cells CSR."""
    out = np.empty(umi.shape[0], dtype=np.float64)
    for i in range(umi.shape[0]):
        row = umi.getrow(i).toarray().ravel()
        # One row × n_cells — mean axis=1 matches R rowMeans(log(x + eps)).
        out[i] = _row_gmean(row.reshape(1, -1), gmean_eps)[0]
    return out


def _dds(genes_log10_gmean: pd.Series, grid_points: int = 512) -> np.ndarray:
    """
    Inverse density weights for gene subsampling.

    Matches R ``vst()``: ``density(genes_log_gmean_step1, bw='nrd', adjust=1)``
    then ``approx(..., xout=genes_log_gmean_step1)`` for sampling weights.
    """
    values = np.asarray(genes_log10_gmean, dtype=np.float64)
    if values.size == 0:
        return np.array([], dtype=np.float64)
    density_at = r_density_interp(values, values, bw="nrd", adjust=1.0, n=grid_points)
    # R passes unnormalized ``1/(density + eps)`` to ``sample(..., prob=)``.
    return 1.0 / (density_at + np.finfo(float).eps)


@numba.jit(cache=True, forceobj=True, nogil=True)
def _qpois_reg(
    x: np.ndarray,
    y: np.ndarray,
    tol: float = 1e-9,
    maxiters: int = 100,
    minphi: float = 1.0001,
    returnfit: bool = True,
) -> dict[str, Any]:
    n, pcols = x.shape
    b_old = np.zeros(pcols, dtype=np.float64)
    b_new = np.zeros(pcols, dtype=np.float64)
    y = y.reshape(-1)

    for i in range(pcols):
        unique_vals = np.unique(x[:, i])
        if unique_vals.shape[0] == 1:
            b_old[i] = np.log(np.mean(y))
            break
        if unique_vals.shape[0] == 2 and (unique_vals[0] == 0 or unique_vals[1] == 0):
            b_old[i] = np.log(max(1e-9, np.dot(y, x[:, i]) / np.sum(x[:, i])))

    dif = 1.0
    it = 0
    while dif > tol and it < maxiters:
        yhat = np.clip(np.dot(x, b_old), -708, 709)
        m = np.exp(yhat)
        phi = y - m
        l1 = np.dot(x.T, phi)
        l2 = np.dot(x.T, x * m[:, None])
        det = np.linalg.det(l2)
        l2_inv = np.linalg.pinv(l2) if abs(det) < 1e-15 else np.linalg.inv(l2)
        b_new = b_old + np.dot(l2_inv, l1)
        dif = np.sum(np.abs(b_new - b_old))
        b_old = b_new
        it += 1

    phi[np.abs(phi) < 1.64487933e-154] = 0.0
    p = np.sum(phi * phi / m) / (n - pcols)
    return {
        "coefficients": b_new,
        "phi": p,
        "theta.guesstimate": np.mean(m) / (max(p, minphi) - 1),
        "fitted": m if returnfit else None,
    }


@numba.jit(cache=True, forceobj=True, nogil=True)
def _trigamma(x: float) -> float:
    return polygamma(1, x)


@numba.jit(cache=True, forceobj=True, nogil=True)
def _score(n: float, th: float, mu: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    a = th + y
    b = th + mu
    return float(np.sum(w * (digamma(a) - digamma(th) + np.log(th) + 1 - np.log(b) - a / b)))


@numba.jit(cache=True, forceobj=True, nogil=True)
def _info(n: float, th: float, mu: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    a = th + y
    b = th + mu
    return float(
        np.sum(w * (-_trigamma(a) + _trigamma(th) - 1 / th + 2 / b - a / np.power(b, 2)))
    )


def _theta_ml(y: np.ndarray, mu: np.ndarray, limit: int = 10, eps: float = 0.0001220703) -> float:
    weights = np.ones(len(y))
    n = float(np.sum(weights))
    t0 = n / np.sum(weights * (y / mu - 1) ** 2)
    it = 1
    delta = 1.0
    while it < limit and abs(delta) > eps:
        t0 = abs(t0)
        i_val = _info(n, t0, mu, y, weights)
        delta = _score(n, t0, mu, y, weights) / i_val
        t0 = t0 + delta
        it += 1
    return max(t0, 0.0)


def _one_row_fit_poisson(
    regressor_data: pd.DataFrame,
    y: np.ndarray,
) -> Optional[list[float]]:
    try:
        fit = _qpois_reg(regressor_data.to_numpy(), y, returnfit=True)
        theta = _theta_ml(y=y, mu=fit["fitted"])
        return [float(theta), *[float(v) for v in fit["coefficients"]]]
    except (np.linalg.LinAlgError, ValueError, FloatingPointError):
        return None


def _fit_poisson(
    umi: csr_matrix,
    regressor_data: pd.DataFrame,
    gene_index: pd.Index,
) -> pd.DataFrame:
    columns = ["theta", *regressor_data.columns.tolist()]
    results = Parallel(n_jobs=cpu_count(), backend="threading")(
        delayed(_one_row_fit_poisson)(regressor_data, row.toarray().ravel())
        for row in umi
    )
    failed = sum(r is None for r in results)
    if failed:
        warnings.warn(
            f"{failed} gene(s) failed Poisson regression and will be excluded.",
            stacklevel=3,
        )
    n_coef = len(columns)
    valid = [r if r is not None else [np.nan] * n_coef for r in results]
    return pd.DataFrame(valid, columns=columns, index=gene_index)


def _mad(x: np.ndarray) -> float:
    return float(np.median(np.abs(x - np.median(x))) * 1.4826)


def _robust_scale(x: np.ndarray) -> np.ndarray:
    return (x - np.median(x)) / (_mad(x) + np.finfo(float).eps)


def _robust_scale_binned(y: pd.Series, x: pd.Series, breaks: np.ndarray) -> np.ndarray:
    x_aligned = x.reindex(y.index)
    bins = pd.cut(x_aligned, bins=breaks, ordered=True)
    scores = np.zeros(len(y), dtype=np.float64)
    for bin_cat, idx in y.groupby(bins, observed=False).groups.items():
        if pd.isna(bin_cat):
            continue
        ys = y.loc[idx].to_numpy(dtype=np.float64)
        if len(ys) == 0:
            continue
        pos = y.index.get_indexer(idx)
        scores[pos] = _robust_scale(ys)
    return scores


def _seq(start: float, stop: float, step: float) -> np.ndarray:
    return np.arange(start, stop, step)


def _is_outlier(y: pd.Series, x: pd.Series, th: float = 10.0) -> np.ndarray:
    x_aligned = x.reindex(y.index)
    if y.isna().all() or x_aligned.isna().all():
        return np.zeros(len(y), dtype=bool)
    finite = np.isfinite(y.to_numpy()) & np.isfinite(x_aligned.to_numpy())
    if finite.sum() < 3:
        return np.zeros(len(y), dtype=bool)
    eps = np.finfo(float).eps * 10
    x_max, x_min = float(x_aligned.max()), float(x_aligned.min())
    if not np.isfinite(x_max - x_min) or (x_max - x_min) <= eps:
        return np.zeros(len(y), dtype=bool)
    bin_width = (x_max - x_min) * bw_sj(x_aligned.to_numpy()) / 2
    breaks1 = _seq(x_min - eps, x_max + bin_width, bin_width)
    breaks2 = _seq(x_min - eps - bin_width / 2, x_max + bin_width, bin_width)
    score1 = _robust_scale_binned(y, x_aligned, breaks1)
    score2 = _robust_scale_binned(y, x_aligned, breaks2)
    return np.minimum(np.abs(score1), np.abs(score2)) > th


def bw_sj(x: np.ndarray, nb: int = 1000) -> float:
    """Sheather-Jones bandwidth selector (R stats::bw.SJ)."""
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 3:
        return 0.5
    if np.allclose(x, x[0]):
        return 0.5

    d, cnt = _bw_pair_cnts(x, nb, n > nb // 2)

    def sdh(h: float) -> float:
        return _bw_phi4(n, d, cnt, h)

    def tdh(h: float) -> float:
        return _bw_phi6(n, d, cnt, h)

    q75, q25 = np.percentile(x, [75, 25])
    x_iqr = q75 - q25
    scale = min(float(np.std(x, ddof=1)), x_iqr / 1.349)
    if not np.isfinite(scale) or scale <= 0:
        scale = max(float(np.std(x, ddof=1)), 1e-6)
    silverman = 0.9 * scale * np.power(n, -0.2)
    a = 1.24 * scale * np.power(n, -1 / 7)
    b = 1.23 * scale * np.power(n, -1 / 9)
    c1 = 1 / (2 * np.sqrt(np.pi) * n)
    td = -tdh(b)
    if not np.isfinite(td) or td <= 0:
        return float(silverman)
    hmax = 1.144 * scale * np.power(n, -1 / 5)
    lower, upper = 0.1 * hmax, hmax
    sd_a = sdh(a)
    if not np.isfinite(sd_a) or sd_a <= 0:
        return float(silverman)
    alph2 = 1.357 * np.power(sd_a / td, 1 / 7)

    def fsd(h: float) -> float:
        return np.power(c1 / sdh(alph2 * np.power(h, 5 / 7)), 1 / 5) - h

    itry = 1
    while fsd(lower) * fsd(upper) > 0:
        if itry >= 99:
            return float(silverman)
        if itry % 2:
            upper *= 1.2
        else:
            lower /= 1.2
        itry += 1
    try:
        return float(optimize.brentq(fsd, lower, upper, xtol=0.1 * lower))
    except ValueError:
        return float(silverman)


def _bw_pair_cnts(x: np.ndarray, nb: int, binned: bool) -> tuple[float, np.ndarray]:
    if binned:
        r = (float(np.min(x)), float(np.max(x)))
        d = (r[1] - r[0]) * 1.01 / nb
        xx = np.trunc(np.abs(x) / d) * np.sign(x)
        xx = xx - np.min(xx) + 1
        xxx = np.bincount(xx.astype(np.int64), minlength=nb + 1)[1:]
        return d, _bw_den_binned(xxx)
    return _bw_den(nb, x)


def _bw_den(nb: int, sx: np.ndarray) -> tuple[float, np.ndarray]:
    sx = sx[np.isfinite(sx)]
    n = len(sx)
    if n == 0:
        return 1.0, np.zeros(nb, dtype=float)
    cnt = np.zeros(nb, dtype=float)
    xmin = xmax = float(sx[0])
    for val in sx:
        xmin = min(xmin, float(val))
        xmax = max(xmax, float(val))
    dd = max((xmax - xmin) * 1.01 / nb, np.finfo(float).eps)
    for i in range(n):
        ii = int(sx[i] / dd)
        for j in range(i):
            jj = int(sx[j] / dd)
            cnt[abs(ii - jj)] += 1
    return dd, cnt


def _bw_den_binned(sx: np.ndarray) -> np.ndarray:
    nb = len(sx)
    cnt = np.zeros(nb, dtype=float)
    for ii in range(nb):
        w = sx[ii]
        cnt[0] += w * (w - 1.0)
        for jj in range(ii):
            cnt[ii - jj] += w * sx[jj]
    cnt[0] *= 0.5
    return cnt


def _bw_phi4(sn: int, sd: float, cnt: np.ndarray, sh: float) -> float:
    total = 0.0
    for i, val in enumerate(cnt):
        delta = (i * sd / sh) ** 2
        if delta >= 1000:
            break
        term = np.exp(-delta / 2) * (delta * delta - 6 * delta + 3)
        total += term * val
    total = 2 * total + sn * 3
    return total / (sn * (sn - 1) * sh**5 * np.sqrt(2 * np.pi))


def _bw_phi6(sn: int, sd: float, cnt: np.ndarray, sh: float) -> float:
    total = 0.0
    for i, val in enumerate(cnt):
        delta = (i * sd / sh) ** 2
        if delta >= 1000:
            break
        term = np.exp(-delta / 2) * (delta**3 - 15 * delta**2 + 45 * delta - 15)
        total += term * val
    total = 2 * total - 15 * sn
    return total / (sn * (sn - 1) * sh**7 * np.sqrt(2 * np.pi))


@numba.jit(nopython=True)
def _bdr_ksmooth(
    x: np.ndarray,
    y: np.ndarray,
    xp: np.ndarray,
    kern: int,
    bw: float,
) -> np.ndarray:
    n = len(x)
    nxp = len(xp)
    yp = np.zeros(nxp, dtype=np.float64)
    imin = 0
    cutoff = 0.0
    if kern == 1:
        bw *= 0.5
        cutoff = bw
    elif kern == 2:
        bw *= 0.3706506
        cutoff = 4 * bw

    while imin < n and x[imin] < xp[0] - cutoff:
        imin += 1

    for j in range(nxp):
        num = den = 0.0
        x0 = xp[j]
        for i in range(imin, n):
            if x[i] < x0 - cutoff:
                imin = i
            elif x[i] > x0 + cutoff:
                break
            else:
                delta = abs(x[i] - x0) / bw
                w = 1.0 if kern == 1 else np.exp(-0.5 * delta * delta)
                num += w * y[i]
                den += w
        yp[j] = num / den if den > 0 else 0.0
    return yp


def _ksmooth(
    x: pd.Series,
    y: pd.Series,
    xp: np.ndarray,
    kern: int = 2,
    bw: float = 1.0,
) -> np.ndarray:
    aligned = pd.DataFrame({"x": x, "y": y}).dropna()
    order = np.argsort(aligned["x"].to_numpy())
    xs = aligned["x"].to_numpy()[order]
    ys = aligned["y"].to_numpy()[order]
    xp_sorted = np.sort(xp.copy())
    return _bdr_ksmooth(xs, ys, xp_sorted, kern, bw)


@numba.jit(cache=True, forceobj=True, nogil=True)
def _deviance_residual(
    y: np.ndarray,
    mu: np.ndarray,
    theta: np.ndarray,
) -> np.ndarray:
    y_safe = np.maximum(y, 0.0)
    y_log = np.zeros_like(y_safe)
    mask = y_safe > 0
    y_log[mask] = y_safe[mask] * np.log(y_safe[mask] / np.maximum(mu[mask], 1e-9))
    theta_row = theta.reshape(-1, 1)
    r = 2.0 * (y_log - (y_safe + theta_row) * np.log((y_safe + theta_row) / (mu + theta_row)))
    return np.sqrt(np.maximum(r, 0.0)) * np.sign(y - mu)


@numba.jit(cache=True, forceobj=True, nogil=True)
def _pearson_residual(
    y: np.ndarray,
    mu: np.ndarray,
    theta: np.ndarray,
    min_var: float = -np.inf,
) -> np.ndarray:
    theta_row = theta.reshape(-1, 1)
    poisson = ~np.isfinite(theta_row)
    variance = mu + np.divide(mu**2, theta_row, where=~poisson, out=np.zeros_like(mu))
    variance = np.where(poisson, mu, variance)
    variance[variance < min_var] = min_var
    return np.divide(y - mu, np.sqrt(variance))


def _compute_residuals(
    y: np.ndarray,
    mu: np.ndarray,
    theta: np.ndarray,
    residual_type: str,
    min_var: float = -np.inf,
) -> np.ndarray:
    if residual_type == "pearson":
        return _pearson_residual(y, mu, theta, min_var)
    if residual_type == "deviance":
        return _deviance_residual(y, mu, theta)
    raise ValueError(f"Unsupported residual_type: {residual_type}")


def _multi_compute_residual(
    bin_idx: int,
    model_pars_fit: pd.DataFrame,
    regressor_data: pd.DataFrame,
    umi: csr_matrix,
    genes: np.ndarray,
    bin_ind: np.ndarray,
    residual_type: str,
    min_variance: float,
) -> np.ndarray:
    genes_bin = genes[bin_ind == bin_idx]
    mu = _predict_mu(model_pars_fit, genes_bin, regressor_data)
    theta = model_pars_fit.loc[genes_bin, "theta"]
    mask = np.isin(genes, genes_bin)
    y = umi[mask, :].toarray()
    return _compute_residuals(y, mu, theta.to_numpy(), residual_type, min_variance)


def _multi_pearson_residual(
    bin_idx: int,
    model_pars_fit: pd.DataFrame,
    regressor_data: pd.DataFrame,
    umi: csr_matrix,
    genes: np.ndarray,
    bin_ind: np.ndarray,
    min_variance: float,
) -> np.ndarray:
    return _multi_compute_residual(
        bin_idx,
        model_pars_fit,
        regressor_data,
        umi,
        genes,
        bin_ind,
        "pearson",
        min_variance,
    )


def _reg_model_pars_pre_ksmooth(
    model_pars: pd.DataFrame,
    genes_log_gmean_step1: pd.Series,
    genes_log_gmean: pd.Series,
    *,
    theta_regularization: str = "od_factor",
    umi: Optional[csr_matrix] = None,
    genes: Optional[pd.Index] = None,
    genes_step1: Optional[pd.Index] = None,
    exclude_poisson: bool = False,
    bw_adjust: float = 3.0,
) -> dict[str, object]:
    """Mirror ``reg_model_pars`` steps before ksmooth (for R parity debugging)."""
    if genes is None:
        genes = genes_log_gmean.index

    all_poisson_genes = pd.Index([])
    if exclude_poisson and umi is not None and genes_step1 is not None:
        all_poisson_genes, _ = build_poisson_gene_tables(umi, genes, genes_step1, model_pars)

    model_pars_all = model_pars.copy()
    gmean_step1 = genes_log_gmean_step1.reindex(model_pars.index)
    theta_vals = model_pars["theta"].to_numpy(dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        if theta_regularization == "theta":
            dispersion = np.log10(theta_vals)
        elif theta_regularization == "od_factor":
            dispersion = np.log10(1 + np.power(10, gmean_step1.to_numpy()) / theta_vals)
        else:
            raise ValueError(f"Unknown theta_regularization: {theta_regularization}")
    dispersion[~np.isfinite(theta_vals)] = 0.0

    model_pars_pre_outlier = model_pars.drop(columns=["theta"]).copy()
    model_pars_pre_outlier["dispersion_par"] = dispersion

    outliers = np.zeros(len(model_pars_pre_outlier), dtype=bool)
    gmean_for_outlier = genes_log_gmean_step1.reindex(model_pars_pre_outlier.index)
    for col in model_pars_pre_outlier.columns:
        outliers |= _is_outlier(model_pars_pre_outlier[col], gmean_for_outlier)
    if exclude_poisson:
        outliers |= ~np.isfinite(model_pars_all["theta"].reindex(model_pars_pre_outlier.index).to_numpy())
    outliers_series = pd.Series(outliers, index=model_pars_pre_outlier.index, dtype=bool)

    model_pars_work = model_pars_pre_outlier
    gmean_work = genes_log_gmean_step1.reindex(model_pars_work.index)
    if outliers.any():
        model_pars_work = model_pars_work.loc[~outliers]
        gmean_work = gmean_work.reindex(model_pars_work.index)

    overdispersed = model_pars_work.index
    if exclude_poisson and len(all_poisson_genes):
        overdispersed = model_pars_work.index.difference(all_poisson_genes)
        model_pars_work = model_pars_work.loc[overdispersed]
        gmean_work = gmean_work.reindex(model_pars_work.index)

    bw = bw_sj(gmean_work.to_numpy()) * bw_adjust
    return {
        "all_poisson_genes": all_poisson_genes,
        "outliers": outliers_series,
        "model_pars_pre_outlier": model_pars_pre_outlier,
        "model_pars_pre_ksmooth": model_pars_work,
        "genes_log_gmean_step1_pre_ksmooth": gmean_work,
        "overdispersed_genes": overdispersed,
        "bw": bw,
    }


def _reg_model_pars(
    model_pars: pd.DataFrame,
    genes_log_gmean_step1: pd.Series,
    genes_log_gmean: pd.Series,
    bw_adjust: float,
    theta_regularization: str,
    *,
    cell_attr: Optional[pd.DataFrame] = None,
    batch_var: Optional[str] = None,
    cells_step1: Optional[pd.Index] = None,
    genes_step1: Optional[pd.Index] = None,
    umi: Optional[csr_matrix] = None,
    gmean_eps: float = 1.0,
    exclude_poisson: bool = False,
    fix_slope: bool = False,
    genes: Optional[pd.Index] = None,
) -> tuple[pd.DataFrame, pd.Series]:
    if genes is None:
        genes = genes_log_gmean.index
    if batch_var is not None:
        exclude_poisson = False
        fix_slope = False

    all_poisson_genes = pd.Index([])
    vst_out_offset: Optional[pd.DataFrame] = None
    if (exclude_poisson or fix_slope) and umi is not None and genes_step1 is not None:
        all_poisson_genes, vst_out_offset = build_poisson_gene_tables(
            umi, genes, genes_step1, model_pars
        )

    model_pars_all = model_pars.copy()
    gmean_step1 = genes_log_gmean_step1.reindex(model_pars.index)
    theta_vals = model_pars["theta"].to_numpy(dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        if theta_regularization == "theta":
            dispersion = np.log10(theta_vals)
        elif theta_regularization == "od_factor":
            dispersion = np.log10(1 + np.power(10, gmean_step1.to_numpy()) / theta_vals)
        else:
            raise ValueError(f"Unknown theta_regularization: {theta_regularization}")
    dispersion[~np.isfinite(theta_vals)] = 0.0

    model_pars = model_pars.drop(columns=["theta"]).copy()
    model_pars["dispersion_par"] = dispersion

    outliers = np.zeros(len(model_pars), dtype=bool)
    gmean_for_outlier = genes_log_gmean_step1.reindex(model_pars.index)
    for col in model_pars.columns:
        outliers |= _is_outlier(model_pars[col], gmean_for_outlier)
    if exclude_poisson:
        outliers |= ~np.isfinite(model_pars_all["theta"].reindex(model_pars.index).to_numpy())
    if outliers.any():
        outliers_series = pd.Series(outliers, index=model_pars.index, dtype=bool)
        model_pars = model_pars.loc[~outliers]
        genes_log_gmean_step1 = genes_log_gmean_step1.reindex(model_pars.index)
    else:
        outliers_series = pd.Series(outliers, index=model_pars.index, dtype=bool)

    if exclude_poisson and len(all_poisson_genes):
        overdispersed = model_pars.index.difference(all_poisson_genes)
        model_pars = model_pars.loc[overdispersed]
        genes_log_gmean_step1 = genes_log_gmean_step1.reindex(model_pars.index)

    genes_arr = genes.to_numpy()
    min_step1, max_step1 = genes_log_gmean_step1.min(), genes_log_gmean_step1.max()
    x_points = genes_log_gmean.clip(lower=min_step1, upper=max_step1).to_numpy()
    order = np.argsort(x_points)
    model_pars_fit = pd.DataFrame(
        np.zeros((len(genes_arr), len(model_pars.columns))),
        index=genes,
        columns=model_pars.columns,
    )
    bw = bw_sj(genes_log_gmean_step1.to_numpy()) * bw_adjust
    gmean_step1_fit = genes_log_gmean_step1.reindex(model_pars.index)

    model_pars_fit.iloc[order, model_pars_fit.columns.get_loc("dispersion_par")] = _ksmooth(
        gmean_step1_fit,
        model_pars["dispersion_par"],
        x_points.copy(),
        kern=2,
        bw=bw,
    )

    coef_cols = [col for col in model_pars.columns if col != "dispersion_par"]
    if batch_var is None or cell_attr is None or umi is None:
        for col in coef_cols:
            model_pars_fit.iloc[order, model_pars_fit.columns.get_loc(col)] = _ksmooth(
                gmean_step1_fit,
                model_pars[col],
                x_points.copy(),
                kern=2,
                bw=bw,
            )
    else:
        if cells_step1 is None or genes_step1 is None:
            raise ValueError("cells_step1 and genes_step1 are required when batch_var is set.")
        genes_step1_arr = np.asarray(genes_step1)
        genes_step1_idx = np.array(
            [i for i, g in enumerate(genes_arr) if g in set(genes_step1_arr)]
        )
        batches = cell_attr[batch_var].cat.categories
        for batch_level in batches:
            sel_step1 = cell_attr.index[
                (cell_attr[batch_var] == batch_level) & cell_attr.index.isin(cells_step1)
            ]
            cells_step1_mask = cell_attr.index.isin(sel_step1)
            umi_batch_step1 = umi[genes_step1_idx][:, cells_step1_mask]
            batch_gmean_step1 = pd.Series(
                np.log10(np.maximum(_row_gmean_sparse(umi_batch_step1, gmean_eps), 1e-9)),
                index=genes_step1_arr,
            )
            if np.isinf(batch_gmean_step1.to_numpy()).any():
                finite = batch_gmean_step1[np.isfinite(batch_gmean_step1)]
                floor = finite.min() if len(finite) else min_step1
                batch_gmean_step1 = batch_gmean_step1.replace(-np.inf, floor)

            sel = cell_attr.index[cell_attr[batch_var] == batch_level]
            cells_mask = cell_attr.index.isin(sel)
            batch_gmean = pd.Series(
                np.log10(np.maximum(_row_gmean_sparse(umi[:, cells_mask], gmean_eps), 1e-9)),
                index=genes,
            )
            batch_gmean = np.maximum(
                batch_gmean,
                batch_gmean_step1.reindex(genes, fill_value=min_step1),
            )
            batch_order = np.argsort(batch_gmean.to_numpy())
            batch_x = batch_gmean.to_numpy()

            for col in _batch_coef_columns(model_pars.columns, batch_var, batch_level):
                smoothed = _ksmooth(
                    batch_gmean_step1.reindex(model_pars.index),
                    model_pars[col],
                    batch_x.copy(),
                    kern=2,
                    bw=bw,
                )
                col_idx = model_pars_fit.columns.get_loc(col)
                values = model_pars_fit[col].to_numpy()
                values[batch_order] = smoothed
                model_pars_fit[col] = values

    if exclude_poisson and len(all_poisson_genes):
        if theta_regularization == "theta":
            poisson_disp = np.full(len(all_poisson_genes), np.inf)
        else:
            poisson_disp = np.zeros(len(all_poisson_genes))
        model_pars_fit.loc[all_poisson_genes, "dispersion_par"] = poisson_disp

    if theta_regularization == "theta":
        theta = np.power(10, model_pars_fit["dispersion_par"])
    else:
        theta = np.power(10, genes_log_gmean.reindex(model_pars_fit.index)) / (
            np.power(10, model_pars_fit["dispersion_par"]) - 1
        )

    model_pars_fit = model_pars_fit.drop(columns=["dispersion_par"])
    model_pars_fit["theta"] = theta.to_numpy()

    if exclude_poisson and vst_out_offset is not None and len(all_poisson_genes):
        shared_cols = [c for c in model_pars_fit.columns if c in vst_out_offset.columns]
        model_pars_fit.loc[all_poisson_genes, shared_cols] = vst_out_offset.loc[
            all_poisson_genes, shared_cols
        ]

    if fix_slope and "log_umi" in model_pars_fit.columns:
        model_pars_fit["log_umi"] = LOG10_SLOPE

    return model_pars_fit, outliers_series


def _get_model_pars_nonreg(
    genes: pd.Index,
    bin_size: int,
    model_pars_fit: pd.DataFrame,
    regressor_data: pd.DataFrame,
    umi: csr_matrix,
    model_str_nonreg: str,
    cell_attr: pd.DataFrame,
) -> pd.DataFrame:
    formula = _get_model_formula(model_str_nonreg)
    design = dmatrix(formula, cell_attr, return_type="dataframe")
    bin_ind = np.ceil(np.arange(1, len(genes) + 1) / bin_size)
    max_bin = int(np.max(bin_ind))
    genes_arr = genes.to_numpy()
    blocks: list[pd.DataFrame] = []

    for bin_idx in range(1, max_bin + 1):
        genes_bin = genes_arr[bin_ind == bin_idx]
        mu = _predict_mu(model_pars_fit, genes_bin, regressor_data)
        umi_bin = umi[np.isin(genes_arr, genes_bin), :].toarray()
        rows: list[dict[str, float]] = []
        for row_idx, gene in enumerate(genes_bin):
            y = umi_bin[row_idx]
            theta = float(model_pars_fit.loc[gene, "theta"])
            offset = np.log(np.maximum(mu[row_idx], 1e-9))
            family = sm.families.NegativeBinomial(alpha=1.0 / max(theta, 1e-7), link=sm.families.links.Log())
            try:
                fit = sm.GLM(y, design, family=family, offset=offset).fit()
                coefs = {design.columns[i]: float(fit.params[i]) for i in range(len(design.columns))}
            except (np.linalg.LinAlgError, ValueError, FloatingPointError):
                coefs = {col: np.nan for col in design.columns}
            rows.append(coefs)
        blocks.append(pd.DataFrame(rows, index=genes_bin, columns=design.columns))

    return pd.concat(blocks)


def _prepare_regressor_data(
    vst_out: dict[str, Any],
    cell_attr: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    cell_attr = cell_attr if cell_attr is not None else vst_out["cell_attr"]
    regressor_data = _build_design_matrix(cell_attr, vst_out["model_str"])
    model_pars_nonreg = vst_out.get("model_pars_nonreg")
    model_str_nonreg = vst_out.get("model_str_nonreg")
    if isinstance(model_pars_nonreg, pd.DataFrame) and model_str_nonreg:
        regressor_data = pd.concat(
            [regressor_data, _build_design_matrix(cell_attr, model_str_nonreg)],
            axis=1,
        )
    return regressor_data


def _combined_model_pars(vst_out: dict[str, Any]) -> pd.DataFrame:
    model_pars = vst_out["model_pars_fit"].copy()
    model_pars_nonreg = vst_out.get("model_pars_nonreg")
    if isinstance(model_pars_nonreg, pd.DataFrame) and not model_pars_nonreg.empty:
        model_pars = pd.concat([model_pars, model_pars_nonreg], axis=1)
    return model_pars


def _get_correct_data(
    coefs: pd.DataFrame,
    regressor_data: pd.DataFrame,
    theta: pd.Series,
    pearson_residuals: pd.DataFrame,
) -> np.ndarray:
    mu = np.exp(np.dot(coefs.to_numpy(dtype=np.float64), regressor_data.T.to_numpy(dtype=np.float64)))
    variance = mu + np.divide(mu**2, theta.to_numpy().reshape(-1, 1))
    return mu + pearson_residuals.to_numpy() * np.sqrt(variance)


def _multi_correct_data(
    vst_out: dict[str, Any],
    genes: np.ndarray,
    bin_ind: np.ndarray,
    residuals: pd.DataFrame,
    bin_idx: int,
    regressor_data: pd.DataFrame,
) -> pd.DataFrame:
    genes_bin = genes[bin_ind == bin_idx]
    model_pars = _combined_model_pars(vst_out)
    coef_cols = _coef_columns(model_pars, regressor_data)
    coefs = model_pars.loc[genes_bin, coef_cols]
    theta = model_pars.loc[genes_bin, "theta"]
    return pd.DataFrame(
        _get_correct_data(coefs, regressor_data, theta, residuals.loc[genes_bin]),
        index=genes_bin,
        columns=regressor_data.index,
    )


def _correct_counts(
    vst_out: dict[str, Any],
    genes: np.ndarray,
    bin_size: int = 500,
    do_round: bool = True,
) -> csr_matrix:
    cell_attr = vst_out["cell_attr"].copy()
    if "log_umi" in cell_attr.columns:
        cell_attr["log_umi"] = float(np.median(cell_attr["log_umi"]))
    regressor_data = _prepare_regressor_data(vst_out, cell_attr)
    bin_ind = np.ceil(np.arange(1, len(genes) + 1) / bin_size)
    max_bin = int(np.max(bin_ind))
    corrected = pd.concat(
        Parallel(n_jobs=cpu_count(), backend="threading")(
            delayed(_multi_correct_data)(
                vst_out, genes, bin_ind, vst_out["y"], i, regressor_data
            )
            for i in range(1, max_bin + 1)
        )
    )
    if do_round:
        corrected = np.round(corrected.to_numpy(), 0)
    corrected[corrected < 0] = 0
    return csr_matrix(corrected)


def _clip_matrix_values(mat: np.ndarray, clip_range: tuple[float, float]) -> np.ndarray:
    out = mat.copy()
    out[out < clip_range[0]] = clip_range[0]
    out[out > clip_range[1]] = clip_range[1]
    return out


def get_residuals(
    vst_out: dict[str, Any],
    umi: csr_matrix,
    genes: pd.Index,
    umi_genes: pd.Index,
    *,
    residual_type: str = "pearson",
    res_clip_range: Optional[tuple[float, float]] = None,
    min_variance: float = -np.inf,
    cell_attr: Optional[pd.DataFrame] = None,
    bin_size: int = 256,
) -> pd.DataFrame:
    """Compute Pearson or deviance residuals for genes using a fitted VST model.

    By default returns UNCLIPPED residuals (matching R ``sctransform::get_residuals``).
    Pass ``res_clip_range`` to clip (e.g. for scale.data storage).
    """
    if residual_type not in RESIDUAL_TYPES - {"none"}:
        raise ValueError(f"Unsupported residual_type: {residual_type}")

    regressor_data = _prepare_regressor_data(vst_out, cell_attr)
    model_pars = _combined_model_pars(vst_out)
    target_genes = pd.Index(
        [g for g in genes if g in model_pars.index and g in umi_genes]
    )
    umi_sub = umi[np.isin(umi_genes, target_genes), :]
    genes_arr = umi_genes[np.isin(umi_genes, target_genes)].to_numpy()
    bin_ind = np.ceil(np.arange(1, len(genes_arr) + 1) / bin_size)
    max_bin = int(np.max(bin_ind))
    blocks = Parallel(n_jobs=cpu_count(), backend="threading")(
        delayed(_multi_compute_residual)(
            i,
            model_pars,
            regressor_data,
            umi_sub,
            genes_arr,
            bin_ind,
            residual_type,
            min_variance,
        )
        for i in range(1, max_bin + 1)
    )
    residuals = pd.DataFrame(np.vstack(blocks), index=target_genes, columns=regressor_data.index)
    if res_clip_range is not None:
        return pd.DataFrame(
            _clip_matrix_values(residuals.to_numpy(), res_clip_range),
            index=target_genes,
            columns=regressor_data.index,
        )
    return residuals


def get_residual_var(
    vst_out: dict[str, Any],
    umi: csr_matrix,
    genes: pd.Index,
    umi_genes: pd.Index,
    *,
    residual_type: str = "pearson",
    res_clip_range: Optional[tuple[float, float]] = None,
    min_variance: float = -np.inf,
    cell_attr: Optional[pd.DataFrame] = None,
    bin_size: int = 256,
) -> pd.Series:
    """Return per-gene residual variance without materializing the full matrix."""
    residuals = get_residuals(
        vst_out,
        umi,
        genes=genes,
        umi_genes=umi_genes,
        residual_type=residual_type,
        res_clip_range=res_clip_range,
        min_variance=min_variance,
        cell_attr=cell_attr,
        bin_size=bin_size,
    )
    return residuals.var(axis=1)


def _load_reference_model(reference_sct_model: dict[str, Any]) -> dict[str, Any]:
    model = dict(reference_sct_model)
    pars = model.get("model_pars_fit")
    if not isinstance(pars, pd.DataFrame):
        if isinstance(pars, dict) and "index" in pars:
            model["model_pars_fit"] = pd.DataFrame(pars["data"], index=pars["index"], columns=pars["columns"])
        else:
            model["model_pars_fit"] = pd.DataFrame(pars)
    gene_attr = model.get("gene_attr")
    if not isinstance(gene_attr, pd.DataFrame):
        model["gene_attr"] = pd.DataFrame.from_dict(gene_attr, orient="index")
    return model


def pack_sct_model(vst_out: dict[str, Any], clip_range: Optional[tuple[float, float]] = None) -> dict[str, Any]:
    """Serialize a VST result for reuse as a reference SCT model."""
    packed = {
        "model_str": vst_out["model_str"],
        "model_pars_fit": vst_out["model_pars_fit"].copy(),
        "gene_attr": vst_out["gene_attr"].copy(),
        "arguments": {
            "sct.clip.range": clip_range,
            "sct.method": "default",
            "min_variance": vst_out.get("arguments", {}).get("min_variance", -np.inf),
        },
    }
    if vst_out.get("model_str_nonreg"):
        packed["model_str_nonreg"] = vst_out["model_str_nonreg"]
    if isinstance(vst_out.get("model_pars_nonreg"), pd.DataFrame):
        packed["model_pars_nonreg"] = vst_out["model_pars_nonreg"].copy()
    return packed


def _regress_out_matrix(
    data_expr: pd.DataFrame,
    latent_data: pd.DataFrame,
    *,
    model_use: str = "linear",
) -> pd.DataFrame:
    """Regression residuals (Seurat ScaleData / RegressOutMatrix)."""
    if model_use not in SCALE_MODELS:
        raise ValueError(f"model_use must be one of {sorted(SCALE_MODELS)}.")
    cells = data_expr.columns
    latent_data = latent_data.reindex(cells)
    vars_to_regress = latent_data.columns.tolist()
    out = np.empty(data_expr.shape, dtype=np.float64)

    if model_use == "linear":
        x = np.column_stack([np.ones(len(cells)), latent_data.to_numpy(dtype=np.float64)])
        for i, gene in enumerate(data_expr.index):
            y = data_expr.loc[gene].to_numpy(dtype=np.float64)
            coeffs = np.linalg.lstsq(x, y, rcond=None)[0]
            out[i] = y - x @ coeffs
        return pd.DataFrame(out, index=data_expr.index, columns=cells)

    fmla_rhs = " + ".join(vars_to_regress)
    for i, gene in enumerate(data_expr.index):
        y = data_expr.loc[gene].to_numpy(dtype=np.float64)
        regression_mat = pd.DataFrame(latent_data.to_numpy(dtype=np.float64), columns=vars_to_regress)
        regression_mat["GENE"] = y
        if model_use == "poisson":
            try:
                fit = sm.GLM(
                    y,
                    sm.add_constant(latent_data.to_numpy(dtype=np.float64), has_constant="add"),
                    family=sm.families.Poisson(),
                ).fit()
                out[i] = fit.resid_pearson
            except (np.linalg.LinAlgError, ValueError, FloatingPointError):
                scaled = np.log(y + 1.0)
                out[i] = (scaled - scaled.mean()) / (scaled.std() or 1.0)
        else:
            out[i] = _nb_regress_out_gene(fmla_rhs, regression_mat, gene)
    return pd.DataFrame(out, index=data_expr.index, columns=cells)


def _nb_regress_out_gene(formula_rhs: str, regression_mat: pd.DataFrame, gene: Any) -> np.ndarray:
    y = regression_mat["GENE"].to_numpy(dtype=np.float64)
    x = sm.add_constant(regression_mat.drop(columns=["GENE"]).to_numpy(dtype=np.float64), has_constant="add")
    try:
        fit = sm.NegativeBinomial(y, x).fit(disp=0, maxiter=100)
        return fit.resid_pearson
    except (np.linalg.LinAlgError, ValueError, FloatingPointError):
        scaled = np.log(y + 1.0)
        return (scaled - scaled.mean()) / (scaled.std() or 1.0)


def _scale_data(
    residuals: pd.DataFrame,
    features: Optional[np.ndarray] = None,
    vars_to_regress: Optional[Sequence[str]] = None,
    latent_data: Optional[pd.DataFrame] = None,
    do_center: bool = True,
    do_scale: bool = False,
    scale_max: float = np.inf,
    block_size: int = 750,
    model_use: str = "linear",
    count_data: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    if features is not None:
        residuals = residuals.loc[features]
    genes = residuals.index.to_numpy()
    mat = residuals.copy()
    if vars_to_regress:
        if latent_data is None:
            raise ValueError("latent_data is required when vars_to_regress is set.")
        missing = [v for v in vars_to_regress if v not in latent_data.columns]
        if missing:
            raise ValueError(f"vars_to_regress not found in latent_data: {missing}")
        regress_expr = mat
        if model_use != "linear":
            if count_data is None:
                raise ValueError("count_data is required when model_use is poisson or negbinom.")
            regress_expr = count_data.loc[mat.index, mat.columns]
        mat = _regress_out_matrix(regress_expr, latent_data[vars_to_regress], model_use=model_use)

    scaled = pd.DataFrame(
        np.zeros(mat.shape, dtype=np.float64),
        index=genes,
        columns=mat.columns,
    )
    max_block = int(np.ceil(len(genes) / block_size))
    for block in range(max_block):
        start = block * block_size
        stop = min((block + 1) * block_size, len(genes))
        block_genes = genes[start:stop]
        block_mat = mat.loc[block_genes].to_numpy()
        if do_center:
            block_mat = block_mat - block_mat.mean(axis=1, keepdims=True)
        if do_scale:
            std = block_mat.std(axis=1, keepdims=True)
            std[std == 0] = 1.0
            block_mat = block_mat / std
        if np.isfinite(scale_max):
            block_mat = np.clip(block_mat, -scale_max, scale_max)
        scaled.loc[block_genes] = block_mat
    return scaled


def vst(
    umi: csr_matrix,
    genes: pd.Index,
    cells: pd.Index,
    *,
    cell_attr: Optional[pd.DataFrame] = None,
    latent_var: Optional[Sequence[str]] = None,
    batch_var: Optional[str] = None,
    latent_var_nonreg: Optional[Sequence[str]] = None,
    n_genes: int = 2000,
    n_cells: Optional[int] = 5000,
    min_cells: int = 5,
    bin_size: int = 500,
    residual_type: str = "pearson",
    do_regularize: bool = True,
    theta_regularization: str = "od_factor",
    return_corrected_umi: bool = True,
    res_clip_range: Optional[tuple[float, float]] = None,
    bw_adjust: float = 3.0,
    gmean_eps: float = 1.0,
    seed: int = 1448145,
    min_variance: Union[float, str] = -np.inf,
    vst_flavor: Optional[str] = None,
    method: str = "poisson",
    exclude_poisson: bool = False,
    cells_step1: Optional[Sequence[str]] = None,
    genes_step1: Optional[Sequence[str]] = None,
    model_pars_step1: Optional[pd.DataFrame] = None,
    model_pars_fit_override: Optional[pd.DataFrame] = None,
    genes_log_gmean: Optional[pd.Series] = None,
    genes_log_gmean_step1: Optional[pd.Series] = None,
) -> dict[str, Any]:
    """
    Variance stabilizing transform (sctransform::vst default path).

    Parameters
    ----------
    umi
        UMI count matrix, genes x cells.
    genes, cells
        Index labels for rows and columns.
    cell_attr
        Optional cell metadata aligned to ``cells``.
    latent_var
        Covariates for the regularized NB regression (default ``log_umi``).
    batch_var
        Batch column in ``cell_attr``; fits batch-specific coefficients.
    latent_var_nonreg
        Additional covariates fit without regularization (GLM offset step).
    vst_flavor
        When ``"v2"``, use offset NB model (``nb_offset``), ``exclude_poisson=True``,
        and ``min_variance='umi_median'`` (Seurat/SCT v2 defaults).
    method
        Initial parameter estimation method: ``poisson`` or ``nb_offset``.
    exclude_poisson
        Exclude low-dispersion genes from regularization (SCT v2).
    model_pars_step1
        Precomputed step-1 model parameters (e.g. R ``vst`` export). Skips offset
        fitting; runs Python ksmooth + residuals only.
    model_pars_fit_override
        Precomputed regularized parameters (e.g. R ``model_pars_fit`` export).
        Skips step-1 fitting and ksmooth; computes residuals only.
    genes_log_gmean, genes_log_gmean_step1
        Optional precomputed log10 geometric means (e.g. R ``vst`` export).
        When ``model_pars_step1`` is set, pass R values for ksmooth parity.
    """
    if model_pars_step1 is not None and model_pars_fit_override is not None:
        raise ValueError("Pass only one of model_pars_step1 or model_pars_fit_override.")

    if model_pars_step1 is not None and genes_step1 is None:
        genes_step1 = pd.Index(model_pars_step1.index.astype(str))

    flavor = apply_vst_flavor_settings(
        vst_flavor,
        method=method,
        exclude_poisson=exclude_poisson,
        min_variance=min_variance,
        n_cells=n_cells,
    )
    method = flavor["method"]
    exclude_poisson = flavor["exclude_poisson"]
    fix_slope = flavor["fix_slope"]
    n_cells = flavor["n_cells"]
    min_variance = resolve_min_variance(flavor["min_variance"], umi)

    if residual_type not in RESIDUAL_TYPES:
        raise ValueError(f"Unsupported residual_type: {residual_type}")
    if res_clip_range is None:
        # Seurat SCTransform default: sct.clip.range = c(-sqrt(ncol/30), sqrt(ncol/30))
        res_clip_range = (-np.sqrt(umi.shape[1] / 30), np.sqrt(umi.shape[1] / 30))

    model_str = _build_model_str(latent_var or ("log_umi",), batch_var)
    cell_attr = _make_cell_attr(
        umi,
        cells,
        cell_attr=cell_attr,
        latent_var=latent_var,
        batch_var=batch_var,
        latent_var_nonreg=latent_var_nonreg,
    )
    genes_cell_count = np.asarray((umi >= 0.01).sum(axis=1)).ravel()
    keep_genes = genes_cell_count >= min_cells
    umi = umi[keep_genes]
    genes = genes[keep_genes]

    if genes_log_gmean is not None:
        genes_log_gmean = genes_log_gmean.reindex(genes)
        keep = genes_log_gmean.notna().to_numpy()
        if not keep.all():
            umi = umi[keep]
            genes = genes[keep]
            genes_log_gmean = genes_log_gmean.loc[genes]
    else:
        genes_log_gmean = pd.Series(
            np.log10(np.maximum(_row_gmean_sparse(umi, gmean_eps), 1e-9)),
            index=genes,
        )
    rng = make_vst_rng(seed)

    if cells_step1 is not None:
        cells_step1 = pd.Index(cells_step1).intersection(cells)
        if cells_step1.empty:
            raise ValueError("cells_step1 has no overlap with cells.")
    elif n_cells and n_cells < umi.shape[1]:
        cells_step1 = pd.Index(
            rng.choice(cells.to_numpy(), size=n_cells, replace=False)
        )
    else:
        cells_step1 = cells

    if genes_step1 is not None:
        genes_step1 = pd.Index(genes_step1).intersection(genes)
        if genes_step1.empty:
            raise ValueError("genes_step1 has no overlap with genes.")
        genes_log_gmean_step1 = genes_log_gmean.loc[genes_step1]
    else:
        genes_step1_mask = (
            np.asarray((umi[:, cells.isin(cells_step1)] > 0).sum(axis=1)).ravel() >= min_cells
        )
        genes_step1 = genes[genes_step1_mask]

        if exclude_poisson:
            amean = np.asarray(umi.mean(axis=1)).ravel()
            var = np.asarray(umi.power(2).mean(axis=1)).ravel() - amean**2
            overdispersed = genes_step1[
                (var[genes.get_indexer(genes_step1)] - amean[genes.get_indexer(genes_step1)]) > 0
            ]
            genes_step1 = pd.Index(overdispersed)

        if n_genes and n_genes < len(genes_step1):
            genes_log_gmean_step1 = genes_log_gmean.loc[genes_step1]
            sampling_prob = _dds(genes_log_gmean_step1)
            genes_step1 = pd.Index(
                rng.choice(
                    genes_step1.to_numpy(),
                    size=n_genes,
                    replace=False,
                    p=sampling_prob,
                )
            )
            genes_log_gmean_step1 = genes_log_gmean.loc[genes_step1]
        else:
            genes_log_gmean_step1 = genes_log_gmean.loc[genes_step1]

    if model_pars_fit_override is not None:
        from ._r_sctransform import normalize_r_model_pars

        model_pars_fit = normalize_r_model_pars(model_pars_fit_override)
        model_pars_fit = model_pars_fit.reindex(genes)
        fit_genes = model_pars_fit.index[model_pars_fit.notna().any(axis=1)]
        keep_mask = np.asarray(genes.isin(fit_genes), dtype=bool)
        genes = genes[keep_mask]
        umi = umi[keep_mask]
        model_pars_fit = model_pars_fit.loc[genes]
        model_pars = None
        outliers = None
    elif model_pars_step1 is not None:
        from ._r_sctransform import align_r_gene_list, align_r_model_pars, normalize_r_model_pars

        if cells_step1 is None or genes_step1 is None:
            raise ValueError(
                "cells_step1 and genes_step1 are required when model_pars_step1 is set."
            )
        genes_step1 = pd.Index(align_r_gene_list(list(genes_step1), genes)).intersection(genes)
        model_pars = align_r_model_pars(normalize_r_model_pars(model_pars_step1), genes)
        model_pars = model_pars.reindex(genes_step1.intersection(model_pars.index))
        model_pars = model_pars.loc[model_pars["theta"].notna()]
        genes_step1 = model_pars.index
        if genes_log_gmean_step1 is not None:
            genes_log_gmean_step1 = genes_log_gmean_step1.reindex(genes_step1)
            if genes_log_gmean_step1.isna().any():
                raise ValueError(
                    "genes_log_gmean_step1 missing values for some step1 genes; "
                    "pass R-exported genes_log_gmean_step1.csv for parity."
                )
        else:
            genes_log_gmean_step1 = genes_log_gmean.loc[genes_step1]
        model_pars["theta"] = model_pars["theta"].clip(lower=1e-7)
        model_pars.loc[~np.isfinite(model_pars["theta"]), "theta"] = np.inf
        if do_regularize:
            model_pars_fit, outliers = _reg_model_pars(
                model_pars,
                genes_log_gmean_step1,
                genes_log_gmean,
                bw_adjust=bw_adjust,
                theta_regularization=theta_regularization,
                cell_attr=cell_attr,
                batch_var=batch_var,
                cells_step1=cells_step1,
                genes_step1=genes_step1,
                umi=umi,
                gmean_eps=gmean_eps,
                exclude_poisson=exclude_poisson,
                fix_slope=fix_slope,
                genes=genes,
            )
        else:
            model_pars_fit, outliers = model_pars.reindex(genes), None
    else:
        data_step1 = cell_attr.loc[cells_step1]
        regressor_data_step1 = _build_design_matrix(data_step1, model_str)
        cells_step1_idx = cells.get_indexer(cells_step1)
        if (cells_step1_idx < 0).any():
            missing = cells_step1[cells_step1_idx < 0]
            raise ValueError(f"cells_step1 not found in cells: {missing[:5].tolist()}")
        umi_step1 = umi[genes.get_indexer(genes_step1)][:, cells_step1_idx]
        if method == "poisson":
            model_pars = _fit_poisson(umi_step1, regressor_data_step1, genes_step1)
        elif method in {"nb_offset", "glmGamPoi_offset"}:
            model_pars = fit_offset_model(
                method,
                umi_step1,
                regressor_data_step1,
                genes_step1,
                allow_inf_theta=exclude_poisson,
            )
        else:
            raise ValueError(f"Unsupported VST method: {method}")
        if exclude_poisson:
            model_pars = mark_suspicious_theta(model_pars, umi, genes)

        coef_valid = model_pars.drop(columns=["theta"]).notna().all(axis=1)
        valid = coef_valid & (model_pars["theta"].notna())
        if (~valid).any():
            warnings.warn(
                f"Dropping {(~valid).sum()} gene(s) with failed Poisson regression.",
                stacklevel=2,
            )
            model_pars = model_pars.loc[valid]
            genes_step1 = model_pars.index
            genes_log_gmean_step1 = genes_log_gmean_step1.loc[genes_step1]

        model_pars["theta"] = model_pars["theta"].clip(lower=1e-7)
        model_pars.loc[~np.isfinite(model_pars["theta"]), "theta"] = np.inf

        if do_regularize:
            model_pars_fit, outliers = _reg_model_pars(
                model_pars,
                genes_log_gmean_step1,
                genes_log_gmean,
                bw_adjust=bw_adjust,
                theta_regularization=theta_regularization,
                cell_attr=cell_attr,
                batch_var=batch_var,
                cells_step1=cells_step1,
                genes_step1=genes_step1,
                umi=umi,
                gmean_eps=gmean_eps,
                exclude_poisson=exclude_poisson,
                fix_slope=fix_slope,
                genes=genes,
            )
        else:
            model_pars_fit, outliers = model_pars.reindex(genes), None

    regressor_data = _build_design_matrix(cell_attr, model_str)
    model_str_nonreg = ""
    model_pars_nonreg: Optional[pd.DataFrame] = None
    if latent_var_nonreg:
        model_str_nonreg = _build_model_str(latent_var_nonreg, batch_var)
        model_pars_nonreg = _get_model_pars_nonreg(
            genes,
            bin_size,
            model_pars_fit,
            regressor_data,
            umi,
            model_str_nonreg,
            cell_attr,
        )
        regressor_data = pd.concat(
            [regressor_data, _build_design_matrix(cell_attr, model_str_nonreg)],
            axis=1,
        )

    model_pars_final = (
        pd.concat([model_pars_fit, model_pars_nonreg], axis=1)
        if isinstance(model_pars_nonreg, pd.DataFrame)
        else model_pars_fit
    )

    residuals: Optional[pd.DataFrame] = None
    if residual_type and residual_type != "none":
        bin_ind = np.ceil(np.arange(1, len(genes) + 1) / bin_size)
        max_bin = int(np.max(bin_ind))
        res_blocks = Parallel(n_jobs=cpu_count(), backend="threading")(
            delayed(_multi_compute_residual)(
                i,
                model_pars_final,
                regressor_data,
                umi,
                genes.to_numpy(),
                bin_ind,
                residual_type,
                min_variance,
            )
            for i in range(1, max_bin + 1)
        )
        residuals = pd.DataFrame(np.vstack(res_blocks), index=genes, columns=regressor_data.index)
        residuals = residuals.clip(lower=res_clip_range[0], upper=res_clip_range[1])

    vst_out: dict[str, Any] = {
        "y": residuals,
        "model_str": model_str,
        "model_pars": model_pars,
        "model_outlier": outliers,
        "model_pars_fit": model_pars_fit,
        "model_str_nonreg": model_str_nonreg,
        "model_pars_nonreg": model_pars_nonreg,
        "genes_log_gmean_step1": genes_log_gmean_step1,
        "cells_step1": cells_step1,
        "cell_attr": cell_attr,
        "umi_genes": genes,
        "umi_cells": regressor_data.index,
        "batch_var": batch_var,
        "latent_var": _as_list(latent_var, ("log_umi",)),
        "latent_var_nonreg": list(latent_var_nonreg) if latent_var_nonreg else None,
        "arguments": {
            "min_variance": min_variance,
            "vst_flavor": vst_flavor,
            "method": method,
            "exclude_poisson": exclude_poisson,
        },
    }

    if return_corrected_umi and residual_type == "pearson":
        vst_out["umi_corrected"] = _correct_counts(vst_out, genes.to_numpy(), bin_size=bin_size)

    gene_attr = pd.DataFrame(index=genes)
    gene_attr["detection_rate"] = np.asarray((umi >= 0.01).sum(axis=1)).ravel() / umi.shape[1]
    gene_attr["gmean"] = np.power(10, genes_log_gmean)
    umi_mean = np.asarray(umi.mean(axis=1)).ravel()
    gene_attr["amean"] = umi_mean
    gene_attr["variance"] = np.asarray(umi.power(2).mean(axis=1)).ravel() - umi_mean**2
    if residuals is not None:
        gene_attr["residual_mean"] = residuals.mean(axis=1)
        gene_attr["residual_variance"] = residuals.var(axis=1)
    else:
        gene_attr["residual_mean"] = np.nan
        gene_attr["residual_variance"] = np.nan
    vst_out["gene_attr"] = gene_attr
    return vst_out


def _select_top_features(
    gene_attr: pd.DataFrame,
    variable_features_n: Optional[int],
    variable_features_rv_th: float,
) -> pd.Index:
    feature_variance = gene_attr["residual_variance"].sort_values(ascending=False)
    feature_variance = feature_variance.dropna()
    if variable_features_n is not None:
        return feature_variance.index[: min(variable_features_n, len(feature_variance))]
    return feature_variance.index[feature_variance >= variable_features_rv_th]


def _resolve_sct_method(
    reference_sct_model: Optional[dict[str, Any]],
    residual_features: Optional[Sequence[str]],
    conserve_memory: bool,
) -> str:
    if reference_sct_model is not None:
        return "reference.model"
    if residual_features is not None:
        return "residual.features"
    if conserve_memory:
        return "conserve.memory"
    return "default"


def run_sctransform(
    umi: ArrayLike,
    genes: Sequence[str],
    cells: Sequence[str],
    *,
    variable_features_n: Optional[int] = 3000,
    variable_features_rv_th: float = 1.3,
    n_cells: int = 5000,
    n_genes: int = 2000,
    do_correct_umi: bool = True,
    do_center: bool = True,
    do_scale: bool = False,
    clip_range: Optional[tuple[float, float]] = None,
    res_clip_range: Optional[tuple[float, float]] = None,
    return_only_var_genes: bool = True,
    residual_features: Optional[Sequence[str]] = None,
    conserve_memory: bool = False,
    reference_sct_model: Optional[dict[str, Any]] = None,
    vars_to_regress: Optional[Sequence[str]] = None,
    latent_data: Optional[pd.DataFrame] = None,
    cell_attr: Optional[pd.DataFrame] = None,
    latent_var: Optional[Sequence[str]] = None,
    batch_var: Optional[str] = None,
    latent_var_nonreg: Optional[Sequence[str]] = None,
    residual_type: str = "pearson",
    model_use: str = "linear",
    seed: int = 1448145,
    vst_flavor: Optional[str] = None,
    method: str = "poisson",
    exclude_poisson: bool = False,
    min_variance: Union[float, str] = -np.inf,
    cells_step1: Optional[Sequence[str]] = None,
    genes_step1: Optional[Sequence[str]] = None,
    backend: str = "python",
    r_vst_export_dir: Optional[Union[str, Path]] = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """
    Run the Seurat SCTransform workflow on a count matrix.

    Supports default, conserve.memory, residual.features, and reference.model
    modes from Seurat 4.4.0.
    """
    gene_index = pd.Index(genes)
    full_gene_index = gene_index.copy()
    cell_index = pd.Index(cells)
    umi_gbc = _as_genes_by_cells(umi, genes, cells)

    if clip_range is None:
        clip_range = (-np.sqrt(umi_gbc.shape[1] / 30), np.sqrt(umi_gbc.shape[1] / 30))
    if res_clip_range is None:
        # Seurat SCTransform default: sct.clip.range = c(-sqrt(ncol/30), sqrt(ncol/30))
        res_clip_range = (-np.sqrt(umi_gbc.shape[1] / 30), np.sqrt(umi_gbc.shape[1] / 30))

    if model_use not in SCALE_MODELS:
        raise ValueError(f"model_use must be one of {sorted(SCALE_MODELS)}.")
    if reference_sct_model is not None and (latent_var or batch_var or latent_var_nonreg):
        raise ValueError("Custom latent variables are not supported with reference_sct_model.")

    sct_method = _resolve_sct_method(reference_sct_model, residual_features, conserve_memory)
    if backend not in {"python", "r", "r_step1", "r_fit"}:
        raise ValueError("backend must be 'python', 'r', 'r_step1', or 'r_fit'.")
    r_export: Optional[dict[str, Any]] = None
    if backend in {"r_step1", "r_fit"}:
        if r_vst_export_dir is None:
            raise ValueError(f"backend='{backend}' requires r_vst_export_dir.")
        from ._r_sctransform import (
            align_r_gene_list,
            align_r_gmean_series,
            align_r_model_pars,
            load_r_vst_export,
            validate_r_vst_export_artifacts,
        )

        validate_r_vst_export_artifacts(r_vst_export_dir, backend=backend)
        r_export = load_r_vst_export(r_vst_export_dir)
        if backend == "r_step1":
            if r_export.get("genes_log_gmean") is None or r_export.get("genes_log_gmean_step1") is None:
                raise ValueError(
                    f"backend='r_step1' requires genes_log_gmean.csv and "
                    f"genes_log_gmean_step1.csv in {r_vst_export_dir}."
                )
        if backend == "r_step1" and r_export.get("genes_vst"):
            genes_vst = pd.Index(align_r_gene_list(r_export["genes_vst"], gene_index))
            gene_pos = gene_index.get_indexer(genes_vst)
            keep = gene_pos >= 0
            gene_pos = gene_pos[keep]
            umi_gbc = umi_gbc[gene_pos]
            gene_index = gene_index[gene_pos]
        cells_step1 = r_export["cells_step1"]
        genes_step1 = align_r_gene_list(r_export["genes_step1"], gene_index)
        if r_export["meta"].get("min_variance") is not None:
            min_variance = r_export["meta"]["min_variance"]
        if vst_flavor is None:
            vst_flavor = "v2"
        if r_export.get("model_pars_step1") is not None:
            r_export["model_pars_step1"] = align_r_model_pars(
                r_export["model_pars_step1"], gene_index
            )
        if r_export.get("model_pars_fit") is not None:
            r_export["model_pars_fit"] = align_r_model_pars(
                r_export["model_pars_fit"], gene_index
            )
        if backend == "r_step1":
            if r_export.get("genes_log_gmean") is not None:
                r_export["genes_log_gmean"] = align_r_gmean_series(
                    r_export["genes_log_gmean"], gene_index
                )
            if r_export.get("genes_log_gmean_step1") is not None:
                r_export["genes_log_gmean_step1"] = align_r_gmean_series(
                    r_export["genes_log_gmean_step1"], gene_index
                )
    if backend == "r":
        if sct_method != "default":
            raise ValueError(f"backend='r' does not support sct_method '{sct_method}'.")
        if vars_to_regress or batch_var or latent_var_nonreg:
            raise ValueError("backend='r' does not support vars_to_regress, batch_var, or latent_var_nonreg.")
        if cells_step1 is not None or genes_step1 is not None:
            raise ValueError(
                "backend='r' uses R subsampling; do not pass cells_step1/genes_step1."
            )
        from ._r_sctransform import run_vst_r, vst_r_available

        if not vst_r_available():
            raise RuntimeError(
                "backend='r' requires R with sctransform. Set TRACKCELL_RSCRIPT or install sctransform."
            )
        cell_attr_r = _make_cell_attr(
            umi_gbc,
            cell_index,
            cell_attr=cell_attr,
            latent_var=latent_var,
        )
        return run_vst_r(
            umi_gbc,
            gene_index,
            cell_index,
            cell_attr=cell_attr_r,
            n_cells=n_cells,
            n_genes=n_genes,
            variable_features_n=variable_features_n,
            variable_features_rv_th=variable_features_rv_th,
            do_correct_umi=do_correct_umi,
            return_only_var_genes=return_only_var_genes,
            clip_range=clip_range,
            residual_type=residual_type,
            vst_flavor=vst_flavor,
            seed=seed,
        )

    return_only_var = return_only_var_genes
    correct_umi = do_correct_umi
    center = do_center
    vst_common = dict(
        cell_attr=cell_attr,
        latent_var=latent_var,
        batch_var=batch_var,
        latent_var_nonreg=latent_var_nonreg,
        n_genes=n_genes,
        n_cells=n_cells,
        seed=seed,
        vst_flavor=vst_flavor,
        method=method,
        exclude_poisson=exclude_poisson,
        min_variance=min_variance,
        cells_step1=cells_step1,
        genes_step1=genes_step1,
    )

    if sct_method == "reference.model":
        if reference_sct_model is None:
            raise ValueError("reference_sct_model is required for reference.model mode.")
        reference_sct_model = _load_reference_model(reference_sct_model)
        if reference_sct_model.get("model_str") != "y ~ log_umi":
            raise ValueError("reference_sct_model must use model_str 'y ~ log_umi'.")
        center = False
        correct_umi = False
        clip_range = reference_sct_model.get("arguments", {}).get("sct.clip.range") or clip_range
        vst_out = {
            "model_str": reference_sct_model["model_str"],
            "model_pars_fit": reference_sct_model["model_pars_fit"].copy(),
            "gene_attr": reference_sct_model["gene_attr"].copy(),
            "cell_attr": _make_cell_attr(umi_gbc, cell_index, cell_attr=cell_attr),
            "y": None,
            "model_str_nonreg": reference_sct_model.get("model_str_nonreg", ""),
            "model_pars_nonreg": reference_sct_model.get("model_pars_nonreg"),
            "arguments": reference_sct_model.get("arguments", {"min_variance": -np.inf}),
        }
        common_genes = gene_index.intersection(vst_out["model_pars_fit"].index)
        vst_out["model_pars_fit"] = vst_out["model_pars_fit"].loc[common_genes]
        vst_out["gene_attr"] = vst_out["gene_attr"].loc[common_genes]
    elif sct_method == "residual.features":
        return_only_var = True
        correct_umi = False
        vst_out = vst(
            umi_gbc,
            gene_index,
            cell_index,
            **vst_common,
            return_corrected_umi=False,
            residual_type="none",
        )
        vst_out["gene_attr"]["residual_variance"] = np.nan
    elif sct_method == "conserve.memory":
        return_only_var = True
        vst_out = vst(
            umi_gbc,
            gene_index,
            cell_index,
            **vst_common,
            return_corrected_umi=False,
            residual_type="none",
        )
        feature_var = get_residual_var(
            vst_out,
            umi_gbc,
            gene_index,
            gene_index,
            residual_type=residual_type,
            res_clip_range=res_clip_range,
        )
        vst_out["gene_attr"]["residual_variance"] = np.nan
        vst_out["gene_attr"].loc[feature_var.index, "residual_variance"] = feature_var
    else:
        extra_vst: dict[str, Any] = {}
        if backend == "r_step1" and r_export is not None:
            extra_vst["model_pars_step1"] = r_export["model_pars_step1"]
            if r_export.get("genes_log_gmean") is not None:
                extra_vst["genes_log_gmean"] = r_export["genes_log_gmean"]
            if r_export.get("genes_log_gmean_step1") is not None:
                extra_vst["genes_log_gmean_step1"] = r_export["genes_log_gmean_step1"]
        elif backend == "r_fit" and r_export is not None:
            fit = r_export.get("model_pars_fit")
            if fit is None:
                raise FileNotFoundError(
                    f"backend='r_fit' requires model_pars_fit.csv in {r_vst_export_dir}"
                )
            extra_vst["model_pars_fit_override"] = fit
        vst_out = vst(
            umi_gbc,
            gene_index,
            cell_index,
            **vst_common,
            **extra_vst,
            return_corrected_umi=correct_umi and residual_type == "pearson",
            residual_type=residual_type,
            res_clip_range=res_clip_range,
        )

    top_features = _select_top_features(
        vst_out["gene_attr"],
        variable_features_n,
        variable_features_rv_th,
    )
    residual_genes = (
        pd.Index(residual_features)
        if residual_features is not None
        else top_features
    )

    if sct_method == "reference.model":
        residual_genes = residual_genes.intersection(gene_index).intersection(
            vst_out["model_pars_fit"].index
        )
        residual_mat = get_residuals(
            vst_out,
            umi_gbc,
            residual_genes,
            gene_index,
            cell_attr=vst_out["cell_attr"],
            residual_type=residual_type,
            res_clip_range=res_clip_range,
        )
        ref_means = vst_out["gene_attr"].loc[residual_genes, "residual_mean"]
        vst_out["y"] = residual_mat.sub(ref_means, axis=0)
        vst_out["gene_attr"] = vst_out["gene_attr"].loc[residual_genes]
    elif sct_method == "residual.features":
        residual_genes = pd.Index(residual_features).intersection(gene_index).intersection(
            vst_out["gene_attr"].index
        )
        residual_mat = get_residuals(
            vst_out,
            umi_gbc,
            residual_genes,
            gene_index,
            cell_attr=vst_out["cell_attr"],
            residual_type=residual_type,
            res_clip_range=res_clip_range,
        )
        vst_out["y"] = residual_mat
        vst_out["gene_attr"].loc[residual_genes, "residual_mean"] = residual_mat.mean(axis=1)
        vst_out["gene_attr"].loc[residual_genes, "residual_variance"] = residual_mat.var(axis=1)
    elif sct_method == "conserve.memory":
        residual_mat = get_residuals(
            vst_out,
            umi_gbc,
            top_features,
            gene_index,
            residual_type=residual_type,
            res_clip_range=res_clip_range,
            cell_attr=vst_out["cell_attr"],
        )
        vst_out["y"] = residual_mat
        vst_out["gene_attr"].loc[top_features, "residual_mean"] = residual_mat.mean(axis=1)
        if correct_umi and residual_type == "pearson":
            warnings.warn(
                "do_correct_umi with conserve_memory runs an extra full residual pass.",
                stacklevel=2,
            )
            full_residuals = get_residuals(
                vst_out,
                umi_gbc,
                gene_index,
                gene_index,
                residual_type=residual_type,
                res_clip_range=res_clip_range,
                cell_attr=vst_out["cell_attr"],
            )
            vst_tmp = dict(vst_out)
            vst_tmp["y"] = full_residuals
            vst_out["umi_corrected"] = _correct_counts(
                vst_tmp, gene_index.to_numpy(), bin_size=500
            )
    else:
        if return_only_var:
            vst_out["y"] = vst_out["y"].loc[top_features]

    if vst_out["y"] is None:
        raise RuntimeError("SCTransform failed to compute residuals.")

    if correct_umi and "umi_corrected" in vst_out:
        counts = vst_out["umi_corrected"]
    elif correct_umi and sct_method == "default" and residual_type == "pearson":
        counts = vst_out.get("umi_corrected", umi_gbc)
    else:
        counts = umi_gbc

    data = counts.copy()
    data.data = np.log1p(data.data)

    count_data = None
    if vars_to_regress and model_use != "linear":
        count_data = pd.DataFrame(
            counts.toarray(),
            index=gene_index,
            columns=cell_index,
        )

    scale_data = _scale_data(
        vst_out["y"],
        features=top_features.to_numpy() if return_only_var else None,
        vars_to_regress=vars_to_regress,
        latent_data=latent_data,
        do_center=center,
        do_scale=do_scale,
        model_use=model_use,
        count_data=count_data,
    )
    scale_data = scale_data.clip(lower=clip_range[0], upper=clip_range[1])

    count_genes = pd.Index(vst_out["gene_attr"].index)
    if counts.shape[0] == len(count_genes) and len(count_genes) != len(full_gene_index):
        counts = _expand_genes_gbc(counts, count_genes, full_gene_index)
        data = _expand_genes_gbc(data, count_genes, full_gene_index)

    variable_features = (
        residual_genes.to_numpy() if residual_features is not None else top_features.to_numpy()
    )
    vst_out["arguments"] = {
        "sct.method": sct_method,
        "sct.clip.range": clip_range,
        "res_clip_range": res_clip_range,
    }

    assay_out = {
        "counts": counts,
        "data": data,
        "scale.data": scale_data,
        "variable_features": variable_features,
        "sct_method": sct_method,
        "vst_method": vst_out.get("arguments", {}).get("method"),
        "backend": backend,
    }
    return assay_out, vst_out
