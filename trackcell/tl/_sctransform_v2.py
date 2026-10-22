"""SCTransform v2 / offset-model helpers (sctransform::vst with vst.flavor='v2')."""

from __future__ import annotations

from typing import Any, Literal, Optional, Union

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.special import digamma, polygamma

LOG10_SLOPE = float(np.log(10.0))


def pyglmGamPoi_available() -> bool:
    """Return True when the native pyglmGamPoi extension is importable."""
    try:
        from pyglmGamPoi import is_available

        return bool(is_available())
    except ImportError:
        return False


V2Backend = Literal["python", "pyglmGamPoi"]


def resolve_v2_step1_method(v2_backend: V2Backend) -> str:
    """Resolve an explicitly selected SCT v2 step-1 backend."""
    if v2_backend == "python":
        return "nb_offset"
    if v2_backend == "pyglmGamPoi":
        if not pyglmGamPoi_available():
            raise RuntimeError(
                "v2_backend='pyglmGamPoi' requires the optional pyglmGamPoi package."
            )
        return "glmGamPoi_offset"
    raise ValueError("v2_backend must be 'python' or 'pyglmGamPoi'.")


def apply_vst_flavor_settings(
    vst_flavor: Optional[str],
    *,
    method: str,
    exclude_poisson: bool,
    min_variance: Union[float, str],
    n_cells: Optional[int],
    v2_backend: V2Backend,
) -> dict[str, Any]:
    """Apply Seurat/sctransform vst.flavor presets."""
    settings: dict[str, Any] = {
        "method": method,
        "exclude_poisson": exclude_poisson,
        "min_variance": min_variance,
        "n_cells": n_cells,
        "fix_slope": False,
    }
    if vst_flavor is None:
        return settings
    if vst_flavor != "v2":
        raise ValueError(f"Unsupported vst_flavor: {vst_flavor}. Only 'v2' is implemented.")
    settings.update(
        {
            "method": resolve_v2_step1_method(v2_backend),
            "exclude_poisson": True,
            "fix_slope": True,
            "min_variance": "umi_median" if min_variance == -np.inf else min_variance,
            "n_cells": 2000 if n_cells is None or n_cells == 5000 else n_cells,
        }
    )
    return settings


def _natural_log_umi(log10_umi: np.ndarray) -> np.ndarray:
    return np.log(10.0) * np.asarray(log10_umi, dtype=np.float64)


# ---------------------------------------------------------------------------
# Pure-Python NB offset fit (inspired by shanuz/sctransform.py _fit_nb_offset)
# ---------------------------------------------------------------------------


def _theta_ml_vectorized(
    y: np.ndarray,
    mu: np.ndarray,
    limit: int = 10,
    eps: float = 0.0001220703,
) -> np.ndarray:
    """Vectorized port of ``MASS::theta.ml``; rows are independent genes."""
    n = y.shape[1]
    with np.errstate(divide="ignore", invalid="ignore"):
        theta = n / np.sum((y / mu - 1.0) ** 2, axis=1)
    active = np.isfinite(theta) & (theta > 0)
    theta = np.where(active, theta, np.inf)
    for _ in range(limit - 1):
        finite = np.isfinite(theta)
        if not finite.any():
            break
        th = np.maximum(np.abs(theta[finite]), 1e-10)[:, None]
        yy = y[finite]
        mm = mu[finite]
        score = np.sum(
            digamma(th + yy)
            - digamma(th)
            + np.log(th)
            + 1.0
            - np.log(th + mm)
            - (yy + th) / (mm + th),
            axis=1,
        )
        info = np.sum(
            -polygamma(1, th + yy)
            + polygamma(1, th)
            - 1.0 / th
            + 2.0 / (mm + th)
            - (yy + th) / (mm + th) ** 2,
            axis=1,
        )
        delta = np.divide(score, info, out=np.zeros_like(score), where=info != 0)
        idx = np.flatnonzero(finite)
        theta[idx] = np.maximum(theta[idx] + delta, 0.0)
        if np.all(np.abs(delta) <= eps):
            break
    return theta


def _fit_nb_offset_py(
    y: np.ndarray,
    log10_umi: np.ndarray,
    gene_chunk: int = 500,
    n_outer: int = 10,
    n_inner: int = 25,
    tol: float = 1e-8,
    allow_inf_theta: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Pure-Python NB GLM with fixed depth slope, vectorised over gene chunks.

    Fits ``y ~ 1 + offset(log_umi)`` per gene via alternating IRLS (intercept only)
    and :func:`theta.ml` (MASS::theta.ml). The depth slope is fixed at ``log(10)``
    on the log10-UMI scale, matching Seurat 5's ``vst.flavor="v2"`` and R's
    ``glm.nb(y ~ 1 + offset(log_umi))``.

    Inspired by shanuz's ``_fit_nb_offset``, verified within 1e-4 relative of
    R 4.6.1 glm.nb on PBMC 3k.

    Parameters
    ----------
    y : (G, N_cells) float64 array
        Gene counts for the step-1 subsample.
    log10_umi : (N_cells,) float64 array
        Cell total UMI on log10 scale.
    gene_chunk : int
        Genes processed per chunk to bound memory.

    Returns
    -------
    b0 : (G,) float64 — per-gene intercept (natural-log scale).
    theta : (G,) float64 — per-gene NB dispersion (inf for Poisson genes).
    """
    G = y.shape[0]
    offset = _natural_log_umi(log10_umi)  # log(umi) on natural-log scale
    b0 = np.zeros(G)
    theta = np.full(G, np.inf)

    for start in range(0, G, gene_chunk):
        end = min(start + gene_chunk, G)
        Y = y[start:end]

        # Poisson start: closed form for an intercept-only log-link model.
        mu_offset_sum = np.exp(offset).sum()
        bb = np.log(np.maximum(Y.sum(axis=1), 1e-9)) - np.log(max(mu_offset_sum, 1e-9))
        mu = np.exp(np.clip(bb[:, None] + offset[None, :], -30, 30))
        th = _theta_ml_vectorized(Y, mu)

        for _ in range(n_outer):
            for _ in range(n_inner):
                eta = np.clip(bb[:, None] + offset[None, :], -30, 30)
                mu = np.exp(eta)
                w = mu / (1.0 + mu / np.maximum(th, 1e-10)[:, None])
                z = eta + (Y - mu) / mu
                nb = np.sum(w * (z - offset[None, :]), axis=1) / np.maximum(
                    np.sum(w, axis=1), 1e-10
                )
                if np.max(np.abs(nb - bb)) < tol:
                    break
                bb = nb
            mu = np.exp(np.clip(bb[:, None] + offset[None, :], -30, 30))
            th_new = _theta_ml_vectorized(Y, mu)
            finite = np.isfinite(th_new) & np.isfinite(th)
            if finite.any():
                change = np.max(np.abs(th_new[finite] - th[finite]))
                scale = max(1.0, float(np.max(th[finite])))
                converged = change < tol * scale
            else:
                converged = True
            th = th_new
            if converged:
                break

        if not allow_inf_theta:
            th = np.minimum(th, Y.mean(axis=1) / 1e-4)
        b0[start:end] = bb
        theta[start:end] = th

    return b0, np.maximum(theta, 1e-10)


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------


def _normalize_model_pars_columns(model_pars: pd.DataFrame) -> pd.DataFrame:
    rename = {col: col.strip("()") for col in model_pars.columns if col.startswith("(")}
    if rename:
        model_pars = model_pars.rename(columns=rename)
    if "Intercept" not in model_pars.columns and "(Intercept)" in model_pars.columns:
        model_pars = model_pars.rename(columns={"(Intercept)": "Intercept"})
    return model_pars


def fit_offset_model(
    method: str,
    umi: csr_matrix,
    regressor_data: pd.DataFrame,
    gene_index: pd.Index,
    *,
    allow_inf_theta: bool = True,
) -> pd.DataFrame:
    """Dispatch glmGamPoi_offset (pyglmGamPoi or pure-Python fallback).

    The pure-Python fallback uses :func:`_fit_nb_offset_py`, a vectorised
    NB GLM inspired by shanuz and verified against R's glm.nb.
    """
    if "log_umi" not in regressor_data.columns:
        raise ValueError(f"{method} requires log_umi in regressor_data.")
    if method == "glmGamPoi_offset":
        if not pyglmGamPoi_available():
            raise RuntimeError(
                "v2_backend='pyglmGamPoi' requires the optional pyglmGamPoi package."
            )
        from pyglmGamPoi.trackcell_compat import fit_offset_model as py_fit

        return _normalize_model_pars_columns(
            py_fit(umi, regressor_data, gene_index, allow_inf_theta=allow_inf_theta)
        )
    if method == "nb_offset":
        log10_umi = regressor_data["log_umi"].to_numpy(dtype=np.float64)
        y = umi.toarray().astype(np.float64)
        b0, theta = _fit_nb_offset_py(
            y,
            log10_umi,
            allow_inf_theta=allow_inf_theta,
        )
        return pd.DataFrame(
            {
                "theta": theta,
                "Intercept": b0,
                "log_umi": np.full(len(b0), LOG10_SLOPE),
            },
            index=gene_index,
        )
    raise ValueError(f"Unsupported offset method: {method}")


def get_nz_median2(umi: csr_matrix, gene_indices: Optional[np.ndarray] = None) -> float:
    """Median of non-zero UMIs (sctransform::get_nz_median2)."""
    if gene_indices is None:
        data = umi.data
    else:
        data = umi[gene_indices, :].data
    if data.size == 0:
        return 0.0
    return float(np.median(data))


def resolve_min_variance(min_variance: Union[float, str], umi: csr_matrix) -> float:
    if isinstance(min_variance, str):
        if min_variance == "umi_median":
            return (get_nz_median2(umi) / 5.0) ** 2
        raise ValueError(
            f"Unsupported min_variance mode: {min_variance}. "
            "Only 'umi_median' is implemented."
        )
    return float(min_variance)


def row_var_sparse(umi: csr_matrix) -> np.ndarray:
    mean = np.asarray(umi.mean(axis=1)).ravel()
    sqmean = np.asarray(umi.power(2).mean(axis=1)).ravel()
    return sqmean - mean**2


def mark_suspicious_theta(
    model_pars: pd.DataFrame,
    umi: csr_matrix,
    genes: pd.Index,
) -> pd.DataFrame:
    """Set theta=Inf when MM/MLE ratio is extreme (sctransform v2 exclude_poisson).

    Mirrors ``get_model_pars`` in R ``vst()``: ``rowMeans`` / ``row_var`` on the full
    count matrix, then compare method-of-moments theta to the fitted MLE for each
    step-1 gene in ``model_pars``.
    """
    model_pars = model_pars.copy()
    amean = pd.Series(np.asarray(umi.mean(axis=1)).ravel(), index=genes)
    var = pd.Series(row_var_sparse(umi), index=genes)
    step_genes = model_pars.index.intersection(genes)
    amean_step = amean.reindex(step_genes)
    var_step = var.reindex(step_genes)
    with np.errstate(divide="ignore", invalid="ignore"):
        predicted = amean_step**2 / (var_step - amean_step)
        actual = model_pars.loc[step_genes, "theta"].to_numpy(dtype=np.float64)
        ratio = predicted.to_numpy(dtype=np.float64) / actual
    suspicious = np.isfinite(ratio) & (ratio < 1e-3)
    if suspicious.any():
        flagged = step_genes[suspicious]
        model_pars.loc[flagged, "theta"] = np.inf
    return model_pars


def build_poisson_gene_tables(
    umi: csr_matrix,
    genes: pd.Index,
    genes_step1: pd.Index,
    model_pars: pd.DataFrame,
) -> tuple[pd.Index, pd.DataFrame]:
    """Identify Poisson genes and build offset replacement parameters."""
    amean = pd.Series(np.asarray(umi.mean(axis=1)).ravel(), index=genes)
    var = pd.Series(row_var_sparse(umi), index=genes)
    overdispersion = var - amean
    all_poisson = genes[(overdispersion <= 0) | (amean < 1e-3)]
    step1_poisson = genes_step1.intersection(all_poisson)
    inf_theta = model_pars.index[~np.isfinite(model_pars["theta"])]
    poisson_step1 = genes_step1.intersection(all_poisson.union(inf_theta))

    mean_cell_sum = float(np.asarray(umi.sum(axis=0)).mean())
    offset_rows: list[dict[str, float]] = []
    for gene in all_poisson:
        row = {
            "theta": np.inf,
            "Intercept": float(
                np.log(max(amean[gene], 1e-9)) - np.log(max(mean_cell_sum, 1e-9))
            ),
            "log_umi": LOG10_SLOPE,
            "dispersion_par": 0.0,
        }
        for col in model_pars.columns:
            if col not in row:
                row[col] = 0.0
        offset_rows.append(row)
    vst_out_offset = pd.DataFrame(offset_rows, index=all_poisson)
    for col in model_pars.columns:
        if col not in vst_out_offset.columns:
            vst_out_offset[col] = 0.0
    _ = poisson_step1
    return all_poisson, vst_out_offset
