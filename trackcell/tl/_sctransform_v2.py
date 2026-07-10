"""SCTransform v2 / offset-model helpers (sctransform::vst with vst.flavor='v2')."""

from __future__ import annotations

from typing import Any, Optional, Union

import numpy as np
import pandas as pd
import statsmodels.api as sm
from joblib import Parallel, cpu_count, delayed
from scipy.sparse import csr_matrix

from ._r_sctransform import glmGamPoi_r_available

LOG10_SLOPE = float(np.log(10.0))


def pyglmGamPoi_available() -> bool:
    """Return True when the native pyglmGamPoi extension is importable."""
    try:
        from pyglmGamPoi import is_available

        return bool(is_available())
    except ImportError:
        return False


def glmGamPoi_offset_available() -> bool:
    """Return True when glmGamPoi_offset can run (pyglmGamPoi or R bridge)."""
    return pyglmGamPoi_available() or glmGamPoi_r_available()


def resolve_v2_step1_method() -> str:
    """Return the step-1 offset model for ``vst.flavor='v2'`` (never ``nb_offset``)."""
    if glmGamPoi_offset_available():
        return "glmGamPoi_offset"
    raise RuntimeError(
        "SCT v2 requires glmGamPoi_offset (pyglmGamPoi or R glmGamPoi). "
        "Install pyglmGamPoi, configure R via TRACKCELL_RSCRIPT, or use "
        "sctransform(..., method='r_fit') / method='r' with R exports."
    )


def apply_vst_flavor_settings(
    vst_flavor: Optional[str],
    *,
    method: str,
    exclude_poisson: bool,
    min_variance: Union[float, str],
    n_cells: Optional[int],
) -> dict[str, Any]:
    """Apply Seurat/sctransform vst.flavor presets."""
    settings = {
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
    if method == "nb_offset":
        raise ValueError(
            "nb_offset is not supported with vst_flavor='v2'. "
            "Use glmGamPoi_offset (pyglmGamPoi / R) or method='r_fit' / method='r'."
        )
    v2_method = resolve_v2_step1_method()
    settings.update(
        {
            "method": v2_method,
            "exclude_poisson": True,
            "fix_slope": True,
            "min_variance": "umi_median" if min_variance == -np.inf else min_variance,
            "n_cells": 2000 if n_cells is None or n_cells == 5000 else n_cells,
        }
    )
    return settings


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
    """Dispatch nb_offset or glmGamPoi_offset (pyglmGamPoi / R) fitting."""
    if method == "glmGamPoi_offset":
        if pyglmGamPoi_available():
            from pyglmGamPoi.trackcell_compat import fit_offset_model as py_fit

            return _normalize_model_pars_columns(
                py_fit(
                    umi,
                    regressor_data,
                    gene_index,
                    allow_inf_theta=allow_inf_theta,
                )
            )
        if glmGamPoi_r_available():
            from ._r_sctransform import fit_glmGamPoi_offset_r

            model_pars = fit_glmGamPoi_offset_r(
                umi,
                regressor_data,
                gene_index,
                allow_inf_theta=allow_inf_theta,
            )
            return _normalize_model_pars_columns(model_pars)
        raise RuntimeError(
            "glmGamPoi_offset requested but neither pyglmGamPoi nor R glmGamPoi is available."
        )
    if method == "nb_offset":
        return fit_nb_offset(
            umi,
            regressor_data,
            gene_index,
            allow_inf_theta=allow_inf_theta,
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


def _natural_log_umi(log10_umi: np.ndarray) -> np.ndarray:
    return np.log(np.power(10.0, log10_umi))


def _one_row_nb_offset(
    y: np.ndarray,
    offset: np.ndarray,
    allow_inf_theta: bool,
) -> list[float]:
    from ._sctransform import _theta_ml

    x = np.ones((len(y), 1), dtype=np.float64)
    try:
        fit = sm.GLM(
            y,
            x,
            family=sm.families.NegativeBinomial(alpha=1.0),
            offset=offset,
        ).fit()
        mu = fit.predict(offset=offset)
        theta = _theta_ml(y=y, mu=mu)
        if not allow_inf_theta and np.isfinite(theta):
            theta = min(theta, float(np.mean(y)) / 1e-4)
        return [float(theta), float(fit.params[0]), LOG10_SLOPE]
    except (np.linalg.LinAlgError, ValueError, FloatingPointError):
        fit = sm.GLM(y, x, family=sm.families.Poisson(), offset=offset).fit()
        theta = np.inf
        if not allow_inf_theta:
            theta = float(np.mean(y) / 1e-4)
        return [theta, float(fit.params[0]), LOG10_SLOPE]


def fit_nb_offset(
    umi: csr_matrix,
    regressor_data: pd.DataFrame,
    gene_index: pd.Index,
    *,
    allow_inf_theta: bool = True,
) -> pd.DataFrame:
    """Native nb_offset fallback for glmGamPoi_offset (sctransform::fit_nb_offset)."""
    if "log_umi" not in regressor_data.columns:
        raise ValueError("nb_offset requires log_umi in regressor_data.")
    offset = _natural_log_umi(regressor_data["log_umi"].to_numpy(dtype=np.float64))
    columns = ["theta", "Intercept", "log_umi"]
    results = Parallel(n_jobs=cpu_count(), backend="threading")(
        delayed(_one_row_nb_offset)(row.toarray().ravel(), offset, allow_inf_theta)
        for row in umi
    )
    out = pd.DataFrame(results, columns=columns, index=gene_index)
    return out


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
            "Intercept": float(np.log(max(amean[gene], 1e-9)) - np.log(max(mean_cell_sum, 1e-9))),
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
