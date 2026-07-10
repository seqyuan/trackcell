"""Metrics for SCT benchmark comparisons."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def pearson_corr(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return float("nan")
    return float(np.corrcoef(a[mask], b[mask])[0, 1])


def compare_series(a: pd.Series, b: pd.Series, *, log10: bool = False) -> dict[str, float]:
    shared = a.index.intersection(b.index)
    x = a.loc[shared].to_numpy(dtype=np.float64)
    y = b.loc[shared].to_numpy(dtype=np.float64)
    if log10:
        x = np.log10(np.clip(x, 1e-9, None))
        y = np.log10(np.clip(y, 1e-9, None))
    diff = x - y
    return {
        "n": int(len(shared)),
        "corr": pearson_corr(x, y),
        "median_abs_diff": float(np.nanmedian(np.abs(diff))),
        "max_abs_diff": float(np.nanmax(np.abs(diff))),
    }


def compare_model_pars(py: pd.DataFrame, r: pd.DataFrame) -> dict[str, Any]:
    py = py.rename(columns={"(Intercept)": "Intercept"})
    r = r.rename(columns={"(Intercept)": "Intercept"})
    shared = py.index.intersection(r.index)
    out: dict[str, Any] = {"n_genes": int(len(shared))}
    py_t = py.loc[shared, "theta"]
    r_t = r.loc[shared, "theta"]
    out["n_inf_py"] = int((~np.isfinite(py_t)).sum())
    out["n_inf_r"] = int((~np.isfinite(r_t)).sum())
    out["n_inf_both"] = int((~np.isfinite(py_t) & ~np.isfinite(r_t)).sum())
    for col in ("theta", "Intercept", "log_umi"):
        if col in py.columns and col in r.columns:
            out[col] = compare_series(py.loc[shared, col], r.loc[shared, col], log10=(col == "theta"))
    return out


def compare_residual_matrices(py: pd.DataFrame, r: pd.DataFrame) -> dict[str, float]:
    shared_genes = py.index.intersection(r.index)
    shared_cells = py.columns.intersection(r.columns)
    if not len(shared_genes) or not len(shared_cells):
        return {"n_genes": len(shared_genes), "n_cells": len(shared_cells)}
    pm = py.loc[shared_genes, shared_cells].to_numpy(dtype=np.float64)
    rm = r.loc[shared_genes, shared_cells].to_numpy(dtype=np.float64)
    diff = pm - rm
    return {
        "n_genes": int(len(shared_genes)),
        "n_cells": int(len(shared_cells)),
        "corr": pearson_corr(pm.ravel(), rm.ravel()),
        "median_abs_diff": float(np.nanmedian(np.abs(diff))),
        "max_abs_diff": float(np.nanmax(np.abs(diff))),
    }


def compare_gene_attr(py: pd.DataFrame, r: pd.DataFrame) -> dict[str, Any]:
    out: dict[str, Any] = {}
    if "residual_variance" in py.columns and "residual_variance" in r.columns:
        rv = compare_series(py["residual_variance"], r["residual_variance"])
        rv["spearman"] = float(
            spearmanr(
                py.loc[py.index.intersection(r.index), "residual_variance"],
                r.loc[py.index.intersection(r.index), "residual_variance"],
            ).correlation
        )
        out["residual_variance"] = rv
    return out


def hvg_jaccard(py_gene_attr: pd.DataFrame, r_hvg: set[str], *, n_top: int = 3000) -> float:
    py_hvg = set(
        py_gene_attr["residual_variance"].sort_values(ascending=False).index[:n_top].astype(str)
    )
    return jaccard(py_hvg, r_hvg)


def summarize_distribution(values: list[float]) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    if arr.size == 0:
        return {}
    return {
        "n": int(arr.size),
        "min": float(arr.min()),
        "q25": float(np.quantile(arr, 0.25)),
        "median": float(np.median(arr)),
        "mean": float(arr.mean()),
        "q75": float(np.quantile(arr, 0.75)),
        "max": float(arr.max()),
    }
