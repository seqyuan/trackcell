"""Shared Tier-B benchmark loaders for SCT golden / parity tests.

Primary loader :func:`load_seurat_parity_umi` mirrors R ``export_sct_stepwise_r.R``:

- ``Read10X(..., gene.column = 2)`` → ``var_names='gene_symbols'``
- ``CreateSeuratObject(min.cells = 3, min.features = 0)``
- QC: ``nFeature_RNA`` 200–7500, ``percent.mt`` < 15
- ``RenameCells(paste0(sample, '_', bc))``

Returns genes × cells CSR **before** ``vst()``'s internal ``min_cells`` filter
(``rowSums(umi >= 0.01) >= 5``), matching the matrix passed into R ``vst()``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix, issparse

from trackcell.tl._sctransform import _as_genes_by_cells

BENCHMARK_ROOT = Path("/Volumes/process/tmp/tcl_test")
DEFAULT_SAMPLE = "GSM8779707"
DEFAULT_VAR_NAMES = "gene_symbols"

# export_sct_stepwise_r.R :: filter_params
SEURAT_FILTER_DEFAULTS = {
    "min_features": 200,
    "max_features": 7500,
    "max_percent_mt": 15.0,
    "min_cells_gene": 3,
}


@dataclass(frozen=True)
class SeuratFilterParams:
    min_features: int = 200
    max_features: int = 7500
    max_percent_mt: float = 15.0
    min_cells_gene: int = 3


def verify_seurat_parity_vs_r(
    umi: csr_matrix,
    genes: pd.Index,
    cells: pd.Index,
    cell_attr: pd.DataFrame,
    *,
    sample: str,
    benchmark_root: Path = BENCHMARK_ROOT,
) -> dict[str, float | int | bool]:
    """
    Compare a loader result to R ``export_sct_stepwise_r.R`` artifacts.

    Count matrices match exactly when R rownames are applied by position; scanpy
    ``-N`` duplicate suffixes are not order-compatible with Seurat ``.N`` names.
    """
    step_dir = benchmark_root / "steps/r_sct_stepwise" / sample
    meta = pd.read_json(step_dir / "meta.json", typ="series")
    r_cells = (benchmark_root / "steps/r/01_filter" / f"{sample}_cells.txt").read_text().strip().split("\n")
    r_cells = [c for c in r_cells if c]

    col_sums = np.asarray(umi.sum(axis=0)).ravel()
    log_umi_ref = np.log10(col_sums + 1e-9)

    out: dict[str, float | int | bool] = {
        "n_genes": int(umi.shape[0]),
        "n_cells": int(umi.shape[1]),
        "meta_n_genes": int(meta["n_genes"]),
        "meta_n_cells": int(meta["n_cells"]),
        "shape_matches_meta": umi.shape == (int(meta["n_genes"]), int(meta["n_cells"])),
        "cells_match_r": set(map(str, cells)) == set(r_cells),
        "log_umi_max_diff": float(np.max(np.abs(cell_attr["log_umi"].to_numpy() - log_umi_ref))),
    }
    r_colsums_path = benchmark_root / "steps" / "py_loader_check" / "r_colsums.txt"
    if r_colsums_path.exists():
        r_col = np.loadtxt(r_colsums_path)
        out["colsums_max_diff_vs_r"] = float(np.max(np.abs(col_sums - r_col)))
    r_genes = _load_r_reference_gene_names(benchmark_root, sample, int(umi.shape[0]))
    if r_genes is not None:
        out["gene_names_match_r_reference"] = list(map(str, genes)) == list(map(str, r_genes))
    return out


def _load_r_reference_gene_names(
    benchmark_root: Path,
    sample: str,
    n_genes: int,
) -> pd.Index | None:
    """
    Seurat ``Read10X`` rownames after ``CreateSeuratObject(min.cells=…)``.

    Duplicate symbols use ``.N`` suffixes in an order that does not match scanpy
    ``make_unique`` (``-N``). Matrix row order matches R; only names differ — assign
    R names by position when a sidecar list is available.
    """
    candidates = [
        benchmark_root / "ref" / f"{sample}_genes_min_cells.txt",
        benchmark_root / "steps" / "py_loader_check" / "r_genes.txt",
    ]
    for path in candidates:
        if not path.exists():
            continue
        names = [line for line in path.read_text().splitlines() if line]
        if len(names) == n_genes:
            return pd.Index(names)
    return None


def _percent_mt(gene_names: pd.Index, umi_gbc: csr_matrix) -> np.ndarray:
    """Seurat ``PercentageFeatureSet(pattern = '^MT-')`` on genes × cells counts."""
    mt_mask = np.asarray(gene_names.str.match(r"^MT-", case=False), dtype=bool)
    if not mt_mask.any():
        return np.zeros(umi_gbc.shape[1], dtype=np.float64)
    total = np.asarray(umi_gbc.sum(axis=0)).ravel().astype(np.float64)
    mt = np.asarray(umi_gbc[mt_mask, :].sum(axis=0)).ravel().astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        pct = np.where(total > 0, mt / total * 100.0, 0.0)
    return pct


def load_seurat_parity_umi(
    sample: str = DEFAULT_SAMPLE,
    *,
    benchmark_root: Path = BENCHMARK_ROOT,
    filter_params: SeuratFilterParams | None = None,
    var_names: str = DEFAULT_VAR_NAMES,
) -> tuple[csr_matrix, pd.Index, pd.Index, pd.DataFrame]:
    """
    Load counts aligned with Seurat ``load_filtered()`` in ``export_sct_stepwise_r.R``.

    Returns
    -------
    umi, genes, cells, cell_attr
        ``umi`` is genes × cells. ``cell_attr`` has ``log_umi = log10(nCount_RNA)``.
    """
    fp = filter_params or SeuratFilterParams(**SEURAT_FILTER_DEFAULTS)
    data_dir = benchmark_root / "data" / sample

    adata = sc.read_10x_mtx(
        data_dir,
        gex_only=True,
        make_unique=True,
        var_names=var_names,
    )
    umi = _as_genes_by_cells(adata.X, adata.var_names, adata.obs_names)
    genes = pd.Index(adata.var_names)

    # CreateSeuratObject(min.cells = min_cells_gene, min.features = 0)
    gene_detected = np.asarray((umi > 0).sum(axis=1)).ravel()
    keep_genes = gene_detected >= fp.min_cells_gene
    umi = umi[keep_genes]
    genes = genes[keep_genes]
    r_genes = _load_r_reference_gene_names(benchmark_root, sample, int(umi.shape[0]))
    if r_genes is not None:
        genes = r_genes

    n_count = np.asarray(umi.sum(axis=0)).ravel()
    n_feature = np.asarray((umi > 0).sum(axis=0)).ravel()
    pct_mt = _percent_mt(genes, umi)

    keep_cells = (
        (n_feature >= fp.min_features)
        & (n_feature <= fp.max_features)
        & (pct_mt < fp.max_percent_mt)
    )
    umi = umi[:, keep_cells]
    cells = pd.Index([f"{sample}_{bc}" for bc in np.asarray(adata.obs_names)[keep_cells]])

    cell_attr = pd.DataFrame(
        {
            "log_umi": np.log10(n_count[keep_cells]),
            "umi": n_count[keep_cells],
            "gene": n_feature[keep_cells],
        },
        index=cells,
    )
    return umi, genes, cells, cell_attr


def load_tcl_filtered_umi(
    sample: str = DEFAULT_SAMPLE,
    *,
    benchmark_root: Path = BENCHMARK_ROOT,
    min_cells: Optional[int] = None,
    var_names: str = DEFAULT_VAR_NAMES,
    seurat_parity: bool = True,
) -> tuple[csr_matrix, pd.Index, pd.Index]:
    """
    Load Tier-B benchmark counts.

    When ``seurat_parity=True`` (default), uses :func:`load_seurat_parity_umi`.
    ``min_cells`` is ignored in that mode — ``vst()`` applies the same filter.
    """
    if seurat_parity:
        umi, genes, cells, _cell_attr = load_seurat_parity_umi(
            sample, benchmark_root=benchmark_root, var_names=var_names
        )
        return umi, genes, cells

    # Legacy scanpy path (01_filter cells + vst-style gene prefilter).
    filter_dir = benchmark_root / "steps/r/01_filter"
    data_dir = benchmark_root / "data" / sample
    min_cells = 5 if min_cells is None else min_cells
    keep = set((filter_dir / f"{sample}_cells.txt").read_text().strip().split("\n"))
    keep.discard("")

    adata = sc.read_10x_mtx(
        data_dir,
        gex_only=True,
        make_unique=True,
        var_names=var_names,
    )
    adata.obs_names = [f"{sample}_{bc}" for bc in adata.obs_names]
    adata = adata[adata.obs_names.isin(keep)].copy()

    umi = _as_genes_by_cells(adata.X, adata.var_names, adata.obs_names)
    genes = pd.Index(adata.var_names)
    cells = pd.Index(adata.obs_names)
    keep_genes = np.asarray((umi >= 0.01).sum(axis=1)).ravel() >= min_cells
    return umi[keep_genes], genes[keep_genes], cells


def load_tcl_filtered_adata(
    sample: str = DEFAULT_SAMPLE,
    *,
    benchmark_root: Path = BENCHMARK_ROOT,
    min_cells: Optional[int] = None,
    var_names: str = DEFAULT_VAR_NAMES,
    seurat_parity: bool = True,
) -> sc.AnnData:
    """AnnData (cells × genes) with ``layers['counts']``."""
    umi, genes, cells = load_tcl_filtered_umi(
        sample,
        benchmark_root=benchmark_root,
        min_cells=min_cells,
        var_names=var_names,
        seurat_parity=seurat_parity,
    )
    x = umi.T.tocsr() if issparse(umi) else csr_matrix(umi.T)
    adata = sc.AnnData(X=x, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
    adata.layers["counts"] = adata.X.copy()
    return adata


def load_seurat_parity_with_cell_attr(
    sample: str = DEFAULT_SAMPLE,
    *,
    benchmark_root: Path = BENCHMARK_ROOT,
) -> tuple[csr_matrix, pd.Index, pd.Index, pd.DataFrame]:
    """Convenience wrapper returning ``cell_attr`` for direct ``vst(..., cell_attr=)``."""
    return load_seurat_parity_umi(sample, benchmark_root=benchmark_root)
