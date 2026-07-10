"""Load GSE288946 benchmark counts with strict QC (100 cells/gene, 50 genes/cell)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix, issparse

from trackcell.tl._sctransform import _as_genes_by_cells

BENCHMARK_ROOT = Path("/Volumes/process/tmp/tcl_test")
BENCHMARK_SAMPLES = ("GSM8779707", "GSM8779708", "GSM8779709")
BENCHMARK_SEEDS = (1448145, 1448146, 42, 123, 999, 2024)
REF_SEED = 1448145
DEFAULT_VAR_NAMES = "gene_symbols"


@dataclass(frozen=True)
class BenchmarkFilterParams:
    """Benchmark QC: min 100 cells per gene, min 50 genes per cell."""

    min_features: int = 50
    min_cells_gene: int = 100
    max_features: int = 1_000_000
    max_percent_mt: float = 100.0


BENCHMARK_FILTER = BenchmarkFilterParams()


def _load_r_reference_gene_names(benchmark_root: Path, sample: str, n_genes: int) -> pd.Index | None:
    path = benchmark_root / "ref" / f"{sample}_genes_min_cells.txt"
    if not path.exists():
        return None
    names = [line for line in path.read_text().splitlines() if line]
    if len(names) == n_genes:
        return pd.Index(names)
    return None


def load_benchmark_umi(
    sample: str,
    *,
    benchmark_root: Path = BENCHMARK_ROOT,
    filter_params: BenchmarkFilterParams | None = None,
    var_names: str = DEFAULT_VAR_NAMES,
) -> tuple[csr_matrix, pd.Index, pd.Index, pd.DataFrame]:
    """
    Load counts mirroring R ``export_sct_benchmark_r.R``.

    Returns genes × cells CSR, Seurat rownames (from ``ref/`` sidecar when present).
    """
    fp = filter_params or BENCHMARK_FILTER
    data_dir = benchmark_root / "data" / sample

    adata = sc.read_10x_mtx(data_dir, gex_only=True, make_unique=True, var_names=var_names)
    umi = _as_genes_by_cells(adata.X, adata.var_names, adata.obs_names)
    genes = pd.Index(adata.var_names)

    gene_detected = np.asarray((umi > 0).sum(axis=1)).ravel()
    keep_genes = gene_detected >= fp.min_cells_gene
    umi = umi[keep_genes]
    genes = genes[keep_genes]
    r_genes = _load_r_reference_gene_names(benchmark_root, sample, int(umi.shape[0]))
    if r_genes is not None:
        genes = r_genes

    n_count = np.asarray(umi.sum(axis=0)).ravel()
    n_feature = np.asarray((umi > 0).sum(axis=0)).ravel()

    keep_cells = n_feature >= fp.min_features
    if fp.max_features < 1_000_000:
        keep_cells &= n_feature <= fp.max_features

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


def r_export_dir(
    sample: str,
    seed: int,
    *,
    benchmark_root: Path = BENCHMARK_ROOT,
) -> Path:
    return benchmark_root / "steps" / "sct_benchmark" / "r" / f"seed_{seed}" / sample


def load_r_export(
    sample: str,
    seed: int,
    *,
    benchmark_root: Path = BENCHMARK_ROOT,
) -> dict:
    """Load R benchmark export for one sample × seed."""
    from trackcell.tl._r_sctransform import load_r_vst_export

    return load_r_vst_export(r_export_dir(sample, seed, benchmark_root=benchmark_root))
