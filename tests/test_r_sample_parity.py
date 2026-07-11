"""Tests for R-compatible SCT step-1 subsampling."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from trackcell.benchmark.sct.data import BENCHMARK_ROOT, BENCHMARK_SAMPLES, REF_SEED, load_benchmark_umi
from trackcell.benchmark.sct.metrics import hvg_jaccard
from trackcell.benchmark.sct.data import load_r_export, r_export_dir
from trackcell.tl._r_density import r_density_fft
from trackcell.tl._r_sample import r_sample_available
from trackcell.tl._sctransform import _dds, _row_gmean_sparse, vst

pytestmark = pytest.mark.skipif(not r_sample_available(), reason="rrng not installed")


def _genes_step1_pool(umi, genes, cells, cells_step1: pd.Index):
    """Pre-subsample gene pool matching ``vst()`` step-1 filters."""
    cells_idx = [list(cells).index(c) for c in cells_step1]
    sub = umi[:, cells_idx]
    genes_cell_count = np.asarray((sub > 0).sum(axis=1)).ravel()
    genes_step1 = genes[genes_cell_count >= 5]
    amean = np.asarray(umi.mean(axis=1)).ravel()
    var = np.asarray(umi.power(2).mean(axis=1)).ravel() - amean**2
    over = (var - amean)[genes.get_indexer(genes_step1)] > 0
    return genes_step1[over]


@pytest.mark.parametrize("sample", BENCHMARK_SAMPLES)
def test_cells_step1_matches_r_export(sample: str) -> None:
    """R ``sample(colnames, n_cells)`` parity via rrng."""
    from trackcell.tl._r_sample import _RSampleRNG

    umi, genes, cells, cell_attr = load_benchmark_umi(sample, benchmark_root=BENCHMARK_ROOT)
    r = load_r_export(sample, REF_SEED, benchmark_root=BENCHMARK_ROOT)
    r_cells = pd.Index(r["cells_step1"]).intersection(cells)

    rng = _RSampleRNG(REF_SEED)
    picked = pd.Index(rng.choice(cells.to_numpy(), size=2000, replace=False))
    assert list(picked) == list(r_cells)


def test_r_density_fft_matches_r_golden() -> None:
    """FFT ``density.default`` grid matches R (legacy kernel coords)."""
    gvals_path = Path("/tmp/r_gvals.csv")
    golden_path = Path("/tmp/r_density_golden.csv")
    if not gvals_path.exists() or not golden_path.exists():
        pytest.skip("R golden vectors not exported (/tmp/r_gvals.csv, /tmp/r_density_golden.csv)")
    gvals = pd.read_csv(gvals_path)["g"].to_numpy()
    golden = pd.read_csv(golden_path)
    _, py_y, py_bw = r_density_fft(gvals, bw="nrd", adjust=1, n=512)
    assert abs(py_bw - 0.1008455) < 1e-6
    assert np.max(np.abs(py_y - golden["y"].to_numpy())) < 1e-12


@pytest.mark.parametrize("sample", BENCHMARK_SAMPLES)
def test_genes_step1_jaccard_vs_r(sample: str) -> None:
    """Weighted ``genes_step1`` via R FFT density + rrng ``sample()``."""
    from trackcell.tl._r_sample import _RSampleRNG

    umi, genes, cells, _ = load_benchmark_umi(sample, benchmark_root=BENCHMARK_ROOT)
    r = load_r_export(sample, REF_SEED, benchmark_root=BENCHMARK_ROOT)
    cells_step1 = pd.Index(r["cells_step1"])
    pool = _genes_step1_pool(umi, genes, cells, cells_step1)
    gmean = pd.Series(
        np.log10(np.maximum(_row_gmean_sparse(umi, 1), 1e-9)),
        index=genes,
    )
    rng = _RSampleRNG(REF_SEED)
    _ = rng.choice(cells.to_numpy(), size=2000, replace=False)
    picked = set(
        rng.choice(pool.to_numpy(), size=2000, replace=False, p=_dds(gmean.loc[pool]))
    )
    rset = set(r["genes_step1"])
    j = len(picked & rset) / len(picked | rset)
    assert j >= 0.99, f"{sample}: genes_step1 Jaccard {j:.4f} < 0.99"


@pytest.mark.parametrize("sample", BENCHMARK_SAMPLES)
def test_native_vst_hvg_jaccard_vs_r(sample: str) -> None:
    """Full native ``vst()`` HVG improves with R subsampling + FFT density."""
    umi, genes, cells, cell_attr = load_benchmark_umi(sample, benchmark_root=BENCHMARK_ROOT)
    r_hvg = set(
        (
            r_export_dir(sample, REF_SEED, benchmark_root=BENCHMARK_ROOT) / "hvg_top3000.txt"
        )
        .read_text()
        .strip()
        .split("\n")
    ) - {""}

    out = vst(
        umi,
        genes,
        cells,
        cell_attr=cell_attr,
        vst_flavor="v2",
        seed=REF_SEED,
        return_corrected_umi=False,
    )
    j = hvg_jaccard(out["gene_attr"], r_hvg)
    # Step-1 parity is tight; remaining chain drift is mostly pyglmGamPoi / HVG ranking.
    assert j >= 0.925, f"{sample}: native HVG Jaccard {j:.4f} < 0.925"
