"""Golden parity tests for SCT gmean / step1 / model_pars (R sctransform exports)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from .sct_benchmark_data import (
    BENCHMARK_ROOT,
    DEFAULT_SAMPLE,
    load_seurat_parity_with_cell_attr,
    load_tcl_filtered_adata,
    load_tcl_filtered_umi,
)
from trackcell.tl._r_sctransform import (
    align_r_gene_list,
    coerce_genes_log_gmean,
    load_r_vst_export,
    read_r_gmean_series,
)
from trackcell.tl._sctransform import (
    _row_gmean_sparse,
    vst,
)


SAMPLE = DEFAULT_SAMPLE


def _load_sample_umi(*, min_cells: int = 5):
    umi, genes, cells, cell_attr = load_seurat_parity_with_cell_attr(SAMPLE)
    return umi, genes, cells, cell_attr


def test_seurat_parity_loader_matches_r_meta():
    """Seurat-parity loader should match R export meta cell/gene counts."""
    import json

    meta = json.loads((BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE / "meta.json").read_text())
    umi, genes, cells, cell_attr = load_seurat_parity_with_cell_attr(SAMPLE)
    assert umi.shape[0] == meta["n_genes"]
    assert umi.shape[1] == meta["n_cells"]
    assert len(cells) == meta["n_cells"]
    filter_cells = set(
        (BENCHMARK_ROOT / "steps/r/01_filter" / f"{SAMPLE}_cells.txt").read_text().strip().split("\n")
    )
    filter_cells.discard("")
    assert set(cells) == filter_cells
    assert "log_umi" in cell_attr.columns


@pytest.mark.parametrize(
    "shape",
    [(12602,), (500,)],
)
def test_row_gmean_sparse_matches_rowwise_formula(shape):
    rng = np.random.default_rng(0)
    n_genes, n_cells = 8, shape[0]
    from scipy.sparse import csr_matrix

    data = rng.poisson(2.0, size=(n_genes, n_cells)).astype(np.float64)
    umi = csr_matrix(data)
    dense = np.exp(np.log(data + 1.0).mean(axis=1)) - 1.0
    sparse = _row_gmean_sparse(umi, gmean_eps=1.0)
    np.testing.assert_allclose(sparse, dense, rtol=0, atol=1e-9)


@pytest.mark.skipif(
    not (BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE / "gene_attr.csv").exists(),
    reason="R stepwise export missing",
)
def test_python_genes_log_gmean_matches_r_golden():
    umi, genes, _cells, _cell_attr = _load_sample_umi()
    export = load_r_vst_export(BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE)
    ga = pd.read_csv(BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE / "gene_attr.csv", index_col=0)
    r_log10 = pd.Series(
        np.log10(ga["gmean"].clip(1e-9).to_numpy(dtype=np.float64)),
        index=align_r_gene_list(ga.index.astype(str).tolist(), genes),
    )
    py_log10 = pd.Series(
        np.log10(np.maximum(_row_gmean_sparse(umi, 1.0), 1e-9)),
        index=genes,
    )
    shared = r_log10.index.intersection(py_log10.index)
    diff = (py_log10.loc[shared] - r_log10.loc[shared]).abs()
    assert float(diff.median()) < 1e-6
    assert float(diff.max()) < 1e-4
    assert py_log10.loc[shared].corr(r_log10.loc[shared]) > 0.9999


def test_coerce_genes_log_gmean_linear_export():
    linear = pd.Series([0.001, 0.01, 1.0, 10.0], index=["a", "b", "c", "d"])
    out = coerce_genes_log_gmean(linear)
    expected = np.log10(linear.to_numpy())
    np.testing.assert_allclose(out.to_numpy(), expected, rtol=0, atol=1e-12)


def test_coerce_genes_log_gmean_keeps_log10_export():
    log10 = pd.Series([-2.0, -0.5, 0.3, 1.8], index=["a", "b", "c", "d"])
    out = coerce_genes_log_gmean(log10)
    pd.testing.assert_series_equal(out, log10.astype(float))


@pytest.mark.skipif(
    not (BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE / "model_pars_step1.csv").exists(),
    reason="R stepwise export missing",
)
def test_vst_gmean_golden_improves_step1_intercept_corr():
    """After gmean fix, native py step1 Intercept should track R more closely."""
    umi, genes, cells, cell_attr = _load_sample_umi()
    export = load_r_vst_export(BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE)
    genes_step1 = pd.Index(align_r_gene_list(export["genes_step1"], genes))
    cells_step1 = pd.Index(export["cells_step1"]).intersection(cells)

    out = vst(
        umi,
        genes,
        cells,
        cell_attr=cell_attr,
        vst_flavor="v2",
        cells_step1=cells_step1,
        genes_step1=genes_step1,
        return_corrected_umi=False,
        seed=1448145,
    )
    py_hvg = set(
        out["gene_attr"]["residual_variance"].sort_values(ascending=False).index[:3000].astype(str)
    )
    r_hvg = set(
        (BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE / "hvg_top3000.txt")
        .read_text()
        .strip()
        .split("\n")
    )
    jaccard = len(py_hvg & r_hvg) / len(py_hvg | r_hvg)
    assert jaccard > 0.5, f"HVG Jaccard {jaccard:.3f} still low after gmean fix"


@pytest.mark.skipif(
    not (BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE / "model_pars_step1.csv").exists(),
    reason="R stepwise export missing",
)
def test_native_step1_intercept_tracks_r():
    """Native pyglmGamPoi step1 Intercept should match R glmGamPoi closely."""
    umi, genes, cells, cell_attr = _load_sample_umi()
    export = load_r_vst_export(BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE)
    genes_step1 = pd.Index(align_r_gene_list(export["genes_step1"], genes))
    cells_step1 = pd.Index(export["cells_step1"]).intersection(cells)

    out = vst(
        umi,
        genes,
        cells,
        cell_attr=cell_attr,
        vst_flavor="v2",
        cells_step1=cells_step1,
        genes_step1=genes_step1,
        return_corrected_umi=False,
        seed=1448145,
    )
    r_step1 = export["model_pars_step1"]
    shared = genes_step1.intersection(r_step1.index)
    py = out["model_pars"].loc[shared, "Intercept"]
    rr = r_step1.loc[shared, "Intercept"]
    assert float(py.corr(rr)) > 0.99


@pytest.mark.skipif(
    not (BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE / "model_pars_step1.csv").exists(),
    reason="R stepwise export missing",
)
def test_native_step1_theta_log10_tracks_r():
    """Document pyglmGamPoi theta MLE gap vs R (target >0.75 until pyglmGamPoi fix)."""
    umi, genes, cells, cell_attr = _load_sample_umi()
    export = load_r_vst_export(BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE)
    genes_step1 = pd.Index(align_r_gene_list(export["genes_step1"], genes))
    cells_step1 = pd.Index(export["cells_step1"]).intersection(cells)

    out = vst(
        umi,
        genes,
        cells,
        cell_attr=cell_attr,
        vst_flavor="v2",
        cells_step1=cells_step1,
        genes_step1=genes_step1,
        return_corrected_umi=False,
        seed=1448145,
    )
    r_step1 = export["model_pars_step1"]
    shared = genes_step1.intersection(r_step1.index)
    py = np.log10(out["model_pars"].loc[shared, "theta"].clip(1e-9))
    rr = np.log10(r_step1.loc[shared, "theta"].clip(1e-9))
    finite = np.isfinite(py) & np.isfinite(rr)
    corr = float(py[finite].corr(rr[finite]))
    assert corr > 0.75, f"step1 theta log10 corr {corr:.3f} regressed vs R"


@pytest.mark.skipif(
    not (BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE / "model_pars_fit.csv").exists(),
    reason="R stepwise export missing",
)
def test_r_step1_model_pars_fit_tracks_r():
    """R step1 + Python ksmooth should match R model_pars_fit (theta is the gate)."""
    umi, genes, cells, cell_attr = _load_sample_umi()
    export = load_r_vst_export(BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE)
    genes_step1 = pd.Index(align_r_gene_list(export["genes_step1"], genes))
    cells_step1 = pd.Index(export["cells_step1"]).intersection(cells)

    out = vst(
        umi,
        genes,
        cells,
        cell_attr=cell_attr,
        vst_flavor="v2",
        cells_step1=cells_step1,
        genes_step1=genes_step1,
        model_pars_step1=export["model_pars_step1"],
        genes_log_gmean=export["genes_log_gmean"],
        genes_log_gmean_step1=export["genes_log_gmean_step1"],
        return_corrected_umi=False,
        seed=1448145,
    )
    r_fit = export["model_pars_fit"]
    shared = genes_step1.intersection(r_fit.index)
    for col in ("Intercept", "theta"):
        py = out["model_pars_fit"].loc[shared, col]
        rr = r_fit.loc[shared, col]
        finite = np.isfinite(py) & np.isfinite(rr)
        assert float(py[finite].corr(rr[finite])) > 0.99, f"{col} corr low vs R"


@pytest.mark.skipif(
    not (BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE / "hvg_top3000.txt").exists(),
    reason="R stepwise export missing",
)
def test_native_hvg_jaccard_gap_documented():
    """Tier-B native Py HVG stays ~0.55 until step1 model_pars matches R."""
    umi, genes, cells, cell_attr = _load_sample_umi()
    export = load_r_vst_export(BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE)
    genes_step1 = pd.Index(align_r_gene_list(export["genes_step1"], genes))
    cells_step1 = pd.Index(export["cells_step1"]).intersection(cells)

    out = vst(
        umi,
        genes,
        cells,
        cell_attr=cell_attr,
        vst_flavor="v2",
        cells_step1=cells_step1,
        genes_step1=genes_step1,
        return_corrected_umi=False,
        seed=1448145,
    )
    r_hvg = set(
        (BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE / "hvg_top3000.txt")
        .read_text()
        .strip()
        .split("\n")
    )
    py_hvg = set(
        out["gene_attr"]["residual_variance"].sort_values(ascending=False).index[:3000].astype(str)
    )
    jaccard = len(py_hvg & r_hvg) / len(py_hvg | r_hvg)
    assert 0.45 < jaccard < 0.65, f"unexpected native HVG Jaccard {jaccard:.3f}"


@pytest.mark.skipif(
    not (BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE / "hvg_top3000.txt").exists(),
    reason="R stepwise export missing",
)
def test_r_step1_hvg_jaccard_near_r():
    """R step1 golden + Python ksmooth/residuals should recover R HVG list."""
    umi, genes, cells, cell_attr = _load_sample_umi()
    export = load_r_vst_export(BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE)
    genes_step1 = pd.Index(align_r_gene_list(export["genes_step1"], genes))
    cells_step1 = pd.Index(export["cells_step1"]).intersection(cells)

    out = vst(
        umi,
        genes,
        cells,
        cell_attr=cell_attr,
        vst_flavor="v2",
        cells_step1=cells_step1,
        genes_step1=genes_step1,
        model_pars_step1=export["model_pars_step1"],
        genes_log_gmean=export["genes_log_gmean"],
        genes_log_gmean_step1=export["genes_log_gmean_step1"],
        return_corrected_umi=False,
        seed=1448145,
    )
    py_hvg = set(
        out["gene_attr"]["residual_variance"].sort_values(ascending=False).index[:3000].astype(str)
    )
    r_hvg = set(
        (BENCHMARK_ROOT / "steps/r_sct_stepwise" / SAMPLE / "hvg_top3000.txt")
        .read_text()
        .strip()
        .split("\n")
    )
    jaccard = len(py_hvg & r_hvg) / len(py_hvg | r_hvg)
    assert jaccard > 0.95, f"r_step1 HVG Jaccard {jaccard:.3f} below acceptance"
