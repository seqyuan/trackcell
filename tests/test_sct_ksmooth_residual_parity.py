"""Parity tests for SCT ksmooth / residual helpers vs R sctransform."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sps

from trackcell.tl._r_sctransform import (
    align_r_gene_list,
    align_r_gmean_series,
    align_r_model_pars,
    load_r_vst_export,
    normalize_r_model_pars,
)
from trackcell.tl._sctransform import (
    _as_genes_by_cells,
    _is_outlier,
    _ksmooth,
    _make_cell_attr,
    _pearson_residual,
    _reg_model_pars,
    _reg_model_pars_pre_ksmooth,
    _robust_scale_binned,
    bw_sj,
    get_residuals,
)


def test_is_outlier_uses_minimum_of_binned_scores():
    rng = np.random.default_rng(0)
    x = pd.Series(rng.normal(size=200))
    y = pd.Series(rng.normal(size=200))
    x_max, x_min = float(x.max()), float(x.min())
    eps = np.finfo(float).eps * 10
    bin_width = (x_max - x_min) * bw_sj(x.to_numpy()) / 2
    breaks1 = np.arange(x_min - eps, x_max + bin_width, bin_width)
    breaks2 = np.arange(x_min - eps - bin_width / 2, x_max + bin_width, bin_width)
    score1 = np.abs(_robust_scale_binned(y, x, breaks1))
    score2 = np.abs(_robust_scale_binned(y, x, breaks2))
    expected = np.minimum(score1, score2) > 10.0
    assert np.array_equal(_is_outlier(y, x, th=10.0), expected)


def test_robust_scale_binned_matches_length():
    x = pd.Series(np.linspace(-2, 2, 50))
    y = pd.Series(np.sin(x) + np.random.default_rng(1).normal(0, 0.05, 50))
    eps = np.finfo(float).eps * 10
    bin_width = (x.max() - x.min()) * bw_sj(x.to_numpy()) / 2
    breaks = np.arange(float(x.min()) - eps, float(x.max()) + bin_width, bin_width)
    scores = _robust_scale_binned(y, x, breaks)
    assert scores.shape == (len(y),)


@pytest.mark.skipif(
    not __import__("trackcell.tl._r_sctransform", fromlist=["glmGamPoi_r_available"]).glmGamPoi_r_available(),
    reason="R not available for ksmooth reference",
)
def test_ksmooth_matches_r_reference():
    rng = np.random.default_rng(42)
    x = np.sort(rng.uniform(-2, 2, 40))
    y = np.sin(x) + rng.normal(0, 0.05, 40)
    xp = np.sort(rng.uniform(-2, 2, 60))
    bw = bw_sj(x) * 3
    py = _ksmooth(pd.Series(x), pd.Series(y), xp.copy(), kern=2, bw=bw)

    import subprocess
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        np.savetxt(tmp_path / "x.csv", x, delimiter=",")
        np.savetxt(tmp_path / "y.csv", y, delimiter=",")
        np.savetxt(tmp_path / "xp.csv", xp, delimiter=",")
        r_code = f"""
        x <- scan("{tmp_path}/x.csv", quiet=TRUE, sep=",")
        y <- scan("{tmp_path}/y.csv", quiet=TRUE, sep=",")
        xp <- scan("{tmp_path}/xp.csv", quiet=TRUE, sep=",")
        bw <- {bw}
        write.table(ksmooth(x, y, x.points=xp, bandwidth=bw, kernel="normal")$y,
                    "{tmp_path}/r_y.csv", row.names=FALSE, col.names=FALSE)
        """
        subprocess.run(
            ["conda", "run", "-n", "st", "Rscript", "-e", r_code],
            check=True,
            capture_output=True,
            text=True,
        )
        r_y = np.loadtxt(tmp_path / "r_y.csv")
    np.testing.assert_allclose(py, r_y, rtol=0, atol=1e-12)


def test_pearson_residual_theta_inf_is_poisson():
    y = np.array([[0.0, 2.0, 4.0]])
    mu = np.array([[1.0, 2.0, 3.0]])
    theta = np.array([np.inf])
    res = _pearson_residual(y, mu, theta, min_var=0.0)
    expected = (y - mu) / np.sqrt(mu)
    np.testing.assert_allclose(res, expected, rtol=1e-12, atol=1e-12)


def test_reg_model_pars_marks_poisson_theta_inf():
    genes = pd.Index([f"g{i}" for i in range(6)])
    cells = pd.Index([f"c{i}" for i in range(20)])
    rng = np.random.default_rng(3)
    counts = rng.poisson(1.5, size=(len(genes), len(cells))).astype(np.float64)
    umi = sps.csr_matrix(counts)
    genes_step1 = genes[:4]
    model_pars = pd.DataFrame(
        {
            "theta": [5.0, 8.0, np.inf, 3.0],
            "Intercept": [-2.0, -1.5, -3.0, -2.5],
            "log_umi": [np.log(10.0)] * 4,
        },
        index=genes_step1,
    )
    genes_log_gmean = pd.Series(np.log10(np.maximum(counts.mean(1), 1e-9)), index=genes)
    genes_log_gmean_step1 = genes_log_gmean.loc[genes_step1]
    fit, _ = _reg_model_pars(
        model_pars,
        genes_log_gmean_step1,
        genes_log_gmean,
        bw_adjust=3.0,
        theta_regularization="od_factor",
        umi=umi,
        exclude_poisson=True,
        genes=genes,
        genes_step1=genes_step1,
    )
    poisson_mask = fit["theta"].to_numpy() == np.inf
    assert poisson_mask.any()
    assert np.isfinite(fit.loc[~poisson_mask, "theta"]).all()


@pytest.mark.skipif(
    not Path("/Volumes/process/tmp/tcl_test/steps/r_sct_stepwise/GSM8779707/model_pars_fit.csv").exists(),
    reason="R stepwise export missing",
)
def test_residuals_match_r_model_pars_fit():
    """Pearson residuals from R model_pars_fit should match R vst residuals."""
    import scanpy as sc

    from trackcell.tl._sctransform import _as_genes_by_cells, _make_cell_attr

    sample = "GSM8779707"
    root = Path("/Volumes/process/tmp/tcl_test")
    r_dir = root / "steps/r_sct_stepwise" / sample
    meta = json.loads((r_dir / "meta.json").read_text()) if (r_dir / "meta.json").exists() else {}
    keep = set((root / "steps/r/01_filter" / f"{sample}_cells.txt").read_text().strip().split("\n"))
    keep.discard("")
    adata = sc.read_10x_mtx(root / "data" / sample, gex_only=True, make_unique=True)
    adata.obs_names = [f"{sample}_{bc}" for bc in adata.obs_names]
    adata = adata[adata.obs_names.isin(keep)].copy()
    umi = _as_genes_by_cells(adata.X, adata.var_names, adata.obs_names)
    genes = pd.Index(adata.var_names)
    keep_genes = np.asarray((umi >= 0.01).sum(axis=1)).ravel() >= 5
    umi = umi[keep_genes]
    genes = genes[keep_genes]
    cells = pd.Index(adata.obs_names)
    cell_attr = _make_cell_attr(umi, cells)

    r_fit = pd.read_csv(r_dir / "model_pars_fit.csv", index_col=0).rename(
        columns={"(Intercept)": "Intercept"}
    )
    r_res = pd.read_csv(r_dir / "residuals.csv", index_col=0)
    min_var = meta.get("min_variance")
    if min_var is None:
        from trackcell.tl._sctransform_v2 import resolve_min_variance

        min_var = resolve_min_variance("umi_median", umi)

    vst_stub = {
        "model_str": "y ~ log_umi",
        "model_pars_fit": r_fit,
        "model_pars_nonreg": None,
        "model_str_nonreg": "",
        "cell_attr": cell_attr,
    }
    py_res = get_residuals(
        vst_stub,
        umi,
        genes=genes,
        umi_genes=genes,
        min_variance=min_var,
    )
    shared_genes = py_res.index.intersection(r_res.index)
    shared_cells = py_res.columns.intersection(r_res.columns)
    py_mat = py_res.loc[shared_genes, shared_cells].to_numpy()
    r_mat = r_res.loc[shared_genes, shared_cells].to_numpy()
    corr = np.corrcoef(py_mat.ravel(), r_mat.ravel())[0, 1]
    assert corr > 0.9999
    assert np.nanmedian(np.abs(py_mat - r_mat)) < 1e-10


@pytest.mark.skipif(
    not Path("/Volumes/process/tmp/tcl_test/steps/r_sct_stepwise/GSM8779707/model_pars_step1.csv").exists(),
    reason="R stepwise export missing",
)
def test_load_r_vst_export_and_r_fit_hvg():
    """R model_pars_fit export -> Python residuals should match R HVG (~1.0 Jaccard)."""
    import scanpy as sc

    from trackcell.tl.sctransform import sctransform

    sample = "GSM8779707"
    root = Path("/Volumes/process/tmp/tcl_test")
    export_dir = root / "steps/r_sct_stepwise" / sample
    exported = load_r_vst_export(export_dir)
    assert exported["model_pars_step1"] is not None
    assert exported["model_pars_fit"] is not None
    assert len(exported["cells_step1"]) > 0
    assert len(exported["genes_step1"]) > 0
    assert exported.get("genes_log_gmean") is not None
    assert exported.get("genes_log_gmean_step1") is not None

    keep = set((root / "steps/r/01_filter" / f"{sample}_cells.txt").read_text().strip().split("\n"))
    keep.discard("")
    adata = sc.read_10x_mtx(root / "data" / sample, gex_only=True, make_unique=True)
    adata.obs_names = [f"{sample}_{bc}" for bc in adata.obs_names]
    adata = adata[adata.obs_names.isin(keep)].copy()
    adata.layers["counts"] = adata.X.copy()

    sctransform(
        adata,
        layer="counts",
        n_top_genes=3000,
        method="r_fit",
        r_vst_export_dir=export_dir,
        vst_flavor="v2",
        seed=1448145,
        do_correct_umi=False,
    )
    py_hvg = set(adata.uns["sct"]["variable_features"])
    r_hvg = set((export_dir / "hvg_top3000.txt").read_text().strip().split("\n"))
    jaccard = len(py_hvg & r_hvg) / len(py_hvg | r_hvg)
    assert jaccard > 0.99


@pytest.mark.skipif(
    not Path("/Volumes/process/tmp/tcl_test/steps/r_sct_reg_debug/GSM8779707/summary.json").exists(),
    reason="R reg_model_pars debug export missing",
)
def test_reg_model_pars_pre_ksmooth_matches_r_export():
    """Pre-ksmooth outlier/poisson filtering matches R when using R gmean exports."""
    import scanpy as sc

    sample = "GSM8779707"
    root = Path("/Volumes/process/tmp/tcl_test")
    step_dir = root / "steps/r_sct_stepwise" / sample
    dbg_dir = root / "steps/r_sct_reg_debug" / sample
    exported = load_r_vst_export(step_dir)

    keep = set((root / "steps/r/01_filter" / f"{sample}_cells.txt").read_text().strip().split("\n"))
    keep.discard("")
    adata = sc.read_10x_mtx(root / "data" / sample, gex_only=True, make_unique=True)
    adata.obs_names = [f"{sample}_{bc}" for bc in adata.obs_names]
    adata = adata[adata.obs_names.isin(keep)].copy()
    umi = _as_genes_by_cells(adata.X, adata.var_names, adata.obs_names)
    genes = pd.Index(adata.var_names)
    keep_genes = np.asarray((umi >= 0.01).sum(axis=1)).ravel() >= 5
    umi = umi[keep_genes]
    genes = genes[keep_genes]

    genes_vst = pd.Index(align_r_gene_list(exported["genes_vst"], genes))
    gene_pos = genes.get_indexer(genes_vst)
    keep = gene_pos >= 0
    umi = umi[gene_pos[keep]]
    genes = genes[gene_pos[keep]]

    py = _reg_model_pars_pre_ksmooth(
        align_r_model_pars(exported["model_pars_step1"], genes),
        align_r_gmean_series(exported["genes_log_gmean_step1"], genes),
        align_r_gmean_series(exported["genes_log_gmean"], genes),
        umi=umi,
        genes=genes,
        genes_step1=pd.Index(align_r_gene_list(exported["genes_step1"], genes)),
        exclude_poisson=True,
    )

    r_out = pd.read_csv(dbg_dir / "outliers_step1.csv", index_col=0)["outlier"]
    r_out.index = align_r_gene_list(r_out.index.astype(str).tolist(), genes)
    shared = py["outliers"].index.intersection(r_out.index)
    assert bool((py["outliers"].loc[shared].to_numpy() == r_out.loc[shared].to_numpy()).all())

    r_over = set(align_r_gene_list((dbg_dir / "overdispersed_genes.txt").read_text().strip().split("\n"), genes))
    py_over = set(map(str, py["overdispersed_genes"]))
    assert py_over == r_over


@pytest.mark.skipif(
    not Path("/Volumes/process/tmp/tcl_test/steps/r_sct_stepwise/GSM8779707/genes_log_gmean_step1.csv").exists(),
    reason="R gmean export missing",
)
def test_r_step1_hvg_with_r_gmean_export():
    """R step1 + Python ksmooth/residuals with R gmean exports should beat naive Py gmean."""
    import scanpy as sc

    from trackcell.tl._sctransform import vst

    sample = "GSM8779707"
    root = Path("/Volumes/process/tmp/tcl_test")
    export_dir = root / "steps/r_sct_stepwise" / sample
    exported = load_r_vst_export(export_dir)

    keep = set((root / "steps/r/01_filter" / f"{sample}_cells.txt").read_text().strip().split("\n"))
    keep.discard("")
    adata = sc.read_10x_mtx(root / "data" / sample, gex_only=True, make_unique=True)
    adata.obs_names = [f"{sample}_{bc}" for bc in adata.obs_names]
    adata = adata[adata.obs_names.isin(keep)].copy()
    umi = _as_genes_by_cells(adata.X, adata.var_names, adata.obs_names)
    genes = pd.Index(adata.var_names)
    cells = pd.Index(adata.obs_names)
    keep_genes = np.asarray((umi >= 0.01).sum(axis=1)).ravel() >= 5
    umi = umi[keep_genes]
    genes = genes[keep_genes]
    cell_attr = _make_cell_attr(umi, cells)

    genes_vst = pd.Index(align_r_gene_list(exported["genes_vst"], genes))
    gene_pos = genes.get_indexer(genes_vst)
    keep = gene_pos >= 0
    umi = umi[gene_pos[keep]]
    genes = genes[gene_pos[keep]]

    out = vst(
        umi,
        genes,
        cells,
        cell_attr=cell_attr,
        vst_flavor="v2",
        cells_step1=pd.Index(exported["cells_step1"]).intersection(cells),
        genes_step1=pd.Index(align_r_gene_list(exported["genes_step1"], genes)),
        model_pars_step1=align_r_model_pars(exported["model_pars_step1"], genes),
        genes_log_gmean=align_r_gmean_series(exported["genes_log_gmean"], genes),
        genes_log_gmean_step1=align_r_gmean_series(exported["genes_log_gmean_step1"], genes),
        min_variance=exported["meta"].get("min_variance", -np.inf),
        return_corrected_umi=False,
        seed=1448145,
    )
    py_hvg = set(
        out["gene_attr"]["residual_variance"].sort_values(ascending=False).index[:3000].astype(str)
    )
    r_hvg = set((export_dir / "hvg_top3000.txt").read_text().strip().split("\n"))
    jaccard = len(py_hvg & r_hvg) / len(py_hvg | r_hvg)
    assert jaccard > 0.75
