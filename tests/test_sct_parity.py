"""Tests for SCT integration parity helpers."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from trackcell.benchmark.rpca.data import RPCA_BENCHMARK_ROOT, RPCA_SAMPLES
from trackcell.tl.integration.sct_parity import (
    load_rpca_subsample_indices,
    materialize_rpca_r_step1_export_dir,
    patch_r_sct_scale_data,
    resolve_r_vst_exports,
    rpca_exports_support_r_step1,
    rpca_scale_data_available,
    rpca_step1_available,
    seurat_sct_clip_range,
)

RPCA_SEED = 1448145


def test_seurat_sct_clip_range():
    lo, hi = seurat_sct_clip_range(3000)
    bound = math.sqrt(3000 / 30.0)
    assert lo == pytest.approx(-bound)
    assert hi == pytest.approx(bound)


@pytest.mark.skipif(
    not RPCA_BENCHMARK_ROOT.is_dir(),
    reason="RPCA benchmark root not available",
)
def test_rpca_step1_exports_available():
    for batch in RPCA_SAMPLES:
        assert rpca_step1_available(RPCA_BENCHMARK_ROOT, RPCA_SEED, batch)


@pytest.mark.skipif(
    not rpca_exports_support_r_step1(RPCA_BENCHMARK_ROOT, RPCA_SEED, RPCA_SAMPLES),
    reason="RPCA step-1 exports missing",
)
def test_load_rpca_subsample_indices():
    subs = load_rpca_subsample_indices(RPCA_BENCHMARK_ROOT, RPCA_SEED, RPCA_SAMPLES)
    assert set(subs) == set(RPCA_SAMPLES)
    for batch in RPCA_SAMPLES:
        assert len(subs[batch]["cells_step1"]) > 0
        assert len(subs[batch]["genes_step1"]) > 0


@pytest.mark.skipif(
    not rpca_exports_support_r_step1(RPCA_BENCHMARK_ROOT, RPCA_SEED, (RPCA_SAMPLES[0],)),
    reason="RPCA step-1 exports missing",
)
def test_materialize_and_resolve_r_vst_exports(tmp_path: Path):
    cache = materialize_rpca_r_step1_export_dir(
        RPCA_BENCHMARK_ROOT,
        RPCA_SEED,
        RPCA_SAMPLES[0],
        force_refresh=True,
    )
    assert (cache / "cells_step1.txt").exists()
    assert (cache / "genes_step1.txt").exists()
    assert (cache / "model_pars_step1.csv").exists()
    assert (cache / "genes_log_gmean.csv").exists()

    exports, layout = resolve_r_vst_exports(
        (RPCA_SAMPLES[0],),
        export_root=RPCA_BENCHMARK_ROOT,
        export_seed=RPCA_SEED,
        layout="rpca_benchmark",
    )
    assert layout == "rpca_benchmark"
    assert RPCA_SAMPLES[0] in exports
    assert Path(exports[RPCA_SAMPLES[0]]).is_dir()


@pytest.mark.skipif(
    not rpca_scale_data_available(RPCA_BENCHMARK_ROOT, RPCA_SEED, RPCA_SAMPLES),
    reason="RPCA residuals_hvg exports missing",
)
def test_patch_r_sct_scale_data():
    from trackcell.benchmark.rpca.data import BATCH_KEY, load_merged_adata_native
    from trackcell.tl.sctransform import sctransform

    adata = load_merged_adata_native(benchmark_root=RPCA_BENCHMARK_ROOT)
    sctransform(
        adata,
        layer="counts",
        batch_key=BATCH_KEY,
        n_top_genes=50,
        n_cells=200,
        n_genes=200,
        seed=RPCA_SEED,
        key_added="sct",
    )
    assert patch_r_sct_scale_data(
        adata,
        sct_key="sct",
        batch_key=BATCH_KEY,
        benchmark_root=RPCA_BENCHMARK_ROOT,
        seed=RPCA_SEED,
        batches=RPCA_SAMPLES,
    )
    entry = adata.uns["sct"]["batch_models"][RPCA_SAMPLES[0]]
    assert "scale_data" in entry
