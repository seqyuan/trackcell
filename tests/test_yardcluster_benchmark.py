"""Smoke tests for YardCluster benchmark suite."""

import numpy as np
import pytest

from trackcell.benchmark.yardcluster.data import make_patch_adata
from trackcell.benchmark.yardcluster.metrics import evaluate_clustering, spatial_coherence
from trackcell.benchmark.yardcluster.runner import run_yardcluster_benchmark


def test_make_patch_adata_shape():
    adata = make_patch_adata(n_per_patch=10, n_patches=2, n_genes=20, seed=1)
    assert adata.n_obs == 20
    assert adata.n_vars == 20
    assert "spatial" in adata.obsm
    assert adata.obs["patch"].nunique() == 2


def test_run_yardcluster_benchmark_smoke():
    pytest.importorskip("leidenalg")
    summary = run_yardcluster_benchmark(seed=0, n_per_patch=25, n_patches=3)
    assert len(summary["configs"]) == 3
    for row in summary["configs"]:
        assert "ari" in row
        assert row["ari"] >= 0.0
        assert row["spatial_coherence"] >= 0.0


def test_evaluate_clustering_perfect_labels():
    from trackcell.tl.spatial_cluster import spatial_neighbors

    adata = make_patch_adata(n_per_patch=8, n_patches=2, seed=2)
    adata.obs["pred"] = adata.obs["patch"].astype(str)
    spatial_neighbors(adata, k=5)
    out = evaluate_clustering(adata, "pred", "patch")
    assert out["ari"] == pytest.approx(1.0)
    assert out["spatial_coherence"] == pytest.approx(1.0, abs=0.05)
