"""Tests for DBSCAN spatial slice separation."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from trackcell.tl.slice_separation import (
    dbscan_slice_labels,
    distance_to_nearest_centroids,
    mark_colony_centroids,
    slice_cluster_summary,
    spatial_colony_cluster,
    spatial_slice_cluster,
    split_by_slice,
    write_slice_annotation,
)


def _make_tma_adata(n_cores: int = 4, cells_per_core: int = 200, gap: float = 500.0):
    """Simulate TMA cores on a grid with large gaps."""
    rng = np.random.default_rng(0)
    coords_list = []
    for row in range(2):
        for col in range(n_cores // 2):
            cx = col * gap + 1500
            cy = row * gap + 1500
            pts = rng.normal(scale=40, size=(cells_per_core, 2)) + [cx, cy]
            coords_list.append(pts)
    coords = np.vstack(coords_list)
    adata = AnnData(X=rng.normal(size=(len(coords), 10)).astype(np.float32))
    adata.obsm["spatial"] = coords.astype(np.float32)
    adata.obs["x_centroid"] = coords[:, 0]
    adata.obs["y_centroid"] = coords[:, 1]
    adata.var_names = [f"g{i}" for i in range(10)]
    adata.obs_names = [f"c{i}" for i in range(len(coords))]
    return adata


def test_dbscan_slice_labels_finds_cores():
    pytest.importorskip("sklearn")
    adata = _make_tma_adata(n_cores=4, cells_per_core=150)
    clusters, slice_ids = dbscan_slice_labels(
        adata.obsm["spatial"], eps=120, min_samples=5, min_cells=50
    )
    valid = {s for s in slice_ids if s != "debris"}
    assert len(valid) == 4
    assert (clusters >= -1).all()
    assert set(slice_ids) - {"debris"} == {"S001", "S002", "S003", "S004"}


def test_spatial_slice_cluster_inplace():
    pytest.importorskip("sklearn")
    adata = _make_tma_adata(n_cores=2, cells_per_core=120)
    spatial_slice_cluster(adata, eps=120, min_samples=5, min_cells=50)
    assert "slice_id" in adata.obs
    assert "slice_cluster" in adata.obs
    assert adata.uns["slice_id_summary"]["n_slices"] == 2


def test_slice_cluster_summary():
    pytest.importorskip("sklearn")
    adata = _make_tma_adata(n_cores=2, cells_per_core=120)
    spatial_slice_cluster(adata, eps=120, min_samples=5, min_cells=50)
    summary = slice_cluster_summary(adata)
    assert len(summary) == 2
    assert {"slice_id", "n_cells", "x_span", "y_span"}.issubset(summary.columns)
    assert (summary["n_cells"] >= 50).all()


def test_split_by_slice():
    pytest.importorskip("sklearn")
    adata = _make_tma_adata(n_cores=2, cells_per_core=80)
    spatial_slice_cluster(adata, eps=120, min_samples=5, min_cells=30)
    adata.uns["cell_boundaries"] = {
        "vertex_x": np.arange(adata.n_obs * 4, dtype=np.float32),
        "vertex_y": np.arange(adata.n_obs * 4, dtype=np.float32),
        "n_vertices": np.full(adata.n_obs, 4, dtype=np.int32),
        "cell_idx": (np.arange(adata.n_obs) * 4).astype(np.int32),
    }
    parts = split_by_slice(adata, reindex_boundaries=True)
    assert len(parts) == 2
    total = sum(p.n_obs for p in parts.values())
    assert total == adata.n_obs
    for sub in parts.values():
        cb = sub.uns["cell_boundaries"]
        assert len(cb["cell_idx"]) == sub.n_obs
        assert cb["cell_idx"][0] == 0


def test_write_slice_annotation(tmp_path):
    pytest.importorskip("sklearn")
    adata = _make_tma_adata(n_cores=2, cells_per_core=80)
    spatial_slice_cluster(adata, eps=120, min_samples=5, min_cells=30)
    out = tmp_path / "slice_annotation.parquet"
    write_slice_annotation(adata, out)
    df = pd.read_parquet(out)
    assert list(df.columns) == [
        "cell_id",
        "x_centroid",
        "y_centroid",
        "slice_cluster",
        "slice_id",
    ]
    assert len(df) == adata.n_obs


def _make_colony_adata(n_colonies: int = 3, cells_per: int = 60, gap: float = 200.0):
    """Simulate dispersed colonies within one section."""
    rng = np.random.default_rng(1)
    coords_list = []
    for i in range(n_colonies):
        cx = (i % 3) * gap + 500
        cy = (i // 3) * gap + 500
        pts = rng.normal(scale=15, size=(cells_per, 2)) + [cx, cy]
        coords_list.append(pts)
    coords = np.vstack(coords_list)
    adata = AnnData(X=rng.normal(size=(len(coords), 5)).astype(np.float32))
    adata.obsm["spatial"] = coords.astype(np.float32)
    adata.obs_names = [f"c{i}" for i in range(len(coords))]
    adata.var_names = [f"g{i}" for i in range(5)]
    return adata


def test_spatial_colony_cluster():
    pytest.importorskip("sklearn")
    adata = _make_colony_adata(n_colonies=3, cells_per=50)
    spatial_colony_cluster(adata, eps=40, min_samples=5, min_cluster_size=30)
    assert "spatial_colony" in adata.obs
    assigned = adata.obs["spatial_colony"].dropna()
    assert len(assigned.unique()) == 3
    assert adata.uns["spatial_colony_summary"]["n_colonies"] == 3


def test_mark_colony_centroids_and_distance():
    pytest.importorskip("sklearn")
    adata = _make_colony_adata(n_colonies=2, cells_per=40)
    spatial_colony_cluster(adata, eps=40, min_samples=5, min_cluster_size=20)
    mark_colony_centroids(adata, centroid_label="GCC")
    assert (adata.obs["cell_centroid_type"] == "GCC").sum() == 2
    distance_to_nearest_centroids(adata, distance_key="distance_to_nearest_gcc")
    assert "distance_to_nearest_gcc" in adata.obs
    assert (adata.obs["distance_to_nearest_gcc"] >= 0).all()
