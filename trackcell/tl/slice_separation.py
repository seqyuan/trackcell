"""
DBSCAN-based spatial grouping for Xenium / Visium workflows.

Two CCHD use cases share the same sklearn DBSCAN backend:

1. **TMA slice separation** — split an entire cassette region into physical
   cores (``spatial_slice_cluster``, eps≈80 µm, min_cells≈1000).
2. **Micro-region colonies** — within a BANKSY-filtered cell type (GC B cells,
   tumor nests), split dispersed spatial colonies (``spatial_colony_cluster``,
   eps≈50 µm on Xenium or ≈14 on Visium HD bins, min_cluster_size≈40).
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
from anndata import AnnData

try:
    from sklearn.cluster import DBSCAN
except ImportError:  # pragma: no cover - sklearn is a scanpy dependency
    DBSCAN = None  # type: ignore[misc, assignment]


def _require_dbscan():
    if DBSCAN is None:
        raise ImportError(
            "spatial_slice_cluster requires scikit-learn. "
            "Install with: pip install scikit-learn"
        )


def _spatial_coords(adata: AnnData, spatial_key: str = "spatial") -> np.ndarray:
    if spatial_key in adata.obsm:
        coords = np.asarray(adata.obsm[spatial_key], dtype=float)
    elif {"x_centroid", "y_centroid"}.issubset(adata.obs.columns):
        coords = adata.obs[["x_centroid", "y_centroid"]].to_numpy(dtype=float)
    else:
        raise ValueError(
            f"Spatial coordinates not found. Provide obsm['{spatial_key}'] "
            "or obs columns 'x_centroid' and 'y_centroid'."
        )
    if coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError(f"Expected N×2 spatial coordinates, got shape {coords.shape}")
    return coords


def dbscan_slice_labels(
    coords: np.ndarray,
    eps: float = 80.0,
    min_samples: int = 10,
    min_cells: int = 1000,
    slice_start: int = 1,
    slice_prefix: str = "S",
    debris_label: str = "debris",
    sort_spatial: bool = True,
    n_jobs: int = -1,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Run DBSCAN on 2D coordinates and assign human-readable slice IDs.

    Clusters with fewer than ``min_cells`` observations are labeled ``debris``.
    Valid clusters receive sequential IDs ``{slice_prefix}{NNN}`` (e.g. S001).

    Parameters
    ----------
    coords
        N×2 array of spatial coordinates (typically µm).
    eps
        DBSCAN neighborhood radius (µm for Xenium centroids).
    min_samples
        Minimum neighbors to form a dense region (DBSCAN ``min_samples``).
    min_cells
        Minimum cells per cluster to keep as a valid slice.
    slice_start
        Starting index for numbered slice IDs.
    slice_prefix
        Prefix for slice IDs (default ``S`` → S001, S002, …).
    debris_label
        Label for noise and undersized clusters.
    sort_spatial
        If True, order slices by mean Y then mean X before numbering.
    n_jobs
        Parallel jobs passed to sklearn DBSCAN (-1 = all cores).

    Returns
    -------
    slice_cluster
        Integer DBSCAN labels (-1 = noise).
    slice_id
        String slice IDs including ``debris_label`` for discarded cells.
    """
    _require_dbscan()
    coords = np.asarray(coords, dtype=float)
    if coords.shape[0] == 0:
        return (
            np.array([], dtype=int),
            np.array([], dtype=object),
        )

    clustering = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=n_jobs)
    labels = clustering.fit_predict(coords)

    cluster_sizes = pd.Series(labels[labels >= 0]).value_counts()
    valid_clusters = cluster_sizes[cluster_sizes >= min_cells].index.tolist()

    if sort_spatial and valid_clusters:
        cluster_info = []
        for cid in valid_clusters:
            mask = labels == cid
            cl_coords = coords[mask]
            cluster_info.append(
                {
                    "cluster": cid,
                    "y_mean": float(cl_coords[:, 1].mean()),
                    "x_mean": float(cl_coords[:, 0].mean()),
                }
            )
        cluster_info.sort(key=lambda d: (d["y_mean"], d["x_mean"]))
        ordered_clusters = [info["cluster"] for info in cluster_info]
    else:
        ordered_clusters = sorted(valid_clusters)

    slice_map = {
        cid: f"{slice_prefix}{slice_start + i:03d}"
        for i, cid in enumerate(ordered_clusters)
    }

    slice_ids = np.full(len(coords), debris_label, dtype=object)
    for cid, sid in slice_map.items():
        slice_ids[labels == cid] = sid

    return labels.astype(int), slice_ids


def spatial_slice_cluster(
    adata: AnnData,
    eps: float = 80.0,
    min_samples: int = 10,
    min_cells: int = 1000,
    slice_start: int = 1,
    slice_prefix: str = "S",
    debris_label: str = "debris",
    spatial_key: str = "spatial",
    key_added: str = "slice_id",
    cluster_key: str = "slice_cluster",
    sort_spatial: bool = True,
    n_jobs: int = -1,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Annotate an AnnData object with DBSCAN-based tissue slice labels.

    Intended for Xenium TMA regions where multiple tissue cores share one
    coordinate system. Operates on ``obsm[spatial_key]`` or
    ``obs[['x_centroid', 'y_centroid']]``.

    Parameters
    ----------
    adata
        AnnData with spatial coordinates.
    eps, min_samples, min_cells, slice_start, slice_prefix, debris_label
        Passed to :func:`dbscan_slice_labels`.
    spatial_key
        Key in ``adata.obsm`` for coordinates.
    key_added
        ``obs`` column for string slice IDs.
    cluster_key
        ``obs`` column for raw DBSCAN cluster integers.
    sort_spatial
        Sort slices by spatial position before numbering.
    n_jobs
        DBSCAN parallel jobs.
    copy
        If True, return a modified copy instead of updating in place.

    Returns
    -------
    AnnData or None
        Copy when ``copy=True``; otherwise None (in-place update).
    """
    target = adata.copy() if copy else adata
    coords = _spatial_coords(target, spatial_key=spatial_key)

    clusters, slice_ids = dbscan_slice_labels(
        coords,
        eps=eps,
        min_samples=min_samples,
        min_cells=min_cells,
        slice_start=slice_start,
        slice_prefix=slice_prefix,
        debris_label=debris_label,
        sort_spatial=sort_spatial,
        n_jobs=n_jobs,
    )

    target.obs[cluster_key] = clusters
    target.obs[key_added] = pd.Categorical(slice_ids)

    params = {
        "eps": eps,
        "min_samples": min_samples,
        "min_cells": min_cells,
        "slice_start": slice_start,
        "slice_prefix": slice_prefix,
        "debris_label": debris_label,
        "sort_spatial": sort_spatial,
    }
    target.uns[f"{key_added}_params"] = params

    n_slices = len([s for s in pd.unique(slice_ids) if s != debris_label])
    target.uns[f"{key_added}_summary"] = {
        "n_slices": n_slices,
        "n_debris": int((slice_ids == debris_label).sum()),
        "slice_ids": sorted(s for s in pd.unique(slice_ids) if s != debris_label),
    }

    return target if copy else None


def slice_cluster_summary(
    adata: AnnData,
    slice_key: str = "slice_id",
    spatial_key: str = "spatial",
    debris_label: str = "debris",
) -> pd.DataFrame:
    """
    Summarize slice assignments with cell counts and spatial extents.

    Returns
    -------
    DataFrame
        One row per slice with n_cells, x/y min/max/span, and center.
    """
    if slice_key not in adata.obs.columns:
        raise ValueError(f"`{slice_key}` not found in adata.obs")
    coords = _spatial_coords(adata, spatial_key=spatial_key)

    rows = []
    for sid in sorted(adata.obs[slice_key].unique(), key=str):
        if sid == debris_label:
            continue
        mask = (adata.obs[slice_key] == sid).to_numpy()
        xy = coords[mask]
        rows.append(
            {
                "slice_id": sid,
                "n_cells": int(mask.sum()),
                "x_min": round(float(xy[:, 0].min()), 1),
                "x_max": round(float(xy[:, 0].max()), 1),
                "y_min": round(float(xy[:, 1].min()), 1),
                "y_max": round(float(xy[:, 1].max()), 1),
                "x_span": round(float(xy[:, 0].max() - xy[:, 0].min()), 1),
                "y_span": round(float(xy[:, 1].max() - xy[:, 1].min()), 1),
                "x_center": round(float(xy[:, 0].mean()), 1),
                "y_center": round(float(xy[:, 1].mean()), 1),
            }
        )
    return pd.DataFrame(rows)


def _reindex_cell_boundaries(sub: AnnData) -> None:
    """Rebuild compact cell_boundaries vertex offsets after subsetting."""
    if "cell_boundaries" not in sub.uns:
        return
    cb = sub.uns["cell_boundaries"]
    old_idx = cb["cell_idx"]
    old_nv = cb["n_vertices"]
    old_vx = cb["vertex_x"]
    old_vy = cb["vertex_y"]
    desc = cb.get("description", "")

    new_nv_list: list[int] = []
    new_vx_list: list[np.ndarray] = []
    new_vy_list: list[np.ndarray] = []
    new_cell_idx: list[int] = []
    cum = 0

    for i in range(sub.n_obs):
        start = int(old_idx[i])
        if start >= 0:
            nv = int(old_nv[i])
            new_cell_idx.append(cum)
            new_nv_list.append(nv)
            cum += nv
            new_vx_list.append(old_vx[start : start + nv])
            new_vy_list.append(old_vy[start : start + nv])
        else:
            new_cell_idx.append(-1)
            new_nv_list.append(0)

    if new_vx_list:
        sub.uns["cell_boundaries"] = {
            "vertex_x": np.concatenate(new_vx_list).astype(np.float32),
            "vertex_y": np.concatenate(new_vy_list).astype(np.float32),
            "n_vertices": np.array(new_nv_list, dtype=np.int32),
            "cell_idx": np.array(new_cell_idx, dtype=np.int32),
            "description": desc,
        }
    else:
        del sub.uns["cell_boundaries"]


def split_by_slice(
    adata: AnnData,
    slice_key: str = "slice_id",
    exclude: Sequence[str] = ("debris",),
    copy: bool = True,
    reindex_boundaries: bool = True,
) -> dict[str, AnnData]:
    """
    Split an AnnData into one object per slice ID.

    Parameters
    ----------
    adata
        Annotated AnnData (e.g. after :func:`spatial_slice_cluster`).
    slice_key
        Column in ``adata.obs`` with slice labels.
    exclude
        Slice labels to skip (default: ``debris``).
    copy
        If True, each slice is a deep copy; if False, views where possible.
    reindex_boundaries
        Rebuild ``uns['cell_boundaries']`` vertex offsets for each subset.

    Returns
    -------
    dict
        Mapping ``slice_id`` → AnnData subset.
    """
    if slice_key not in adata.obs.columns:
        raise ValueError(f"`{slice_key}` not found in adata.obs")

    exclude_set = set(exclude)
    out: dict[str, AnnData] = {}
    for sid in sorted(adata.obs[slice_key].unique(), key=str):
        if sid in exclude_set:
            continue
        mask = adata.obs[slice_key] == sid
        sub = adata[mask].copy() if copy else adata[mask]
        if reindex_boundaries:
            _reindex_cell_boundaries(sub)
        out[str(sid)] = sub
    return out


def write_slice_annotation(
    adata: AnnData,
    path: Union[str, Path],
    slice_key: str = "slice_id",
    cluster_key: str = "slice_cluster",
    spatial_key: str = "spatial",
) -> Path:
    """
    Write slice labels to parquet (CCHD-compatible ``slice_annotation.parquet``).

    Columns: cell_id, x_centroid, y_centroid, slice_cluster, slice_id.
    """
    path = Path(path)
    if slice_key not in adata.obs.columns:
        raise ValueError(f"Run spatial_slice_cluster first; missing obs['{slice_key}']")

    coords = _spatial_coords(adata, spatial_key=spatial_key)
    df = pd.DataFrame(
        {
            "cell_id": adata.obs_names,
            "x_centroid": coords[:, 0],
            "y_centroid": coords[:, 1],
            "slice_cluster": adata.obs[cluster_key].values,
            "slice_id": adata.obs[slice_key].astype(str).values,
        }
    )
    df.to_parquet(path, index=False)
    return path


def spatial_colony_cluster(
    adata: AnnData,
    eps: float = 50.0,
    min_samples: int = 5,
    min_cluster_size: int = 40,
    spatial_key: str = "spatial",
    cluster_key: str = "spatial_colony",
    metric: str = "euclidean",
    n_jobs: int = -1,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    DBSCAN colony labels for dispersed micro-regions within one tissue section.

    Intended for **pre-filtered** subsets (e.g. BANKSY ``Germinal Center B
    Cells`` or tumor ROI clusters). Clusters smaller than ``min_cluster_size``
    and DBSCAN noise (``-1``) are stored as ``NaN`` in ``obs[cluster_key]``.

    Default parameters match CCHD ``54_1_GC_dist.ipynb`` (Xenium GC colonies).

    Parameters
    ----------
    adata
        AnnData subset with spatial coordinates (typically after BANKSY).
    eps, min_samples, min_cluster_size
        DBSCAN radius and density filters; small clusters → NaN.
    spatial_key
        Key in ``adata.obsm`` for coordinates.
    cluster_key
        ``obs`` column for colony IDs (float categorical with NaN).
    metric
        Distance metric for DBSCAN.
    n_jobs
        DBSCAN parallel jobs.
    copy
        If True, return a modified copy instead of updating in place.

    Returns
    -------
    AnnData or None
        Copy when ``copy=True``; otherwise None (in-place update).
    """
    _require_dbscan()
    target = adata.copy() if copy else adata
    coords = _spatial_coords(target, spatial_key=spatial_key)

    clustering = DBSCAN(eps=eps, min_samples=min_samples, metric=metric, n_jobs=n_jobs)
    labels = clustering.fit_predict(coords)

    filtered = labels.astype(object).copy()
    if labels.max() >= 0:
        unique, counts = np.unique(labels[labels >= 0], return_counts=True)
        for cid, count in zip(unique, counts):
            if count < min_cluster_size:
                filtered[filtered == cid] = np.nan
    filtered[filtered == -1] = np.nan

    target.obs[cluster_key] = pd.Categorical(filtered)
    target.uns[f"{cluster_key}_params"] = {
        "eps": eps,
        "min_samples": min_samples,
        "min_cluster_size": min_cluster_size,
        "metric": metric,
    }
    valid = filtered[~pd.isna(filtered)]
    target.uns[f"{cluster_key}_summary"] = {
        "n_colonies": len(np.unique(valid)) if len(valid) else 0,
        "n_assigned": int((~pd.isna(filtered)).sum()),
        "n_noise_or_small": int(pd.isna(filtered).sum()),
    }
    return target if copy else None


def mark_colony_centroids(
    adata: AnnData,
    colony_key: str = "spatial_colony",
    centroid_key: str = "cell_centroid_type",
    centroid_label: str = "GCC",
    spatial_key: str = "spatial",
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Mark one representative cell per colony (nearest to the geometric centroid).

    Stores ``obs[centroid_key]`` and colony metadata in
    ``uns['colony_centroids']``. Default label ``GCC`` follows CCHD GC analysis.

    Parameters
    ----------
    adata
        AnnData after :func:`spatial_colony_cluster`.
    colony_key
        Column with colony IDs (NaN = unassigned).
    centroid_key
        New column for centroid cell labels.
    centroid_label
        Label assigned to centroid cells (e.g. ``GCC`` or ``TNC``).
    spatial_key
        Key in ``adata.obsm`` for coordinates.
    copy
        If True, return a modified copy instead of updating in place.

    Returns
    -------
    AnnData or None
        Copy when ``copy=True``; otherwise None (in-place update).
    """
    target = adata.copy() if copy else adata
    if colony_key not in target.obs.columns:
        raise ValueError(f"Column '{colony_key}' not found in adata.obs")

    coords = _spatial_coords(target, spatial_key=spatial_key)
    colonies = target.obs[colony_key]
    centroid_labels = np.full(target.n_obs, np.nan, dtype=object)
    centroid_info: list[dict] = []

    for colony_id in colonies.dropna().unique():
        mask = colonies == colony_id
        indices = np.where(mask.to_numpy())[0]
        colony_coords = coords[mask.to_numpy()]
        if len(indices) == 0:
            continue
        center = colony_coords.mean(axis=0)
        dists = np.linalg.norm(colony_coords - center, axis=1)
        nearest = indices[int(np.argmin(dists))]
        centroid_labels[nearest] = centroid_label
        centroid_info.append(
            {
                "colony_id": colony_id,
                "centroid_cell_index": int(nearest),
                "colony_size": len(indices),
                "centroid_coord": center.tolist(),
                "min_distance": float(dists.min()),
            }
        )

    target.obs[centroid_key] = pd.Categorical(centroid_labels)
    target.uns["colony_centroids"] = centroid_info
    return target if copy else None


def distance_to_nearest_centroids(
    adata: AnnData,
    centroid_key: str = "cell_centroid_type",
    distance_key: str = "distance_to_nearest_centroid",
    centroid_label: str = "GCC",
    spatial_key: str = "spatial",
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Per-cell distance to the nearest marked colony centroid (cKDTree).

    Parameters
    ----------
    adata
        AnnData after :func:`mark_colony_centroids` (full section or subset).
    centroid_key
        Column with centroid labels.
    distance_key
        New ``obs`` column for distances (same units as ``spatial_key`` coords).
    centroid_label
        Label identifying centroid cells in ``centroid_key``.
    spatial_key
        Key in ``adata.obsm`` for coordinates.
    copy
        If True, return a modified copy instead of updating in place.

    Returns
    -------
    AnnData or None
        Copy when ``copy=True``; otherwise None (in-place update).
    """
    from scipy.spatial import cKDTree

    target = adata.copy() if copy else adata
    if centroid_key not in target.obs.columns:
        raise ValueError(f"Column '{centroid_key}' not found in adata.obs")

    coords = _spatial_coords(target, spatial_key=spatial_key)
    gcc_mask = target.obs[centroid_key] == centroid_label
    gcc_coords = coords[gcc_mask.to_numpy()]
    if len(gcc_coords) == 0:
        raise ValueError(f"No cells labeled '{centroid_label}' in obs['{centroid_key}']")

    tree = cKDTree(gcc_coords)
    distances, _ = tree.query(coords, k=1)
    target.obs[distance_key] = distances
    return target if copy else None
