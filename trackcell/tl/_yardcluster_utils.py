"""Helpers for YardCluster: presets, caching, sketch sampling, hybrid graphs."""

from __future__ import annotations

from typing import Any, Callable, Literal, Optional

import numpy as np
import scipy.sparse as sps
from anndata import AnnData

PlatformType = Literal["visium_hd", "sto", "xenium", "xenium_tma"]
GraphModeType = Literal["expr", "spatial", "hybrid"]

_GLOBAL_DEFAULTS: dict[str, Any] = {
    "k_spatial": 15,
    "k_expr": 50,
    "n_pcs": 20,
    "resolution": 1.0,
    "n_top_genes": 2000,
}

PLATFORM_PRESETS: dict[str, dict[str, Any]] = {
    "visium_hd": {
        "k_spatial": 15,
        "k_expr": 50,
        "n_pcs": 20,
        "resolution": 1.0,
        "n_top_genes": 2000,
    },
    "sto": {
        "k_spatial": 20,
        "k_expr": 50,
        "n_pcs": 25,
        "resolution": 0.8,
        "n_top_genes": 2000,
    },
    "xenium": {
        "k_spatial": 12,
        "k_expr": 30,
        "n_pcs": 20,
        "resolution": 1.2,
        "n_top_genes": 2000,
    },
    "xenium_tma": {
        # Multi-core TMA (e.g. CCHD cervical cancer Xenium): many slices, ~0.21 µm/px
        "k_spatial": 10,
        "k_expr": 30,
        "n_pcs": 25,
        "resolution": 0.9,
        "n_top_genes": 3000,
    },
}


def resolve_platform_params(
    platform: Optional[str],
    *,
    k_spatial: Optional[int],
    k_expr: Optional[int],
    n_pcs: Optional[int],
    resolution: Optional[float],
    n_top_genes: Optional[int],
) -> dict[str, Any]:
    """Fill ``None`` params from platform preset or global defaults."""
    base = dict(_GLOBAL_DEFAULTS)
    if platform is not None:
        if platform not in PLATFORM_PRESETS:
            raise ValueError(
                f"Unknown platform '{platform}'. "
                f"Choose from: {', '.join(PLATFORM_PRESETS)}"
            )
        base.update(PLATFORM_PRESETS[platform])

    return {
        "k_spatial": k_spatial if k_spatial is not None else base["k_spatial"],
        "k_expr": k_expr if k_expr is not None else base["k_expr"],
        "n_pcs": n_pcs if n_pcs is not None else base["n_pcs"],
        "resolution": resolution if resolution is not None else base["resolution"],
        "n_top_genes": n_top_genes if n_top_genes is not None else base["n_top_genes"],
    }


def spatial_neighbors_cached(
    adata: AnnData,
    *,
    k: int,
    batch_key: Optional[str],
    key_added: str,
    build_fn: Callable[..., None],
) -> bool:
    """Return True if an existing spatial graph was reused."""
    conn_key = f"{key_added}_connectivities"
    params_key = f"{key_added}_neighbors_params"
    existing = adata.uns.get(params_key)
    if (
        conn_key in adata.obsp
        and existing
        and existing.get("k") == k
        and existing.get("batch_key") == batch_key
    ):
        return True
    build_fn(adata, k=k, batch_key=batch_key, key_added=key_added)
    return False


def neighborhood_features_cached(
    adata: AnnData,
    *,
    expected_params: dict[str, Any],
    key_added: str,
    compute_fn: Callable[[], None],
) -> bool:
    """Return True if neighborhood M/G were reused."""
    params_key = f"{key_added}_features_params"
    m_key = f"{key_added}_M"
    existing = adata.uns.get(params_key)
    if m_key not in adata.obsm or not existing:
        compute_fn()
        return False

    for key, val in expected_params.items():
        if existing.get(key) != val:
            compute_fn()
            return False
    return True


def spatial_stratified_sketch_indices(
    coords: np.ndarray,
    n_sketch: int,
    *,
    n_bins: int = 50,
    seed: int = 0,
) -> np.ndarray:
    """
    Sample cells with approximate spatial coverage using grid stratification.

    Each occupied spatial bin contributes roughly equal quota before filling remainder.
    """
    coords = np.asarray(coords, dtype=np.float64)
    n_obs = coords.shape[0]
    if n_sketch >= n_obs:
        return np.arange(n_obs, dtype=np.int64)

    rng = np.random.default_rng(seed)
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    span = np.maximum(maxs - mins, 1e-9)
    norm = (coords - mins) / span
    bins = np.clip((norm * n_bins).astype(int), 0, n_bins - 1)
    bin_ids = bins[:, 0] * n_bins + bins[:, 1]

    selected: list[int] = []
    unique_bins = np.unique(bin_ids)
    per_bin = max(1, n_sketch // max(len(unique_bins), 1))

    for b in unique_bins:
        idx = np.where(bin_ids == b)[0]
        take = min(per_bin, len(idx))
        if take > 0:
            selected.extend(rng.choice(idx, size=take, replace=False).tolist())

    selected = list(dict.fromkeys(selected))
    if len(selected) < n_sketch:
        remaining = np.setdiff1d(np.arange(n_obs), np.array(selected, dtype=int), assume_unique=False)
        extra = min(n_sketch - len(selected), len(remaining))
        if extra > 0:
            selected.extend(rng.choice(remaining, size=extra, replace=False).tolist())

    if len(selected) > n_sketch:
        selected = rng.choice(np.array(selected, dtype=int), size=n_sketch, replace=False).tolist()

    return np.asarray(selected, dtype=np.int64)


def _row_normalize_csr(mat: sps.csr_matrix) -> sps.csr_matrix:
    mat = mat.tocsr().astype(np.float64)
    row_sums = np.asarray(mat.sum(axis=1)).ravel()
    for u in range(mat.shape[0]):
        start, end = mat.indptr[u], mat.indptr[u + 1]
        denom = row_sums[u]
        if denom > 0:
            mat.data[start:end] /= denom
    return mat


def mix_spatial_expr_connectivities(
    spatial: sps.csr_matrix,
    expr: sps.csr_matrix,
    alpha: float,
) -> sps.csr_matrix:
    """Convex mix of row-normalized spatial and expression connectivities."""
    alpha = float(alpha)
    if not 0.0 <= alpha <= 1.0:
        raise ValueError("hybrid_alpha must be in [0, 1].")
    if spatial.shape != expr.shape:
        raise ValueError("spatial and expr connectivities must share shape.")
    s = _row_normalize_csr(spatial)
    e = _row_normalize_csr(expr)
    mixed = ((1.0 - alpha) * e + alpha * s).tocsr()
    return _row_normalize_csr(mixed)
