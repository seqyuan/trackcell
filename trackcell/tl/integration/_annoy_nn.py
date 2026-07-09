"""Annoy nearest-neighbor search (Seurat ``NNHelper`` with method='annoy')."""

from __future__ import annotations

from typing import Optional

import numpy as np

try:
    from annoy import AnnoyIndex
except ImportError:  # pragma: no cover - optional at install time
    AnnoyIndex = None


def annoy_available() -> bool:
    return AnnoyIndex is not None


def nn_helper(
    data: np.ndarray,
    query: Optional[np.ndarray] = None,
    *,
    k: int,
    n_trees: int = 50,
    search_k: int = -1,
    metric: str = "euclidean",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Find k nearest neighbors using Annoy (Seurat ``AnnoyNN``).

    Parameters
    ----------
    data
        Reference points (n_ref x dims).
    query
        Query points (n_query x dims). Defaults to ``data``.
    """
    if AnnoyIndex is None:
        raise ImportError(
            "annoy is required for Annoy-based integration. "
            "Install with: pip install annoy"
        )
    if query is None:
        query = data
    data = np.asarray(data, dtype=np.float32)
    query = np.asarray(query, dtype=np.float32)
    if data.shape[0] == 0 or query.shape[0] == 0:
        return (
            np.zeros((query.shape[0], 0), dtype=int),
            np.zeros((query.shape[0], 0), dtype=np.float32),
        )

    dim = data.shape[1]
    k_eff = min(k, data.shape[0])
    if search_k < 0:
        search_k = n_trees * k_eff
    index = AnnoyIndex(dim, metric)
    for i, vec in enumerate(data):
        index.add_item(i, vec.tolist())
    index.build(n_trees)

    idx = np.zeros((query.shape[0], k_eff), dtype=np.int64)
    dists = np.zeros((query.shape[0], k_eff), dtype=np.float32)
    for i, vec in enumerate(query):
        neighbors, distances = index.get_nns_by_vector(
            vec.tolist(),
            k_eff,
            search_k=search_k,
            include_distances=True,
        )
        idx[i] = neighbors
        dists[i] = distances
    return idx, dists
