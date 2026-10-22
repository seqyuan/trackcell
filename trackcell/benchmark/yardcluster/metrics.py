"""Metrics for YardCluster benchmark."""

from __future__ import annotations

from typing import Any

import numpy as np
from anndata import AnnData
from sklearn.metrics import adjusted_rand_score


def spatial_coherence(
    adata: AnnData,
    cluster_key: str,
    *,
    spatial_key: str = "spatial_connectivities",
    k: int = 15,
) -> float:
    """
    Fraction of each cell's top spatial neighbors sharing the same cluster label.

    Uses ``obsp[spatial_key]`` when present; otherwise returns NaN.
    """
    if spatial_key not in adata.obsp:
        return float("nan")
    if cluster_key not in adata.obs:
        return float("nan")

    weights = adata.obsp[spatial_key].tocsr()
    labels = adata.obs[cluster_key].astype(str).to_numpy()
    n_obs = adata.n_obs
    scores = np.zeros(n_obs, dtype=np.float64)

    for u in range(n_obs):
        start, end = weights.indptr[u], weights.indptr[u + 1]
        if start == end:
            scores[u] = 1.0
            continue
        nbr = weights.indices[start:end]
        nbr = nbr[nbr != u]
        if nbr.size == 0:
            scores[u] = 1.0
            continue
        if nbr.size > k:
            w = weights.data[start:end]
            nbr_local = weights.indices[start:end]
            mask = nbr_local != u
            w = w[mask]
            nbr_local = nbr_local[mask]
            top = np.argsort(w)[-k:]
            nbr = nbr_local[top]
        same = labels[nbr] == labels[u]
        scores[u] = float(same.mean())

    return float(scores.mean())


def evaluate_clustering(
    adata: AnnData,
    pred_key: str,
    truth_key: str,
    *,
    spatial_key: str = "spatial_connectivities",
) -> dict[str, Any]:
    """ARI vs ground truth and spatial coherence of predicted labels."""
    pred = adata.obs[pred_key].astype(str)
    truth = adata.obs[truth_key].astype(str)
    return {
        "pred_key": pred_key,
        "truth_key": truth_key,
        "n_clusters": int(pred.nunique()),
        "ari": float(adjusted_rand_score(truth, pred)),
        "spatial_coherence": spatial_coherence(
            adata, pred_key, spatial_key=spatial_key
        ),
    }
