"""Integrate batch pairs using anchor vectors (Seurat ``IntegrateData`` SCT path)."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from ._annoy_nn import nn_helper


def integrate_pair(
    data_ref: np.ndarray,
    data_query: np.ndarray,
    anchors: pd.DataFrame,
    *,
    k_weight: int = 100,
    sd_weight: float = 1.0,
    dims: int = 30,
    n_trees: int = 50,
) -> np.ndarray:
    """
    Integrate query batch toward reference using anchors (``RunIntegration`` for one pair).

    ``anchors`` must use local indices: ``cell1`` in ref, ``cell2`` in query.
    """
    if anchors.empty:
        return data_query.copy()

    cell1 = anchors["cell1"].astype(int).to_numpy()
    cell2 = anchors["cell2"].astype(int).to_numpy()
    integration_matrix = data_query[cell2] - data_ref[cell1]

    merged = np.vstack([data_ref, data_query])
    merged_centered = merged - merged.mean(axis=0, keepdims=True)
    n_comp = min(dims, merged.shape[0] - 1, merged.shape[1])
    n_comp = max(n_comp, 1)
    dr = PCA(n_components=n_comp, random_state=0).fit_transform(merged_centered)
    dr_query = dr[len(data_ref) :]

    anchor_query_cells = np.unique(cell2)
    if len(anchor_query_cells) == 0:
        return data_query.copy()

    k = min(k_weight, len(anchor_query_cells))
    anchor_dr = dr_query[anchor_query_cells]
    search_k = n_trees * k
    nbr_idx, dists = nn_helper(
        anchor_dr,
        query=dr_query,
        k=k,
        n_trees=n_trees,
        search_k=search_k,
    )

    corrected = data_query.copy()
    scores = anchors["score"].to_numpy() if "score" in anchors.columns else np.ones(len(anchors))
    for q_idx in range(data_query.shape[0]):
        d_row = dists[q_idx]
        max_d = float(d_row.max()) if d_row.size else 1.0
        weights = 1.0 - d_row / max(max_d, 1e-9)
        weights = np.exp(-np.square(weights) / max(sd_weight, 1e-9))
        if weights.sum() <= 0:
            continue
        weights = weights / weights.sum()
        delta = np.zeros(data_query.shape[1], dtype=np.float64)
        for w, anchor_local in zip(weights, nbr_idx[q_idx]):
            anchor_cell = int(anchor_query_cells[anchor_local])
            rows = np.where(cell2 == anchor_cell)[0]
            if rows.size == 0:
                continue
            best = rows[np.argmax(scores[rows])]
            delta += w * integration_matrix[best]
        corrected[q_idx] += delta
    return corrected


def stack_integrated(
    integrated: dict[str, np.ndarray],
    batch_labels: list[str],
    cell_offsets: dict[str, int],
) -> np.ndarray:
    n_cells = sum(integrated[b].shape[0] for b in batch_labels)
    n_feat = integrated[batch_labels[0]].shape[1]
    out = np.zeros((n_cells, n_feat), dtype=np.float32)
    for label in batch_labels:
        off = cell_offsets[label]
        block = integrated[label]
        out[off : off + block.shape[0]] = block
    return out
