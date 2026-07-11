"""Sample tree construction and recursive batch merging (Seurat ``IntegrateData``)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

from .integrate_data import DEFAULT_WEIGHT_QUERY_CHUNK_SIZE, integrate_pair
from ._xgboost_integration import integrate_pair_xgb


@dataclass
class _IntegrateCluster:
    data: np.ndarray
    batch_indices: set[int] = field(default_factory=set)
    global_indices: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))


def count_anchor_similarity(
    anchors: pd.DataFrame,
    n_batches: int,
    n_cells: list[int],
) -> np.ndarray:
    """Seurat ``CountAnchors`` similarity matrix (anchor count / min n_cells)."""
    sim = np.zeros((n_batches, n_batches), dtype=np.float64)
    for i in range(n_batches):
        for j in range(i + 1, n_batches):
            mask = (
                ((anchors["dataset1"] == i) & (anchors["dataset2"] == j))
                | ((anchors["dataset1"] == j) & (anchors["dataset2"] == i))
            )
            score = float(mask.sum()) / min(n_cells[i], n_cells[j])
            sim[i, j] = score
            sim[j, i] = score
    return sim


def build_sample_tree(similarity: np.ndarray) -> np.ndarray:
    """Seurat ``BuildSampleTree`` via hierarchical clustering on 1/similarity."""
    n = similarity.shape[0]
    if n <= 1:
        return np.zeros((0, 2), dtype=int)
    if n == 2:
        return np.array([[0, 1]], dtype=int)

    dist = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            s = max(similarity[i, j], 1e-9)
            dist[i, j] = 1.0 / s
            dist[j, i] = dist[i, j]
    condensed = squareform(dist, checks=False)
    z = linkage(condensed, method="average")
    return z[:, :2].astype(int)


def _filter_pair_anchors(
    anchors: pd.DataFrame,
    ref_obj: _IntegrateCluster,
    query_obj: _IntegrateCluster,
) -> pd.DataFrame:
    """Select anchors between two composite clusters and map to local indices."""
    batches_a = ref_obj.batch_indices
    batches_b = query_obj.batch_indices
    global_to_local_a = {int(g): i for i, g in enumerate(ref_obj.global_indices)}
    global_to_local_b = {int(g): i for i, g in enumerate(query_obj.global_indices)}
    rows: list[dict[str, float]] = []
    scores = anchors["score"].to_numpy() if "score" in anchors.columns else np.ones(len(anchors))

    for k, (_, row) in enumerate(anchors.iterrows()):
        d1, d2 = int(row["dataset1"]), int(row["dataset2"])
        c1, c2 = int(row["cell1"]), int(row["cell2"])
        score = float(scores[k])
        if d1 in batches_a and d2 in batches_b and c1 in global_to_local_a and c2 in global_to_local_b:
            rows.append(
                {
                    "cell1": global_to_local_a[c1],
                    "cell2": global_to_local_b[c2],
                    "score": score,
                }
            )
        elif d2 in batches_a and d1 in batches_b and c2 in global_to_local_a and c1 in global_to_local_b:
            rows.append(
                {
                    "cell1": global_to_local_a[c2],
                    "cell2": global_to_local_b[c1],
                    "score": score,
                }
            )
    return pd.DataFrame(rows)


def integrate_batches_sample_tree(
    batch_data: dict[str, np.ndarray],
    batch_labels: list[str],
    anchors: pd.DataFrame,
    global_indices: dict[str, np.ndarray],
    *,
    k_weight: int = 100,
    sd_weight: float = 1.0,
    dims: int = 30,
    n_trees: int = 50,
    n_trees_weight: Optional[int] = None,
    weight_query_chunk_size: int = DEFAULT_WEIGHT_QUERY_CHUNK_SIZE,
    integration_dtype: Literal["float32", "float64"] = "float32",
    preserve_order: bool = False,
    correction_method: str = "seurat",
    xgb_n_delta_pcs: int = 20,
    xgb_n_estimators: int = 150,
    xgb_max_depth: int = 4,
    xgb_seed: int = 0,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """
    Recursively merge batches using a sample tree (Seurat ``PairwiseIntegrateReference``).

    Returns integrated per-original-batch data and the scipy linkage merge matrix.
    """
    n_batches = len(batch_labels)
    n_cells = [batch_data[label].shape[0] for label in batch_labels]
    similarity = count_anchor_similarity(anchors, n_batches, n_cells)
    merge_order = build_sample_tree(similarity)

    clusters: dict[int, _IntegrateCluster] = {}
    for i, label in enumerate(batch_labels):
        clusters[i] = _IntegrateCluster(
            data=batch_data[label].copy(),
            batch_indices={i},
            global_indices=global_indices[label].copy(),
        )

    next_id = n_batches
    for step in range(n_batches - 1):
        left_id, right_id = int(merge_order[step, 0]), int(merge_order[step, 1])
        obj_a = clusters[left_id]
        obj_b = clusters[right_id]

        ref_obj, query_obj = obj_a, obj_b
        if (not preserve_order) and query_obj.data.shape[0] > ref_obj.data.shape[0]:
            ref_obj, query_obj = query_obj, ref_obj

        local_anchors = _filter_pair_anchors(anchors, ref_obj, query_obj)
        if correction_method == "xgboost":
            corrected_query = integrate_pair_xgb(
                ref_obj.data,
                query_obj.data,
                local_anchors,
                n_delta_pcs=xgb_n_delta_pcs,
                n_estimators=xgb_n_estimators,
                max_depth=xgb_max_depth,
                seed=xgb_seed,
            )
        else:
            corrected_query = integrate_pair(
                ref_obj.data,
                query_obj.data,
                local_anchors,
                k_weight=k_weight,
                sd_weight=sd_weight,
                dims=dims,
                n_trees=n_trees,
                n_trees_weight=n_trees_weight,
                weight_query_chunk_size=weight_query_chunk_size,
                integration_dtype=integration_dtype,
            )
        merged = _IntegrateCluster(
            data=np.vstack([ref_obj.data, corrected_query]),
            batch_indices=ref_obj.batch_indices | query_obj.batch_indices,
            global_indices=np.concatenate([ref_obj.global_indices, query_obj.global_indices]),
        )
        clusters[next_id] = merged
        del clusters[left_id], clusters[right_id]
        next_id += 1

    final = clusters[next_id - 1]
    global_to_row = {int(g): i for i, g in enumerate(final.global_indices)}
    integrated: dict[str, np.ndarray] = {}
    for label in batch_labels:
        gidx = global_indices[label]
        rows = [global_to_row[int(g)] for g in gidx]
        integrated[label] = final.data[rows]
    return integrated, merge_order
