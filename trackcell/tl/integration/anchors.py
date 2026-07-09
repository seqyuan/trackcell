"""Find integration anchors (Seurat RPCA + Annoy MNN)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ._annoy_nn import nn_helper


def _normalize_anchor_scores(scores: np.ndarray) -> np.ndarray:
    """Seurat ``ScoreAnchors`` quantile normalization."""
    if scores.size == 0:
        return scores
    lo = float(np.quantile(scores, 0.01))
    hi = float(np.quantile(scores, 0.9))
    if hi <= lo:
        return np.ones_like(scores, dtype=np.float64)
    out = (scores - lo) / (hi - lo)
    return np.clip(out, 0.0, 1.0)


def score_anchors_seurat(
    anchors: pd.DataFrame,
    *,
    nnab_idx: np.ndarray,
    nnba_idx: np.ndarray,
    nnaa_idx: np.ndarray | None,
    nnbb_idx: np.ndarray | None,
    n_ref: int,
    k_score: int,
) -> np.ndarray:
    """
    Score anchors by shared neighbor count (Seurat ``ScoreAnchors``).

    ``anchors`` must use stacked indices: batch1 ``0..n_ref-1``, batch2 ``n_ref..``.
    """
    if anchors.empty or k_score <= 0:
        return np.ones(len(anchors), dtype=np.float64)

    scores = np.zeros(len(anchors), dtype=np.float64)
    for idx, row in enumerate(anchors.itertuples(index=False)):
        c1, c2 = int(row.cell1), int(row.cell2)
        local1 = c1
        local2 = c2 - n_ref
        nbr_a = set(int(x) for x in nnab_idx[local1, :k_score])
        nbr_b = set(int(x) for x in nnba_idx[local2, :k_score])
        if nnaa_idx is not None:
            nbr_a.update(int(x) for x in nnaa_idx[local1, :k_score])
        if nnbb_idx is not None:
            nbr_b.update(int(x) + n_ref for x in nnbb_idx[local2, :k_score])
        nbr_a.update(int(x) + n_ref for x in nnab_idx[local1, :k_score])
        scores[idx] = float(len(nbr_a & nbr_b))
    return _normalize_anchor_scores(scores)


def find_pairwise_anchors(
    ref_embedding: np.ndarray,
    query_embedding: np.ndarray,
    *,
    n_ref: int,
    pca1_embedding: np.ndarray | None = None,
    pca2_embedding: np.ndarray | None = None,
    k_anchor: int = 5,
    k_score: int = 30,
    n_trees: int = 50,
) -> pd.DataFrame:
    """
    Find cross-batch anchors with Annoy MNN (Seurat ``FindAnchorPairs`` RPCA path).

    ``ref_embedding`` = projectedpca.ref, ``query_embedding`` = projectedpca.query.
    Layout: rows ``0..n_ref-1`` = batch1, ``n_ref..`` = batch2.
    """
    n_total = ref_embedding.shape[0]
    n_query = n_total - n_ref
    if n_query <= 0:
        return pd.DataFrame(columns=["cell1", "cell2", "score"])

    dims = ref_embedding.shape[1]
    ref_emb = ref_embedding[:, :dims]
    query_emb = query_embedding[:, :dims]

    cells1 = np.arange(n_ref, dtype=int)
    cells2 = np.arange(n_ref, n_total, dtype=int)
    emb1_ref = ref_emb[cells1]
    emb2_ref = ref_emb[cells2]
    emb1_query = query_emb[cells1]
    emb2_query = query_emb[cells2]

    k_neighbor = max(k_anchor, k_score) if k_score > 0 else k_anchor
    search_k = n_trees * k_neighbor

    # nnab / nnba (Seurat FindNN with reduction.2 set)
    nnab_idx, _ = nn_helper(
        emb2_query,
        query=emb1_query,
        k=k_neighbor,
        n_trees=n_trees,
        search_k=search_k,
    )
    nnba_idx, _ = nn_helper(
        emb1_ref,
        query=emb2_ref,
        k=k_neighbor,
        n_trees=n_trees,
        search_k=search_k,
    )

    # within-dataset neighbors for ScoreAnchors (k+1 includes self; drop self)
    nnaa_idx = None
    nnbb_idx = None
    if k_score > 0 and pca1_embedding is not None and pca2_embedding is not None:
        k_internal = min(k_neighbor + 1, pca1_embedding.shape[0])
        nnaa_idx, _ = nn_helper(
            pca1_embedding[:, :dims],
            k=k_internal,
            n_trees=n_trees,
            search_k=search_k,
        )
        k_internal_b = min(k_neighbor + 1, pca2_embedding.shape[0])
        nnbb_idx, _ = nn_helper(
            pca2_embedding[:, :dims],
            k=k_internal_b,
            n_trees=n_trees,
            search_k=search_k,
        )

    pairs: list[tuple[int, int]] = []
    for cell1_idx in range(n_ref):
        neighbors_ab = nnab_idx[cell1_idx, :k_anchor]
        for local2 in neighbors_ab:
            local2 = int(local2)
            ba_neighbors = nnba_idx[local2, :k_anchor]
            if cell1_idx in ba_neighbors:
                pairs.append((cell1_idx, int(cells2[local2])))

    if not pairs:
        return pd.DataFrame(columns=["cell1", "cell2", "score"])

    anchors = pd.DataFrame(sorted(set(pairs)), columns=["cell1", "cell2"])
    anchors["score"] = score_anchors_seurat(
        anchors,
        nnab_idx=nnab_idx,
        nnba_idx=nnba_idx,
        nnaa_idx=nnaa_idx,
        nnbb_idx=nnbb_idx,
        n_ref=n_ref,
        k_score=k_score,
    )
    return anchors.reset_index(drop=True)


def mirror_anchors(anchors: pd.DataFrame) -> pd.DataFrame:
    """Duplicate anchors with swapped cell order (Seurat ``FindIntegrationAnchors``)."""
    if anchors.empty:
        return anchors.copy()
    swapped = anchors.rename(columns={"cell1": "cell2", "cell2": "cell1"})
    if "dataset1" in anchors.columns and "dataset2" in anchors.columns:
        swapped = swapped.rename(columns={"dataset1": "dataset2", "dataset2": "dataset1"})
    if "batch1" in anchors.columns and "batch2" in anchors.columns:
        swapped = swapped.rename(columns={"batch1": "batch2", "batch2": "batch1"})
    return pd.concat([anchors, swapped], ignore_index=True)


def map_anchors_to_global(
    anchors: pd.DataFrame,
    offset1: int,
    offset2: int,
    n_local1: int,
    dataset1: int,
    dataset2: int,
) -> pd.DataFrame:
    """Map local stacked indices to global cell indices with dataset ids."""

    def _map_idx(c: int) -> int:
        c = int(c)
        if c < n_local1:
            return offset1 + c
        return offset2 + (c - n_local1)

    out = anchors.copy()
    out["cell1"] = out["cell1"].map(_map_idx)
    out["cell2"] = out["cell2"].map(_map_idx)
    out["dataset1"] = dataset1
    out["dataset2"] = dataset2
    return out
