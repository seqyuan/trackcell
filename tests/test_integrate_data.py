"""Tests for IntegrateData weight building and memory-friendly paths."""

import numpy as np
import pandas as pd
import pytest

from trackcell.tl.integration.integrate_data import (
    DEFAULT_WEIGHT_QUERY_CHUNK_SIZE,
    _build_weights_seurat,
    _build_weights_seurat_legacy,
    integrate_pair,
)


def _synthetic_pair(n_ref=200, n_query=500, n_anchors=800, n_feat=50, seed=0):
    rng = np.random.default_rng(seed)
    data_ref = rng.standard_normal((n_ref, n_feat))
    data_q = rng.standard_normal((n_query, n_feat))
    cell1 = rng.integers(0, n_ref, n_anchors)
    cell2 = rng.integers(0, n_query, n_anchors)
    scores = rng.random(n_anchors)
    anchors = pd.DataFrame({"cell1": cell1, "cell2": cell2, "score": scores})
    return data_ref, data_q, anchors


def test_build_weights_chunked_matches_legacy():
    data_ref, data_q, anchors = _synthetic_pair()
    cell1 = anchors["cell1"].to_numpy()
    cell2 = anchors["cell2"].to_numpy()
    scores = anchors["score"].to_numpy()
    n_ref = data_ref.shape[0]
    n_query = data_q.shape[0]

    from sklearn.decomposition import PCA

    merged = np.vstack([data_ref, data_q])
    dr = PCA(n_components=10, random_state=0).fit_transform(merged - merged.mean(0, keepdims=True))
    dr_q = dr[n_ref:]
    aqc = np.unique(cell2)
    k_nn = min(20, len(aqc))
    from trackcell.tl.integration._annoy_nn import nn_helper

    nbr_idx, dists = nn_helper(dr_q[aqc], query=dr_q, k=k_nn, n_trees=10)

    w_legacy = _build_weights_seurat_legacy(
        len(cell1), n_query, cell2, scores, aqc, nbr_idx, dists, k_weight=20, sd_weight=1.0
    )
    for chunk in (1, 7, DEFAULT_WEIGHT_QUERY_CHUNK_SIZE):
        w_new = _build_weights_seurat(
            len(cell1),
            n_query,
            cell2,
            scores,
            aqc,
            nbr_idx,
            dists,
            k_weight=20,
            sd_weight=1.0,
            query_chunk_size=chunk,
        )
        np.testing.assert_allclose(w_legacy.toarray(), w_new.toarray(), rtol=0, atol=1e-12)


def test_integrate_pair_dtypes_and_trees_weight():
    data_ref, data_q, anchors = _synthetic_pair()
    out64 = integrate_pair(
        data_ref,
        data_q,
        anchors,
        n_trees=50,
        n_trees_weight=50,
        integration_dtype="float64",
        weight_query_chunk_size=64,
    )
    out32 = integrate_pair(
        data_ref,
        data_q,
        anchors,
        n_trees=50,
        n_trees_weight=10,
        integration_dtype="float32",
        weight_query_chunk_size=64,
    )
    assert out64.shape == out32.shape == data_q.shape
    assert np.isfinite(out64).all()
    assert np.isfinite(out32).all()
    corr = np.corrcoef(out64.ravel(), out32.ravel())[0, 1]
    assert corr > 0.99
