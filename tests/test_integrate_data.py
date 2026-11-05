"""Tests for IntegrateData weight building and memory-friendly paths."""

import numpy as np
import pandas as pd
import pytest

from trackcell.tl.integration.integrate_data import (
    DEFAULT_WEIGHT_QUERY_CHUNK_SIZE,
    _build_weights_seurat,
    _build_weights_seurat_legacy,
    _distance_to_similarity,
    integrate_pair,
)

pytest.importorskip("annoy")


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
    dr = PCA(n_components=10, random_state=0).fit_transform(
        merged - merged.mean(0, keepdims=True)
    )
    dr_q = dr[n_ref:]
    aqc = np.unique(cell2)
    k_nn = min(20, len(aqc))
    from trackcell.tl.integration._annoy_nn import nn_helper

    nbr_idx, dists = nn_helper(dr_q[aqc], query=dr_q, k=k_nn, n_trees=10)
    dists = _distance_to_similarity(dists)

    w_legacy = _build_weights_seurat_legacy(
        len(cell1),
        n_query,
        cell2,
        scores,
        aqc,
        nbr_idx,
        dists,
        k_weight=20,
        sd_weight=1.0,
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
        np.testing.assert_allclose(
            w_legacy.toarray(), w_new.toarray(), rtol=0, atol=1e-12
        )


def test_distance_to_similarity_is_seurat_normalised():
    """``1 - d/d_k``: closest neighbour -> 1, k-th neighbour -> 0."""
    dists = np.array(
        [[0.0, 0.5, 2.0], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], dtype=np.float64
    )
    sim = _distance_to_similarity(dists)

    np.testing.assert_allclose(sim[0], [1.0, 0.75, 0.0])
    # All neighbours coincide (self included): d_k falls back to _MIN_DIST so the
    # row is 1 rather than 1 - 0/0 = NaN.
    np.testing.assert_allclose(sim[1], [1.0, 1.0, 1.0])
    # Equidistant but non-zero neighbours: d_k equals d, so every sim is 0. This
    # is Seurat's own behaviour; the column normaliser then zeroes the
    # correction for that cell. Unreachable in practice because a kNN row always
    # contains the query point itself at distance 0.
    np.testing.assert_allclose(sim[2], [0.0, 0.0, 0.0])
    assert np.isfinite(sim).all()
    # Monotonically decreasing in distance.
    assert (np.diff(sim[0]) <= 0).all()
    # Empty input round-trips.
    assert _distance_to_similarity(np.empty((0, 0))).shape == (0, 0)


def test_weights_decrease_with_distance():
    """Closer anchors must receive *more* weight than farther ones.

    Regression guard for the inverted-kernel bug: the factories expect
    Seurat-style similarities, so feeding them raw Euclidean distances must
    not silently produce a backwards kernel.
    """
    near, far = 0.1, 10.0
    nbr_idx = np.array([[0, 1], [1, 0]])
    dists = _distance_to_similarity(
        np.array([[near, far], [near, far]], dtype=np.float64)
    )
    weights = _build_weights_seurat(
        n_anchors=2,
        n_query=2,
        cell2=np.array([0, 1]),
        scores=np.ones(2),
        anchor_query_cells=np.array([0, 1]),
        nbr_idx=nbr_idx,
        dists=dists,
        k_weight=2,
        sd_weight=1.0,
    ).toarray()

    assert weights[0, 0] > weights[1, 0], (
        f"near anchor weight {weights[0, 0]:.4f} should exceed "
        f"far anchor weight {weights[1, 0]:.4f}"
    )
    # Column-normalised, so each query cell's weights sum to 1.
    np.testing.assert_allclose(weights.sum(axis=0), [1.0, 1.0])


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
