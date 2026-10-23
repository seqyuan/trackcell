"""Tests for Seurat-style anchor scoring (combined index space)."""

import numpy as np
import pandas as pd

from trackcell.tl.integration.anchors import (
    _normalize_anchor_scores,
    _raw_anchor_scores,
    score_anchors_seurat,
)


def test_raw_anchor_scores_does_not_mix_index_spaces():
    """Batch1/batch2 indices must not collide: nnab=1 is not nnba=1."""
    n_ref = 3
    anchors = pd.DataFrame({"cell1": [0], "cell2": [n_ref + 0]})

    nnab = np.array([[1, 1]])       # batch2 local neighbor of cell1
    nnba = np.array([[1, 1]])       # batch1 local neighbor of cell2
    nnaa = np.array([[0, 0]])       # batch1 neighbor of cell1
    nnbb = np.array([[0, 0]])       # batch2 neighbor of cell2

    raw = _raw_anchor_scores(
        anchors,
        nnab_idx=nnab,
        nnba_idx=nnba,
        nnaa_idx=nnaa,
        nnbb_idx=nnbb,
        n_ref=n_ref,
        k_score=2,
    )
    # cell1 neighbours = {0, n_ref + 1}; cell2 neighbours = {1, n_ref + 0}
    assert raw.tolist() == [0.0]


def test_raw_anchor_scores_counts_shared_neighbors_per_dataset():
    n_ref = 3
    anchors = pd.DataFrame({"cell1": [0], "cell2": [n_ref + 0]})

    nnab = np.array([[0, 1]])
    nnba = np.array([[1, 2]])
    nnaa = np.array([[1, 2]])
    nnbb = np.array([[0, 1]])

    raw = _raw_anchor_scores(
        anchors,
        nnab_idx=nnab,
        nnba_idx=nnba,
        nnaa_idx=nnaa,
        nnbb_idx=nnbb,
        n_ref=n_ref,
        k_score=2,
    )
    # {1, 2, 3, 4} ∩ {1, 2, 3, 4}
    assert raw.tolist() == [4.0]


def test_score_anchors_seurat_normalizes_raw_counts():
    n_ref = 3
    anchors = pd.DataFrame(
        {"cell1": [0, 1], "cell2": [n_ref + 0, n_ref + 1]}
    )
    # Anchor 0 shares everything; anchor 1 shares nothing.
    nnab = np.array([[0, 1], [2, 2]])
    nnba = np.array([[1, 2], [1, 1]])
    nnaa = np.array([[1, 2], [0, 0]])
    nnbb = np.array([[0, 1], [0, 0]])

    raw = _raw_anchor_scores(
        anchors,
        nnab_idx=nnab,
        nnba_idx=nnba,
        nnaa_idx=nnaa,
        nnbb_idx=nnbb,
        n_ref=n_ref,
        k_score=2,
    )
    assert raw.tolist() == [4.0, 0.0]

    scores = score_anchors_seurat(
        anchors,
        nnab_idx=nnab,
        nnba_idx=nnba,
        nnaa_idx=nnaa,
        nnbb_idx=nnbb,
        n_ref=n_ref,
        k_score=2,
    )
    assert np.allclose(scores, _normalize_anchor_scores(raw))
    assert scores[0] > scores[1]
