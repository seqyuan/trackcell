"""Tests for the matplotlib ROI selector (non-interactive backend)."""

import matplotlib

matplotlib.use("Agg")

import warnings
import numpy as np
from anndata import AnnData
from matplotlib.figure import Figure
from shapely.geometry import Polygon

from trackcell.pl.select import RegionSelector


def _make_adata() -> AnnData:
    adata = AnnData(X=np.ones((6, 1), dtype=np.float32))
    adata.obs_names = [f"c{i}" for i in range(6)]
    adata.obsm["spatial"] = np.column_stack([np.arange(6, dtype=float), np.zeros(6)])
    return adata


def test_save_with_key_added_writes_that_column():
    adata = _make_adata()
    fig = Figure()
    ax = fig.add_subplot(111)

    selector = RegionSelector(
        adata,
        ax=ax,
        mode="squarebin",
        basis="spatial",
        valid_ids=adata.obs_names,
        key_added="ROI",
        inplace=False,
    )
    selector.add_roi(np.array([[0.5, -1], [3.5, -1], [3.5, 1], [0.5, 1]]))
    selector.save(key_added="NEW_ROI")

    assert "NEW_ROI" in adata.obs.columns
    assert "ROI" not in adata.obs.columns
    assert (adata.obs["NEW_ROI"] == "ROI_1").sum() == 3
    selector.disconnect()


def test_add_roi_without_name_autonumbers():
    adata = _make_adata()
    fig = Figure()
    ax = fig.add_subplot(111)
    selector = RegionSelector(
        adata,
        ax=ax,
        mode="squarebin",
        basis="spatial",
        valid_ids=adata.obs_names,
        inplace=True,
    )
    name = selector.add_roi(np.array([[0.5, -1], [3.5, -1], [3.5, 1], [0.5, 1]]))
    assert name == "ROI_1"
    assert (adata.obs["ROI"] == "ROI_1").sum() == 3
    selector.disconnect()


def _make_selector_adata():
    adata = AnnData(X=np.ones((4, 1), dtype="float32"))
    adata.obs_names = [f"c{i}" for i in range(4)]
    adata.var_names = ["g0"]
    adata.obsm["spatial"] = np.array([[0.5, 0.5], [1.5, 0.5], [5.0, 5.0], [6.0, 6.0]])
    return adata


def _new_selector(adata, **kw):
    fig = Figure()
    ax = fig.add_subplot(111)
    return RegionSelector(
        adata,
        ax=ax,
        mode="squarebin",
        basis="spatial",
        valid_ids=adata.obs_names,
        inplace=True,
        **kw,
    )


def test_save_to_adata_preserves_existing_roi_labels():
    """Creating a new ROI must not wipe labels restored from a previous h5ad."""
    adata = _make_selector_adata()
    adata.obs["roi"] = ["old_a", "old_a", "old_b", None]

    sel = _new_selector(adata, key_added="roi")
    sel._save_to_adata()

    # Cells outside the new (empty) ROI set keep their previous labels.
    assert adata.obs["roi"].tolist()[:3] == ["old_a", "old_a", "old_b"]
    assert adata.obs["roi"].isna().sum() == 1


def test_to_adata_preserves_existing_roi_labels():
    adata = _make_selector_adata()
    adata.obs["roi"] = ["old_a", "old_a", "old_b", None]
    fig = Figure()
    ax = fig.add_subplot(111)
    sel = RegionSelector(
        adata,
        ax=ax,
        mode="squarebin",
        basis="spatial",
        valid_ids=adata.obs_names,
        key_added="roi",
        inplace=False,
    )
    sel.add_roi(np.array([[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]]))
    copy = sel.to_adata()

    assert copy.obs["roi"].tolist()[:3] == ["ROI_1", "ROI_1", "old_b"]
    assert copy.obs["roi"].isna().sum() == 1
    # ``to_adata`` returns a copy; the original column is untouched.
    assert adata.obs["roi"].tolist()[:2] == ["old_a", "old_a"]


def test_extract_ids_repairs_self_intersecting_polygon():
    """A self-intersecting lasso must not produce arbitrary membership."""
    adata = _make_selector_adata()
    sel = _new_selector(adata, key_added="roi")

    # A lasso that traces the square and then retraces a diagonal — a normal
    # mouse artefact. The ring is invalid; contains() on the unrepaired polygon
    # returned only c1, silently dropping c0 from a ROI the user clearly drew.
    lasso = np.array(
        [
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
            [0.0, 0.0],
            [2.0, 2.0],
            [0.0, 0.0],
        ]
    )
    assert not Polygon(lasso).is_valid

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ids = sel._extract_ids(lasso)

    assert set(ids).issubset(set(adata.obs_names))
    assert set(ids) == {"c0", "c1"}
