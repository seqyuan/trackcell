"""Tests for the matplotlib ROI selector (non-interactive backend)."""

import matplotlib

matplotlib.use("Agg")

import numpy as np
from anndata import AnnData
from matplotlib.figure import Figure

from trackcell.pl.select import RegionSelector


def _make_adata() -> AnnData:
    adata = AnnData(X=np.ones((6, 1), dtype=np.float32))
    adata.obs_names = [f"c{i}" for i in range(6)]
    adata.obsm["spatial"] = np.column_stack(
        [np.arange(6, dtype=float), np.zeros(6)]
    )
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
    name = selector.add_roi(
        np.array([[0.5, -1], [3.5, -1], [3.5, 1], [0.5, 1]])
    )
    assert name == "ROI_1"
    assert (adata.obs["ROI"] == "ROI_1").sum() == 3
    selector.disconnect()
