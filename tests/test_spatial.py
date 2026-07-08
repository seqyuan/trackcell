"""Tests for spatial analysis helpers."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from trackcell.tl.spatial import hd_labeldist


def _make_labeldist_adata(include_hires_scale: bool = True) -> AnnData:
    adata = AnnData(
        X=np.ones((3, 1), dtype=np.float32),
        obs=pd.DataFrame({"celltype": ["target", "other", "other"]}, index=["c0", "c1", "c2"]),
        var=pd.DataFrame(index=["g0"]),
    )
    adata.obsm["spatial"] = np.array([[0.0, 0.0], [10.0, 0.0], [0.0, 20.0]])
    scalefactors = {"microns_per_pixel": 0.5}
    if include_hires_scale:
        scalefactors["tissue_hires_scalef"] = 0.25
    adata.uns["spatial"] = {"sample": {"scalefactors": scalefactors}}
    return adata


def test_hd_labeldist_explicit_coordinate_systems():
    adata = _make_labeldist_adata()

    fullres = hd_labeldist(
        adata,
        groupby="celltype",
        label="target",
        inplace=False,
        coordinate_system="fullres",
    )
    hires = hd_labeldist(
        adata,
        groupby="celltype",
        label="target",
        inplace=False,
        coordinate_system="hires",
    )

    np.testing.assert_allclose(fullres["target_dist"].values, [0.0, 5.0, 10.0])
    np.testing.assert_allclose(hires["target_dist"].values, [0.0, 20.0, 40.0])
    assert fullres.attrs["params"]["coordinate_system"] == "fullres"
    assert hires.attrs["params"]["coordinate_system"] == "hires"


def test_hd_labeldist_auto_warns_when_ambiguous():
    adata = _make_labeldist_adata()

    with pytest.warns(UserWarning, match="Could not confidently infer"):
        result = hd_labeldist(
            adata,
            groupby="celltype",
            label="target",
            inplace=False,
            coordinate_system="auto",
        )

    assert result.attrs["params"]["coordinate_system"] == "hires"


def test_hd_labeldist_hires_requires_scale():
    adata = _make_labeldist_adata(include_hires_scale=False)

    with pytest.raises(ValueError, match="requires `hires_scale`"):
        hd_labeldist(
            adata,
            groupby="celltype",
            label="target",
            inplace=False,
            coordinate_system="hires",
        )
