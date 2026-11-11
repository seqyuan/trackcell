"""Tests for spatial analysis helpers."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from trackcell.tl.spatial import hd_labeldist


def _make_labeldist_adata(include_hires_scale: bool = True) -> AnnData:
    adata = AnnData(
        X=np.ones((3, 1), dtype=np.float32),
        obs=pd.DataFrame(
            {"celltype": ["target", "other", "other"]}, index=["c0", "c1", "c2"]
        ),
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


def test_hd_labeldist_physical_coords_ignore_metadata_pixel_size():
    """Xenium/Atera centres are µm; ``metadata['pixel_size']`` is the image
    pixel size and must not be applied (README: do not multiply by it)."""
    adata = _make_labeldist_adata()
    scalefactors = adata.uns["spatial"]["sample"]["scalefactors"]
    scalefactors.pop("microns_per_pixel")
    scalefactors.pop("tissue_hires_scalef", None)
    adata.uns["spatial"]["sample"]["metadata"] = {
        "pixel_size": 0.2125,
        "coord_unit": "um",
    }

    with pytest.warns(UserWarning, match="ignored"):
        result = hd_labeldist(
            adata,
            groupby="celltype",
            label="target",
            inplace=False,
            coordinate_system="fullres",
        )

    # 10 and 20 coordinate units are 10 and 20 µm, not 2.1 and 4.3 µm.
    np.testing.assert_allclose(result["target_dist"].values, [0.0, 10.0, 20.0])
    assert result.attrs["params"]["microns_per_pixel"] == 1.0
    assert result.attrs["params"]["unit_family"] == "um"
    assert result.attrs["params"]["coordinate_system"] == "data"


def test_hd_labeldist_physical_scale_from_stomics_metadata():
    """STOmics coordinates are resolution units (500 nm by default)."""
    adata = _make_labeldist_adata()
    scalefactors = adata.uns["spatial"]["sample"]["scalefactors"]
    scalefactors.pop("microns_per_pixel")
    scalefactors.pop("tissue_hires_scalef", None)
    adata.uns["spatial"]["sample"]["metadata"] = {"pixel_size_um": 0.5}

    result = hd_labeldist(
        adata,
        groupby="celltype",
        label="target",
        inplace=False,
        coordinate_system="auto",
    )

    np.testing.assert_allclose(result["target_dist"].values, [0.0, 5.0, 10.0])
    assert result.attrs["params"]["microns_per_pixel"] == 0.5
    assert result.attrs["params"]["unit_family"] == "um"


def test_hd_labeldist_hires_rejected_for_physical_coordinates():
    adata = _make_labeldist_adata()
    adata.uns["spatial"]["sample"]["metadata"] = {"coord_unit": "um"}

    with pytest.raises(ValueError, match="only meaningful"):
        hd_labeldist(
            adata,
            groupby="celltype",
            label="target",
            inplace=False,
            coordinate_system="hires",
        )


def test_hd_labeldist_unit_hint_beats_scalefactor_pixel_size():
    """``coord_unit='um'`` wins over a stray ``microns_per_pixel`` entry."""
    adata = _make_labeldist_adata()  # mpp=0.5 + hires scale
    adata.uns["spatial"]["sample"]["metadata"] = {
        "coord_unit": "um",
        "pixel_size": 0.2125,
    }

    result = hd_labeldist(
        adata,
        groupby="celltype",
        label="target",
        inplace=False,
        coordinate_system="auto",
    )

    np.testing.assert_allclose(result["target_dist"].values, [0.0, 10.0, 20.0])
    assert result.attrs["params"]["microns_per_pixel"] == 1.0
