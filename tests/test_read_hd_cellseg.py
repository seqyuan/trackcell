"""Tests for ``read_hd_cellseg`` (SpaceRanger HD segmented output).

Bundles come from ``tests/conftest.py``: a 10x cell-feature h5 plus a
``graphclust_annotated_cell_segmentations.geojson`` of square polygons.
"""

import geopandas as gpd
import numpy as np
import pytest
from shapely import wkt

from trackcell.io.read_data import read_hd_cellseg
from tests.conftest import (
    DEFAULT_GENES,
    SEG_CELLID_NAMES,
    SEG_CELL_IDS,
    write_10x_h5,
    write_cell_segmentations,
)


def test_read_hd_cellseg_happy_path(hd_cellseg_bundle):
    adata = read_hd_cellseg(hd_cellseg_bundle)

    assert adata.shape == (3, len(DEFAULT_GENES))
    assert list(adata.obs_names) == list(SEG_CELLID_NAMES)

    # Centroid coordinates of the three 8x8 squares.
    np.testing.assert_allclose(
        adata.obsm["spatial"], [[4.0, 4.0], [14.0, 14.0], [24.0, 24.0]]
    )

    # WKT strings in obs, GeoDataFrame in uns.
    assert all(
        isinstance(g, str) and g.startswith("POLYGON") for g in adata.obs["geometry"]
    )
    parsed = [wkt.loads(g) for g in adata.obs["geometry"]]
    assert [poly.area for poly in parsed] == [64.0, 64.0, 64.0]

    spatial = adata.uns["spatial"]["Cse1"]
    geometries = spatial["geometries"]
    assert isinstance(geometries, gpd.GeoDataFrame)
    assert len(geometries) == adata.n_obs == 3
    # GeoJSON's default WGS84 must not leak onto pixel coordinates.
    assert geometries.crs is None
    assert spatial["images"]["hires"].shape == (12, 16, 3)
    # spot_diameter_fullres was absent from the scalefactors JSON, so the
    # reader estimates it as fiducial_diameter_fullres / 40.
    assert spatial["scalefactors"]["spot_diameter_fullres"] == pytest.approx(
        145.8916900543866 / 40.0
    )


def test_read_hd_cellseg_missing_images_warns(hd_cellseg_bundle):
    (hd_cellseg_bundle / "spatial" / "tissue_hires_image.png").unlink()
    (hd_cellseg_bundle / "spatial" / "tissue_lowres_image.png").unlink()

    with pytest.warns(RuntimeWarning, match="Could not load tissue images"):
        adata = read_hd_cellseg(hd_cellseg_bundle, sample="Cse1")

    assert "images" not in adata.uns["spatial"]["Cse1"]
    assert len(adata.uns["spatial"]["Cse1"]["geometries"]) == 3


def test_read_hd_cellseg_extra_segmentations_warn_and_are_dropped(hd_cellseg_bundle):
    # Segmentation is a superset of the expression matrix (cell 4 has no counts).
    write_cell_segmentations(
        hd_cellseg_bundle / "graphclust_annotated_cell_segmentations.geojson",
        cell_ids=(1, 2, 3, 4),
    )

    with pytest.warns(
        UserWarning, match="segmented cells have no expression matrix row"
    ):
        adata = read_hd_cellseg(hd_cellseg_bundle, sample="Cse1")

    assert adata.n_obs == 3
    assert len(adata.uns["spatial"]["Cse1"]["geometries"]) == 3
    assert list(adata.obs_names) == list(SEG_CELLID_NAMES)


def test_read_hd_cellseg_extra_expression_cells_warn_and_are_dropped(hd_cellseg_bundle):
    # Expression is a superset of the segmentation table.
    extra_cell = "cellid_000000004-1"
    write_10x_h5(
        hd_cellseg_bundle / "filtered_feature_cell_matrix.h5",
        cells=SEG_CELLID_NAMES + (extra_cell,),
    )

    with pytest.warns(UserWarning, match="have no segmentation geometry"):
        adata = read_hd_cellseg(hd_cellseg_bundle, sample="Cse1")

    assert adata.n_obs == 3
    assert extra_cell not in set(adata.obs_names)
    assert list(adata.obs_names) == list(SEG_CELLID_NAMES)


def test_read_hd_cellseg_no_intersection_raise(hd_cellseg_bundle):
    write_cell_segmentations(
        hd_cellseg_bundle / "graphclust_annotated_cell_segmentations.geojson",
        cell_ids=(91, 92, 93),
    )

    with pytest.raises(ValueError, match="No expression barcodes match"):
        read_hd_cellseg(hd_cellseg_bundle, sample="Cse1")


def test_read_hd_cellseg_duplicate_cell_ids_raise(hd_cellseg_bundle):
    write_cell_segmentations(
        hd_cellseg_bundle / "graphclust_annotated_cell_segmentations.geojson",
        cell_ids=(1, 1, 2),
    )

    with pytest.raises(ValueError, match="duplicate cell IDs"):
        read_hd_cellseg(hd_cellseg_bundle, sample="Cse1")


def test_read_hd_cellseg_missing_cell_id_column_raise(hd_cellseg_bundle):
    write_cell_segmentations(
        hd_cellseg_bundle / "graphclust_annotated_cell_segmentations.geojson",
        cell_id_column="segment_id",
    )

    with pytest.raises(ValueError, match="must contain a 'cell_id' column"):
        read_hd_cellseg(hd_cellseg_bundle, sample="Cse1")


def test_read_hd_cellseg_alternative_filename_warns(hd_cellseg_bundle):
    default = hd_cellseg_bundle / "graphclust_annotated_cell_segmentations.geojson"
    default.rename(hd_cellseg_bundle / "cell_segmentations.geojson")

    with pytest.warns(UserWarning, match="Using alternative file"):
        adata = read_hd_cellseg(hd_cellseg_bundle, sample="Cse1")

    assert adata.n_obs == 3
    assert len(adata.uns["spatial"]["Cse1"]["geometries"]) == 3


def test_read_hd_cellseg_missing_segmentations_file_raise(hd_cellseg_bundle):
    (hd_cellseg_bundle / "graphclust_annotated_cell_segmentations.geojson").unlink()

    with pytest.raises(FileNotFoundError, match="Cell segmentations file not found"):
        read_hd_cellseg(hd_cellseg_bundle, sample="Cse1")


def test_read_hd_cellseg_classification_metadata(hd_cellseg_bundle):
    write_cell_segmentations(
        hd_cellseg_bundle / "graphclust_annotated_cell_segmentations.geojson",
        cell_ids=SEG_CELL_IDS,
        classification=[
            '{"name": "B cell", "color": [255, 0, 0]}',
            '{"name": "T cell", "color": [0, 128, 255]}',
            '{"name": "B cell", "color": [255, 0, 0]}',
        ],
    )

    adata = read_hd_cellseg(hd_cellseg_bundle, sample="Cse1")

    assert adata.obs["classification"].tolist() == ["B cell", "T cell", "B cell"]
    assert adata.uns["classification_colors"] == {
        "B cell": "#ff0000",
        "T cell": "#0080ff",
    }


def test_read_hd_cellseg_without_vendor_colors_publishes_no_palette(
    hd_cellseg_bundle,
):
    """A classification without colours must not publish an empty palette.

    ``pl.spatial_cell`` treated ``classification_colors={}`` as a real palette,
    warned that every category was missing, and painted them all grey.
    """
    write_cell_segmentations(
        hd_cellseg_bundle / "graphclust_annotated_cell_segmentations.geojson",
        cell_ids=SEG_CELL_IDS,
        classification=[
            '{"name": "B cell"}',
            '{"name": "T cell"}',
            '{"name": "B cell"}',
        ],
    )

    adata = read_hd_cellseg(hd_cellseg_bundle, sample="Cse1")

    assert adata.obs["classification"].tolist() == ["B cell", "T cell", "B cell"]
    assert "classification_colors" not in adata.uns
