"""Tests for ``read_hd_bin`` (SpaceRanger HD bin-level output).

All bundles are synthetic: see ``tests/conftest.py`` for the 10x h5 / parquet /
image writers that reproduce the SpaceRanger on-disk layout.
"""

import numpy as np
import pytest

from trackcell.io.read_data import read_hd_bin
from tests.conftest import (
    DEFAULT_BARCODES,
    DEFAULT_GENES,
    dense_counts,
    write_10x_mtx_dir,
    write_tissue_positions,
)


def test_read_hd_bin_happy_path(hd_bin_bundle):
    adata = read_hd_bin(hd_bin_bundle, binsize=16)

    assert adata.shape == (len(DEFAULT_BARCODES), len(DEFAULT_GENES))
    assert list(adata.var_names) == list(DEFAULT_GENES)
    assert list(adata.obs_names) == list(DEFAULT_BARCODES)
    assert int(adata.X.sum()) == int(dense_counts(len(DEFAULT_GENES), len(DEFAULT_BARCODES)).sum())

    # Positions were merged into obs and spatial coordinates come from
    # the full-resolution pixel columns, ordered (col, row).
    for column in ("in_tissue", "array_row", "array_col",
                   "pxl_row_in_fullres", "pxl_col_in_fullres"):
        assert column in adata.obs.columns
    expected = np.column_stack(
        [
            np.array([60.0, 180.0, 300.0]),
            np.array([50.0, 150.0, 250.0]),
        ]
    )
    np.testing.assert_allclose(adata.obsm["spatial"], expected)

    spatial = adata.uns["spatial"]["Cse1"]
    assert spatial["binsize"] == 16
    assert spatial["images"]["hires"].shape == (12, 16, 3)
    assert spatial["images"]["lowres"].shape == (6, 8, 3)
    assert spatial["scalefactors"]["tissue_hires_scalef"] == pytest.approx(
        0.1502397209405899
    )


def test_read_hd_bin_infers_sample_from_spaceranger_layout(hd_bin_bundle):
    # datapath is <tmp>/Cse1/outs/binned_outputs/square_016um
    assert hd_bin_bundle.name == "square_016um"

    adata = read_hd_bin(hd_bin_bundle)

    assert list(adata.uns["spatial"]) == ["Cse1"]


def test_read_hd_bin_csv_positions_fallback(hd_bin_bundle):
    (hd_bin_bundle / "spatial" / "tissue_positions.parquet").unlink()
    write_tissue_positions(
        hd_bin_bundle / "spatial" / "tissue_positions.csv", writer="csv"
    )
    assert not (hd_bin_bundle / "spatial" / "tissue_positions.parquet").exists()

    adata = read_hd_bin(hd_bin_bundle, sample="Cse1")

    assert adata.shape == (3, 3)
    np.testing.assert_allclose(adata.obsm["spatial"][:, 0], [60.0, 180.0, 300.0])
    assert list(adata.obs["array_row"]) == [0, 2, 4]


def test_read_hd_bin_duplicate_barcodes_raise(hd_bin_bundle):
    write_tissue_positions(
        hd_bin_bundle / "spatial" / "tissue_positions.parquet",
        barcodes=DEFAULT_BARCODES,
        extra_rows=(DEFAULT_BARCODES[0],),
    )

    with pytest.raises(ValueError, match="duplicate barcodes"):
        read_hd_bin(hd_bin_bundle, sample="Cse1")


def test_read_hd_bin_missing_barcode_raise(hd_bin_bundle):
    # Only the first two barcodes have positions.
    write_tissue_positions(
        hd_bin_bundle / "spatial" / "tissue_positions.parquet",
        barcodes=DEFAULT_BARCODES[:2],
    )

    with pytest.raises(ValueError, match="missing from tissue positions"):
        read_hd_bin(hd_bin_bundle, sample="Cse1")


def test_read_hd_bin_missing_coordinate_columns_raise(hd_bin_bundle):
    write_tissue_positions(
        hd_bin_bundle / "spatial" / "tissue_positions.parquet",
        include_coords=False,
        include_array=False,
    )

    with pytest.raises(ValueError, match="Could not find spatial coordinate columns"):
        read_hd_bin(hd_bin_bundle, sample="Cse1")


def test_read_hd_bin_without_images_warns_but_returns_adata(hd_bin_bundle):
    (hd_bin_bundle / "spatial" / "tissue_hires_image.png").unlink()
    (hd_bin_bundle / "spatial" / "tissue_lowres_image.png").unlink()

    with pytest.warns(RuntimeWarning, match="Could not load tissue images"):
        adata = read_hd_bin(hd_bin_bundle, sample="Cse1")

    assert adata.shape == (3, 3)
    assert adata.obsm["spatial"].shape == (3, 2)
    assert "images" not in adata.uns["spatial"]["Cse1"]
    assert adata.uns["spatial"]["Cse1"]["scalefactors"]["tissue_hires_scalef"] > 0


def test_read_hd_bin_without_scalefactors_warns(hd_bin_bundle):
    (hd_bin_bundle / "spatial" / "scalefactors_json.json").unlink()

    with pytest.warns(RuntimeWarning, match="Could not load scalefactors"):
        adata = read_hd_bin(hd_bin_bundle, sample="Cse1")

    assert adata.uns["spatial"]["Cse1"]["scalefactors"] == {}
    # Images still load, so the overlay keys survive.
    assert "images" in adata.uns["spatial"]["Cse1"]


def test_read_hd_bin_directory_matrix_form(hd_bin_bundle):
    # SpaceRanger also emits the plain mtx directory; no h5 at all here.
    (hd_bin_bundle / "filtered_feature_bc_matrix.h5").unlink()
    write_10x_mtx_dir(hd_bin_bundle / "filtered_feature_bc_matrix")

    adata = read_hd_bin(hd_bin_bundle, sample="Cse1")

    assert adata.shape == (3, 3)
    assert list(adata.var_names) == list(DEFAULT_GENES)
    assert adata.var["feature_types"].unique().tolist() == ["Gene Expression"]
    expected_total = int(
        dense_counts(len(DEFAULT_GENES), len(DEFAULT_BARCODES)).sum()
    )
    assert int(adata.X.sum()) == expected_total


def test_read_hd_bin_corrupt_h5_falls_back_to_directory(hd_bin_bundle):
    # Truncated/not-an-HDF5 file: sc.read_10x_h5 raises, reader must fall back.
    (hd_bin_bundle / "filtered_feature_bc_matrix.h5").write_bytes(b"not an h5 file")
    write_10x_mtx_dir(hd_bin_bundle / "filtered_feature_bc_matrix")

    adata = read_hd_bin(hd_bin_bundle, sample="Cse1")

    assert adata.shape == (3, 3)
    expected_total = int(
        dense_counts(len(DEFAULT_GENES), len(DEFAULT_BARCODES)).sum()
    )
    assert int(adata.X.sum()) == expected_total


def test_read_hd_bin_missing_matrix_raise(tmp_path):
    datapath = tmp_path / "Cse1" / "outs"
    datapath.mkdir(parents=True)

    with pytest.raises(FileNotFoundError, match="filtered_feature_bc_matrix"):
        read_hd_bin(datapath, sample="Cse1")
