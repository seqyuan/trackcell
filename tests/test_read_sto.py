"""Tests for the STOmics (Stereo-seq) GEF readers.

``read_sto_cellbin`` is the primary target; ``read_sto_bin`` and the image
auto-discovery helper are covered as well because they share the GEF layout.
Every bundle is synthesised in ``tests/conftest.py`` with the exact
``cellBin/{cell,cellBorder,gene,cellExp}`` and ``geneExp/bin1/{expression,gene}``
group/dataset/dtype structure the readers consume.
"""

import geopandas as gpd
import numpy as np
import pytest
from shapely import wkt

from trackcell.io.read_sto import read_sto_bin, read_sto_cellbin
from tests.conftest import (
    STO_GENES,
    sentinel_border,
    write_cellbin_gef,
    write_gem,
    write_image,
    write_tif_image,
)

# ---------------------------------------------------------------------------
# read_sto_cellbin
# ---------------------------------------------------------------------------


def test_read_sto_cellbin_happy_path(sto_cellbin_bundle):
    adata = read_sto_cellbin(sto_cellbin_bundle)

    assert adata.shape == (3, len(STO_GENES))
    assert list(adata.var_names) == list(STO_GENES)
    assert list(adata.obs_names) == ["cell_1", "cell_2", "cell_3"]
    assert adata.obs.index.name == "cell_id"

    # cellExp offsets: cell1 -> GeneA:5, GeneB:2 ; cell2 -> GeneC:3 ; cell3 -> GeneB:4, GeneA:1
    np.testing.assert_array_equal(
        adata.X.toarray(),
        np.array([[5.0, 2.0, 0.0], [0.0, 0.0, 3.0], [1.0, 4.0, 0.0]]),
    )
    assert adata.obs["gene_count"].tolist() == [2, 1, 2]
    assert adata.obs["exp_count"].tolist() == [7, 3, 5]
    assert adata.obs["dnb_count"].tolist() == [10, 20, 30]
    assert adata.obs["area"].tolist() == [64, 64, 64]

    # Centroids come straight from the GEF cell table.
    np.testing.assert_allclose(
        adata.obsm["spatial"], [[100.0, 200.0], [300.0, 400.0], [500.0, 600.0]]
    )

    # Sample name was inferred from the GEF filename.
    spatial = adata.uns["spatial"]["B000001"]
    assert spatial["metadata"]["sample_id"] == "B000001"
    assert spatial["metadata"]["resolution_nm"] == 500
    assert spatial["metadata"]["pixel_size_um"] == pytest.approx(0.5)
    assert spatial["metadata"]["coordinate_range"] == {"x": [100, 500], "y": [200, 600]}
    assert spatial["scalefactors"]["spot_diameter_fullres"] == pytest.approx(
        2.0 * np.sqrt(64 / np.pi)
    )

    # Border offsets are relative to the centroid, so 8x8 squares of area 64.
    geometries = spatial["geometries"]
    assert isinstance(geometries, gpd.GeoDataFrame)
    assert len(geometries) == adata.n_obs == 3
    assert [poly.area for poly in geometries.geometry] == [64.0, 64.0, 64.0]
    assert [poly.centroid.x for poly in geometries.geometry] == [100.0, 300.0, 500.0]

    assert all(g.startswith("POLYGON") for g in adata.obs["geometry"])
    assert wkt.loads(adata.obs["geometry"].iloc[0]).bounds == (
        96.0,
        196.0,
        104.0,
        204.0,
    )
    assert "images" not in spatial  # no 03.register/ sibling in this bundle


def test_read_sto_cellbin_integer_index_var_names(sto_cellbin_bundle):
    adata = read_sto_cellbin(sto_cellbin_bundle, gene_name_index=False)

    assert list(adata.var_names) == ["0", "1", "2"]
    assert adata.shape == (3, 3)
    # Expression is unchanged, only the var index differs.
    assert adata.X[0, 0] == 5.0


def test_read_sto_cellbin_accepts_directory(sto_cellbin_bundle):
    # The GEF lives in <tmp>/C57_7/041.cellcut/ ; pass the directory instead.
    assert sto_cellbin_bundle.parent.name == "041.cellcut"

    from_directory = read_sto_cellbin(sto_cellbin_bundle.parent, sample="C57_7")
    from_file = read_sto_cellbin(sto_cellbin_bundle, sample="C57_7")

    assert from_directory.shape == from_file.shape == (3, 3)
    assert list(from_directory.uns["spatial"]) == ["C57_7"]
    np.testing.assert_array_equal(from_directory.X.toarray(), from_file.X.toarray())


def test_read_sto_cellbin_without_cellbin_group_raises(tmp_path):
    gef = write_cellbin_gef(
        tmp_path / "B000009.cellbin.gef", include_cellbin_group=False
    )

    with pytest.raises(ValueError, match="does not contain a 'cellBin' group"):
        read_sto_cellbin(gef)


def test_read_sto_cellbin_all_sentinel_borders_use_fallback_squares(tmp_path):
    n_cells = 3
    gef = write_cellbin_gef(
        tmp_path / "B000002.cellbin.gef",
        borders=[sentinel_border() for _ in range(n_cells)],
    )

    adata = read_sto_cellbin(gef, sample="sentinel")

    geometries = adata.uns["spatial"]["sentinel"]["geometries"]
    assert len(geometries) == 3
    # Every cell falls back to a 10x10 square centred on the centroid.
    assert [poly.area for poly in geometries.geometry] == [100.0, 100.0, 100.0]
    assert [round(poly.centroid.x) for poly in geometries.geometry] == [100, 300, 500]
    assert all(not poly.is_empty for poly in geometries.geometry)
    assert all("POLYGON" in g for g in adata.obs["geometry"])


def test_read_sto_cellbin_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="GEF file not found"):
        read_sto_cellbin(tmp_path / "nope.cellbin.gef")


def test_read_sto_cellbin_directory_without_gef_raises(tmp_path):
    empty = tmp_path / "041.cellcut"
    empty.mkdir()

    with pytest.raises(FileNotFoundError, match="No cellbin GEF file found"):
        read_sto_cellbin(empty)


def test_read_sto_cellbin_autodiscovers_registered_image(sto_cellbin_bundle):
    # STOmics keeps the ssDNA image in a sibling 03.register/ directory.
    write_tif_image(
        sto_cellbin_bundle.parents[1] / "03.register" / "ssDNA_1_regist.tif"
    )

    adata = read_sto_cellbin(sto_cellbin_bundle.parent)

    images = adata.uns["spatial"]["B000001"]["images"]
    # Grayscale TIFF is converted to RGB so sc.pl.spatial renders it correctly.
    assert images["hires"].shape == (10, 12, 3)
    assert images["hires"].dtype == np.uint8
    np.testing.assert_array_equal(images["hires"], images["lowres"])
    assert adata.uns["spatial"]["B000001"]["scalefactors"]["tissue_hires_scalef"] == 1.0


def test_read_sto_cellbin_explicit_image_path(sto_cellbin_bundle, tmp_path):
    image = write_tif_image(tmp_path / "custom_image.tif", height=4, width=5)

    adata = read_sto_cellbin(sto_cellbin_bundle, image_path=image)

    assert adata.uns["spatial"]["B000001"]["images"]["hires"].shape == (4, 5, 3)


def test_read_sto_cellbin_unsupported_image_format_warns(sto_cellbin_bundle, tmp_path):
    image = write_image(tmp_path / "not_supported.png")

    with pytest.warns(UserWarning, match="Unsupported image format"):
        adata = read_sto_cellbin(sto_cellbin_bundle, image_path=image)

    assert "images" not in adata.uns["spatial"]["B000001"]
    assert adata.n_obs == 3


# ---------------------------------------------------------------------------
# read_sto_bin (squarebin GEF / GEM)
# ---------------------------------------------------------------------------


def test_read_sto_bin_gef_happy_path(sto_tissue_bundle):
    adata = read_sto_bin(sto_tissue_bundle, bin_size=2)

    # Entries (0,0) (1,1) (2,0) (3,3) (4,1) collapse into 4 bins of size 2.
    assert adata.shape == (4, len(STO_GENES))
    assert adata.obs.index.name == "bin_id"
    assert list(adata.var_names) == list(STO_GENES)
    # Bin ids are the packed (bin_x << 32) | bin_y keys.
    assert set(adata.obs_names) == {"0", "4294967296", "4294967297", "8589934592"}
    np.testing.assert_allclose(
        adata.obsm["spatial"],
        [[1.0, 1.0], [3.0, 1.0], [3.0, 3.0], [5.0, 1.0]],
    )

    spatial = adata.uns["spatial"]["B000001"]
    assert spatial["metadata"]["bin_size"] == 2
    assert spatial["metadata"]["bin_level"] == 1
    assert spatial["metadata"]["bin_type"] == "square_bins"
    assert spatial["scalefactors"]["spot_diameter_fullres"] == 2.0  # == bin_size
    assert adata.X.sum() == 15.0


def test_read_sto_bin_without_gene_name_index(sto_tissue_bundle):
    adata = read_sto_bin(sto_tissue_bundle, bin_size=2, gene_name_index=False)

    assert list(adata.var_names) == ["0", "1", "2"]
    assert adata.n_vars == 3


def test_read_sto_bin_accepts_directory(sto_tissue_bundle, tmp_path):
    # A cellbin GEF next to the tissue GEF must be skipped by the tissue
    # auto-discovery glob.
    write_cellbin_gef(sto_tissue_bundle.parent / "B000001.cellbin.gef")

    adata = read_sto_bin(sto_tissue_bundle.parent, bin_size=2)

    assert adata.shape == (4, 3)
    assert list(adata.uns["spatial"]) == ["B000001"]


def test_read_sto_bin_directory_without_tissue_gef_raises(tmp_path):
    cut_dir = tmp_path / "041.cellcut"
    write_cellbin_gef(cut_dir / "B000001.cellbin.gef")

    with pytest.raises(FileNotFoundError, match="No tissue/squarebin GEF file found"):
        read_sto_bin(cut_dir)


def test_read_sto_bin_gem_text_fallback(tmp_path):
    gem = write_gem(
        tmp_path / "B000003.gem",
        rows=(
            ("GeneA", 0, 0, 1),
            ("GeneB", 1, 1, 2),
            ("GeneA", 2, 2, 3),
            ("GeneB", 3, 3, 4),
        ),
    )

    with pytest.warns(UserWarning, match="Reading large GEM files can be slow"):
        adata = read_sto_bin(gem, bin_size=2)

    assert adata.shape == (2, 2)
    assert list(adata.var_names) == ["GeneA", "GeneB"]
    assert list(adata.obs_names) == ["0_0", "1_1"]
    np.testing.assert_allclose(adata.obsm["spatial"], [[1.0, 1.0], [3.0, 3.0]])
    np.testing.assert_array_equal(adata.X.toarray(), [[1.0, 2.0], [3.0, 4.0]])
    assert adata.uns["spatial"]["B000003"]["metadata"]["bin_size"] == 2


def test_read_sto_bin_unsupported_format_raises(tmp_path):
    bogus = tmp_path / "B000004.txt"
    bogus.write_text("not a spatial matrix\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Unsupported file format"):
        read_sto_bin(bogus)


def test_read_sto_cellbin_records_gef_offsets(tmp_path):
    # GEF offsets position the crop within the chip; they are recorded in
    # metadata (not added to coordinates, which are already in image frame).
    gef = write_cellbin_gef(
        tmp_path / "B000010.cellbin.gef",
        cells=(("A", 100, 200), ("B", 300, 400), ("C", 500, 600)),
        resolution=500,
        offset_x=1000,
        offset_y=2000,
    )

    adata = read_sto_cellbin(gef, sample="off")
    meta = adata.uns["spatial"]["off"]["metadata"]

    assert meta["offset_x"] == 1000
    assert meta["offset_y"] == 2000
    # resolution 500 nm/unit -> 0.5 um/unit; offset um = offset_units * 0.5
    assert meta["offset_um"] == pytest.approx([500.0, 1000.0])
    # Coordinates stay in the GEF/image frame (offsets are NOT applied).
    np.testing.assert_allclose(
        adata.obsm["spatial"], [[100.0, 200.0], [300.0, 400.0], [500.0, 600.0]]
    )
