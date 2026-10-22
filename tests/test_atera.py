"""Tests for the Atera readers (synthetic fixtures, no external data)."""

import json
import zipfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import zarr

import trackcell as tcl

# zarr 3 renamed ``Group.create_dataset`` to ``Group.create_array``; provide the
# newer name on zarr 2 so the fixtures run against either major version.
if not hasattr(zarr.Group, "create_array"):
    zarr.Group.create_array = zarr.Group.create_dataset  # type: ignore[attr-defined]


# ---------------------------------------------------------------------------
# synthetic official bundle
# ---------------------------------------------------------------------------
PREFIXES = np.array([100, 101, 102, 103, 104, 105], dtype=np.uint32)
CENTROIDS = np.array(
    [[10.0, 10.0], [20.0, 20.0], [30.0, 30.0],
     [40.0, 40.0], [999.0, 999.0], [50.0, 50.0]],
    dtype=np.float32,
)
GENES = ["GENE1", "GENE2", "NegCtrl1", "Unassigned1"]
GENE_IDS = ["ENSG0001", "ENSG0002", "NC0001", "UN0001"]
FEATURE_TYPES = [
    "Gene Expression", "Gene Expression",
    "Negative Control Codeword", "Unassigned Codeword",
]


def _write_matrix_zarr(root: Path, with_de: bool = True) -> None:
    n_cells, n_genes = len(CENTROIDS), len(GENES)
    dense = np.array(
        [[2, 1, 0, 0],
         [0, 3, 0, 0],
         [1, 0, 1, 0],
         [0, 0, 0, 1],
         [4, 0, 0, 0],
         [0, 1, 0, 0]],
        dtype=np.int32,
    )
    X = sp.csr_matrix(dense)
    barcodes = pd.Index([f"cell{i}-1" for i in range(n_cells)], name="barcode")
    obs = pd.DataFrame({
        "barcode": barcodes,
        "cell_id": (1 << 32) + PREFIXES.astype(np.int64),
        "centroid_column": CENTROIDS[:, 0],
        "centroid_row": CENTROIDS[:, 1],
        "cell_area": np.full(n_cells, 20.0),
        "nucleus_area": np.full(n_cells, 9.0),
        "transcript_counts": np.array([3, 3, 2, 1, 4, 1], dtype=np.int32),
        "control_probe_counts": np.zeros(n_cells, dtype=np.int32),
        "control_codeword_counts": np.array([0, 0, 1, 0, 0, 0], dtype=np.int32),
        "total_counts": np.array([3, 3, 3, 1, 4, 1], dtype=np.int32),
    }, index=barcodes)
    var = pd.DataFrame({
        "feature_name": GENES,
        "feature_id": GENE_IDS,
        "feature_type": FEATURE_TYPES,
        "genome": "hg38",
    }, index=pd.Index(GENE_IDS, name="feature_id"))
    if with_de:
        var["de_leiden_c0_logfc"] = np.arange(n_genes, dtype=float)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm["X_pca"] = np.zeros((n_cells, 2), dtype=np.float32)
    adata.write_zarr(str(root / "cell_feature_matrix.zarr"))


def _write_cells_zarr(root: Path) -> None:
    group = zarr.open_group(str(root / "cells.zarr"), mode="w")
    cell_ids = np.column_stack([PREFIXES, np.ones_like(PREFIXES)])
    group.create_array("cell_id", data=cell_ids.astype(np.uint32))

    summary = np.zeros((len(CENTROIDS), 8), dtype=np.float64)
    summary[:, 0] = CENTROIDS[:, 0]
    summary[:, 1] = CENTROIDS[:, 1]
    summary[:, 2] = 20.0            # cell_area
    summary[:, 3] = CENTROIDS[:, 0]
    summary[:, 4] = CENTROIDS[:, 1]
    summary[:, 5] = 9.0             # nucleus_area
    node = group.create_array("cell_summary", data=summary)
    node.attrs["column_names"] = [
        "cell_centroid_x", "cell_centroid_y", "cell_area",
        "nucleus_centroid_x", "nucleus_centroid_y", "nucleus_area",
        "z_level", "nucleus_count",
    ]

    def write_polygons(cell_indices, offset, path):
        n = len(cell_indices)
        vertices = np.zeros((n, 50), dtype=np.float32)
        num_vertices = np.zeros(n, dtype=np.int32)
        for row, cell in enumerate(cell_indices):
            x, y = CENTROIDS[cell]
            ring = np.array(
                [[x - 2, y - 2], [x + 2, y - 2], [x + 2, y + 2], [x - 2, y + 2]],
                dtype=np.float32,
            )
            vertices[row, :8] = ring.ravel()
            num_vertices[row] = 4
        poly = group.require_group(path)
        poly.create_array("cell_index", data=np.asarray(cell_indices, dtype=np.uint32))
        poly.create_array("num_vertices", data=num_vertices)
        poly.create_array("vertices", data=vertices)
        poly.create_array("method", data=np.full(n, offset, dtype=np.uint32))

    write_polygons(np.arange(len(CENTROIDS)), 0, "polygon_sets/1")
    write_polygons(np.arange(len(CENTROIDS) - 1), 3, "polygon_sets/0")  # last cell missing


TX_ROWS = [
    # gene, qv, cell (None = unassigned), x, y, overlaps_nucleus
    ("GENE1", 30.0, 0, 10.0, 10.0, True),
    ("GENE1", 10.0, 0, 11.0, 10.0, False),
    ("GENE1", 20.0, 4, 999.0, 999.0, False),
    ("GENE2", 40.0, 1, 20.0, 20.0, False),
    ("GENE2", 25.0, None, 21.0, 20.0, False),
    ("NegCtrl1", 30.0, 2, 30.0, 30.0, True),
    ("Unassigned1", 20.0, 3, 40.0, 40.0, False),
]


def _write_transcripts_zarr(root: Path) -> None:
    group = zarr.open_group(str(root / "transcripts.zarr"), mode="w")
    group.attrs["gene_names"] = list(GENES)
    group.attrs["number_genes"] = len(GENES)
    group.attrs["number_rnas"] = len(TX_ROWS)
    group.attrs["gene_index_map"] = {gene: i for i, gene in enumerate(GENES)}
    group.attrs["codeword_gene_mapping"] = list(range(len(GENES)))
    group.attrs["spatial_units"] = "micron"
    group.attrs["coordinate_space"] = "refined-final_global_micron"
    group.attrs["major_version"] = 1
    group.attrs["minor_version"] = 1

    grid = group.create_group("grid")
    grid.attrs["grid_keys"] = ["0,0", "-1,0"]
    grid.attrs["grid_number_objects"] = [len(TX_ROWS), 1]
    grid.attrs["grid_size"] = 2000.0

    gene_to_index = {gene: i for i, gene in enumerate(GENES)}
    n = len(TX_ROWS)
    location = np.array([[r[3], r[4], 0.0] for r in TX_ROWS], dtype=np.float32)
    quality = np.array([[r[1]] for r in TX_ROWS], dtype=np.float16)
    cell_ids = np.array([
        [PREFIXES[r[2]], 1] if r[2] is not None else [4294967295, 4294967295]
        for r in TX_ROWS
    ], dtype=np.uint32)
    codeword = np.array(
        [[gene_to_index[r[0]]] for r in TX_ROWS], dtype=np.uint32
    )
    nucleus = np.array([[1 if r[5] else 0] for r in TX_ROWS], dtype=np.uint8)
    gene_offset = np.array(
        [[0, 3], [3, 5], [5, 6], [6, 7]], dtype=np.uint32
    )

    tile = grid.create_group("0,0")
    tile.create_array("location", data=location)
    tile.create_array("quality_score", data=quality)
    tile.create_array("cell_id", data=cell_ids)
    tile.create_array("codeword_identity", data=codeword)
    tile.create_array("overlaps_nucleus", data=nucleus)
    tile.create_array("gene_offset", data=gene_offset)

    # Atera writes a single placeholder row into background tiles.
    empty = grid.create_group("-1,0")
    empty.create_array(
        "location", data=np.array([[-0.4, 10.0, 0.0]], dtype=np.float32)
    )
    empty.create_array("quality_score", data=np.array([[40.0]], dtype=np.float16))
    empty.create_array(
        "cell_id",
        data=np.array([[4294967295, 4294967295]], dtype=np.uint32),
    )
    empty.create_array(
        "codeword_identity", data=np.array([[0]], dtype=np.uint32)
    )
    empty.create_array(
        "overlaps_nucleus", data=np.array([[0]], dtype=np.uint8)
    )
    empty.create_array("gene_offset", data=gene_offset)


def make_official_bundle(root: Path, with_de: bool = True) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    _write_matrix_zarr(root, with_de=with_de)
    _write_cells_zarr(root)
    _write_transcripts_zarr(root)
    (root / "experiment.spatial").write_text(json.dumps({
        "major_version": 1,
        "minor_version": 2,
        "run_name": "RUN1",
        "slide_name": "SLIDE1",
        "region_name": "ROI1",
        "pixel_size": 0.2125,
        "z_step_size": 3.0,
        "chemistry_version": "Atera v1",
        "analysis_sw_version": "atera-test",
    }))
    (root / "panel_config.json").write_text(json.dumps({
        "spec_version": "1",
        "metadata": {"customer": "SHOULD_NOT_BE_LOADED_BY_DEFAULT"},
        "payloads": [],
    }))
    return root


def zip_bundle(src: Path, dst: Path) -> Path:
    """Zip each zarr directory so the reader must use ZipStore."""
    dst.mkdir(parents=True, exist_ok=True)
    for member, name in (
        ("cell_feature_matrix.zarr", "cell_feature_matrix.zarr.zip"),
        ("cells.zarr", "cells.zarr.zip"),
        ("transcripts.zarr", "transcripts.zarr.zip"),
    ):
        with zipfile.ZipFile(dst / name, "w", zipfile.ZIP_DEFLATED) as handle:
            for path in (src / member).rglob("*"):
                if path.is_file():
                    handle.write(path, path.relative_to(src / member))
    for name in ("experiment.spatial", "panel_config.json"):
        (dst / name).write_bytes((src / name).read_bytes())
    return dst


# ---------------------------------------------------------------------------
# synthetic preview bundle (Xenium-v4 compatible)
# ---------------------------------------------------------------------------
def make_preview_bundle(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    n_cells, n_genes = len(CENTROIDS), len(GENES)
    dense = np.array(
        [[2, 1, 0, 0],
         [0, 3, 0, 0],
         [1, 0, 1, 0],
         [0, 0, 0, 1],
         [4, 0, 0, 0],
         [0, 1, 0, 0]],
        dtype=np.int32,
    )
    csc = sp.csc_matrix(dense.T)  # genes x cells
    barcodes = [f"cell{i}-1" for i in range(n_cells)]

    with h5py.File(root / "cell_feature_matrix.h5", "w") as handle:
        matrix = handle.create_group("matrix")
        matrix.create_dataset("shape", data=np.array([n_genes, n_cells]))
        matrix.create_dataset("data", data=csc.data)
        matrix.create_dataset("indices", data=csc.indices)
        matrix.create_dataset("indptr", data=csc.indptr)
        matrix.create_dataset(
            "barcodes", data=np.array([b.encode() for b in barcodes])
        )
        features = matrix.create_group("features")
        features.create_dataset("name", data=np.array([g.encode() for g in GENES]))
        features.create_dataset("id", data=np.array([g.encode() for g in GENE_IDS]))
        features.create_dataset(
            "feature_type", data=np.array([t.encode() for t in FEATURE_TYPES])
        )

    pd.DataFrame({
        "cell_id": barcodes,
        "x_centroid": CENTROIDS[:, 0],
        "y_centroid": CENTROIDS[:, 1],
        "cell_area": 20.0,
        "nucleus_area": 9.0,
        "transcript_counts": np.array([3, 3, 2, 1, 4, 1]),
        "total_counts": np.array([3, 3, 3, 1, 4, 1]),
        "segmentation_method": "test",
    }).to_parquet(root / "cells.parquet", index=False)

    rows = []
    for i, (x, y) in enumerate(CENTROIDS):
        for dx, dy in ((-2, -2), (2, -2), (2, 2), (-2, 2)):
            rows.append({"cell_id": barcodes[i],
                         "vertex_x": x + dx, "vertex_y": y + dy})
    pd.DataFrame(rows).to_parquet(root / "cell_boundaries.parquet", index=False)
    pd.DataFrame(rows[: len(rows) - 4]).to_parquet(
        root / "nucleus_boundaries.parquet", index=False
    )
    (root / "experiment.xenium").write_text(json.dumps({
        "run_name": "RUN1",
        "slide_id": "SLIDE1",
        "region_name": "PREVIEW1",
        "pixel_size": 0.2125,
        "z_step_size": 0.75,
    }))
    return root


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------
@pytest.fixture()
def official_bundle(tmp_path):
    return make_official_bundle(tmp_path / "official")


def test_read_official_schema(official_bundle):
    adata = tcl.io.read_atera(official_bundle)

    assert adata.shape == (6, 4)
    assert adata.var_names.tolist() == GENES
    assert adata.var.index.name == "gene_symbol"
    assert adata.var["gene_ids"].tolist() == GENE_IDS
    assert adata.var["feature_types"].tolist() == FEATURE_TYPES
    assert int(adata.var["is_gene_expression"].sum()) == 2
    assert "de_leiden_c0_logfc" not in adata.var.columns

    np.testing.assert_allclose(adata.obsm["spatial"], CENTROIDS)
    assert adata.obs["barcode"].tolist() == [
        f"cell{i}-1" for i in range(6)
    ]
    assert adata.obs["cell_prefix"].tolist() == PREFIXES.tolist()
    assert (adata.obs["dataset_suffix"] == 1).all()
    assert adata.obs["cell_id"].tolist() == ((1 << 32) + PREFIXES).tolist()

    assert adata.uns["atera"]["sample"] == "ROI1"
    assert adata.uns["atera"]["roi"] == "ROI1"
    assert adata.uns["atera"]["coord_unit"] == "um"
    assert adata.uns["atera"]["pixel_size"] == 0.2125
    assert adata.uns["atera"]["roi"] == adata.uns["experiment"]["region_name"]

    sample = adata.uns["atera"]["sample"]
    geometries = adata.uns["spatial"][sample]["geometries"]
    assert len(geometries) == 6
    assert set(geometries.index).issubset(set(adata.obs_names))
    assert len(adata.uns["spatial"][sample]["nucleus_geometries"]) == 5
    assert set(adata.uns["cell_boundaries"]) == {
        "vertex_x", "vertex_y", "n_vertices", "cell_idx"
    }
    assert adata.uns["cell_boundaries"]["cell_idx"].shape == (6,)
    # store_wkt is off by default
    assert "geometry" not in adata.obs.columns
    # panel is not loaded by default
    assert "gene_panel" not in adata.uns
    assert adata.uns["atera"]["panel_path"].endswith("panel_config.json")


def test_read_official_options(official_bundle):
    # keep_de
    adata = tcl.io.read_atera(official_bundle, keep_de=True)
    assert "de_leiden_c0_logfc" in adata.var.columns

    # store_wkt / add_wkt_geometry
    adata = tcl.io.read_atera(official_bundle, store_wkt=True)
    assert adata.obs["geometry"].notna().sum() == 6
    adata2 = tcl.io.read_atera(official_bundle)
    tcl.io.add_wkt_geometry(adata2)
    assert adata2.obs["geometry"].notna().sum() == 6

    # load_panel
    adata = tcl.io.read_atera(official_bundle, load_panel=True)
    assert adata.uns["gene_panel"]["metadata"]["customer"] == (
        "SHOULD_NOT_BE_LOADED_BY_DEFAULT"
    )

    # no polygons (nucleus defaults to off when polygons are off)
    adata = tcl.io.read_atera(official_bundle, build_polygons=False)
    assert "nucleus_geometries" not in adata.uns["spatial"]["ROI1"]
    with pytest.raises(ValueError, match="include_nucleus"):
        tcl.io.read_atera(
            official_bundle, build_polygons=False, include_nucleus=True
        )
    assert "geometries" not in adata.uns["spatial"]["ROI1"]
    assert "cell_boundaries" not in adata.uns

    # no nucleus
    adata = tcl.io.read_atera(official_bundle, include_nucleus=False)
    assert "nucleus_geometries" not in adata.uns["spatial"]["ROI1"]
    assert "nucleus_boundaries" not in adata.uns


def test_region_and_cell_subset(official_bundle):
    adata = tcl.io.read_atera(official_bundle, region=(0, 0, 45, 45))
    assert adata.n_obs == 4
    geometries = adata.uns["spatial"]["ROI1"]["geometries"]
    assert set(geometries.index) == set(adata.obs_names)
    assert adata.obsm["spatial"][:, 0].max() <= 45

    adata = tcl.io.read_atera(official_bundle, n_cells=3, random_state=0)
    assert adata.n_obs == 3
    assert len(adata.uns["spatial"]["ROI1"]["geometries"]) <= 3

    # shapely geometry is accepted too
    from shapely.geometry import box
    adata = tcl.io.read_atera(official_bundle, region=box(0, 0, 25, 25))
    assert adata.n_obs == 2

    with pytest.raises(ValueError, match="xmin"):
        tcl.io.read_atera(official_bundle, region=(45, 0, 0, 45))


def test_filter_genes_and_cell_qc(official_bundle):
    adata = tcl.io.read_atera(official_bundle)
    filtered = tcl.io.filter_atera_genes(adata)
    assert filtered.shape == (6, 2)
    assert filtered.var_names.tolist() == ["GENE1", "GENE2"]
    assert adata.shape == (6, 4)

    inplace = tcl.io.read_atera(official_bundle)
    tcl.io.filter_atera_genes(inplace, inplace=True)
    assert inplace.shape == (6, 2)

    qc = tcl.io.atera_cell_qc(adata)
    assert qc["n_genes"].tolist() == [2, 1, 2, 1, 1, 1]
    assert qc["n_genes_expression"].tolist() == [2, 1, 1, 0, 1, 1]
    assert qc["transcript_counts"].tolist() == [3, 3, 2, 1, 4, 1]
    assert qc.loc["cell2-1", "control_counts"] == 1


def test_atera_info(official_bundle):
    info = tcl.io.atera_info(official_bundle, verbose=False)
    assert info["format"] == "cellseg"
    assert info["n_cells"] == 6
    assert info["n_features"] == 4
    assert info["n_gene_expression"] == 2
    assert info["n_cell_polygons"] == 6
    assert info["n_nucleus_polygons"] == 5
    assert info["number_rnas"] == len(TX_ROWS)
    assert info["n_tiles"] == 2
    assert info["roi"] == "ROI1"


def test_transcript_streaming(official_bundle):
    df = tcl.io.read_atera_transcripts(
        official_bundle, genes=["GENE1"], min_qv=20
    )
    assert len(df) == 2  # Q10 row removed, out-of-region row kept
    assert df["gene"].tolist() == ["GENE1", "GENE1"]
    assert str(df["gene"].dtype) == "category"
    assert df["cell_index"].tolist() == [0, 4]

    df = tcl.io.read_atera_transcripts(
        official_bundle, genes=["GENE1"], min_qv=20, region=(0, 0, 50, 50)
    )
    assert len(df) == 1
    assert df.iloc[0]["x"] == 10.0

    # gene_offset fast path must not read other genes
    df = tcl.io.read_atera_transcripts(
        official_bundle, genes=["GENE2"], min_qv=None, include_unassigned=True
    )
    assert set(df["gene_index"]) == {1}
    assert len(df) == 2

    # unassigned molecules
    df = tcl.io.read_atera_transcripts(
        official_bundle, min_qv=20, include_unassigned=True, force=True
    )
    unassigned = df[df["cell_index"] < 0]
    assert len(unassigned) == 1
    assert unassigned.iloc[0]["gene"] == "GENE2"

    # background placeholder tile is skipped even with include_unassigned
    assert (df["x"] > 0).all()

    # full-scene guard
    with pytest.raises(ValueError, match="Refusing"):
        tcl.io.read_atera_transcripts(official_bundle)

    # bounded read
    df = tcl.io.read_atera_transcripts(
        official_bundle, min_qv=None, include_unassigned=True, max_rows=2
    )
    assert len(df) == 2

    # cell_index honours a subsetted AnnData
    adata = tcl.io.read_atera(official_bundle, region=(0, 0, 25, 25))
    assert adata.n_obs == 2
    df = tcl.io.read_atera_transcripts(
        official_bundle, genes=["GENE1"], min_qv=20, adata=adata
    )
    assert df["cell_index"].tolist() == [0, -1]


def test_iter_transcripts_closes_store(official_bundle):
    iterator = tcl.io.iter_atera_transcripts(official_bundle, min_qv=None)
    first = next(iterator)
    assert len(first) > 0
    iterator.close()  # must not raise / leak


def test_transcript_qc(official_bundle):
    adata = tcl.io.read_atera(official_bundle)
    summary = tcl.io.atera_transcript_qc(
        official_bundle, adata=adata, min_qv=20, bins=10, verbose=False
    )
    assert summary["read"] == len(TX_ROWS)
    assert summary["q20"] == 6
    assert summary["assigned_q20"] == 5
    assert summary["q20_fraction"] == pytest.approx(6 / 7)
    assert summary["assigned_fraction_q20"] == pytest.approx(5 / 6)
    assert summary["cells_with_transcripts"] == 5

    assert adata.obs["tx_q20_counts"].tolist() == [1, 1, 1, 1, 1, 0]
    assert adata.obs["tx_nuclear_counts"].tolist() == [1, 0, 1, 0, 0, 0]
    assert np.isnan(adata.obs.loc["cell5-1", "nuclear_fraction"])
    assert adata.obs.loc["cell0-1", "nuclear_fraction"] == pytest.approx(1.0)

    atera_uns = adata.uns["atera"]
    assert atera_uns["density"]["H"].sum() == 6  # Q>=20 molecules
    assert set(atera_uns["gene_qc"]["gene"]) == set(GENES)
    assert atera_uns["qscore"]["hq"].sum() == 6

    # plotting functions are re-callable from adata only
    import matplotlib
    matplotlib.use("Agg")
    axes = tcl.pl.atera_qc_summary(adata, show=False)
    assert np.asarray(axes).size == 4
    ax = tcl.pl.atera_transcript_density(adata, show=False)
    assert ax is not None
    axes = tcl.pl.atera_nuclear_fraction(
        adata, genes=["GENE1", "GENE2"], show=False
    )
    assert axes.shape == (3,)


def test_concat_rois(official_bundle):
    a = tcl.io.read_atera(official_bundle, sample="ROI1")
    b = tcl.io.read_atera(official_bundle, sample="ROI2")

    with pytest.raises(ValueError, match="offsets"):
        tcl.io.concat_atera_rois([a, b])
    with pytest.raises(ValueError, match="missing ROIs"):
        tcl.io.concat_atera_rois([a, b], offsets={"ROI1": (0, 0)})

    original = np.asarray(a.obsm["spatial"]).copy()
    merged = tcl.io.concat_atera_rois(
        {"ROI1": a, "ROI2": b},
        offsets={"ROI1": (0, 0), "ROI2": (1000, 500)},
        index_unique="-",
    )
    assert merged.obs_names.is_unique
    assert merged.n_obs == 12
    assert set(merged.obs["roi"]) == {"ROI1", "ROI2"}
    np.testing.assert_allclose(merged.obsm["spatial"][:6], original)
    np.testing.assert_allclose(
        merged.obsm["spatial"][6:], original + np.array([1000, 500])
    )
    # geometries were translated too
    gdf = merged.uns["spatial"]["ROI2"]["geometries"]
    assert gdf.geometry.bounds["maxx"].max() > 1000
    assert set(gdf.index).issubset(set(merged.obs_names))


def test_list_and_iter_rois(tmp_path, official_bundle):
    run = tmp_path / "run"
    run.mkdir()
    (run / "ROI_A").symlink_to(official_bundle)
    second = make_official_bundle(run / "ROI_B")
    (second / "experiment.spatial").write_text(json.dumps({
        "region_name": "ROI_B",
        "slide_name": "SLIDE1",
        "pixel_size": 0.2125,
    }))

    table = tcl.io.list_atera_rois(run)
    assert set(table["roi"]) == {"ROI1", "ROI_B"}
    assert table["n_cells"].tolist() == [6, 6]
    assert table["has_transcripts"].all()

    names = []
    for name, adata in tcl.io.iter_atera_rois(run, rois=["ROI_B"],
                                              build_polygons=False):
        names.append(name)
        assert adata.n_obs == 6
    assert names == ["ROI_B"]


def test_zipped_bundle(tmp_path, official_bundle):
    zipped = zip_bundle(official_bundle, tmp_path / "zipped")
    adata = tcl.io.read_atera(zipped)
    assert adata.shape == (6, 4)
    assert len(adata.uns["spatial"]["ROI1"]["geometries"]) == 6
    df = tcl.io.read_atera_transcripts(zipped, genes=["GENE1"], min_qv=20)
    assert len(df) == 2


def test_preview_bundle(tmp_path):
    preview = make_preview_bundle(tmp_path / "preview")
    adata = tcl.io.read_atera(preview)

    assert adata.shape == (6, 4)
    assert adata.uns["atera"]["format"] == "preview"
    assert adata.uns["atera"]["sample"] == "PREVIEW1"
    assert adata.obs_names.name == "barcode"
    assert adata.obs["barcode"].tolist() == [f"cell{i}-1" for i in range(6)]
    assert adata.var_names.tolist() == GENES
    assert adata.var["gene_ids"].tolist() == GENE_IDS
    assert adata.var["feature_types"].tolist() == FEATURE_TYPES
    assert int(adata.var["is_gene_expression"].sum()) == 2
    assert adata.uns["atera"]["coord_unit"] == "um"
    assert len(adata.uns["spatial"]["PREVIEW1"]["geometries"]) == 6
    assert "geometry" not in adata.obs.columns
    assert "nucleus_boundaries" in adata.uns

    # options behave like the official reader
    adata = tcl.io.read_atera(preview, build_polygons=False)
    assert "geometries" not in adata.uns["spatial"]["PREVIEW1"]
    assert "cell_boundaries" not in adata.uns

    adata = tcl.io.read_atera(preview, include_nucleus=False)
    assert "nucleus_boundaries" not in adata.uns

    adata = tcl.io.read_atera(preview, store_wkt=True)
    assert adata.obs["geometry"].notna().sum() == 6

    adata = tcl.io.read_atera(preview, region=(0, 0, 45, 45))
    assert adata.n_obs == 4
    assert set(adata.uns["spatial"]["PREVIEW1"]["geometries"].index) == set(
        adata.obs_names
    )

    # explicit preview entry point
    adata = tcl.io.read_atera_preview(preview)
    assert adata.n_obs == 6


def test_preview_schema_matches_official(tmp_path, official_bundle):
    official = tcl.io.read_atera(official_bundle)
    preview = tcl.io.read_atera(make_preview_bundle(tmp_path / "preview"))

    for key in ("gene_ids", "feature_types", "is_gene_expression"):
        assert key in official.var.columns
        assert key in preview.var.columns
    assert official.var.index.name == preview.var.index.name == "gene_symbol"
    assert "barcode" in official.obs.columns
    assert "barcode" in preview.obs.columns
    assert official.uns["atera"]["coord_unit"] == preview.uns["atera"]["coord_unit"]
    assert np.asarray(official.obsm["spatial"]).dtype == np.float32
    assert np.asarray(preview.obsm["spatial"]).dtype == np.float32


def test_gene_name_index_false(official_bundle, tmp_path):
    official = tcl.io.read_atera(official_bundle, gene_name_index=False)
    assert official.var_names.tolist() == GENE_IDS
    assert official.var.index.name == "feature_id"
    assert official.var["gene_ids"].tolist() == GENE_IDS

    preview = tcl.io.read_atera(
        make_preview_bundle(tmp_path / "preview"), gene_name_index=False
    )
    assert preview.var_names.tolist() == GENE_IDS
    assert preview.var.index.name == "feature_id"


def test_unzip_support_and_bad_paths(tmp_path, official_bundle):
    # directory-backed zarr is supported (useful for HPC staging)
    adata = tcl.io.read_atera(official_bundle)
    assert adata.n_obs == 6

    with pytest.raises(FileNotFoundError):
        tcl.io.read_atera(tmp_path / "does_not_exist")

    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileNotFoundError, match="No Atera matrix"):
        tcl.io.read_atera(empty)

    run = tmp_path / "run"
    run.mkdir()
    (run / "ROI_A").symlink_to(official_bundle)
    with pytest.raises(FileNotFoundError, match="run directory"):
        tcl.io.read_atera(run)
