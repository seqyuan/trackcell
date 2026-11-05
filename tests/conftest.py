"""Shared synthetic vendor-bundle fixtures for the IO readers.

The builders in this module produce *minimal but structurally faithful*
imitations of real vendor output layouts (10x SpaceRanger HD, STOmics GEF).
Nothing here reaches the network and every bundle is written into
``tmp_path`` so tests stay hermetic and fast.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse as sp_sparse

# ---------------------------------------------------------------------------
# Canned identifiers
# ---------------------------------------------------------------------------

DEFAULT_GENES = ("GeneA", "GeneB", "GeneC")
DEFAULT_BARCODES = ("AAAC-1", "AAAG-1", "AAAT-1")

# SpaceRanger cell-segmentation ids follow ``cellid_{9 digits}-1``.
SEG_CELL_IDS = (1, 2, 3)
SEG_CELLID_NAMES = tuple(f"cellid_{cid:09d}-1" for cid in SEG_CELL_IDS)


def dense_counts(n_genes: int, n_cells: int) -> np.ndarray:
    """Deterministic (genes x cells) non-negative count matrix."""
    counts = np.zeros((n_genes, n_cells), dtype=np.int32)
    for cell in range(n_cells):
        counts[cell % n_genes, cell] += cell + 1
        counts[0, cell] += 1
    return counts


def _as_counts(genes, cells, counts) -> np.ndarray:
    if counts is None:
        return dense_counts(len(genes), len(cells))
    counts = np.asarray(counts, dtype=np.int32)
    assert counts.shape == (len(genes), len(cells)), counts.shape
    return counts


# ---------------------------------------------------------------------------
# 10x feature-barcode matrix writers
# ---------------------------------------------------------------------------


def write_10x_h5(
    path: Path, genes=DEFAULT_GENES, cells=DEFAULT_BARCODES, counts=None
) -> Path:
    """Write a 10x Genomics v3 ``filtered_feature_bc_matrix.h5``.

    Layout (as consumed by ``scanpy.read_10x_h5``)::

        matrix/{data,indices,indptr,shape,barcodes}
        matrix/features/{name,id,feature_type,genome}

    ``data``/``indices``/``indptr`` describe the **CSC** (genes x cells)
    matrix, i.e. ``indptr`` has ``n_cells + 1`` entries.
    """
    genes = [str(g) for g in genes]
    cells = [str(c) for c in cells]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    csc = sp_sparse.csc_matrix(_as_counts(genes, cells, counts))

    with h5py.File(path, "w") as handle:
        matrix = handle.create_group("matrix")
        matrix.create_dataset(
            "shape", data=np.array([len(genes), len(cells)], dtype=np.int32)
        )
        matrix.create_dataset("data", data=csc.data.astype(np.int32))
        matrix.create_dataset("indices", data=csc.indices.astype(np.int32))
        matrix.create_dataset("indptr", data=csc.indptr.astype(np.int32))
        matrix.create_dataset("barcodes", data=np.array([c.encode() for c in cells]))

        features = matrix.create_group("features")
        features.create_dataset("name", data=np.array([g.encode() for g in genes]))
        features.create_dataset(
            "id", data=np.array([f"ENSG{i:011d}".encode() for i in range(len(genes))])
        )
        features.create_dataset(
            "feature_type",
            data=np.array([b"Gene Expression"] * len(genes)),
        )
        features.create_dataset("genome", data=np.array([b"GRCh38"] * len(genes)))
    return path


def write_10x_mtx_dir(
    directory: Path, genes=DEFAULT_GENES, cells=DEFAULT_BARCODES, counts=None
) -> Path:
    """Write the directory form of a 10x v3 matrix (gzipped mtx + tsvs).

    ``matrix.mtx.gz`` is stored genes x cells, matching what Cell Ranger
    writes; ``scanpy.read_10x_mtx`` transposes it on read.
    """
    genes = [str(g) for g in genes]
    cells = [str(c) for c in cells]
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    counts = _as_counts(genes, cells, counts)
    coo = sp_sparse.coo_matrix(counts)

    with gzip.open(directory / "matrix.mtx.gz", "wt") as handle:
        handle.write("%%MatrixMarket matrix coordinate integer general\n")
        handle.write("%\n")
        handle.write(f"{len(genes)} {len(cells)} {coo.nnz}\n")
        for row, col, value in zip(coo.row, coo.col, coo.data):
            handle.write(f"{row + 1} {col + 1} {int(value)}\n")

    with gzip.open(directory / "features.tsv.gz", "wt") as handle:
        for i, gene in enumerate(genes):
            handle.write(f"ENSG{i:011d}\t{gene}\tGene Expression\n")

    with gzip.open(directory / "barcodes.tsv.gz", "wt") as handle:
        for cell in cells:
            handle.write(f"{cell}\n")

    return directory


# ---------------------------------------------------------------------------
# SpaceRanger spatial side-car files
# ---------------------------------------------------------------------------


def write_tissue_positions(
    path: Path,
    barcodes=DEFAULT_BARCODES,
    *,
    include_coords: bool = True,
    include_array: bool = True,
    extra_rows=(),
    writer=None,
) -> Path:
    """Write ``spatial/tissue_positions.parquet`` (or csv when ``writer`` set).

    ``include_coords`` controls the ``pxl_*_in_fullres`` columns and
    ``include_array`` the ``array_row``/``array_col`` columns, so tests can
    exercise the reader's coordinate-column fallback chain.
    """
    barcodes = list(barcodes) + list(extra_rows)
    n = len(barcodes)
    data = {
        "barcode": barcodes,
        "in_tissue": np.ones(n, dtype=np.int8),
    }
    if include_array:
        data["array_row"] = np.arange(n, dtype=np.int64) * 2
        data["array_col"] = np.arange(n, dtype=np.int64) * 3
    if include_coords:
        data["pxl_row_in_fullres"] = np.arange(n, dtype=np.float64) * 100.0 + 50.0
        data["pxl_col_in_fullres"] = np.arange(n, dtype=np.float64) * 120.0 + 60.0
    df = pd.DataFrame(data)

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if writer == "csv":
        df.to_csv(path, index=False)
    else:
        df.to_parquet(path, index=False)
    return path


def write_scalefactors(path: Path, include_spot_diameter: bool = True) -> Path:
    scalefactors = {
        "tissue_hires_scalef": 0.1502397209405899,
        "tissue_lowres_scalef": 0.0450727992124723,
        "fiducial_diameter_fullres": 145.8916900543866,
    }
    if include_spot_diameter:
        scalefactors["spot_diameter_fullres"] = 50.400000000000006
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(scalefactors), encoding="utf-8")
    return path


def write_image(
    path: Path, height: int = 12, width: int = 16, value: int = 128
) -> Path:
    """Write a tiny RGB PNG that imageio can read back."""
    import imageio.v3 as iio

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    img = np.full((height, width, 3), value, dtype=np.uint8)
    iio.imwrite(path, img)
    return path


def write_tif_image(
    path: Path, height: int = 10, width: int = 12, value: int = 42
) -> Path:
    """Write a single-channel (grayscale) TIFF, like STOmics ssDNA captures."""
    import imageio.v3 as iio

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    iio.imwrite(path, np.full((height, width), value, dtype=np.uint8))
    return path


# ---------------------------------------------------------------------------
# SpaceRanger bundles
# ---------------------------------------------------------------------------


@pytest.fixture
def hd_bin_bundle(tmp_path: Path) -> Path:
    """``<tmp>/Cse1/outs/binned_outputs/square_016um`` with h5 + images."""
    square = tmp_path / "Cse1" / "outs" / "binned_outputs" / "square_016um"
    square.mkdir(parents=True)
    write_10x_h5(square / "filtered_feature_bc_matrix.h5")
    write_tissue_positions(square / "spatial" / "tissue_positions.parquet")
    write_scalefactors(square / "spatial" / "scalefactors_json.json")
    write_image(square / "spatial" / "tissue_hires_image.png")
    write_image(square / "spatial" / "tissue_lowres_image.png", height=6, width=8)
    return square


@pytest.fixture
def hd_cellseg_bundle(tmp_path: Path) -> Path:
    """``<tmp>/Cse1/outs/segmented_outputs`` with geojson + cell matrix."""
    seg_dir = tmp_path / "Cse1" / "outs" / "segmented_outputs"
    seg_dir.mkdir(parents=True)
    write_10x_h5(seg_dir / "filtered_feature_cell_matrix.h5", cells=SEG_CELLID_NAMES)
    write_cell_segmentations(
        seg_dir / "graphclust_annotated_cell_segmentations.geojson"
    )
    # Deliberately omit `spot_diameter_fullres` so the reader has to estimate it.
    write_scalefactors(
        seg_dir / "spatial" / "scalefactors_json.json", include_spot_diameter=False
    )
    write_image(seg_dir / "spatial" / "tissue_hires_image.png")
    write_image(seg_dir / "spatial" / "tissue_lowres_image.png", height=6, width=8)
    return seg_dir


def write_cell_segmentations(
    path: Path,
    cell_ids=SEG_CELL_IDS,
    *,
    classification=None,
    cell_id_column: str = "cell_id",
):
    """Write a ``graphclust_annotated_cell_segmentations.geojson``."""
    import geopandas as gpd
    from shapely.geometry import Polygon

    records = {
        cell_id_column: list(cell_ids),
        "geometry": [],
    }
    for i, _ in enumerate(cell_ids):
        x0, y0 = 10.0 * i, 10.0 * i
        records["geometry"].append(
            Polygon([(x0, y0), (x0 + 8.0, y0), (x0 + 8.0, y0 + 8.0), (x0, y0 + 8.0)])
        )
    if classification is not None:
        records["classification"] = list(classification)

    gdf = gpd.GeoDataFrame(pd.DataFrame(records))
    gdf = gdf.set_geometry("geometry")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_file(path, driver="GeoJSON")
    return path


# ---------------------------------------------------------------------------
# STOmics GEF builders
# ---------------------------------------------------------------------------

BORDER_SENTINEL = np.int16(32767)
N_BORDER_SLOTS = 32

CELL_DTYPE = np.dtype(
    [
        ("id", "<u4"),
        ("x", "<i4"),
        ("y", "<i4"),
        ("area", "<u4"),
        ("dnbCount", "<u4"),
        ("geneCount", "<u4"),
        ("expCount", "<u4"),
        ("offset", "<u8"),
        ("clusterID", "<u2"),
    ]
)
GENE_DTYPE = np.dtype([("geneName", "S64")])
CELL_EXP_DTYPE = np.dtype([("geneID", "<u4"), ("count", "<u4")])
BIN_EXPR_DTYPE = np.dtype([("x", "<i4"), ("y", "<i4"), ("count", "<u4")])
BIN_GENE_DTYPE = np.dtype([("gene", "S64"), ("offset", "<u8"), ("count", "<u8")])

STO_GENES = ("GeneA", "GeneB", "GeneC")


def square_border(half: int = 4) -> np.ndarray:
    """A 4-vertex relative-offset square padded to 32 (dx, dy) slots."""
    offsets = np.full((N_BORDER_SLOTS, 2), BORDER_SENTINEL, dtype=np.int16)
    offsets[0] = (-half, -half)
    offsets[1] = (half, -half)
    offsets[2] = (half, half)
    offsets[3] = (-half, half)
    return offsets


def sentinel_border() -> np.ndarray:
    return np.full((N_BORDER_SLOTS, 2), BORDER_SENTINEL, dtype=np.int16)


def write_cellbin_gef(
    path: Path,
    *,
    cells=(("A", 100, 200), ("B", 300, 400), ("C", 500, 600)),
    genes=STO_GENES,
    expressions=(((0, 5), (1, 2)), ((2, 3),), ((1, 4), (0, 1))),
    borders=None,
    resolution: int = 500,
    include_cellbin_group: bool = True,
) -> Path:
    """Write a STOmics ``*.cellbin.gef`` (HDF5) file.

    Topology mirrors the real GEF: ``cellBin/{cell,cellBorder,gene,cellExp}``
    where ``cell['offset']`` indexes into the flattened ``cellExp`` table.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    genes = [str(g) for g in genes]
    n_cells = len(cells)

    if borders is None:
        borders = [square_border() for _ in range(n_cells)]

    cell = np.zeros(n_cells, dtype=CELL_DTYPE)
    cell_exp = []
    offset = 0
    for i, ((cell_id, x, y), gene_counts) in enumerate(zip(cells, expressions)):
        cell["id"][i] = cell_id if isinstance(cell_id, int) else i + 1
        cell["x"][i] = x
        cell["y"][i] = y
        cell["area"][i] = 64
        cell["dnbCount"][i] = 10 * (i + 1)
        cell["geneCount"][i] = len(gene_counts)
        cell["expCount"][i] = sum(count for _, count in gene_counts)
        cell["offset"][i] = offset
        offset += len(gene_counts)
        cell_exp.extend(gene_counts)

    cell_exp_array = np.zeros(len(cell_exp), dtype=CELL_EXP_DTYPE)
    for i, (gene_id, count) in enumerate(cell_exp):
        cell_exp_array["geneID"][i] = gene_id
        cell_exp_array["count"][i] = count

    gene_array = np.zeros(len(genes), dtype=GENE_DTYPE)
    for i, gene in enumerate(genes):
        gene_array["geneName"][i] = gene.encode()

    border_array = np.stack(borders).astype(np.int16)

    with h5py.File(path, "w") as handle:
        handle.attrs["resolution"] = resolution
        handle.attrs["offsetX"] = 0
        handle.attrs["offsetY"] = 0
        if include_cellbin_group:
            cellbin = handle.create_group("cellBin")
            cell_ds = cellbin.create_dataset("cell", data=cell)
            cell_ds.attrs["minX"] = int(cell["x"].min())
            cell_ds.attrs["maxX"] = int(cell["x"].max())
            cell_ds.attrs["minY"] = int(cell["y"].min())
            cell_ds.attrs["maxY"] = int(cell["y"].max())
            cellbin.create_dataset("cellBorder", data=border_array)
            cellbin.create_dataset("gene", data=gene_array)
            cellbin.create_dataset("cellExp", data=cell_exp_array)
    return path


def write_tissue_gef(
    path: Path,
    *,
    entries=((0, 0, 1), (1, 1, 2), (2, 0, 3), (3, 3, 4), (4, 1, 5)),
    genes=STO_GENES,
    bin_key: str = "bin1",
    resolution: int = 500,
) -> Path:
    """Write a STOmics squarebin ``*.tissue.gef`` (``geneExp/bin1``) file.

    ``entries`` is ``(x, y, count)`` with the gene implied by position: the
    expression table is gene-ordered, so gene *g* owns entries
    ``[offset[g], offset[g] + count[g])``.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    genes = [str(g) for g in genes]

    expr = np.zeros(len(entries), dtype=BIN_EXPR_DTYPE)
    gene_dtype = np.zeros(len(genes), dtype=BIN_GENE_DTYPE)
    # Spread entries evenly over genes, gene-ordered.
    per_gene = np.full(len(genes), len(entries) // len(genes), dtype=np.int64)
    per_gene[: len(entries) % len(genes)] += 1
    running = 0
    for g, n in enumerate(per_gene):
        gene_dtype["gene"][g] = genes[g].encode()
        gene_dtype["offset"][g] = running
        gene_dtype["count"][g] = n
        running += n
    for i, (x, y, count) in enumerate(entries):
        expr["x"][i] = x
        expr["y"][i] = y
        expr["count"][i] = count

    with h5py.File(path, "w") as handle:
        gene_exp = handle.create_group("geneExp")
        bin_group = gene_exp.create_group(bin_key)
        expr_ds = bin_group.create_dataset("expression", data=expr)
        expr_ds.attrs["resolution"] = resolution
        expr_ds.attrs["minX"] = int(expr["x"].min())
        expr_ds.attrs["maxX"] = int(expr["x"].max())
        expr_ds.attrs["minY"] = int(expr["y"].min())
        expr_ds.attrs["maxY"] = int(expr["y"].max())
        bin_group.create_dataset("gene", data=gene_dtype)
    return path


def write_gem(
    path: Path, rows=(("GeneA", 0, 0, 1), ("GeneB", 1, 1, 2), ("GeneA", 2, 0, 3))
) -> Path:
    """Write a tiny STOmics ``.gem`` text file (geneID, x, y, MIDCount)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# GEM written by tests.conftest", "geneID\tx\ty\tMIDCount"]
    lines += [f"{gene}\t{x}\t{y}\t{count}" for gene, x, y, count in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


@pytest.fixture
def sto_cellbin_bundle(tmp_path: Path) -> Path:
    """``<tmp>/C57_7/041.cellcut/B000001.cellbin.gef`` in the real nest."""
    cut_dir = tmp_path / "C57_7" / "041.cellcut"
    cut_dir.mkdir(parents=True)
    return write_cellbin_gef(cut_dir / "B000001.cellbin.gef")


@pytest.fixture
def sto_tissue_bundle(tmp_path: Path) -> Path:
    """``<tmp>/C57_7/04.tissuecut/B000001.tissue.gef`` in the real nest."""
    cut_dir = tmp_path / "C57_7" / "04.tissuecut"
    cut_dir.mkdir(parents=True)
    return write_tissue_gef(cut_dir / "B000001.tissue.gef")
