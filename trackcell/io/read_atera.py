"""
10x Atera readers for TrackCell.

Atera produces two bundle flavours:

``preview``
    Xenium-v4 compatible files (``cell_feature_matrix.h5`` +
    ``cells.parquet`` + ``cell_boundaries.parquet``).
``cellseg`` (official, Atera 1.x)
    AnnData Zarr (``cell_feature_matrix.zarr.zip``), segmentation Zarr
    (``cells.zarr.zip``) and single-molecule point cloud
    (``transcripts.zarr.zip``).

Design notes
------------
* ``obsm['spatial']`` is always in micrometres.  Official
  ``centroid_column`` / ``centroid_row`` are already µm; never multiply
  them by ``pixel_size`` (that only applies to morphology images).
* ``var_names`` are gene symbols by default (``gene_name_index=True``),
  matching :func:`trackcell.io.read_xenium_cellseg`.
* ``var['gene_ids']`` / ``var['feature_types']`` / ``var['is_gene_expression']``
  are available for both flavours.
* ``transcripts.zarr.zip`` may contain >1e9 molecules; use
  :func:`iter_atera_transcripts` (streaming) rather than
  :func:`read_atera_transcripts` unless the query is small.
* Atera is gene-level: there is **no** isoform, splice-junction,
  spliced/unspliced or RNA-velocity information.  ``overlaps_nucleus`` is
  a nuclear-localisation flag, not an unspliced count.
* Atera has no customer barcodes.  The delivery unit is an ROI; use
  :func:`list_atera_rois` / :func:`iter_atera_rois` and
  :func:`concat_atera_rois` for the same customer's ROIs.
"""

from __future__ import annotations

import json
import time
import warnings
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Mapping, Optional, Sequence, Tuple, Union

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.utils import make_index_unique
from shapely import wkt
from shapely.geometry import Polygon
from shapely.geometry.base import BaseGeometry

try:  # zarr is a required dependency of anndata, but keep a clear error.
    import zarr
except ImportError:  # pragma: no cover - defensive
    zarr = None

try:
    # zarr 3 moved stores out of the top-level namespace. Import through
    # ``zarr.storage`` so both zarr 2.18 and 3.x work.
    from zarr.storage import ZipStore as _ZipStore
except ImportError:  # pragma: no cover - defensive
    _ZipStore = getattr(zarr, "ZipStore", None)

try:
    import h5py
except ImportError:  # pragma: no cover - defensive
    h5py = None

from .read_data import read_xenium_cellseg

_MATRIX_ZIP = "cell_feature_matrix.zarr.zip"
_CELLS_ZIP = "cells.zarr.zip"
_TX_ZIP = "transcripts.zarr.zip"
_MATRIX_H5 = "cell_feature_matrix.h5"
_EXPERIMENT_SPATIAL = "experiment.spatial"
_EXPERIMENT_XENIUM = "experiment.xenium"
_PANEL_CONFIG = "panel_config.json"
_GENE_PANEL = "gene_panel.json"

_UNASSIGNED = np.uint32(4294967295)
_UM = "um"

__all__ = [
    "atera_info",
    "read_atera",
    "read_atera_cellseg",
    "read_atera_preview",
    "read_atera_panel",
    "list_atera_rois",
    "iter_atera_rois",
    "concat_atera_rois",
    "iter_atera_transcripts",
    "read_atera_transcripts",
    "atera_cell_qc",
    "atera_transcript_qc",
    "add_wkt_geometry",
    "filter_atera_genes",
]


# ---------------------------------------------------------------------------
# low level helpers
# ---------------------------------------------------------------------------
def _require_zarr():
    if zarr is None:  # pragma: no cover - defensive
        raise ImportError(
            "The `zarr` package is required to read Atera Zarr bundles. "
            "Install it with `pip install zarr`."
        )


def _log(verbose: bool):
    return print if verbose else (lambda *args, **kwargs: None)


def _member(root: Path, zipped: str, unzipped: Optional[str] = None) -> Optional[Path]:
    """Return the zipped member if present, else an unzipped directory."""
    p = root / zipped
    if p.exists():
        return p
    if unzipped is not None:
        q = root / unzipped
        if q.exists():
            return q
    return None


@contextmanager
def _zarr_group(path: Union[str, Path]):
    """Open a zarr group from a ``.zip`` store or an unzipped directory."""
    _require_zarr()
    path = Path(path)
    store = None
    if path.is_dir():
        group = zarr.open_group(str(path), mode="r")
    else:
        store = _ZipStore(str(path), mode="r")
        group = zarr.open_group(store, mode="r")
    try:
        yield group
    finally:
        if store is not None:
            store.close()


def _read_adata_zarr(path: Path) -> ad.AnnData:
    _require_zarr()
    if path.is_dir():
        return ad.read_zarr(str(path))
    store = _ZipStore(str(path), mode="r")
    try:
        return ad.read_zarr(store)
    finally:
        store.close()


def _read_zarr_column(parent, name: str) -> np.ndarray:
    """Read a (possibly categorical) zarr dataframe column as strings."""
    node = parent[name]
    attrs = getattr(node, "attrs", {})
    encoding = attrs.get("encoding-type") if attrs else None
    if encoding == "categorical" or (
        hasattr(node, "keys") and "codes" in node and "categories" in node
    ):
        codes = np.asarray(node["codes"][:])
        cats = np.asarray(node["categories"][:])
        cats = np.array(
            [c.decode("utf-8") if isinstance(c, bytes) else str(c) for c in cats],
            dtype=object,
        )
        return cats[codes]
    arr = np.asarray(node[:])
    if arr.dtype.kind in ("S", "O"):
        arr = np.array(
            [a.decode("utf-8") if isinstance(a, bytes) else a for a in arr],
            dtype=object,
        )
    return arr


def _detect_format(path: Union[str, Path], fmt: str = "auto") -> Tuple[str, Path]:
    """Return ``(format, root)`` with format in ``{'cellseg', 'preview'}``."""
    root = Path(path).expanduser()
    if not root.exists():
        raise FileNotFoundError(f"Atera path does not exist: {root}")
    if root.is_file():
        raise ValueError(
            f"Expected an Atera output directory, got a file: {root}. "
            "Pass the directory containing cell_feature_matrix.* instead."
        )

    has_official = (
        _member(root, _MATRIX_ZIP, "cell_feature_matrix.zarr") is not None
    )
    has_preview = (root / _MATRIX_H5).exists()

    if fmt in ("auto", None):
        if has_official and has_preview:
            warnings.warn(
                "Both official Zarr and preview HDF5 matrices found; "
                "reading the official Zarr bundle. Pass format='preview' "
                "to read the preview bundle.",
                stacklevel=3,
            )
            return "cellseg", root
        if has_official:
            return "cellseg", root
        if has_preview:
            return "preview", root
        children = [p for p in sorted(root.iterdir()) if p.is_dir()]
        hint = ""
        if any(
            (child / _MATRIX_ZIP).exists()
            or (child / _MATRIX_H5).exists()
            or (child / "cell_feature_matrix.zarr").exists()
            for child in children
        ):
            hint = (
                " This looks like a run directory. Use `tcl.io.list_atera_rois()` "
                "or `tcl.io.iter_atera_rois()` instead."
            )
        raise FileNotFoundError(
            f"No Atera matrix found in {root}. Expected "
            f"'{_MATRIX_ZIP}' / 'cell_feature_matrix.zarr' (official) or "
            f"'{_MATRIX_H5}' (preview).{hint}"
        )

    fmt = fmt.lower()
    if fmt in ("official", "zarr"):
        fmt = "cellseg"
    if fmt == "cellseg" and not has_official:
        raise FileNotFoundError(
            f"Official Atera bundle not found in {root}: '{_MATRIX_ZIP}' missing."
        )
    if fmt == "preview" and not has_preview:
        raise FileNotFoundError(
            f"Preview Atera bundle not found in {root}: '{_MATRIX_H5}' missing."
        )
    if fmt not in ("cellseg", "preview"):
        raise ValueError(
            f"Unknown Atera format {fmt!r}. Use 'auto', 'cellseg' or 'preview'."
        )
    return fmt, root


def _read_experiment(root: Path) -> dict:
    for name in (_EXPERIMENT_SPATIAL, _EXPERIMENT_XENIUM):
        p = root / name
        if p.exists():
            try:
                return json.loads(p.read_text())
            except (OSError, json.JSONDecodeError) as exc:
                warnings.warn(f"Could not parse {p}: {exc}")
    return {}


def _default_sample(path: Path, experiment: Mapping) -> str:
    for key in ("region_name", "sample_id", "slide_name"):
        value = experiment.get(key)
        if value:
            return str(value)
    return path.name


# ---------------------------------------------------------------------------
# metadata helpers
# ---------------------------------------------------------------------------
def _var_feature_types(adata: ad.AnnData) -> Optional[np.ndarray]:
    for key in ("feature_types", "feature_type"):
        if key in adata.var.columns:
            return adata.var[key].astype(str).to_numpy()
    return None


def _normalize_var(adata: ad.AnnData, gene_name_index: bool = True,
                   keep_de: bool = False) -> ad.AnnData:
    """Normalise ``var`` to the trackcell schema.

    Result
    ------
    ``var_names``      gene symbols (``gene_name_index=True``) or feature IDs
    ``var['gene_ids']`` Ensembl / feature IDs
    ``var['feature_types']`` one of ``Gene Expression``, ``Negative Control
        Codeword``, ...
    ``var['is_gene_expression']`` boolean mask
    """
    var = adata.var

    # 1. Always expose the feature IDs in a stable column.
    if "gene_ids" not in var.columns:
        if "feature_id" in var.columns:
            var["gene_ids"] = var["feature_id"].astype(str)
        else:
            var["gene_ids"] = adata.var_names.astype(str)

    # 2. Choose var_names.  Preview bundles already use symbols as the index;
    #    official bundles index by feature_id but carry `feature_name`.
    if gene_name_index:
        if "feature_name" in var.columns:
            names = var["feature_name"].astype(str).to_numpy()
        else:
            names = adata.var_names.astype(str)
        index = pd.Index(names, name="gene_symbol")
        if not index.is_unique:
            warnings.warn(
                "Duplicate gene symbols in `feature_name`; making var_names "
                "unique with suffixes.",
                stacklevel=3,
            )
            index = make_index_unique(index)
        adata.var_names = index
    else:
        id_col = "feature_id" if "feature_id" in var.columns else (
            "gene_ids" if "gene_ids" in var.columns else None
        )
        if id_col is not None:
            index = pd.Index(var[id_col].astype(str).to_numpy(),
                             name="feature_id")
            if not index.is_unique:
                index = make_index_unique(index)
            adata.var_names = index
        else:
            adata.var.index = adata.var.index.astype(str)
            adata.var.index.name = adata.var.index.name or "feature_id"

    var = adata.var
    if "feature_types" not in var.columns:
        if "feature_type" in var.columns:
            var["feature_types"] = var["feature_type"].astype(str)
        else:
            var["feature_types"] = "Gene Expression"

    var["is_gene_expression"] = (
        var["feature_types"].astype(str).str.strip().str.lower()
        == "gene expression"
    )

    if not keep_de:
        drop = [c for c in var.columns if str(c).startswith("de_")]
        if drop:
            var.drop(columns=drop, inplace=True)

    return adata


def _add_cell_ids(adata: ad.AnnData, cells_path: Optional[Path],
                  verbose: bool = True) -> Optional[np.ndarray]:
    """Attach ``cell_prefix`` / ``dataset_suffix`` and return the
    ``cells.zarr`` row -> ``adata.obs`` row map (``None`` when unavailable)."""
    log = _log(verbose)
    if cells_path is None:
        return None

    with _zarr_group(cells_path) as group:
        cell_ids = np.asarray(group["cell_id"][:], dtype=np.uint32)
        summary = np.asarray(group["cell_summary"][:]) if "cell_summary" in group else None
        summary_names = (
            list(group["cell_summary"].attrs.get("column_names", []))
            if summary is not None else []
        )

    n_cells = len(cell_ids)
    prefix = cell_ids[:, 0]
    suffix = cell_ids[:, 1]

    if n_cells != adata.n_obs:
        warnings.warn(
            f"`cells.zarr` has {n_cells:,} cells but the matrix has "
            f"{adata.n_obs:,}; polygon/cell-id mapping will be skipped.",
            stacklevel=3,
        )
        return None

    if "cell_id" not in adata.obs.columns:
        adata.obs["cell_id"] = (
            (suffix.astype(np.uint64) << np.uint64(32)) | prefix.astype(np.uint64)
        ).astype(np.int64)
    adata.obs["cell_prefix"] = prefix
    adata.obs["dataset_suffix"] = suffix

    # Verify that cells.zarr rows are aligned with matrix obs rows.
    if summary is not None and summary_names:
        try:
            ix = summary_names.index("cell_centroid_x")
            iy = summary_names.index("cell_centroid_y")
        except ValueError:
            ix = iy = None
        if ix is not None:
            obs_x = adata.obs.get("centroid_column")
            obs_y = adata.obs.get("centroid_row")
            if obs_x is not None and obs_y is not None:
                obs_xy = np.column_stack(
                    [np.asarray(obs_x, dtype=np.float64),
                     np.asarray(obs_y, dtype=np.float64)]
                )
                cell_xy = summary[:, [ix, iy]].astype(np.float64)
                diff = np.abs(obs_xy - cell_xy).max()
                if diff > 0.5:
                    mapping = _match_rows_by_xy(obs_xy, cell_xy)
                    if mapping is not None:
                        return mapping
                    warnings.warn(
                        "`cells.zarr` row order does not match the matrix obs "
                        f"order (max centroid diff {diff:.3f} µm); polygons will "
                        "use positional mapping. Please report this dataset.",
                        stacklevel=3,
                    )
    return np.arange(n_cells, dtype=np.int64)


def _match_rows_by_xy(obs_xy: np.ndarray, cell_xy: np.ndarray) -> Optional[np.ndarray]:
    """Map cells.zarr rows to obs rows using rounded centroids."""
    keys = {(round(float(x), 3), round(float(y), 3)): i
            for i, (x, y) in enumerate(obs_xy)}
    mapping = np.empty(len(cell_xy), dtype=np.int64)
    seen = set()
    for i, (x, y) in enumerate(cell_xy):
        key = (round(float(x), 3), round(float(y), 3))
        j = keys.get(key)
        if j is None or j in seen:
            return None
        mapping[i] = j
        seen.add(j)
    return mapping


def _parse_region(
    region,
) -> Optional[Tuple[Tuple[float, float, float, float], Optional[BaseGeometry]]]:
    """Return ``((xmin, ymin, xmax, ymax), geometry)``.

    Accepts ``(xmin, ymin, xmax, ymax)``, a shapely geometry, or a mapping
    with xmin/ymin/xmax/ymax keys.
    """
    if region is None:
        return None
    if isinstance(region, BaseGeometry):
        geom = region
        return tuple(float(v) for v in geom.bounds), geom
    if isinstance(region, Mapping):
        bounds = (
            float(region["xmin"]), float(region["ymin"]),
            float(region["xmax"]), float(region["ymax"]),
        )
        return bounds, None
    values = tuple(region)
    if len(values) != 4:
        raise ValueError(
            "`region` must be (xmin, ymin, xmax, ymax), a shapely geometry, "
            f"or a mapping with xmin/ymin/xmax/ymax; got {region!r}."
        )
    xmin, ymin, xmax, ymax = (float(v) for v in values)
    if xmin > xmax or ymin > ymax:
        raise ValueError(
            "`region` must be (xmin, ymin, xmax, ymax) in µm; "
            f"got {region!r} (xmin > xmax or ymin > ymax)."
        )
    return (xmin, ymin, xmax, ymax), None


def _region_mask(x: np.ndarray, y: np.ndarray, parsed_region) -> np.ndarray:
    if parsed_region is None:
        return np.ones(len(x), dtype=bool)
    (xmin, ymin, xmax, ymax), geom = parsed_region
    mask = (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)
    if geom is not None:
        try:
            from shapely import contains_xy
            mask &= contains_xy(geom, x, y)
        except ImportError:  # pragma: no cover - shapely < 2.0
            from shapely.geometry import Point
            mask &= np.array([geom.contains(Point(xi, yi)) for xi, yi in zip(x, y)])
    return mask


def _cell_mask(adata: ad.AnnData, region=None, n_cells: Optional[int] = None,
               random_state: int = 0) -> np.ndarray:
    parsed = _parse_region(region)
    n = adata.n_obs
    mask = np.ones(n, dtype=bool)
    if parsed is not None:
        if "spatial" not in adata.obsm:
            warnings.warn(
                "`region` was provided but `obsm['spatial']` is missing; "
                "no cells were selected by region.",
                stacklevel=3,
            )
        else:
            coords = np.asarray(adata.obsm["spatial"])
            mask &= _region_mask(coords[:, 0], coords[:, 1], parsed)
    if n_cells is not None and n_cells < int(mask.sum()):
        rng = np.random.default_rng(random_state)
        candidates = np.flatnonzero(mask)
        chosen = np.sort(rng.choice(candidates, size=int(n_cells), replace=False))
        mask = np.zeros(n, dtype=bool)
        mask[chosen] = True
    return mask


def _subset_geometries(adata: ad.AnnData) -> None:
    """Drop geometries that are no longer in ``adata.obs_names``."""
    for sample, info in adata.uns.get("spatial", {}).items():
        if not isinstance(info, Mapping):
            continue
        for key in ("geometries", "nucleus_geometries"):
            gdf = info.get(key)
            if isinstance(gdf, gpd.GeoDataFrame):
                info[key] = gdf.loc[gdf.index.intersection(adata.obs_names)]


def _subset_cells(adata: ad.AnnData, mask: np.ndarray) -> ad.AnnData:
    if mask.all():
        return adata
    adata = adata[mask].copy()
    _subset_geometries(adata)
    return adata


# ---------------------------------------------------------------------------
# polygons
# ---------------------------------------------------------------------------
def _polygons_from_group(group) -> Tuple[list, np.ndarray]:
    """Return ``(polygons, cell_index)`` from a ``polygon_sets/{0,1}`` group."""
    num_vertices = np.asarray(group["num_vertices"][:], dtype=np.int64)
    vertices = np.asarray(group["vertices"][:], dtype=np.float32)
    max_vertices = vertices.shape[1] // 2
    vertices = vertices.reshape(len(num_vertices), max_vertices, 2)
    cell_index = np.asarray(group["cell_index"][:], dtype=np.int64)

    polygons, records = [], []
    for i, n in enumerate(num_vertices):
        if n < 3:
            continue
        poly = Polygon(vertices[i, : int(n)])
        if poly.is_valid and not poly.is_empty:
            polygons.append(poly)
            records.append(i)
    if not polygons:
        return [], np.array([], dtype=np.int64)
    return polygons, cell_index[np.asarray(records, dtype=np.int64)]


def _attach_polygons(adata: ad.AnnData, group, key: str, sample: str,
                     row_map: np.ndarray, original_to_subset: np.ndarray,
                     verbose: bool = True) -> None:
    """Build a GeoDataFrame for ``key`` and attach it to ``uns['spatial']``."""
    log = _log(verbose)
    polygons, cell_index = _polygons_from_group(group)
    if not polygons:
        warnings.warn(f"No valid polygons found in `{key}`.")
        return

    valid = (cell_index >= 0) & (cell_index < len(row_map))
    obs_rows = np.full(len(polygons), -1, dtype=np.int64)
    obs_rows[valid] = original_to_subset[row_map[cell_index[valid]]]
    keep = obs_rows >= 0
    if not keep.any():
        warnings.warn(f"No `{key}` polygons map to current cells.")
        return

    polygons = [p for p, k in zip(polygons, keep) if k]
    obs_rows = obs_rows[keep]

    gdf = gpd.GeoDataFrame(
        {"geometry": polygons}, index=adata.obs_names[obs_rows], geometry="geometry"
    )
    gdf.index.name = "cell_id"
    adata.uns["spatial"][sample][key] = gdf

    atera_uns = adata.uns.setdefault("atera", {})
    atera_uns[f"n_{key}"] = len(polygons)
    log(f"[Atera] {key}: {len(polygons):,} polygons "
        f"({len(polygons) / adata.n_obs * 100:.1f}% of cells)")


def _attach_compact_boundaries(adata: ad.AnnData, sample: str, key: str,
                               uns_key: str) -> None:
    """Store compact polygon vertices (per-cell offsets) in ``uns``."""
    gdf = adata.uns.get("spatial", {}).get(sample, {}).get(key)
    if not isinstance(gdf, gpd.GeoDataFrame) or not len(gdf):
        return
    rings = [np.asarray(p.exterior.coords)[:-1] for p in gdf.geometry]
    n_vertices = np.array([len(r) for r in rings], dtype=np.int32)
    vertex_x = np.concatenate([r[:, 0] for r in rings]).astype(np.float32)
    vertex_y = np.concatenate([r[:, 1] for r in rings]).astype(np.float32)

    cell_idx = np.full(adata.n_obs, -1, dtype=np.int32)
    positions = adata.obs_names.get_indexer(gdf.index)
    valid = positions >= 0
    cell_idx[positions[valid]] = np.arange(len(gdf), dtype=np.int32)[valid]
    adata.uns[uns_key] = {
        "vertex_x": vertex_x,
        "vertex_y": vertex_y,
        "n_vertices": n_vertices,
        "cell_idx": cell_idx,
    }


# ---------------------------------------------------------------------------
# public readers
# ---------------------------------------------------------------------------
def read_atera(
    path: Union[str, Path],
    sample: Optional[str] = None,
    *,
    format: str = "auto",
    gene_name_index: bool = True,
    build_polygons: bool = True,
    include_nucleus: Optional[bool] = None,
    store_wkt: bool = False,
    keep_de: bool = False,
    load_panel: bool = False,
    region=None,
    n_cells: Optional[int] = None,
    random_state: int = 0,
    verbose: bool = True,
) -> ad.AnnData:
    """Read an Atera output bundle (auto-detecting preview vs official).

    Parameters
    ----------
    path
        Atera output directory.  This is usually one ROI bundle, not the run
        directory (use :func:`list_atera_rois` / :func:`iter_atera_rois` for
        a run).
    sample
        Key used for ``uns['spatial'][sample]``.  If ``None``, uses
        ``experiment.spatial['region_name']`` (official) or
        ``experiment.xenium['region_name']`` (preview), falling back to the
        directory name.
    format
        ``'auto'`` (default), ``'cellseg'`` (official Zarr) or ``'preview'``.
    gene_name_index
        If True, use gene symbols as ``var_names`` and store Ensembl IDs in
        ``var['gene_ids']``.  If False, keep the Ensembl IDs as ``var_names``.
    build_polygons
        Build cell polygons from ``cells.zarr`` / ``cell_boundaries.parquet``.
        Has a real cost at ~400k cells (about one minute for official data).
    include_nucleus
        Also build nucleus polygons.  ``None`` (default) follows
        ``build_polygons``: nucleus polygons are built when cell polygons are
        built.  Explicitly passing ``True`` with ``build_polygons=False``
        raises an error.
    store_wkt
        Store WKT strings in ``obs['geometry']``.  Off by default: it costs
        memory and is only needed for serialization.  Use
        :func:`add_wkt_geometry` to add it later.
    keep_de
        Keep the official per-cluster ``var['de_*']`` columns (official data
        has hundreds of them and they are rarely needed at read time).
    load_panel
        Parse ``panel_config.json`` (official, ~11 MB and may contain customer
        information) into ``uns['gene_panel']``.  Off by default; only the
        path is recorded.
    region
        ``(xmin, ymin, xmax, ymax)`` in µm, a shapely geometry, or a mapping
        with ``xmin/ymin/xmax/ymax``.  Cells are subset by centroid.
    n_cells
        Randomly subset to at most this many cells (``random_state``).
    random_state
        Seed for ``n_cells``.
    verbose
        Print progress information.

    Returns
    -------
    AnnData
        With ``obsm['spatial']`` (µm), ``obs['barcode']``, polygons in
        ``uns['spatial'][sample]['geometries']`` and normalised ``var``
        columns.  Official data additionally has ``obs['cell_id']`` (uint32),
        ``obs['cell_prefix']`` and ``obs['dataset_suffix']``.

    Examples
    --------
    >>> adata = tcl.io.read_atera("/path/to/ROI_bundle")
    >>> adata = tcl.io.read_atera("/path/to/preview", format="preview")
    """
    if include_nucleus is None:
        include_nucleus = build_polygons
    elif include_nucleus and not build_polygons:
        raise ValueError(
            "`include_nucleus=True` requires `build_polygons=True`."
        )
    fmt, root = _detect_format(path, format)
    kwargs = dict(
        sample=sample,
        gene_name_index=gene_name_index,
        build_polygons=build_polygons,
        include_nucleus=include_nucleus,
        store_wkt=store_wkt,
        keep_de=keep_de,
        load_panel=load_panel,
        region=region,
        n_cells=n_cells,
        random_state=random_state,
        verbose=verbose,
    )
    if fmt == "preview":
        return read_atera_preview(root, **kwargs)
    return read_atera_cellseg(root, **kwargs)


def read_atera_cellseg(
    path: Union[str, Path],
    sample: Optional[str] = None,
    *,
    gene_name_index: bool = True,
    build_polygons: bool = True,
    include_nucleus: Optional[bool] = None,
    store_wkt: bool = False,
    keep_de: bool = False,
    load_panel: bool = False,
    region=None,
    n_cells: Optional[int] = None,
    random_state: int = 0,
    verbose: bool = True,
) -> ad.AnnData:
    """Read an official Atera 1.x bundle (AnnData Zarr + cells Zarr).

    See :func:`read_atera` for the parameter documentation.
    """
    if include_nucleus is None:
        include_nucleus = build_polygons
    elif include_nucleus and not build_polygons:
        raise ValueError(
            "`include_nucleus=True` requires `build_polygons=True`."
        )
    log = _log(verbose)
    t_total = time.time()
    root = Path(path).expanduser().resolve()
    experiment = _read_experiment(root)
    sample = sample or _default_sample(root, experiment)

    matrix_path = _member(root, _MATRIX_ZIP, "cell_feature_matrix.zarr")
    if matrix_path is None:
        raise FileNotFoundError(
            f"Official Atera matrix not found in {root}: '{_MATRIX_ZIP}'."
        )
    cells_path = _member(root, _CELLS_ZIP, "cells.zarr")
    tx_path = _member(root, _TX_ZIP, "transcripts.zarr")

    t0 = time.time()
    adata = _read_adata_zarr(matrix_path)
    nnz = getattr(adata.X, "nnz", None)
    nnz_msg = f"{nnz / 1e6:.1f}M non-zero" if nnz is not None else "dense"
    log(f"[Atera] Matrix: {adata.n_obs:,} cells × {adata.n_vars:,} features "
        f"({nnz_msg}) [{time.time() - t0:.1f}s]")

    # ── spatial coordinates ──────────────────────────────────────────
    if {"centroid_column", "centroid_row"}.issubset(adata.obs.columns):
        adata.obsm["spatial"] = adata.obs[
            ["centroid_column", "centroid_row"]
        ].to_numpy(np.float32)
    else:
        warnings.warn(
            "Official bundle has no `centroid_column` / `centroid_row`; "
            "`obsm['spatial']` was not set."
        )

    adata.obs["barcode"] = adata.obs_names.astype(str)
    adata.obs_names.name = "barcode"

    # ── metadata ─────────────────────────────────────────────────────
    if experiment:
        adata.uns["experiment"] = experiment
    adata.uns["atera"] = {
        "platform": "atera",
        "format": "cellseg",
        "sample": sample,
        "roi": experiment.get("region_name", sample),
        "slide": experiment.get("slide_name"),
        "run": experiment.get("run_name"),
        "coord_unit": _UM,
        "spatial_units": "micron",
        "pixel_size": experiment.get("pixel_size"),
        "z_step_size": experiment.get("z_step_size"),
        "chemistry_version": experiment.get("chemistry_version"),
        "analysis_version": experiment.get("analysis_sw_version"),
        "major_version": experiment.get("major_version"),
        "minor_version": experiment.get("minor_version"),
        "panel_path": str(root / _PANEL_CONFIG) if (root / _PANEL_CONFIG).exists() else None,
        "transcripts_path": str(tx_path) if tx_path is not None else None,
    }

    # ── cell ids / alignment ─────────────────────────────────────────
    row_map = _add_cell_ids(adata, cells_path, verbose=verbose)

    # ── var schema ───────────────────────────────────────────────────
    _normalize_var(adata, gene_name_index=gene_name_index, keep_de=keep_de)
    n_expr = int(adata.var["is_gene_expression"].sum())
    log(f"[Atera] Features: {adata.n_vars:,} total, {n_expr:,} Gene Expression "
        f"(use `tcl.io.filter_atera_genes(adata)` before HVG/PCA)")

    # ── cell subsetting (before polygon construction) ────────────────
    mask = _cell_mask(adata, region=region, n_cells=n_cells,
                      random_state=random_state)
    if not mask.all():
        original_to_subset = np.full(adata.n_obs, -1, dtype=np.int64)
        original_to_subset[np.flatnonzero(mask)] = np.arange(int(mask.sum()))
        adata = _subset_cells(adata, mask)
        log(f"[Atera] Subset to {adata.n_obs:,} cells")
    else:
        original_to_subset = np.arange(adata.n_obs, dtype=np.int64)

    # ── polygons ─────────────────────────────────────────────────────
    if build_polygons and cells_path is not None and row_map is not None:
        t0 = time.time()
        adata.uns.setdefault("spatial", {})[sample] = {
            "metadata": {
                "sample_id": sample,
                "pixel_size": experiment.get("pixel_size", 0.2125),
                "coord_unit": _UM,
            },
            "scalefactors": {
                "spot_diameter_fullres": 1.0,
                "tissue_hires_scalef": 1.0,
                "tissue_lowres_scalef": 1.0,
            },
        }
        with _zarr_group(cells_path) as group:
            if "polygon_sets/1" in group:
                _attach_polygons(adata, group["polygon_sets/1"],
                                 "geometries", sample, row_map,
                                 original_to_subset, verbose=verbose)
                _attach_compact_boundaries(
                    adata, sample, "geometries", "cell_boundaries"
                )
            else:
                warnings.warn("`polygon_sets/1` (cell polygons) not found.")
            if include_nucleus and "polygon_sets/0" in group:
                _attach_polygons(adata, group["polygon_sets/0"],
                                 "nucleus_geometries", sample, row_map,
                                 original_to_subset, verbose=verbose)
                _attach_compact_boundaries(
                    adata, sample, "nucleus_geometries", "nucleus_boundaries"
                )
        if store_wkt:
            add_wkt_geometry(adata, sample=sample)
        log(f"[Atera] Polygons built [{time.time() - t0:.1f}s]")
    else:
        if build_polygons:
            warnings.warn(
                "Polygons requested but `cells.zarr` is missing or could not be "
                "aligned; `spatial_cell()` will not be able to draw cells."
            )
        adata.uns.setdefault("spatial", {})[sample] = {
            "metadata": {
                "sample_id": sample,
                "pixel_size": experiment.get("pixel_size", 0.2125),
                "coord_unit": _UM,
            },
            "scalefactors": {
                "spot_diameter_fullres": 1.0,
                "tissue_hires_scalef": 1.0,
                "tissue_lowres_scalef": 1.0,
            },
        }

    # ── optional panel ───────────────────────────────────────────────
    panel_path = root / _PANEL_CONFIG
    if load_panel and panel_path.exists():
        adata.uns["gene_panel"] = read_atera_panel(root)

    log(f"[Atera] Total: read_atera_cellseg done [{time.time() - t_total:.1f}s]")
    return adata


def read_atera_preview(
    path: Union[str, Path],
    sample: Optional[str] = None,
    *,
    gene_name_index: bool = True,
    build_polygons: bool = True,
    include_nucleus: Optional[bool] = None,
    store_wkt: bool = False,
    keep_de: bool = False,
    load_panel: bool = False,
    region=None,
    n_cells: Optional[int] = None,
    random_state: int = 0,
    verbose: bool = True,
) -> ad.AnnData:
    """Read an Atera preview bundle (Xenium-v4 compatible formats).

    See :func:`read_atera` for the parameter documentation.  The preview
    bundle uses ``cell_feature_matrix.h5`` + ``cells.parquet`` +
    ``cell_boundaries.parquet``.
    """
    if include_nucleus is None:
        include_nucleus = build_polygons
    elif include_nucleus and not build_polygons:
        raise ValueError(
            "`include_nucleus=True` requires `build_polygons=True`."
        )
    log = _log(verbose)
    root = Path(path).expanduser().resolve()
    experiment = _read_experiment(root)
    sample = sample or _default_sample(root, experiment)

    adata = read_xenium_cellseg(
        str(root),
        sample=sample,
        build_polygons=build_polygons,
        include_nucleus=include_nucleus,
        store_wkt=store_wkt,
        verbose=verbose,
    )

    # unify the obs/var schema with the official reader
    adata.obs["barcode"] = adata.obs_names.astype(str)
    adata.obs_names.name = "barcode"
    _normalize_var(adata, gene_name_index=gene_name_index, keep_de=keep_de)
    adata.uns["atera"] = {
        "platform": "atera",
        "format": "preview",
        "sample": sample,
        "roi": experiment.get("region_name", sample),
        "slide": experiment.get("slide_name") or experiment.get("slide_id"),
        "run": experiment.get("run_name"),
        "coord_unit": _UM,
        "pixel_size": experiment.get("pixel_size"),
        "z_step_size": experiment.get("z_step_size"),
        "chemistry_version": experiment.get("chemistry_version"),
        "analysis_version": experiment.get("analysis_sw_version"),
        "panel_path": str(root / _GENE_PANEL) if (root / _GENE_PANEL).exists() else None,
        "transcripts_path": None,
    }
    adata.uns["spatial"][sample].setdefault("metadata", {})["coord_unit"] = _UM

    if not build_polygons:
        for key in ("geometries", "nucleus_geometries"):
            adata.uns["spatial"][sample].pop(key, None)
        adata.uns.pop("cell_boundaries", None)
        adata.uns.pop("nucleus_boundaries", None)
    elif not include_nucleus:
        adata.uns.pop("nucleus_boundaries", None)

    if not store_wkt:
        adata.obs.drop(columns=["geometry"], inplace=True, errors="ignore")
    if not load_panel:
        adata.uns.pop("gene_panel", None)

    mask = _cell_mask(adata, region=region, n_cells=n_cells,
                      random_state=random_state)
    adata = _subset_cells(adata, mask)
    if not mask.all():
        log(f"[Atera] Subset to {adata.n_obs:,} cells")
    return adata


def read_atera_panel(path: Union[str, Path]) -> dict:
    """Read ``panel_config.json`` from an Atera bundle.

    The file is ~11 MB and may contain customer / designer information, so
    :func:`read_atera` does not load it by default.
    """
    root = Path(path).expanduser()
    if root.is_file() and root.name == _PANEL_CONFIG:
        panel_path = root
    else:
        panel_path = root / _PANEL_CONFIG
    if not panel_path.exists():
        raise FileNotFoundError(f"Atera panel config not found: {panel_path}")
    return json.loads(panel_path.read_text())


def filter_atera_genes(adata: ad.AnnData, inplace: bool = False) -> ad.AnnData:
    """Keep only ``feature_type == 'Gene Expression'`` features.

    Atera panels include thousands of negative-control codewords.  They must
    be removed before HVG selection / PCA / clustering.

    Examples
    --------
    >>> adata = tcl.io.filter_atera_genes(adata)
    """
    if "is_gene_expression" in adata.var.columns:
        mask = adata.var["is_gene_expression"].to_numpy(dtype=bool)
    elif "feature_types" in adata.var.columns:
        mask = (
            adata.var["feature_types"].astype(str).str.strip().str.lower()
            == "gene expression"
        ).to_numpy()
    else:
        raise ValueError(
            "No `var['is_gene_expression']` / `var['feature_types']` column found. "
            "This function only applies to Atera / Xenium data."
        )
    if mask.all():
        return adata
    if inplace:
        adata._inplace_subset_var(mask)
        return adata
    return adata[:, mask].copy()


def add_wkt_geometry(adata: ad.AnnData, sample: Optional[str] = None,
                     key: str = "geometry", nucleus: bool = False) -> ad.AnnData:
    """Add WKT polygon strings to ``obs[key]``.

    This is the opt-in replacement for ``store_wkt=True``.
    """
    sample = sample or next(iter(adata.uns.get("spatial", {})), None)
    if sample is None:
        raise ValueError("No `adata.uns['spatial']` found.")
    geom_key = "nucleus_geometries" if nucleus else "geometries"
    gdf = adata.uns["spatial"][sample].get(geom_key)
    if not isinstance(gdf, gpd.GeoDataFrame):
        raise ValueError(
            f"No `uns['spatial']['{sample}']['{geom_key}']` GeoDataFrame found."
        )
    values = np.empty(adata.n_obs, dtype=object)
    for cid, poly in zip(gdf.index, gdf.geometry):
        values[adata.obs_names.get_loc(cid)] = wkt.dumps(poly)
    adata.obs[key] = pd.Series(values, index=adata.obs_names)
    return adata


# ---------------------------------------------------------------------------
# run / ROI helpers
# ---------------------------------------------------------------------------
def _is_bundle(path: Path) -> bool:
    return (
        _member(path, _MATRIX_ZIP, "cell_feature_matrix.zarr") is not None
        or (path / _MATRIX_H5).exists()
    )


def _quick_counts(path: Path) -> dict:
    """Cheap cell/gene counts without loading the expression matrix."""
    fmt, root = _detect_format(path)
    if fmt == "cellseg":
        matrix_path = _member(root, _MATRIX_ZIP, "cell_feature_matrix.zarr")
        with _zarr_group(matrix_path) as group:
            shape = tuple(int(v) for v in group["X"].attrs["shape"])
            feature_type = _read_zarr_column(group["var"], "feature_type")
        return {
            "format": "cellseg",
            "n_cells": shape[0],
            "n_features": shape[1],
            "n_gene_expression": int(
                (pd.Series(feature_type).astype(str).str.strip().str.lower()
                 == "gene expression").sum()
            ),
        }
    if h5py is None:  # pragma: no cover - defensive
        raise ImportError("`h5py` is required to inspect Atera preview bundles.")
    with h5py.File(root / _MATRIX_H5, "r") as handle:
        shape = tuple(int(v) for v in handle["matrix"]["shape"][:])
        ftype = np.asarray(handle["matrix"]["features"]["feature_type"][:])
    ftype = np.array([f.decode() if isinstance(f, bytes) else f for f in ftype])
    return {
        "format": "preview",
        "n_cells": shape[1],
        "n_features": shape[0],
        "n_gene_expression": int(
            (pd.Series(ftype).astype(str).str.strip().str.lower()
             == "gene expression").sum()
        ),
    }


def atera_info(path: Union[str, Path], verbose: bool = True) -> dict:
    """Return bundle metadata without loading the expression matrix.

    This is the recommended first call on an unfamiliar 10–500 GB Atera
    bundle: it reports cells, genes, polygon coverage, transcript count and
    software version in a couple of seconds.

    Examples
    --------
    >>> tcl.io.atera_info("/path/to/ROI_bundle")
    """
    fmt, root = _detect_format(path)
    info = {"path": str(root), "format": fmt}
    experiment = _read_experiment(root)
    info["version"] = {
        "chemistry": experiment.get("chemistry_version"),
        "analysis": experiment.get("analysis_sw_version"),
        "major_version": experiment.get("major_version"),
        "minor_version": experiment.get("minor_version"),
    }
    info.update({
        "roi": experiment.get("region_name"),
        "slide": experiment.get("slide_name") or experiment.get("slide_id"),
        "run": experiment.get("run_name"),
        "pixel_size": experiment.get("pixel_size"),
        "z_step_size": experiment.get("z_step_size"),
        "region_area_um2": experiment.get("region_area"),
    })
    info.update(_quick_counts(root))

    cells_path = _member(root, _CELLS_ZIP, "cells.zarr")
    if cells_path is not None:
        with _zarr_group(cells_path) as group:
            n_cell_polys = (
                int(group["polygon_sets/1/cell_index"].shape[0])
                if "polygon_sets/1" in group else 0
            )
            n_nuc_polys = (
                int(group["polygon_sets/0/cell_index"].shape[0])
                if "polygon_sets/0" in group else 0
            )
        info.update({
            "n_cell_polygons": n_cell_polys,
            "n_nucleus_polygons": n_nuc_polys,
            "cell_boundary_frac": (
                n_cell_polys / info["n_cells"] if info.get("n_cells") else None
            ),
        })

    tx_path = _member(root, _TX_ZIP, "transcripts.zarr")
    if tx_path is not None:
        with _zarr_group(tx_path) as group:
            info.update({
                "number_rnas": int(group.attrs.get("number_rnas", 0)),
                "n_genes": int(group.attrs.get("number_genes", 0)),
                "n_tiles": len(group["grid"].attrs.get("grid_keys", [])),
                "grid_size_um": float(group["grid"].attrs.get("grid_size", 0.0)),
                "coordinate_space": group.attrs.get("coordinate_space"),
                "spatial_units": group.attrs.get("spatial_units"),
            })

    if verbose:
        print(f"[Atera] {info['path']}")
        print(f"  format          : {info['format']}")
        if info.get("roi"):
            print(f"  roi/slide/run   : {info['roi']} / {info.get('slide')} / {info.get('run')}")
        print(f"  cells × features: {info['n_cells']:,} × {info['n_features']:,}"
              f" ({info['n_gene_expression']:,} Gene Expression)")
        if "n_cell_polygons" in info:
            print(f"  polygons        : {info['n_cell_polygons']:,} cell / "
                  f"{info['n_nucleus_polygons']:,} nucleus")
        if info.get("number_rnas") is not None:
            print(f"  transcripts     : {info['number_rnas']:,} molecules, "
                  f"{info['n_tiles']} grid tiles, {info['grid_size_um']:g} µm/tile")
        print(f"  coord unit      : µm (pixel_size={info.get('pixel_size')}, "
              f"z_step={info.get('z_step_size')})")
        print(f"  software        : chemistry={info['version']['chemistry']}, "
              f"analysis={info['version']['analysis']}")
    return info


def list_atera_rois(run_dir: Union[str, Path]) -> pd.DataFrame:
    """List the ROI bundles contained in an Atera run directory.

    Returns a DataFrame with one row per ROI (``roi``, ``path``, ``format``,
    ``n_cells``, ``n_features``, ``n_gene_expression``, ``has_transcripts``).
    """
    root = Path(run_dir).expanduser().resolve()
    if not root.exists():
        raise FileNotFoundError(f"Atera run directory not found: {root}")
    candidates = [p for p in sorted(root.iterdir()) if p.is_dir() and _is_bundle(p)]
    if _is_bundle(root):
        candidates.insert(0, root)
    if not candidates:
        raise FileNotFoundError(
            f"No Atera ROI bundles found in {root}. Each ROI should contain "
            f"'{_MATRIX_ZIP}' or '{_MATRIX_H5}'."
        )
    rows = []
    for path in candidates:
        try:
            counts = _quick_counts(path)
        except Exception as exc:  # pragma: no cover - malformed bundle
            warnings.warn(f"Skipping {path}: {exc}")
            continue
        experiment = _read_experiment(path)
        rows.append({
            "roi": experiment.get("region_name", path.name),
            "path": str(path),
            **counts,
            "has_transcripts": (
                _member(path, _TX_ZIP, "transcripts.zarr") is not None
            ),
            "number_rnas": _quick_number_rnas(path),
        })
    return pd.DataFrame(rows)


def _quick_number_rnas(path: Path) -> Optional[int]:
    tx_path = _member(path, _TX_ZIP, "transcripts.zarr")
    if tx_path is None:
        return None
    try:
        with _zarr_group(tx_path) as group:
            return int(group.attrs.get("number_rnas", 0))
    except Exception:  # pragma: no cover - defensive
        return None


def iter_atera_rois(run_dir: Union[str, Path],
                    rois: Optional[Sequence[str]] = None,
                    **read_kwargs) -> Iterator[Tuple[str, ad.AnnData]]:
    """Yield ``(roi_name, adata)`` one ROI at a time.

    Unlike a dict-returning helper, this keeps peak memory at one ROI.  Pass
    extra keyword arguments through to :func:`read_atera`.
    """
    table = list_atera_rois(run_dir)
    if rois is not None:
        wanted = set(rois)
        table = table[table["roi"].isin(wanted) | table["path"].isin(wanted)]
        if table.empty:
            raise ValueError(f"None of the requested ROIs were found: {list(rois)}")
    for _, row in table.iterrows():
        kwargs = dict(read_kwargs)
        kwargs.setdefault("sample", row["roi"])
        yield row["roi"], read_atera(row["path"], **kwargs)


def _translate_roi(adata: ad.AnnData, dx: float, dy: float) -> None:
    if "spatial" in adata.obsm:
        adata.obsm["spatial"] = (
            np.asarray(adata.obsm["spatial"], dtype=np.float32)
            + np.array([dx, dy], dtype=np.float32)
        )
    for col, delta in (("x_centroid", dx), ("y_centroid", dy)):
        if col in adata.obs.columns:
            adata.obs[col] = np.asarray(adata.obs[col], dtype=np.float64) + delta
    for info in adata.uns.get("spatial", {}).values():
        if not isinstance(info, Mapping):
            continue
        for key in ("geometries", "nucleus_geometries"):
            gdf = info.get(key)
            if isinstance(gdf, gpd.GeoDataFrame) and len(gdf):
                gdf["geometry"] = gdf.geometry.translate(xoff=dx, yoff=dy)


def concat_atera_rois(adatas: Union[Sequence[ad.AnnData], Mapping[str, ad.AnnData]],
                      offsets=None, *, names: Optional[Sequence[str]] = None,
                      join: str = "inner", label: str = "roi",
                      assume_aligned: bool = False,
                      index_unique: Optional[str] = None,
                      **concat_kwargs) -> ad.AnnData:
    """Concatenate ROIs from the same customer.

    Parameters
    ----------
    adatas
        List of AnnData objects, or a mapping ``{roi_name: adata}``.
    offsets
        ``{roi_name: (dx, dy)}`` in µm, or a sequence matching ``adatas``.
        Required unless ``assume_aligned=True`` (ROIs already share a
        coordinate system).
    names
        ROI names used for ``obs[label]``.  Defaults to the mapping keys or
        ``roi_0``, ``roi_1``, ...
    join
        ``'inner'`` (default) or ``'outer'`` for variable alignment.
    label
        ``obs`` column that records the ROI name.
    assume_aligned
        Set to True only when the ROIs come from the same slide and already
        share coordinates.  There is no reliable way to derive offsets from
        ``overview_scan_fov_locations.json`` (5–10 µm error), so the default
        is to require explicit offsets.
    index_unique
        If given (e.g. ``"-"``), append ``f"{index_unique}{roi_name}"`` to
        ``obs_names`` and geometry indices so that ROIs with identical
        barcodes stay distinguishable.  Atera barcodes can repeat across
        ROIs; pass this when concatenating them.

    Notes
    -----
    The input objects are translated **in place** (copying multi-GB AnnData
    objects is usually worse).
    """
    if isinstance(adatas, Mapping):
        if names is None:
            names = list(adatas.keys())
        adatas = list(adatas.values())
    else:
        adatas = list(adatas)
        if names is None:
            names = [f"roi_{i}" for i in range(len(adatas))]
    if len(names) != len(adatas):
        raise ValueError("`names` length must match `adatas` length.")

    if offsets is None:
        if not assume_aligned:
            raise ValueError(
                "`offsets` is required to place ROIs from different captures. "
                "Pass {roi_name: (dx, dy)} in µm, or set `assume_aligned=True` "
                "if the ROIs already share a coordinate system. "
                "Do not derive offsets from overview_scan_fov_locations.json: "
                "its accuracy is only 5–10 µm."
            )
        offsets = {name: (0.0, 0.0) for name in names}
    elif isinstance(offsets, Mapping):
        missing = [name for name in names if name not in offsets]
        if missing:
            raise ValueError(
                f"`offsets` is missing ROIs: {missing}. Pass offsets for every "
                "ROI (use (0, 0) for the reference ROI)."
            )
        offsets = {name: offsets[name] for name in names}
    else:
        offsets = list(offsets)
        if len(offsets) != len(names):
            raise ValueError(
                "`offsets` must provide one (dx, dy) per ROI, or be a mapping "
                "keyed by ROI name."
            )
        offsets = dict(zip(names, offsets))

    for name, adata in zip(names, adatas):
        dx, dy = offsets[name]
        if dx or dy:
            _translate_roi(adata, dx, dy)

    if index_unique is not None:
        for name, adata in zip(names, adatas):
            suffix = f"{index_unique}{name}"
            adata.obs_names = [f"{value}{suffix}" for value in adata.obs_names]
            for info in adata.uns.get("spatial", {}).values():
                if not isinstance(info, Mapping):
                    continue
                for key in ("geometries", "nucleus_geometries"):
                    gdf = info.get(key)
                    if isinstance(gdf, gpd.GeoDataFrame) and len(gdf):
                        gdf.index = [f"{value}{suffix}" for value in gdf.index]
    elif len({name for adata in adatas for name in adata.obs_names}) < sum(
        adata.n_obs for adata in adatas
    ):
        warnings.warn(
            "Duplicate `obs_names` across ROIs. Pass `index_unique='-'` to "
            "suffix obs names (and geometry indices) per ROI.",
            stacklevel=2,
        )

    merged = ad.concat(
        adatas, join=join, label=label, keys=list(names),
        uns_merge="first", **concat_kwargs
    )

    spatial = {}
    for adata in adatas:
        for key, value in adata.uns.get("spatial", {}).items():
            if key in spatial and spatial[key] is not value:
                warnings.warn(
                    f"Multiple ROIs share the library_id '{key}'; keeping the "
                    "first. Pass distinct `sample=` values to read_atera() to "
                    "avoid this.",
                    stacklevel=2,
                )
                continue
            spatial[key] = value
    if spatial:
        merged.uns["spatial"] = spatial
    merged.uns.setdefault("atera", {})["merged_rois"] = list(names)
    merged.uns["atera"]["offsets"] = {k: tuple(map(float, v)) for k, v in offsets.items()}
    return merged


# ---------------------------------------------------------------------------
# transcripts
# ---------------------------------------------------------------------------
def _transcript_gene_index(group) -> Tuple[list, dict]:
    gene_names = list(group.attrs["gene_names"])
    gene_index_map = group.attrs.get("gene_index_map")
    if gene_index_map is None:
        gene_index_map = {name: i for i, name in enumerate(gene_names)}
    return gene_names, dict(gene_index_map)


def _resolve_genes(genes, gene_names, gene_index_map) -> Optional[np.ndarray]:
    if genes is None:
        return None
    if isinstance(genes, str):
        genes = [genes]
    genes = list(genes)
    missing = [g for g in genes if g not in gene_index_map]
    if missing:
        preview = ", ".join(repr(g) for g in missing[:5])
        more = f" (+{len(missing) - 5} more)" if len(missing) > 5 else ""
        warnings.warn(
            f"{len(missing)} requested gene(s) are not in the Atera panel: "
            f"{preview}{more}. Available names are gene symbols, not Ensembl IDs.",
            stacklevel=4,
        )
    found = sorted({int(gene_index_map[g]) for g in genes if g in gene_index_map})
    if not found:
        raise ValueError(
            "None of the requested genes are present in the Atera panel. "
            "Use `tcl.io.atera_info(path)` / `adata.var_names` to list genes."
        )
    return np.asarray(found, dtype=np.int64)


def _build_key_lookup(keys: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    keys = np.asarray(keys, dtype=np.int64)
    order = np.argsort(keys)
    return keys[order], order.astype(np.int64)


def _lookup_from_adata(adata: ad.AnnData):
    if "cell_id" in adata.obs.columns:
        return _build_key_lookup(adata.obs["cell_id"].to_numpy(np.int64))
    if {"cell_prefix", "dataset_suffix"}.issubset(adata.obs.columns):
        prefix = adata.obs["cell_prefix"].to_numpy(np.uint64)
        suffix = adata.obs["dataset_suffix"].to_numpy(np.uint64)
        keys = ((suffix << np.uint64(32)) | prefix).astype(np.int64)
        return _build_key_lookup(keys)
    return None


def _lookup_from_cells(cells_path: Path):
    with _zarr_group(cells_path) as group:
        cell_ids = np.asarray(group["cell_id"][:], dtype=np.uint32)
    keys = (
        (cell_ids[:, 1].astype(np.uint64) << np.uint64(32))
        | cell_ids[:, 0].astype(np.uint64)
    ).astype(np.int64)
    return _build_key_lookup(keys)


def _map_cell_keys(prefix: np.ndarray, suffix: np.ndarray, lookup) -> np.ndarray:
    """Map transcript ``(prefix, suffix)`` pairs to lookup rows (-1 = none)."""
    rows = np.full(len(prefix), -1, dtype=np.int64)
    if lookup is None:
        return rows
    keys, row_of_key = lookup
    if len(keys) == 0:
        return rows

    query = (
        (suffix.astype(np.uint64) << np.uint64(32))
        | prefix.astype(np.uint64)
    ).astype(np.int64)
    pos = np.searchsorted(keys, query)
    pos_clip = np.minimum(pos, len(keys) - 1)
    ok = (
        (pos < len(keys))
        & (keys[pos_clip] == query)
        & (prefix != _UNASSIGNED)
        & (suffix != _UNASSIGNED)
    )
    rows[ok] = row_of_key[pos_clip[ok]]
    return rows


def _merged_gene_ranges(gene_offset: np.ndarray, gene_indices: np.ndarray):
    ranges = sorted(
        (int(gene_offset[i, 0]), int(gene_offset[i, 1]))
        for i in gene_indices
        if int(gene_offset[i, 1]) > int(gene_offset[i, 0])
    )
    merged = []
    for start, end in ranges:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _iter_slices(start: int, end: int, chunksize: int):
    for s in range(start, end, chunksize):
        yield s, min(s + chunksize, end)


def _open_transcripts(path: Union[str, Path]):
    root = Path(path).expanduser()
    if root.is_file():
        tx_path = root
    else:
        tx_path = _member(root, _TX_ZIP, "transcripts.zarr")
    if tx_path is None:
        raise FileNotFoundError(
            f"No '{_TX_ZIP}' found in {root}. Official Atera bundles include a "
            "transcripts point cloud; preview bundles do not."
        )
    _require_zarr()
    store = None
    try:
        if tx_path.is_dir():
            group = zarr.open_group(str(tx_path), mode="r")
        else:
            store = _ZipStore(str(tx_path), mode="r")
            group = zarr.open_group(store, mode="r")
    except Exception:
        if store is not None:
            store.close()
        raise
    return store, group


def _iter_transcript_arrays(
    path: Union[str, Path],
    *,
    genes=None,
    region=None,
    min_qv: Optional[float] = 20.0,
    include_unassigned: bool = False,
    chunksize: Optional[int] = None,
    cell_lookup=None,
    verbose: bool = False,
):
    """Yield transcript arrays tile-by-tile / chunk-by-chunk."""
    if chunksize is not None and int(chunksize) <= 0:
        raise ValueError("`chunksize` must be a positive integer.")
    log = _log(verbose)
    store, group = _open_transcripts(path)
    try:
        gene_names, gene_index_map = _transcript_gene_index(group)
        gene_indices = _resolve_genes(genes, gene_names, gene_index_map)

        codeword_gene = np.asarray(
            group.attrs["codeword_gene_mapping"], dtype=np.int64
        )
        codeword_mask = None
        if gene_indices is not None:
            codeword_mask = np.isin(codeword_gene, gene_indices)

        if cell_lookup is None:
            root = Path(path).expanduser()
            bundle = root if root.is_dir() else root.parent
            cells_path = _member(bundle, _CELLS_ZIP, "cells.zarr")
            if cells_path is not None:
                cell_lookup = _lookup_from_cells(cells_path)

        parsed_region = _parse_region(region)
        grid_keys = list(group["grid"].attrs.get("grid_keys", []))
        skipped_negative = 0
        n_tiles = 0
        t0 = time.time()

        for tile_key in grid_keys:
            tile = group["grid"][tile_key]
            n_rows = int(tile["location"].shape[0])
            if n_rows == 0:
                continue
            try:
                negative_tile = any(int(part) < 0 for part in str(tile_key).split(","))
            except ValueError:
                negative_tile = False
            if negative_tile:
                # Atera writes a single placeholder row into tiles outside the
                # ROI (e.g. "-1,0").  It looks like a real transcript (Q=40)
                # but has an unassigned cell and invalid coordinates.
                skipped_negative += 1
                continue

            n_tiles += 1
            base_chunk = chunksize
            if base_chunk is None:
                try:
                    base_chunk = max(int(tile["location"].chunks[0]), 1_000_000)
                except (AttributeError, TypeError):
                    base_chunk = 2_000_000

            if gene_indices is not None and "gene_offset" in tile and (
                len(gene_indices) <= 0.5 * len(gene_names)
            ):
                gene_offset = np.asarray(tile["gene_offset"][:])
                slices = [
                    sl
                    for start, end in _merged_gene_ranges(gene_offset, gene_indices)
                    for sl in _iter_slices(start, end, base_chunk)
                ]
            else:
                slices = list(_iter_slices(0, n_rows, base_chunk))

            for start, end in slices:
                location = np.asarray(tile["location"][start:end], dtype=np.float32)
                qv = np.asarray(tile["quality_score"][start:end, 0], dtype=np.float32)
                cell_ids = np.asarray(tile["cell_id"][start:end], dtype=np.uint32)
                codeword = np.asarray(
                    tile["codeword_identity"][start:end, 0], dtype=np.int64
                )
                nucleus = np.asarray(
                    tile["overlaps_nucleus"][start:end, 0], dtype=bool
                )

                mask = codeword < len(codeword_gene)
                if min_qv is not None:
                    mask &= qv >= float(min_qv)
                if codeword_mask is not None:
                    safe_codeword = np.clip(codeword, 0, len(codeword_gene) - 1)
                    mask &= codeword_mask[safe_codeword]
                if not include_unassigned:
                    mask &= (cell_ids[:, 0] != _UNASSIGNED)
                    mask &= (cell_ids[:, 1] != _UNASSIGNED)
                if parsed_region is not None:
                    mask &= _region_mask(location[:, 0], location[:, 1],
                                         parsed_region)
                if not mask.any():
                    continue

                sel = np.flatnonzero(mask)
                prefix = cell_ids[sel, 0]
                suffix = cell_ids[sel, 1]
                rows = _map_cell_keys(prefix, suffix, cell_lookup)

                yield {
                    "x": location[sel, 0],
                    "y": location[sel, 1],
                    "z": location[sel, 2],
                    "qv": qv[sel],
                    "cell_prefix": prefix,
                    "cell_index": rows,
                    "gene_index": codeword_gene[codeword[sel]],
                    "codeword": codeword[sel],
                    "overlaps_nucleus": nucleus[sel],
                }

        if verbose:
            log(f"[Atera] Streamed {n_tiles} transcript tiles "
                f"in {time.time() - t0:.1f}s ({skipped_negative} background "
                f"tile(s) skipped)")
    finally:
        if store is not None:
            store.close()


def iter_atera_transcripts(
    path: Union[str, Path],
    *,
    genes: Optional[Sequence[str]] = None,
    region=None,
    min_qv: float = 20.0,
    include_unassigned: bool = False,
    adata: Optional[ad.AnnData] = None,
    chunksize: Optional[int] = None,
    verbose: bool = False,
) -> Iterator[pd.DataFrame]:
    """Stream Atera transcripts as DataFrames.

    Parameters
    ----------
    path
        ROI bundle directory (or the ``transcripts.zarr.zip`` file itself).
    genes
        Gene symbols to keep.  Uses the per-tile ``gene_offset`` index so only
        the requested gene ranges are read.
    region
        ``(xmin, ymin, xmax, ymax)`` in µm, a shapely geometry, or a mapping.
    min_qv
        Minimum Q-score.  Use ``None`` to keep all molecules.  The official
        cell matrix uses Q ≥ 20.
    include_unassigned
        Keep molecules with no cell assignment (``cell_index = -1``).
    adata
        Optional AnnData whose ``obs`` row order defines ``cell_index``.  Pass
        this after subsetting cells, otherwise ``cell_index`` refers to the
        full bundle (``cells.zarr``) order.
    chunksize
        Molecules per read.  ``None`` picks the zarr chunk size (≥1M).
    verbose
        Print tile progress.

    Yields
    ------
    DataFrame
        Columns: ``x``, ``y``, ``z`` (µm), ``qv``, ``cell_index`` (row in
        ``adata.obs``, -1 = unassigned), ``cell_prefix``, ``gene_index``,
        ``gene`` (categorical), ``codeword``, ``overlaps_nucleus``.
        ``overlaps_nucleus`` is a nuclear-localisation flag, **not** an
        unspliced count.

    Examples
    --------
    >>> for chunk in tcl.io.iter_atera_transcripts(
    ...         roi, genes=["NPHS1", "UMOD"], min_qv=20):
    ...     print(chunk.head())

    Notes
    -----
    The returned ``gene`` column is a ``pandas.Categorical`` to keep memory
    low on billion-row scans.

    Background grid tiles (negative keys such as ``"-1,0"``) are skipped: they
    only contain placeholder rows outside the ROI.
    """
    gene_names = _transcript_gene_names(path)
    dtype = pd.CategoricalDtype(categories=gene_names)
    cell_lookup = _lookup_from_adata(adata) if adata is not None else None

    for arrays in _iter_transcript_arrays(
        path, genes=genes, region=region, min_qv=min_qv,
        include_unassigned=include_unassigned, chunksize=chunksize,
        cell_lookup=cell_lookup, verbose=verbose,
    ):
        yield pd.DataFrame({
            "x": arrays["x"],
            "y": arrays["y"],
            "z": arrays["z"],
            "qv": arrays["qv"],
            "cell_index": arrays["cell_index"],
            "cell_prefix": arrays["cell_prefix"],
            "gene_index": arrays["gene_index"].astype(np.int32),
            "gene": pd.Categorical.from_codes(
                arrays["gene_index"].astype(np.int32), dtype=dtype
            ),
            "codeword": arrays["codeword"].astype(np.uint32),
            "overlaps_nucleus": arrays["overlaps_nucleus"],
        })


def _transcript_gene_names(path: Union[str, Path]) -> list:
    """Read transcript gene names without keeping the store open."""
    store, group = _open_transcripts(path)
    try:
        return list(group.attrs["gene_names"])
    finally:
        if store is not None:
            store.close()


def read_atera_transcripts(
    path: Union[str, Path],
    *,
    genes: Optional[Sequence[str]] = None,
    region=None,
    min_qv: float = 20.0,
    include_unassigned: bool = False,
    adata: Optional[ad.AnnData] = None,
    max_rows: Optional[int] = None,
    force: bool = False,
    verbose: bool = True,
) -> pd.DataFrame:
    """Read Atera transcripts into a single DataFrame.

    This is only safe for small queries.  A full kidney ROI has ~1.19 billion
    molecules, so a query without ``genes`` or ``region`` is rejected unless
    ``max_rows`` or ``force=True`` is provided.

    Use :func:`iter_atera_transcripts` for large scans.
    """
    if genes is None and region is None and max_rows is None and not force:
        raise ValueError(
            "Refusing to load the entire Atera transcript cloud into memory. "
            "Pass `genes=` and/or `region=` for a targeted query, "
            "`max_rows=` to cap the result, or use "
            "`tcl.io.iter_atera_transcripts()` / "
            "`tcl.io.atera_transcript_qc()` for streaming."
        )
    frames = []
    n = 0
    for chunk in iter_atera_transcripts(
        path, genes=genes, region=region, min_qv=min_qv,
        include_unassigned=include_unassigned, adata=adata, verbose=verbose,
    ):
        if max_rows is not None and n + len(chunk) > max_rows:
            chunk = chunk.iloc[: max_rows - n]
        frames.append(chunk)
        n += len(chunk)
        if max_rows is not None and n >= max_rows:
            break
    if not frames:
        return pd.DataFrame({
            "x": pd.Series(dtype=np.float32),
            "y": pd.Series(dtype=np.float32),
            "z": pd.Series(dtype=np.float32),
            "qv": pd.Series(dtype=np.float32),
            "cell_index": pd.Series(dtype=np.int64),
            "cell_prefix": pd.Series(dtype=np.uint32),
            "gene_index": pd.Series(dtype=np.int32),
            "gene": pd.Series(
                dtype=pd.CategoricalDtype(categories=_transcript_gene_names(path))
            ),
            "codeword": pd.Series(dtype=np.uint32),
            "overlaps_nucleus": pd.Series(dtype=bool),
        })
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------------
def atera_cell_qc(adata: ad.AnnData, compute_expression_genes: bool = True
                  ) -> pd.DataFrame:
    """Compute cell-level QC for an Atera AnnData object.

    Atera's ``obs['total_counts']`` includes control probes/codewords; use
    ``obs['transcript_counts']`` for gene-only counts.

    Returns
    -------
    DataFrame
        Indexed by ``adata.obs_names`` with ``total_counts``,
        ``transcript_counts``, ``n_genes`` (matrix nnz per cell, includes
        control features), ``n_genes_expression`` (Gene Expression features
        only, when ``compute_expression_genes=True``), ``cell_area``,
        ``nucleus_area``, ``nuclear_area_ratio`` and control fractions.
    """
    obs = adata.obs
    qc = pd.DataFrame(index=adata.obs_names)

    if "total_counts" in obs.columns:
        qc["total_counts"] = pd.to_numeric(obs["total_counts"], errors="coerce")
    else:
        qc["total_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    if "transcript_counts" in obs.columns:
        qc["transcript_counts"] = pd.to_numeric(
            obs["transcript_counts"], errors="coerce"
        )

    X_matrix = adata.X if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    X_matrix = X_matrix.tocsr()
    qc["n_genes"] = np.diff(X_matrix.indptr).astype(np.int64)

    if compute_expression_genes and "is_gene_expression" in adata.var.columns:
        col_mask = adata.var["is_gene_expression"].to_numpy(dtype=bool)
        X = X_matrix
        counts = np.zeros(adata.n_obs, dtype=np.int64)
        for i in range(adata.n_obs):
            start, end = X.indptr[i], X.indptr[i + 1]
            counts[i] = int(col_mask[X.indices[start:end]].sum())
        qc["n_genes_expression"] = counts

    for col in ("cell_area", "nucleus_area"):
        if col in obs.columns:
            qc[col] = pd.to_numeric(obs[col], errors="coerce")
    if {"cell_area", "nucleus_area"}.issubset(qc.columns):
        qc["nuclear_area_ratio"] = qc["nucleus_area"] / qc["cell_area"].replace(0, np.nan)

    control_cols = [
        c for c in ("control_probe_counts", "control_codeword_counts")
        if c in obs.columns
    ]
    if control_cols:
        control = obs[control_cols].apply(pd.to_numeric, errors="coerce").sum(axis=1)
        qc["control_counts"] = control
        qc["control_fraction"] = control / qc["total_counts"].replace(0, np.nan)

    return qc


def _transcript_bounds(path: Path, adata: Optional[ad.AnnData], region) -> Optional[tuple]:
    parsed = _parse_region(region)
    if parsed is not None:
        return parsed[0]
    xy = None
    if adata is not None and "spatial" in adata.obsm:
        xy = np.asarray(adata.obsm["spatial"])[:, :2]
    else:
        cells_path = _member(path, _CELLS_ZIP, "cells.zarr")
        if cells_path is not None:
            with _zarr_group(cells_path) as group:
                summary = np.asarray(group["cell_summary"][:])
                names = list(group["cell_summary"].attrs.get("column_names", []))
            if {"cell_centroid_x", "cell_centroid_y"}.issubset(names):
                ix = names.index("cell_centroid_x")
                iy = names.index("cell_centroid_y")
                xy = summary[:, [ix, iy]]
    if xy is None or not len(xy):
        return None
    xmin, ymin = float(np.min(xy[:, 0])), float(np.min(xy[:, 1]))
    xmax, ymax = float(np.max(xy[:, 0])), float(np.max(xy[:, 1]))
    if not np.isfinite([xmin, ymin, xmax, ymax]).all() or xmax <= xmin or ymax <= ymin:
        return None
    mx, my = 0.05 * (xmax - xmin), 0.05 * (ymax - ymin)
    return (xmin - mx, ymin - my, xmax + mx, ymax + my)


def atera_transcript_qc(
    path: Union[str, Path],
    adata: Optional[ad.AnnData] = None,
    *,
    min_qv: float = 20.0,
    bins: int = 100,
    region=None,
    verbose: bool = True,
) -> dict:
    """Stream the Atera point cloud and compute Q-score / nuclear / density QC.

    Results are written back to ``adata`` when provided:

    * ``obs['tx_q20_counts']``, ``obs['tx_nuclear_counts']``,
      ``obs['nuclear_fraction']``
    * ``uns['atera']['transcript_qc']`` (summary dict)
    * ``uns['atera']['qscore']`` (histogram)
    * ``uns['atera']['density']`` (``H``, ``xedges``, ``yedges``)
    * ``uns['atera']['gene_qc']`` (per-gene totals and nuclear fractions)

    This makes :func:`trackcell.pl.atera_nuclear_fraction` and
    :func:`trackcell.pl.atera_transcript_density` re-callable, since they only
    read ``adata``.

    Returns
    -------
    dict
        Summary with transcript counts, Q≥min_qv fractions, assignment rate,
        nuclear fraction and runtime.
    """
    log = _log(verbose)
    root = Path(path).expanduser()
    bundle = root if root.is_dir() else root.parent
    t0 = time.time()

    lookup = None
    n_cells = 0
    if adata is not None:
        lookup = _lookup_from_adata(adata)
        n_cells = adata.n_obs
    else:
        cells_path = _member(bundle, _CELLS_ZIP, "cells.zarr")
        if cells_path is not None:
            lookup = _lookup_from_cells(cells_path)
            n_cells = len(lookup[0])

    store, group = _open_transcripts(root)
    try:
        gene_names = list(group.attrs["gene_names"])
        n_genes = len(gene_names)
        number_rnas = int(group.attrs.get("number_rnas", 0))
    finally:
        if store is not None:
            store.close()

    bounds = _transcript_bounds(bundle, adata, region)
    if bounds is None:
        bounds = (0.0, 0.0, 1.0, 1.0)
        warnings.warn(
            "Could not determine a spatial bounding box for the density map; "
            "pass `region=` or `adata=` explicitly."
        )
    xmin, ymin, xmax, ymax = bounds
    xe = np.linspace(xmin, xmax, bins + 1)
    ye = np.linspace(ymin, ymax, bins + 1)
    density = np.zeros((bins, bins), dtype=np.int64)

    # Atera Q-scores are reported in the 10–40 range; keep head-room so no
    # molecule is silently dropped from the histogram.
    q_bins = np.arange(0, 60.5, 0.5)
    q_total = np.zeros(len(q_bins) - 1, dtype=np.int64)
    q_hq = np.zeros(len(q_bins) - 1, dtype=np.int64)

    cell_total = np.zeros(n_cells, dtype=np.int64)
    cell_nuclear = np.zeros(n_cells, dtype=np.int64)
    gene_total = np.zeros(n_genes, dtype=np.int64)
    gene_nuclear = np.zeros(n_genes, dtype=np.int64)

    n_read = 0
    n_q20 = 0
    n_assigned = 0
    n_tiles_done = 0

    for arrays in _iter_transcript_arrays(
        root, min_qv=None, include_unassigned=True, region=region,
        cell_lookup=lookup, verbose=False,
    ):
        qv = arrays["qv"]
        hq = qv >= float(min_qv)
        q_total += np.histogram(qv, bins=q_bins)[0]
        q_hq += np.histogram(qv[hq], bins=q_bins)[0]
        n_read += len(qv)
        n_q20 += int(hq.sum())

        density += np.histogram2d(
            arrays["x"][hq], arrays["y"][hq], bins=[xe, ye]
        )[0].astype(np.int64)

        rows = arrays["cell_index"]
        assigned = hq & (rows >= 0)
        if assigned.any():
            cell_total += np.bincount(rows[assigned], minlength=n_cells)
            nuc = assigned & arrays["overlaps_nucleus"]
            cell_nuclear += np.bincount(rows[nuc], minlength=n_cells)
            n_assigned += int(assigned.sum())

        gene_idx = arrays["gene_index"].astype(np.int64)
        gene_total += np.bincount(gene_idx[hq], minlength=n_genes)
        gene_nuclear += np.bincount(
            gene_idx[hq & arrays["overlaps_nucleus"]], minlength=n_genes
        )

        n_tiles_done += 1
        if verbose and n_tiles_done % 10 == 0:
            log(f"[Atera] QC: {n_read:,} molecules "
                f"({time.time() - t0:.0f}s)")

    elapsed = time.time() - t0
    overall_nuclear = (
        float(cell_nuclear.sum() / cell_total.sum()) if cell_total.sum() else float("nan")
    )
    summary = {
        "number_rnas": number_rnas,
        "read": int(n_read),
        "min_qv": float(min_qv),
        "q20": int(n_q20),
        "q20_fraction": float(n_q20 / n_read) if n_read else float("nan"),
        "assigned_q20": int(n_assigned),
        "assigned_fraction_q20": float(n_assigned / n_q20) if n_q20 else float("nan"),
        "overall_nuclear_fraction": overall_nuclear,
        "cells_with_transcripts": int((cell_total > 0).sum()),
        "n_cells": int(n_cells),
        "n_genes": int(n_genes),
        "seconds": float(elapsed),
    }

    if adata is not None:
        adata.obs["tx_q20_counts"] = cell_total.astype(np.int32)
        adata.obs["tx_nuclear_counts"] = cell_nuclear.astype(np.int32)
        adata.obs["nuclear_fraction"] = np.divide(
            cell_nuclear, cell_total,
            out=np.full(n_cells, np.nan, dtype=np.float64),
            where=cell_total > 0,
        )
        atera_uns = adata.uns.setdefault("atera", {})
        atera_uns["transcript_qc"] = summary
        atera_uns["qscore"] = {
            "centers": (q_bins[:-1] + q_bins[1:]) / 2.0,
            "total": q_total,
            "hq": q_hq,
            "min_qv": float(min_qv),
        }
        atera_uns["density"] = {
            "H": density,
            "xedges": xe,
            "yedges": ye,
            "unit": _UM,
        }
        atera_uns["gene_qc"] = pd.DataFrame({
            "gene": gene_names,
            "total": gene_total,
            "nuclear": gene_nuclear,
            "nuclear_fraction": np.divide(
                gene_nuclear, gene_total,
                out=np.full(n_genes, np.nan, dtype=np.float64),
                where=gene_total > 0,
            ),
        })

    if verbose:
        log(f"[Atera] QC done in {elapsed:.0f}s: "
            f"{n_read:,} molecules, Q≥{min_qv:g} {summary['q20_fraction'] * 100:.1f}%, "
            f"assigned {summary['assigned_fraction_q20'] * 100:.1f}%, "
            f"nuclear {overall_nuclear * 100:.1f}%")
    return summary
