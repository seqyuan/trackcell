"""h5ad serialization helpers for AnnData objects with polygon geometries.

``read_hd_cellseg`` / ``read_xenium_cellseg`` / ``read_atera`` store Shapely
objects inside ``uns['spatial'][sample]['geometries']`` (a GeoDataFrame).
HDF5 cannot serialize those objects, so ``adata.write_h5ad()`` raises
``TypeError: Can't implicitly convert non-string objects to strings``.

Use :func:`save_h5ad` instead of ``adata.write_h5ad``.  Geometries are
temporarily converted to WKT strings for the write and restored in memory
afterwards.  Read the file back with :func:`restore_geometries`.

Two round-trip properties are load-bearing and covered by tests:

* **Every attribute column survives.** Cell type, area, classification and any
  other non-geometry column is preserved; only the active ``geometry`` column
  is stringified.
* **Missing polygons survive.** Cells without a segmentation outline (Xenium
  edge cells, Atera debris) are written as an empty string and read back as
  ``None``; anndata rejects ``None`` inside an object column, so this cannot be
  written directly.
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Optional, Union

import geopandas as gpd
import pandas as pd
from anndata import AnnData
from shapely import wkt

from .._logging import get_logger

logger = get_logger(__name__)

GEOMETRY_KEYS = ("geometries", "nucleus_geometries")

# Sidecar key holding a WKT-encoded CRS next to its geometry table. HDF5 can
# store the string; a ``GeoDataFrame.crs`` object cannot survive ``write_h5ad``.
CRS_SUFFIX = "__crs"
# ``uns`` columns cannot hold ``None`` in an object column (anndata raises
# ``TypeError: Can't implicitly convert non-string objects to strings``), so a
# missing polygon is stored as the empty string and mapped back to ``None``.
EMPTY_GEOMETRY = ""


def _geometry_to_wkt(geom: Any) -> str:
    """Encode one shapely geometry (or a missing value) as a WKT string."""
    if geom is None:
        return EMPTY_GEOMETRY
    try:
        if pd.isna(geom):
            return EMPTY_GEOMETRY
    except (TypeError, ValueError):
        # shapely geometries raise on elementwise isna; that means "present".
        pass
    text = getattr(geom, "wkt", None)
    if not isinstance(text, str) or not text:
        return EMPTY_GEOMETRY
    return text


def _wkt_to_geometry(value: Any, *, where: str = "") -> Any:
    """Decode a WKT string back to a shapely geometry, or ``None`` if absent."""
    if not isinstance(value, str) or not value.strip():
        return None
    try:
        return wkt.loads(value)
    except Exception as exc:  # pragma: no cover - corrupt file
        warnings.warn(
            f"Could not parse WKT geometry{(' in ' + where) if where else ''}: "
            f"{exc}",
            RuntimeWarning,
            stacklevel=3,
        )
        return None


def _to_wkt_frame(obj: Any) -> Any:
    """Return a WKT-backed DataFrame for a GeoDataFrame/GeoSeries, else ``obj``.

    Every non-geometry column (cell type, area, classification, …) is preserved
    — only the active ``geometry`` column is stringified. Dropping the attribute
    columns here used to make them disappear silently on a save/load round trip.
    """
    if isinstance(obj, gpd.GeoSeries):
        obj = gpd.GeoDataFrame({"geometry": obj}, geometry="geometry")
    if not isinstance(obj, gpd.GeoDataFrame):
        return obj

    columns = list(obj.columns)
    frame = pd.DataFrame(obj.drop(columns="geometry"), copy=True)
    frame["geometry"] = [_geometry_to_wkt(g) for g in obj.geometry]
    return frame[columns]


def _from_wkt_frame(frame: pd.DataFrame, *, where: str = "") -> gpd.GeoDataFrame:
    """Rebuild a GeoDataFrame from a WKT-string DataFrame, keeping all columns."""
    columns = list(frame.columns)
    geoms = gpd.GeoSeries(
        [_wkt_to_geometry(v, where=where) for v in frame["geometry"]],
        index=frame.index,
    )
    gdf = gpd.GeoDataFrame(
        pd.DataFrame(frame.drop(columns="geometry"), copy=True),
        geometry=geoms,
        index=frame.index,
    )
    return gdf[columns]


def _spatial_samples(adata: AnnData, sample: Optional[str]) -> list[str]:
    spatial = adata.uns.get("spatial")
    if not isinstance(spatial, Mapping):
        return []
    if sample is not None:
        return [sample] if sample in spatial else []
    return list(spatial.keys())


def flatten_geometries(
    adata: AnnData,
    sample: Optional[str] = None,
    *,
    inplace: bool = False,
) -> AnnData:
    """Convert polygon GeoDataFrames to WKT-string DataFrames for h5ad output.

    Parameters
    ----------
    adata
        AnnData whose ``uns['spatial'][sample]`` may contain
        ``geometries`` / ``nucleus_geometries`` GeoDataFrames.
    sample
        Restrict to one library id.  By default all samples are processed.
    inplace
        If ``True``, modify ``adata`` and return it; otherwise return a copy.

    Returns
    -------
    AnnData
        Object with WKT-backed geometry tables (safe for ``write_h5ad``).
    """
    out = adata if inplace else adata.copy()
    for sid in _spatial_samples(out, sample):
        info = out.uns["spatial"][sid]
        if not isinstance(info, dict):
            continue
        for key in GEOMETRY_KEYS:
            obj = info.get(key)
            if not isinstance(obj, (gpd.GeoSeries, gpd.GeoDataFrame)):
                continue
            info[key] = _to_wkt_frame(obj)
            crs = getattr(obj, "crs", None)
            if crs is not None:
                info[f"{key}{CRS_SUFFIX}"] = crs.to_wkt()
    return out


def save_h5ad(
    adata: AnnData,
    filename: Union[str, Path],
    *,
    compression: Optional[str] = None,
    **kwargs: Any,
) -> Path:
    """Write an AnnData to h5ad while handling GeoDataFrame geometries.

    This is a drop-in replacement for ``adata.write_h5ad()`` for objects
    produced by TrackCell's spatial readers.  Polygon tables are converted to
    WKT strings for the duration of the write and restored in memory
    afterwards, so the caller's ``adata`` is left unchanged.

    Parameters
    ----------
    adata
        AnnData to write.
    filename
        Output path.
    compression
        Passed to :meth:`anndata.AnnData.write_h5ad` (e.g. ``"gzip"``).
    **kwargs
        Additional keyword arguments for :meth:`anndata.AnnData.write_h5ad`.

    Returns
    -------
    pathlib.Path
        The written path.

    Examples
    --------
    >>> import trackcell as tcl
    >>> adata = tcl.io.read_xenium_cellseg("/path/to/xenium")
    >>> tcl.io.save_h5ad(adata, "xenium.h5ad")
    >>> # adata still holds GeoDataFrames in memory
    >>> adata = tcl.io.restore_geometries(__import__("anndata").read_h5ad("xenium.h5ad"))
    """
    path = Path(filename)
    backup: dict[tuple[str, str], Any] = {}
    crs_added: list[tuple[str, str]] = []
    crs_backup: dict[tuple[str, str], Any] = {}

    for sid in _spatial_samples(adata, None):
        info = adata.uns["spatial"][sid]
        if not isinstance(info, dict):
            continue
        for key in GEOMETRY_KEYS:
            obj = info.get(key)
            if isinstance(obj, gpd.GeoSeries) or isinstance(obj, gpd.GeoDataFrame):
                backup[(sid, key)] = obj
                info[key] = _to_wkt_frame(obj)
                crs = getattr(obj, "crs", None)
                if crs is not None:
                    crs_key = f"{key}{CRS_SUFFIX}"
                    # Do not destroy a sidecar the caller already had: the
                    # temporary value is restored below, not popped.
                    if crs_key in info and (sid, crs_key) not in crs_backup:
                        crs_backup[(sid, crs_key)] = info[crs_key]
                    info[crs_key] = crs.to_wkt()
                    crs_added.append((sid, crs_key))

    try:
        adata.write_h5ad(path, compression=compression, **kwargs)
    finally:
        for (sid, key), obj in backup.items():
            adata.uns["spatial"][sid][key] = obj
        for (sid, crs_key), value in crs_backup.items():
            adata.uns["spatial"][sid][crs_key] = value
        for sid, crs_key in crs_added:
            if (sid, crs_key) not in crs_backup:
                adata.uns["spatial"][sid].pop(crs_key, None)

    return path


__all__ = [
    "CRS_SUFFIX",
    "EMPTY_GEOMETRY",
    "GEOMETRY_KEYS",
    "flatten_geometries",
    "restore_geometries",
    "save_h5ad",
]


def restore_geometries(adata: AnnData, sample: Optional[str] = None) -> AnnData:
    """Restore GeoDataFrames from the WKT tables written by :func:`save_h5ad`.

    Inverse of :func:`flatten_geometries`: parses ``geometry`` back into shapely
    objects, keeps every attribute column, and re-attaches the CRS recorded in
    the ``<key>__crs`` sidecar.

    Parameters
    ----------
    adata
        AnnData read back from an h5ad written with :func:`save_h5ad`.
    sample
        Restrict to one library id.  By default all samples are processed.

    Returns
    -------
    AnnData
        The same object, mutated in place, with GeoDataFrames restored.

    Examples
    --------
    >>> import anndata as ad
    >>> import trackcell as tcl
    >>> tcl.io.save_h5ad(adata, "xenium.h5ad")
    >>> adata = tcl.io.restore_geometries(ad.read_h5ad("xenium.h5ad"))
    """
    for sid in _spatial_samples(adata, sample):
        info = adata.uns["spatial"][sid]
        if not isinstance(info, dict):
            continue
        for key in GEOMETRY_KEYS:
            frame = info.get(key)
            if not isinstance(frame, pd.DataFrame) or "geometry" not in frame:
                continue
            gdf = _from_wkt_frame(frame, where=f"uns['spatial']['{sid}']['{key}']")
            crs_wkt = info.pop(f"{key}{CRS_SUFFIX}", None)
            if isinstance(crs_wkt, str) and crs_wkt:
                try:
                    gdf = gdf.set_crs(crs_wkt, allow_override=True)
                except Exception:  # pragma: no cover - unknown CRS string
                    pass
            info[key] = gdf
            logger.info("Restored %s %s for sample '%s'", len(gdf), key, sid)
    return adata
