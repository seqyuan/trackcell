"""h5ad serialization helpers for AnnData objects with polygon geometries.

``read_hd_cellseg`` / ``read_xenium_cellseg`` / ``read_atera`` store Shapely
objects inside ``uns['spatial'][sample]['geometries']`` (a GeoDataFrame).
HDF5 cannot serialize those objects, so ``adata.write_h5ad()`` raises
``TypeError: Can't implicitly convert non-string objects to strings``.

Use :func:`save_h5ad` instead of ``adata.write_h5ad``.  Geometries are
temporarily converted to WKT strings for the write and restored in memory
afterwards.  Read the file back with :func:`~trackcell.io.restore_geometries`.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any, Optional, Union

import geopandas as gpd
import pandas as pd
from anndata import AnnData

GEOMETRY_KEYS = ("geometries", "nucleus_geometries")


def _to_wkt_frame(obj: Any) -> Any:
    """Return a WKT-backed DataFrame for a GeoDataFrame/GeoSeries, else ``obj``."""
    if isinstance(obj, gpd.GeoSeries):
        obj = gpd.GeoDataFrame({"geometry": obj}, geometry="geometry")
    if not isinstance(obj, gpd.GeoDataFrame):
        return obj

    values = []
    for geom in obj.geometry:
        if geom is None:
            values.append(None)
            continue
        try:
            if pd.isna(geom):
                values.append(None)
                continue
        except (TypeError, ValueError):
            pass
        values.append(geom.wkt)
    return pd.DataFrame({"geometry": values}, index=obj.index)


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
        if not isinstance(info, Mapping):
            continue
        for key in GEOMETRY_KEYS:
            if key in info:
                info[key] = _to_wkt_frame(info[key])
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

    for sid in _spatial_samples(adata, None):
        info = adata.uns["spatial"][sid]
        if not isinstance(info, Mapping):
            continue
        for key in GEOMETRY_KEYS:
            obj = info.get(key)
            if isinstance(obj, gpd.GeoSeries) or isinstance(obj, gpd.GeoDataFrame):
                backup[(sid, key)] = obj
                info[key] = _to_wkt_frame(obj)

    try:
        adata.write_h5ad(path, compression=compression, **kwargs)
    finally:
        for (sid, key), obj in backup.items():
            adata.uns["spatial"][sid][key] = obj

    return path


__all__ = ["GEOMETRY_KEYS", "flatten_geometries", "save_h5ad"]
