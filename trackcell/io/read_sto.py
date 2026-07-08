"""
STOmics (Stereo-seq) data reading functions for TrackCell package.

Supports reading STOmics GEF/GEM cellbin and squarebin data formats
and converting them to trackcell-compatible AnnData objects with
spatial coordinates and cell polygon geometries.

GEF format reference:
    https://stereopy.readthedocs.io/en/latest/index.html
"""

import scanpy as sc
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Polygon
from shapely import wkt
from scipy.sparse import csr_matrix
from pathlib import Path
from typing import Optional, Union, List
try:
    import h5py
except ImportError:
    h5py = None
import time
import warnings
import glob
import json


# Sentinel value used to pad cellBorder arrays in GEF cellbin files.
# When a cell has fewer than 32 border vertices, remaining slots are
# filled with int16 max (32767).
_BORDER_SENTINEL = np.int16(32767)


def _extract_polygon_from_relative_offsets(
    border_offsets: np.ndarray,
    cell_x: int,
    cell_y: int,
) -> Optional[Polygon]:
    """
    Convert relative cell border offsets to a Shapely Polygon.

    cellBorder stores vertices as (dx, dy) offsets relative to the cell's
    centroid position (cell_x, cell_y).  Trailing sentinel slots (32767, 32767)
    indicate unused vertices and must be stripped.

    Parameters
    ----------
    border_offsets : np.ndarray
        Array of shape (32, 2) with relative dx, dy offsets.
    cell_x : int
        Absolute x-coordinate of the cell centroid.
    cell_y : int
        Absolute y-coordinate of the cell centroid.

    Returns
    -------
    Polygon or None
        Shapely Polygon in image pixel coordinates, or None if the polygon
        has fewer than 3 valid vertices.
    """
    # Strip trailing sentinel entries
    valid_mask = (border_offsets[:, 0] != _BORDER_SENTINEL) | (border_offsets[:, 1] != _BORDER_SENTINEL)
    valid_offsets = border_offsets[valid_mask]

    if len(valid_offsets) < 3:
        # Not enough points for a polygon; use a small square around centroid.
        half = 5
        return Polygon([
            (cell_x - half, cell_y - half),
            (cell_x + half, cell_y - half),
            (cell_x + half, cell_y + half),
            (cell_x - half, cell_y + half),
        ])

    # Convert offsets + centroid → absolute coordinates
    abs_coords = valid_offsets.astype(np.float64)
    abs_coords[:, 0] += cell_x
    abs_coords[:, 1] += cell_y

    poly = Polygon(abs_coords)
    if not poly.is_valid or poly.is_empty:
        # Fallback: convex hull of points → buffer to ensure polygon
        from shapely.geometry import MultiPoint
        hull = MultiPoint(abs_coords).convex_hull
        if hasattr(hull, 'buffer'):
            return hull.buffer(1.0)
        return None

    return poly


def _build_expression_matrix_from_cell_exp(
    cell: np.ndarray,
    cell_exp: np.ndarray,
    n_genes: int,
    n_cells: int,
) -> csr_matrix:
    """
    Build a CSR expression matrix from the GEF cellExp offset-based storage.

    Parameters
    ----------
    cell : np.ndarray
        Structured array from cellBin/cell with fields including 'offset'.
    cell_exp : np.ndarray
        Structured array from cellBin/cellExp with fields 'geneID' and 'count'.
    n_genes : int
        Total number of genes.
    n_cells : int
        Total number of cells.

    Returns
    -------
    csr_matrix
        Shape (n_cells, n_genes) with dtype np.float32.
    """
    cell_offsets = cell['offset'].astype(np.int64)
    total_entries = len(cell_exp)

    # Determine end offsets for each cell
    end_offsets = np.empty(n_cells, dtype=np.int64)
    end_offsets[:-1] = cell_offsets[1:]
    end_offsets[-1] = total_entries

    # Count non-zero entries per cell
    nz_per_cell = end_offsets - cell_offsets
    total_nz = int(nz_per_cell.sum())

    # Pre-allocate arrays
    rows = np.empty(total_nz, dtype=np.int32)
    cols = np.empty(total_nz, dtype=np.int32)
    data = np.empty(total_nz, dtype=np.float32)
    indptr = np.zeros(n_cells + 1, dtype=np.int32)

    pos = 0
    for i in range(n_cells):
        start = cell_offsets[i]
        end = end_offsets[i]
        n = end - start
        if n > 0:
            indptr[i + 1] = indptr[i] + n
            rows[pos:pos + n] = i
            cols[pos:pos + n] = cell_exp['geneID'][start:end]
            data[pos:pos + n] = cell_exp['count'][start:end].astype(np.float32)
            pos += n
        else:
            indptr[i + 1] = indptr[i]

    return csr_matrix((data, cols, indptr), shape=(n_cells, n_genes), dtype=np.float32)


def _find_register_image(register_dir: Path) -> Optional[Path]:
    """Auto-discover the registered ssDNA image in a register directory.

    Looks for ``ssDNA_*_regist.tif`` first, then falls back to
    ``ssDNA_fov_stitched_transformed.tif`` or any ``*.tif``.
    """
    patterns = [
        '*_regist.tif',
        '*_regist.TIF',
        'ssDNA_fov_stitched_transformed.tif',
        '*.tif',
        '*.TIF',
    ]
    for pat in patterns:
        candidates = list(register_dir.glob(pat))
        if candidates:
            # Prefer files without 'mask' or 'tissue_cut' in name
            for c in candidates:
                if 'mask' not in c.name.lower() and 'tissue_cut' not in c.name.lower():
                    return c
            return candidates[0]
    return None


def _load_tissue_image(
    image_path: Optional[Union[str, Path]],
    gef_dir: Path,
) -> Optional[np.ndarray]:
    """Load a tissue image for STOmics data.

    Parameters
    ----------
    image_path : str, Path, or None
        Explicit image path, or None to auto-discover.
    gef_dir : Path
        The directory containing the GEF file (or its parent),
        used to search for the register directory.

    Returns
    -------
    np.ndarray or None
        The loaded image array, or None if not found.
    """
    if image_path is not None:
        image_path = Path(image_path).resolve()
        if image_path.is_dir():
            found = _find_register_image(image_path)
            if found is not None:
                image_path = found
            else:
                warnings.warn(f"No tissue image found in {image_path}")
                return None
        if image_path.suffix.lower() in ('.tif', '.tiff'):
            import imageio.v3 as iio
            img = iio.imread(image_path)
            print(f"[STO] Loaded tissue image: {image_path.name}"
                  f" ({img.shape[0]}×{img.shape[1]})")
            return img
        elif image_path.suffix.lower() in ('.rpi',):
            warnings.warn(
                "RPI pyramid images are not supported yet. "
                "Use the corresponding .tif file instead."
            )
            return None
        else:
            warnings.warn(f"Unsupported image format: {image_path.suffix}")
            return None

    # Auto-discover from the register directory
    # GEF is typically in .../C57_7/04.tissuecut/ or .../C57_7/041.cellcut/
    # Register dir is .../C57_7/03.register/
    # First check sibling directories
    parent = gef_dir.parent
    for candidate_name in ('03.register', 'register'):
        reg_dir = parent / candidate_name
        if reg_dir.is_dir():
            found = _find_register_image(reg_dir)
            if found is not None:
                import imageio.v3 as iio
                img = iio.imread(found)
                print(f"[STO] Auto-loaded tissue image: {found.name}"
                      f" ({img.shape[0]}×{img.shape[1]})")
                return img

    # Also check two levels up (some pipelines nest deeper)
    grandparent = parent.parent
    for candidate_name in ('03.register', 'register'):
        reg_dir = grandparent / candidate_name
        if reg_dir.is_dir():
            found = _find_register_image(reg_dir)
            if found is not None:
                import imageio.v3 as iio
                img = iio.imread(found)
                print(f"[STO] Auto-loaded tissue image: {found.name}"
                      f" ({img.shape[0]}×{img.shape[1]})")
                return img

    return None


def _find_cellbin_gef(dirpath: Path) -> Path:
    """Auto-discover cellbin GEF file in a directory.

    Looks for files matching ``*cellbin*.gef`` (case-insensitive on macOS).
    """
    candidates = list(dirpath.glob('*cellbin*.gef')) + list(dirpath.glob('*cellbin*.GEF'))
    if not candidates:
        # Also try any .gef file (some STOmics outputs don't have 'cellbin' in name)
        candidates = list(dirpath.glob('*.gef')) + list(dirpath.glob('*.GEF'))
    if not candidates:
        raise FileNotFoundError(
            f"No cellbin GEF file found in {dirpath}. "
            f"Expected a file matching '*cellbin*.gef' or '*.gef'."
        )
    if len(candidates) > 1:
        warnings.warn(
            f"Multiple GEF files found in {dirpath}, using {candidates[0].name}"
        )
    return candidates[0]


def _find_tissue_gef(dirpath: Path) -> Path:
    """Auto-discover tissue/squarebin GEF file in a directory.

    Looks for files matching ``*tissue*.gef`` (case-insensitive on macOS).
    Falls back to any .gef file that is NOT a cellbin file.
    """
    candidates = list(dirpath.glob('*tissue*.gef')) + list(dirpath.glob('*tissue*.GEF'))
    if not candidates:
        # Any .gef that doesn't contain 'cellbin'
        all_gef = list(dirpath.glob('*.gef')) + list(dirpath.glob('*.GEF'))
        candidates = [f for f in all_gef if 'cellbin' not in f.name.lower()]
    if not candidates:
        raise FileNotFoundError(
            f"No tissue/squarebin GEF file found in {dirpath}. "
            f"Expected a file matching '*tissue*.gef' or a .gef without 'cellbin'."
        )
    if len(candidates) > 1:
        warnings.warn(
            f"Multiple tissue GEF files found in {dirpath}, using {candidates[0].name}"
        )
    return candidates[0]


def read_sto_cellbin(
    gef_path: Union[str, Path],
    sample: Optional[str] = None,
    gene_name_index: bool = True,
    image_path: Optional[Union[str, Path]] = None,
) -> sc.AnnData:
    """
    Read STOmics cellbin GEF data and create an AnnData with cell polygons.

    Reads a STOmics cellbin GEF (HDF5) file containing cell segmentation
    boundaries, expression counts, and metadata.  Produces a trackcell-compatible
    AnnData with polygon geometries for each cell.

    Parameters
    ----------
    gef_path : str or Path
        Path to the ``.cellbin.gef`` file, or to the directory containing it
        (e.g. ``.../041.cellcut/``).  When a directory is given the function
        auto-discovers the ``*cellbin*.gef`` file inside.
    sample : str, optional
        Sample name. If None, inferred from the GEF filename (chip name).
    gene_name_index : bool, default True
        If True, use gene names (symbols) as ``.var`` index.
    image_path : str or Path, optional
        Path to the tissue image file (``.tif``), the register directory
        (e.g. ``.../03.register/``), or None to auto-discover.
        When None, the function looks for ``03.register/`` as a sibling or
        parent of the GEF directory.
        The image is stored in ``.uns['spatial'][sample]['images']``.

    Returns
    -------
    sc.AnnData
        AnnData object with:
        - ``.X``: CSR expression matrix (cells × genes)
        - ``.obs``: cell_id, area, dnb_count, gene_count, exp_count, geometry
        - ``.obsm['spatial']``: cell centroid coordinates
        - ``.obs['geometry']``: WKT strings of cell polygons
        - ``.uns['spatial'][sample]['geometries']``: GeoDataFrame with cell polygons
        - ``.uns['spatial'][sample]['images']``: tissue image (if found)
        - ``.uns['spatial'][sample]['scalefactors']`` / ``metadata``

    Examples
    --------
    >>> import trackcell.io as tcio

    >>> # Pass the file directly:
    >>> adata = tcio.read_sto_cellbin(
    ...     "/path/to/C57_7/041.cellcut/B02825B2.cellbin.gef",
    ...     sample="C57_7"
    ... )

    >>> # Or pass the directory — the function finds the right file:
    >>> adata = tcio.read_sto_cellbin(
    ...     "/path/to/C57_7/041.cellcut/",
    ...     sample="C57_7"
    ... )

    >>> # Auto-discover tissue image from register directory:
    >>> adata = tcio.read_sto_cellbin(
    ...     "/path/to/C57_7/041.cellcut/",
    ...     sample="C57_7",
    ...     image_path="/path/to/C57_7/03.register/"
    ... )

    Notes
    -----
    - The GEF file itself does **not** contain the ssDNA tissue image.
      The image is stored separately in the ``03.register/`` directory.
    - The ``cellBorder`` dataset stores up to 32 (dx, dy) offset pairs
      relative to the cell centroid. Trailing slots (sentinel value 32767)
      are stripped during geometry construction.
    - Expression matrix is built directly from the cell offset table
      for efficient CSR construction (no intermediate dense matrix).
    """

    gef_path = Path(gef_path).resolve()
    if gef_path.is_dir():
        gef_path = _find_cellbin_gef(gef_path)
    if not gef_path.exists():
        raise FileNotFoundError(f"GEF file not found: {gef_path}")

    if h5py is None:
        raise ImportError(
            "h5py is required to read GEF files. Install with: pip install h5py"
        )

    t_total = time.time()

    # ── 1. Open GEF and validate ──
    with h5py.File(gef_path, 'r') as f:
        if 'cellBin' not in f:
            raise ValueError(
                f"'{gef_path.name}' does not contain a 'cellBin' group. "
                "Expected a cellbin GEF file."
            )

        cb = f['cellBin']

        # ── 2. Read cell metadata ──
        t0 = time.time()
        cell = cb['cell'][:]
        n_cells = len(cell)
        print(f"[STO cellbin] {n_cells:,} cells [{time.time() - t0:.1f}s]")

        # ── 3. Read cell borders → polygons ──
        t0 = time.time()
        cell_border = cb['cellBorder'][:]  # (n_cells, 32, 2), int16
        border_sentinel = np.int16(_BORDER_SENTINEL)

        cells_x = cell['x'].astype(np.float64)
        cells_y = cell['y'].astype(np.float64)

        polygons = {}
        n_valid = 0
        n_fallback = 0

        for i in range(n_cells):
            offsets = cell_border[i]
            # Quick check: is this cell fully sentinel?
            if np.all(offsets == border_sentinel):
                half = 5.0
                poly = Polygon([
                    (cells_x[i] - half, cells_y[i] - half),
                    (cells_x[i] + half, cells_y[i] - half),
                    (cells_x[i] + half, cells_y[i] + half),
                    (cells_x[i] - half, cells_y[i] + half),
                ])
                n_fallback += 1
            else:
                poly = _extract_polygon_from_relative_offsets(offsets, cells_x[i], cells_y[i])
                if poly is not None:
                    n_valid += 1
                else:
                    n_fallback += 1
            polygons[i] = poly

        print(f"[STO cellbin] Borders: {n_valid:,} valid polygons, "
              f"{n_fallback:,} fallback [{time.time() - t0:.1f}s]")

        # ── 4. Read gene metadata ──
        t0 = time.time()
        gene_data = cb['gene'][:]
        # Some GEF files have gene rows with count==0 (placeholders); keep all.
        gene_names_bytes = gene_data['geneName']
        gene_names = np.array([gn.decode('utf-8') if isinstance(gn, bytes) else str(gn)
                               for gn in gene_names_bytes])
        n_genes = len(gene_names)
        print(f"[STO cellbin] {n_genes:,} genes [{time.time() - t0:.1f}s]")

        # ── 5. Build expression matrix ──
        t0 = time.time()
        cell_exp = cb['cellExp'][:]
        exp_matrix = _build_expression_matrix_from_cell_exp(cell, cell_exp, n_genes, n_cells)
        del cell_exp
        print(f"[STO cellbin] Expression matrix: {exp_matrix.shape} "
              f"({exp_matrix.nnz / 1e6:.1f}M non-zero) [{time.time() - t0:.1f}s]")

        # ── 6. Additional metadata from GEF ──
        # GEF top-level attrs: resolution (nm per unit), offsetX, offsetY
        gef_attrs = dict(f.attrs)
        resolution = int(gef_attrs.get('resolution', [500])[0] if isinstance(gef_attrs.get('resolution'), np.ndarray) else gef_attrs.get('resolution', 500))
        offset_x = int(gef_attrs.get('offsetX', [0])[0] if isinstance(gef_attrs.get('offsetX'), np.ndarray) else 0)
        offset_y = int(gef_attrs.get('offsetY', [0])[0] if isinstance(gef_attrs.get('offsetY'), np.ndarray) else 0)
        # Cell-level attrs give bounding box
        cell_attrs = dict(f['cellBin']['cell'].attrs)
        cell_min_x = int(cell_attrs.get('minX', 0))
        cell_max_x = int(cell_attrs.get('maxX', 0))
        cell_min_y = int(cell_attrs.get('minY', 0))
        cell_max_y = int(cell_attrs.get('maxY', 0))

    # ── 7. Infer sample name ──
    if sample is None:
        # Try to infer from filename: e.g. B02825B2.cellbin.gef -> B02825B2
        sample = gef_path.stem.split('.')[0]

    # ── 8. Build AnnData ──
    t0 = time.time()

    # Cell IDs: use numeric cell IDs from GEF
    cell_ids = np.array([f"cell_{cid}" for cid in cell['id']])

    obs = pd.DataFrame({
        'area': cell['area'].astype(np.int32),
        'dnb_count': cell['dnbCount'].astype(np.int32),
        'gene_count': cell['geneCount'].astype(np.int32),
        'exp_count': cell['expCount'].astype(np.int32),
    }, index=cell_ids)

    if gene_name_index:
        var = pd.DataFrame(index=gene_names)
    else:
        var = pd.DataFrame(index=np.arange(n_genes).astype(str))

    adata = sc.AnnData(X=exp_matrix, obs=obs, var=var)
    adata.obs.index.name = 'cell_id'

    # Spatial coordinates (centroids)
    adata.obsm['spatial'] = np.column_stack([cells_x, cells_y])

    # Store WKT geometry strings in .obs for serialization
    adata.obs['geometry'] = pd.Series(
        {cell_ids[i]: wkt.dumps(poly) for i, poly in polygons.items()},
        name='geometry'
    )

    # ── 9. Setup spatial metadata (trackcell-compatible) ──
    # Physical resolution: 1 coordinate unit = resolution nm
    # STOmics default: 500nm = 0.5μm per unit
    pixel_size_um = resolution / 1000.0  # nm → μm
    # For cellbin, image and GEF coords share the same coordinate system.
    # The registered image covers the full chip (typically 26460×26460 px),
    # while cell coordinates are a subset within that field.
    # scalefactor = 1.0: 1 image pixel = 1 GEF coordinate unit.
    spot_diameter = np.sqrt(np.median(cell['area'])) * 2.0  # approximate

    adata.uns['spatial'] = {
        sample: {
            'metadata': {
                'sample_id': sample,
                'pixel_size_um': pixel_size_um,
                'resolution_nm': resolution,
                'coordinate_range': {
                    'x': [cell_min_x, cell_max_x],
                    'y': [cell_min_y, cell_max_y],
                },
            },
            'scalefactors': {
                'spot_diameter_fullres': float(spot_diameter),
                'tissue_hires_scalef': 1.0,
                'tissue_lowres_scalef': 1.0,
            },
        }
    }

    # Convert polygon dict to GeoDataFrame
    gdf_geometries = [polygons[i] for i in range(n_cells)]
    gdf = gpd.GeoDataFrame(
        {'geometry': gdf_geometries},
        geometry='geometry',
        index=cell_ids,
    )
    gdf.index.name = 'cell_id'
    adata.uns['spatial'][sample]['geometries'] = gdf

    # ── 10. Load tissue image (optional) ──
    tissue_img = _load_tissue_image(image_path, gef_dir=gef_path.parent)
    if tissue_img is not None:
        adata.uns['spatial'][sample]['images'] = {
            'hires': tissue_img,
            'lowres': tissue_img,
        }
        # Update scalefactors based on image dimensions
        # STOmics images are registered to sequencing coordinates,
        # so the scale factor is 1.0 (no scaling needed)
        adata.uns['spatial'][sample]['scalefactors']['tissue_hires_scalef'] = 1.0
        adata.uns['spatial'][sample]['scalefactors']['tissue_lowres_scalef'] = 1.0

    print(f"[STO cellbin] AnnData built [{time.time() - t0:.1f}s]")

    elapsed = time.time() - t_total
    print(f"[STO cellbin] Total: {elapsed:.1f}s")
    print(f"  Output: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    return adata


def read_sto_bin(
    gef_path: Union[str, Path],
    bin_size: int = 100,
    sample: Optional[str] = None,
    gene_name_index: bool = True,
    image_path: Optional[Union[str, Path]] = None,
) -> sc.AnnData:
    """
    Read STOmics squarebin GEF data and create an AnnData object.

    Reads a STOmics tissue/squarebin GEF (HDF5) file and produces a
    bin-level AnnData with spatial coordinates.

    Parameters
    ----------
    gef_path : str or Path
        Path to the ``.tissue.gef`` file, or to the directory containing it
        (e.g. ``.../04.tissuecut/``).  When a directory is given the function
        auto-discovers the ``*tissue*.gef`` file inside.
    bin_size : int, default 100
        Bin size in STOmics coordinate units (~500nm per unit).
        Merges raw 1-unit bins into larger bins.
    sample : str, optional
        Sample name. If None, inferred from the filename.
    gene_name_index : bool, default True
        If True, use gene names as ``.var`` index.
    image_path : str or Path, optional
        Path to the tissue image file (``.tif``), the register directory
        (e.g. ``.../03.register/``), or None to auto-discover.

    Returns
    -------
    sc.AnnData
        AnnData object with bin-level spatial data.

    Examples
    --------
    >>> import trackcell.io as tcio

    >>> # Pass the file directly:
    >>> adata = tcio.read_sto_bin(
    ...     "/path/to/C57_7/04.tissuecut/B02825B2.tissue.gef",
    ...     bin_size=50, sample="C57_7"
    ... )

    >>> # Or pass the directory:
    >>> adata = tcio.read_sto_bin(
    ...     "/path/to/C57_7/04.tissuecut/",
    ...     bin_size=50, sample="C57_7"
    ... )

    >>> # With tissue image auto-discovery:
    >>> adata = tcio.read_sto_bin(
    ...     "/path/to/C57_7/04.tissuecut/",
    ...     bin_size=50,
    ...     image_path="/path/to/C57_7/03.register/"
    ... )

    Notes
    -----
    - GEF is the preferred format (binary HDF5, fast ~15s for 95M entries).
    - GEM text (.gem / .gem.gz) files are also accepted; for large datasets
      the function will try to auto-redirect to the matching .gef file.
    - The GEF file does **not** contain the ssDNA tissue image.
      Images are in ``03.register/`` and can be loaded via ``image_path``.
    """
    path = Path(gef_path).resolve()
    if path.is_dir():
        path = _find_tissue_gef(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    suf = path.suffix.lower()

    if suf in ('.gef',):
        return _read_sto_bin_from_gef(path, bin_size, sample, gene_name_index, image_path)
    elif suf in ('.gem',) or (suf == '.gz' and len(path.suffixes) >= 2 and path.suffixes[-2] == '.gem'):
        # Try to find corresponding GEF
        alt_gef = path.with_suffix('').with_suffix('.gef')
        if alt_gef.exists():
            print(f"[STO bin] Found matching GEF: {alt_gef}, using that instead.")
            return _read_sto_bin_from_gef(alt_gef, bin_size, sample, gene_name_index, image_path)
        # Fallback to GEM reader (may be slow for large datasets)
        warnings.warn(
            "Reading large GEM files can be slow. Consider using the GEF (.gef) format."
        )
        return _read_sto_bin_from_gem(path, bin_size, sample, gene_name_index)
    else:
        raise ValueError(
            f"Unsupported file format: {suf}. Expected .gef, .gem, or .gem.gz"
        )


def _read_sto_bin_from_gef(
    path: Path,
    bin_size: int,
    sample: Optional[str],
    gene_name_index: bool,
    image_path: Optional[Union[str, Path]] = None,
) -> sc.AnnData:
    """Read squarebin data from a GEF (HDF5) file."""
    if h5py is None:
        raise ImportError(
            "h5py is required to read GEF files. Install with: pip install h5py"
        )

    t_total = time.time()

    with h5py.File(path, 'r') as f:
        if 'geneExp' not in f:
            raise ValueError(
                f"'{path.name}' does not contain 'geneExp' group (squarebin format)."
            )

        # GEF structure: geneExp/bin{N}/expression, geneExp/bin{N}/gene
        bin_key = None
        for key in f['geneExp'].keys():
            bin_key = key
            break
        if bin_key is None:
            raise ValueError("No bin group found in geneExp.")

        bin_group = f['geneExp'][bin_key]
        bin_level = int(bin_key.replace('bin', ''))  # e.g. 'bin1' -> 1

        # Extract resolution from expression attrs
        expr_attrs = dict(bin_group['expression'].attrs)
        resolution = int(expr_attrs.get('resolution', 500))
        min_x = int(expr_attrs.get('minX', 0))
        max_x = int(expr_attrs.get('maxX', 0))
        min_y = int(expr_attrs.get('minY', 0))
        max_y = int(expr_attrs.get('maxY', 0))

        t0 = time.time()
        expr = bin_group['expression'][:]
        genes = bin_group['gene'][:]

        xs = expr['x'].astype(np.int64)
        ys = expr['y'].astype(np.int64)
        cnt_vals = expr['count'].astype(np.float32)

        gene_names_bytes = genes['gene']
        gene_names = np.array([
            gn.decode('utf-8') if isinstance(gn, bytes) else str(gn)
            for gn in gene_names_bytes
        ])
        n_genes = len(gene_names)
        n_entries = len(expr)
        print(
            f"[STO bin] {n_entries:,} entries, {n_genes:,} genes "
            f"[{time.time() - t0:.1f}s]"
        )

    # ── Merge bins ──
    t0 = time.time()
    bin_x = (xs // bin_size).astype(np.int64)
    bin_y = (ys // bin_size).astype(np.int64)
    bin_ids_int = (bin_x.astype(np.uint64) << 32) | bin_y.astype(np.uint64)
    unique_bins, inverse = np.unique(bin_ids_int, return_inverse=True)
    n_bins = len(unique_bins)
    print(f"[STO bin] {n_bins:,} unique bins (size={bin_size}) [{time.time() - t0:.1f}s]")

    # ── Gene index vector (the expression array is gene-ordered) ──
    t0 = time.time()
    gene_offsets = genes['offset'].astype(np.int64)
    gene_counts = genes['count'].astype(np.int64)
    col_vec = np.empty(n_entries, dtype=np.int32)
    for g in range(n_genes):
        s = gene_offsets[g]
        e = s + gene_counts[g]
        if e > s:
            col_vec[s:e] = g
    print(f"[STO bin] Gene index vector [{time.time() - t0:.1f}s]")

    # ── Build CSR via COO triples → sum_duplicates → CSR ──
    t0 = time.time()
    coo = csr_matrix(
        (cnt_vals, (inverse.astype(np.int32), col_vec)),
        shape=(n_bins, n_genes),
        dtype=np.float32,
    )
    exp_matrix = coo  # csr_matrix() from COO triplets already handles duplicates
    print(
        f"[STO bin] CSR built: {exp_matrix.shape}, "
        f"nnz={exp_matrix.nnz:,} [{time.time() - t0:.1f}s]"
    )

    # ── Build AnnData ──
    t0 = time.time()
    bin_x_unique = (unique_bins >> 32).astype(np.int64)
    bin_y_unique = (unique_bins & 0xFFFFFFFF).astype(np.int64)

    if sample is None:
        sample = path.stem.split('.')[0]

    obs = pd.DataFrame(index=unique_bins.astype(str))
    adata = sc.AnnData(X=exp_matrix, obs=obs)
    adata.var_names = gene_names if gene_name_index else np.arange(n_genes).astype(str)
    adata.obs.index.name = 'bin_id'

    center_x = (bin_x_unique * bin_size + bin_size // 2).astype(np.float64)
    center_y = (bin_y_unique * bin_size + bin_size // 2).astype(np.float64)
    adata.obsm['spatial'] = np.column_stack([center_x, center_y])

    adata.uns['spatial'] = {
        sample: {
            'metadata': {
                'sample_id': sample,
                'pixel_size_um': resolution / 1000.0,
                'resolution_nm': resolution,
                'bin_size': bin_size,
                'bin_level': bin_level,
                'bin_type': 'square_bins',
                'coordinate_range': {
                    'x': [min_x, max_x],
                    'y': [min_y, max_y],
                },
            },
            'scalefactors': {
                'spot_diameter_fullres': float(bin_size),
                'tissue_hires_scalef': 1.0,
                'tissue_lowres_scalef': 1.0,
            },
        }
    }

    # ── Load tissue image (optional) ──
    tissue_img = _load_tissue_image(image_path, gef_dir=path.parent)
    if tissue_img is not None:
        adata.uns['spatial'][sample]['images'] = {
            'hires': tissue_img,
            'lowres': tissue_img,
        }
        adata.uns['spatial'][sample]['scalefactors']['tissue_hires_scalef'] = 1.0
        adata.uns['spatial'][sample]['scalefactors']['tissue_lowres_scalef'] = 1.0

    print(f"[STO bin] AnnData built [{time.time() - t0:.1f}s]")
    elapsed = time.time() - t_total
    print(f"[STO bin] Total: {elapsed:.1f}s")
    print(f"  Output: {adata.n_obs:,} bins × {adata.n_vars:,} genes")

    return adata


def _read_sto_bin_from_gem(
    path: Path,
    bin_size: int,
    sample: Optional[str],
    gene_name_index: bool,
) -> sc.AnnData:
    """Read squarebin data from a GEM text file (plain or gzip)."""
    import gc
    from collections import defaultdict

    t_total = time.time()
    compression = 'gzip' if path.suffix.lower() == '.gz' else None

    # Streaming groupby using a dict-of-dicts accumulator
    bin_gene_counts: dict = defaultdict(lambda: defaultdict(float))
    gene_set: set = set()
    n_rows = 0

    for chunk in pd.read_csv(
        path, sep='\t', comment='#', header=0,
        dtype={'geneID': str, 'x': np.int64, 'y': np.int64,
               'MIDCount': np.float32, 'MIDCounts': np.float32},
        chunksize=5_000_000, compression=compression,
        engine='python',
    ):
        count_col = 'MIDCounts' if 'MIDCounts' in chunk.columns else 'MIDCount'
        if count_col not in chunk.columns:
            raise ValueError("Missing MIDCount/MIDCounts in GEM file.")

        bx = (chunk['x'].values // bin_size).astype(np.int64)
        by = (chunk['y'].values // bin_size).astype(np.int64)
        gn = chunk['geneID'].values
        cv = chunk[count_col].values.astype(np.float32)

        for i in range(len(bx)):
            key = (int(bx[i]), int(by[i]))
            bin_gene_counts[key][gn[i]] += float(cv[i])
            gene_set.add(str(gn[i]))
            n_rows += 1

        if n_rows % 25_000_000 == 0:
            print(f"[STO bin GEM] Processed {n_rows:,} rows...")

    print(f"[STO bin GEM] {n_rows:,} rows → {len(bin_gene_counts):,} bins, "
          f"{len(gene_set):,} genes [{time.time() - t_total:.1f}s]")

    # ── Build sparse matrix ──
    t0 = time.time()
    genes = sorted(gene_set)
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    n_genes = len(genes)

    bins = sorted(bin_gene_counts.keys())
    bin_to_idx = {b: i for i, b in enumerate(bins)}
    n_bins = len(bins)

    # Build CSR directly
    indptr = np.zeros(n_bins + 1, dtype=np.int32)
    nnz_total = sum(len(d) for d in bin_gene_counts.values())
    data_arr = np.empty(nnz_total, dtype=np.float32)
    cols_arr = np.empty(nnz_total, dtype=np.int32)

    pos = 0
    for bi, bkey in enumerate(bins):
        d = bin_gene_counts[bkey]
        indptr[bi + 1] = indptr[bi] + len(d)
        for gname, cnt in d.items():
            cols_arr[pos] = gene_to_idx[gname]
            data_arr[pos] = cnt
            pos += 1

    exp_matrix = csr_matrix(
        (data_arr, cols_arr, indptr),
        shape=(n_bins, n_genes),
        dtype=np.float32,
    )
    del bin_gene_counts, data_arr, cols_arr; gc.collect()
    print(f"[STO bin GEM] CSR: {exp_matrix.shape}, nnz={exp_matrix.nnz:,} [{time.time() - t0:.1f}s]")

    # ── Build AnnData ──
    t0 = time.time()
    bin_x_arr = np.array([b[0] for b in bins], dtype=np.int64)
    bin_y_arr = np.array([b[1] for b in bins], dtype=np.int64)

    if sample is None:
        sample = path.stem.split('.')[0]
    if sample.endswith('.tissue'):
        sample = sample.rsplit('.', 1)[0]

    center_x = (bin_x_arr * bin_size + bin_size // 2).astype(np.float64)
    center_y = (bin_y_arr * bin_size + bin_size // 2).astype(np.float64)

    obs = pd.DataFrame(
        index=[f"{x}_{y}" for x, y in zip(bin_x_arr, bin_y_arr)]
    )
    adata = sc.AnnData(X=exp_matrix, obs=obs)
    adata.var_names = genes
    adata.obs.index.name = 'bin_id'
    adata.obsm['spatial'] = np.column_stack([center_x, center_y])

    adata.uns['spatial'] = {
        sample: {
            'metadata': {'sample_id': sample, 'bin_size': bin_size, 'bin_type': 'square_bins'},
            'scalefactors': {'spot_diameter_fullres': float(bin_size)},
        }
    }

    elapsed = time.time() - t_total
    print(f"[STO bin GEM] Total: {elapsed:.1f}s")
    print(f"  Output: {adata.n_obs:,} bins × {adata.n_vars:,} genes")
    return adata


def _read_sto_bin_from_gem_gz(
    path: Path,
    bin_size: int,
    sample: Optional[str],
    gene_name_index: bool,
) -> sc.AnnData:
    """Read squarebin data from a gzipped GEM file (delegates to GEM reader)."""
    return _read_sto_bin_from_gem(path, bin_size, sample, gene_name_index)
