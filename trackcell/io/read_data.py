"""
Data reading functions for TrackCell package.

This module provides functions for reading spatial transcriptomics data,
particularly from SpaceRanger output (both bin-level and cell segmentation data).
"""

import scanpy as sc
import geopandas as gpd
import pandas as pd
from shapely import wkt, geometry
import numpy as np
import json
import imageio.v3 as iio
from pathlib import Path
from typing import Optional, Union
import ast
import warnings
import time
import gc
import h5py
from scipy import sparse as sp_sparse

def read_hd_bin(
    datapath: Union[str, Path],
    sample: Optional[str] = None,
    binsize: int = 16,
    matrix_file_h5: str = "filtered_feature_bc_matrix.h5",
    matrix_file_dir: str = "filtered_feature_bc_matrix",
    tissue_positions_file: str = "spatial/tissue_positions.parquet",
    hires_image_file: str = "spatial/tissue_hires_image.png",
    lowres_image_file: str = "spatial/tissue_lowres_image.png",
    scalefactors_file: str = "spatial/scalefactors_json.json"
) -> sc.AnnData:
    """
    Read 10X HD SpaceRanger bin-level output (2um/8um/16um) and create an AnnData object with spatial information.
    
    This function reads the output from SpaceRanger pipeline and creates an AnnData object
    that includes spatial coordinates, tissue images, and scalefactors for bin-level data.
    
    Parameters
    ----------
    datapath : str or Path
        Path to the SpaceRanger output directory containing bin-level outputs.
    sample : str, optional
        Sample name. If None, will be inferred from the path.
    binsize : int, default 16
        Bin size in micrometers. Common values are 2, 8, or 16. This information
        will be stored in adata.uns["spatial"][sample]["binsize"].
    matrix_file_h5 : str, default "filtered_feature_bc_matrix.h5"
        Name of the H5 matrix file. Will be tried first.
    matrix_file_dir : str, default "filtered_feature_bc_matrix"
        Name of the matrix directory. Will be used if H5 file is not available.
    tissue_positions_file : str, default "spatial/tissue_positions.parquet"
        Path to tissue positions file (parquet or csv format).
    hires_image_file : str, default "spatial/tissue_hires_image.png"
        Path to the high-resolution tissue image relative to datapath.
    lowres_image_file : str, default "spatial/tissue_lowres_image.png"
        Path to the low-resolution tissue image relative to datapath.
    scalefactors_file : str, default "spatial/scalefactors_json.json"
        Name of the scalefactors JSON file.
    
    Returns
    -------
    sc.AnnData
        AnnData object containing:
        - Expression matrix in .X
        - Cell metadata in .obs
        - Gene metadata in .var
        - Spatial coordinates in .obsm["spatial"]
        - Tissue images in .uns["spatial"][sample]["images"]
        - Scalefactors in .uns["spatial"][sample]["scalefactors"]
        - Bin size in .uns["spatial"][sample]["binsize"]
    
    Examples
    --------
    >>> import trackcell.io as tcio
    >>> adata = tcio.read_hd_bin("SpaceRanger4.0/Case1/outs", sample="Case1")
    >>> print(adata)
    AnnData object with n_obs × n_vars = 10000 × 2000
        obs: 'barcode'
        obsm: 'spatial'
        uns: 'spatial'
    
    Notes
    -----
    This function expects the SpaceRanger output to have the following structure:
    - filtered_feature_bc_matrix.h5 or filtered_feature_bc_matrix/: Expression matrix
    - spatial/tissue_positions.parquet or spatial/tissue_positions.csv: Spatial coordinates
    - spatial/tissue_hires_image.png: High-resolution tissue image
    - spatial/tissue_lowres_image.png: Low-resolution tissue image
    - spatial/scalefactors_json.json: Image scaling factors
    """
    
    # Convert to Path object for easier handling
    datapath = Path(datapath).resolve()
    
    # If sample is not provided, try to infer from path
    if sample is None:
        sample = datapath.parent.parent.name if datapath.name == "outs" else datapath.name
    
    # Read expression matrix - try H5 first, then directory
    h5_path = datapath / matrix_file_h5
    matrix_dir_path = datapath / matrix_file_dir
    
    if h5_path.exists():
        try:
            adata = sc.read_10x_h5(h5_path)
        except Exception as e:
            print(f"Warning: Could not read H5 file {h5_path}: {e}")
            print(f"Trying to read from directory {matrix_dir_path}")
            if matrix_dir_path.exists():
                adata = sc.read_10x_mtx(matrix_dir_path)
            else:
                raise FileNotFoundError(f"Neither H5 file nor matrix directory found in {datapath}")
    elif matrix_dir_path.exists():
        adata = sc.read_10x_mtx(matrix_dir_path)
    else:
        raise FileNotFoundError(f"Neither {h5_path} nor {matrix_dir_path} found")
    
    # Load the Spatial Coordinates
    tissue_pos_path = datapath / tissue_positions_file
    
    if not tissue_pos_path.exists():
        # Try CSV format as fallback
        tissue_pos_path_csv = datapath / tissue_positions_file.replace('.parquet', '.csv')
        if tissue_pos_path_csv.exists():
            tissue_pos_path = tissue_pos_path_csv
        else:
            raise FileNotFoundError(f"Tissue positions file not found: {tissue_pos_path} or {tissue_pos_path_csv}")
    
    # Read tissue positions file
    try:
        if tissue_pos_path.suffix == '.parquet':
            df_tissue_positions = pd.read_parquet(tissue_pos_path)
        else:
            # Try CSV with different separators
            try:
                df_tissue_positions = pd.read_csv(tissue_pos_path, sep=',')
            except Exception:
                try:
                    df_tissue_positions = pd.read_csv(tissue_pos_path, sep='\t')
                except Exception:
                    df_tissue_positions = pd.read_csv(tissue_pos_path, sep=None, engine='python')
    except Exception as e:
        raise ValueError(f"Could not read tissue positions file {tissue_pos_path}: {e}")
    
    # Set the index of the dataframe to the barcodes
    if 'barcode' in df_tissue_positions.columns:
        df_tissue_positions = df_tissue_positions.set_index('barcode')
    elif df_tissue_positions.index.name == 'barcode' or 'barcode' in str(df_tissue_positions.index.name):
        # Already indexed by barcode
        pass
    else:
        # Try to use first column as barcode if it's not already indexed
        if len(df_tissue_positions.columns) > 0:
            df_tissue_positions = df_tissue_positions.set_index(df_tissue_positions.columns[0])
    
    if df_tissue_positions.index.has_duplicates:
        duplicate_examples = df_tissue_positions.index[df_tissue_positions.index.duplicated()].unique()[:5].tolist()
        raise ValueError(
            "Tissue positions contain duplicate barcodes; cannot align spatial coordinates safely. "
            f"Examples: {duplicate_examples}"
        )

    missing_positions = pd.Index(adata.obs_names).difference(df_tissue_positions.index)
    if len(missing_positions) > 0:
        raise ValueError(
            f"{len(missing_positions)} expression barcodes are missing from tissue positions. "
            f"Examples: {missing_positions[:5].tolist()}"
        )

    # Adding the tissue positions to the metadata
    adata.obs = pd.merge(adata.obs, df_tissue_positions, left_index=True, right_index=True, how='left')
    
    # Extract spatial coordinates
    # Try different possible column names for spatial coordinates
    coord_cols = None
    for col_pair in [
        ["pxl_col_in_fullres", "pxl_row_in_fullres"],
        ["pxl_col", "pxl_row"],
        ["x", "y"],
        ["array_col", "array_row"]
    ]:
        if all(col in adata.obs.columns for col in col_pair):
            coord_cols = col_pair
            break
    
    if coord_cols is None:
        raise ValueError(
            "Could not find spatial coordinate columns. "
            "Expected one of: ['pxl_col_in_fullres', 'pxl_row_in_fullres'], "
            "['pxl_col', 'pxl_row'], ['x', 'y'], or ['array_col', 'array_row']"
        )
    
    coord_values = adata.obs[coord_cols]
    if coord_values.isna().any().any():
        bad_barcodes = coord_values[coord_values.isna().any(axis=1)].index[:5].tolist()
        raise ValueError(
            "Spatial coordinate columns contain missing values after alignment. "
            f"Examples: {bad_barcodes}"
        )
    adata.obsm['spatial'] = coord_values.values
    
    # Read tissue images
    try:
        hires_img = iio.imread(datapath / hires_image_file)
        lowres_img = iio.imread(datapath / lowres_image_file)
    except FileNotFoundError as e:
        print(f"Warning: Could not load tissue images: {e}")
        hires_img = None
        lowres_img = None
    
    # Initialize spatial metadata
    adata.uns["spatial"] = {}
    adata.uns["spatial"][sample] = {}
    
    # Load scalefactors
    try:
        with open(datapath / scalefactors_file, 'r', encoding='utf-8') as file:
            scalefactor = json.load(file)
    except FileNotFoundError as e:
        print(f"Warning: Could not load scalefactors: {e}")
        scalefactor = {}
    
    # Store images and scalefactors
    if hires_img is not None and lowres_img is not None:
        adata.uns["spatial"][sample]["images"] = {
            "hires": hires_img,
            "lowres": lowres_img
        }
    
    adata.uns["spatial"][sample]["scalefactors"] = scalefactor
    adata.uns["spatial"][sample]["binsize"] = binsize
    
    return adata

def read_hd_cellseg(
    datapath: Union[str, Path],
    sample: Optional[str] = None,
    cell_segmentations_file: str = "graphclust_annotated_cell_segmentations.geojson",
    matrix_file: str = "filtered_feature_cell_matrix.h5",
    hires_image_file: str = "spatial/tissue_hires_image.png",
    lowres_image_file: str = "spatial/tissue_lowres_image.png",
    scalefactors_file: str = "spatial/scalefactors_json.json"
) -> sc.AnnData:
    """
    Read 10X HD SpaceRanger cell segmentation output and create an AnnData object with spatial information.
    
    This function reads the output from SpaceRanger pipeline and creates an AnnData object
    that includes spatial coordinates, cell segmentations, and tissue images.
    
    Parameters
    ----------
    datapath : str or Path
        Path to the SpaceRanger output directory containing segmented outputs.
    sample : str, optional
        Sample name. If None, will be inferred from the path.
    cell_segmentations_file : str, default "graphclust_annotated_cell_segmentations.geojson"
        Name of the cell segmentations file.
    matrix_file : str, default "filtered_feature_cell_matrix.h5"
        Name of the filtered feature-cell matrix file.
    hires_image_file : str, default "spatial/tissue_hires_image.png"
        Path to the high-resolution tissue image relative to datapath.
    lowres_image_file : str, default "spatial/tissue_lowres_image.png"
        Path to the low-resolution tissue image relative to datapath.
    scalefactors_file : str, default "spatial/scalefactors_json.json"
        Name of the scalefactors JSON file.
    
    Returns
    -------
    sc.AnnData
        AnnData object containing:
        - Expression matrix in .X
        - Cell metadata in .obs
        - Gene metadata in .var
        - Spatial coordinates in .obsm["spatial"]
        - Tissue images in .uns["spatial"][sample]["images"]
        - Scalefactors in .uns["spatial"][sample]["scalefactors"]
        - Cell geometries in .uns["spatial"][sample]["geometries"] (GeoDataFrame)
        - Cell geometries in .obs["geometry"] (WKT strings for serialization)
    
    Examples
    --------
    >>> import trackcell.io as tcio
    >>> adata = tcio.read_hd_cellseg("SpaceRanger4.0/Case1/outs/segmented_outputs", sample="Case1")
    >>> print(adata)
    AnnData object with n_obs × n_vars = 1000 × 2000
        obs: 'cellid'
        obsm: 'spatial'
        uns: 'spatial'
    
    Notes
    -----
    This function expects the SpaceRanger output to have the following structure:
    - graphclust_annotated_cell_segmentations.geojson: Cell segmentation polygons
    - filtered_feature_cell_matrix.h5: Expression matrix
    - spatial/tissue_hires_image.png: High-resolution tissue image
    - spatial/tissue_lowres_image.png: Low-resolution tissue image
    - spatial/scalefactors_json.json: Image scaling factors
    """
    
    # Convert to Path object for easier handling
    datapath = Path(datapath).resolve()
    
    # If sample is not provided, try to infer from path
    if sample is None:
        sample = datapath.parent.parent.name
    
    # Read cell segmentations
    # Try the specified file first, then try common alternative names if not found
    seg_file_path = datapath / cell_segmentations_file
    if not seg_file_path.exists():
        # Try common alternative filenames
        alternative_names = [
            "cell_segmentations.geojson",
            "cell_segmentations_annotated.geojson",
            "annotated_cell_segmentations.geojson"
        ]
        found_file = None
        for alt_name in alternative_names:
            alt_path = datapath / alt_name
            if alt_path.exists():
                found_file = alt_name
                warnings.warn(
                    f"Specified file '{cell_segmentations_file}' not found. "
                    f"Using alternative file '{alt_name}' instead."
                )
                break
        
        if found_file is None:
            raise FileNotFoundError(
                f"Cell segmentations file not found: {seg_file_path}\n"
                f"Also tried: {alternative_names}\n"
                f"Please specify the correct filename using the `cell_segmentations_file` parameter."
            )
        seg_file_path = datapath / found_file
    
    gdf_seg = gpd.read_file(seg_file_path)
    df = pd.DataFrame(gdf_seg)
    
    if "cell_id" not in df.columns:
        raise ValueError("Cell segmentation file must contain a 'cell_id' column.")
    if "geometry" not in df.columns:
        raise ValueError("Cell segmentation file must contain a geometry column.")

    # Create cellid in the format expected by SpaceRanger
    df['cellid'] = df['cell_id'].apply(lambda x: f"cellid_{str(x).zfill(9)}-1")
    if df['cellid'].duplicated().any():
        duplicate_examples = df.loc[df['cellid'].duplicated(), 'cellid'].head().tolist()
        raise ValueError(
            "Cell segmentation file contains duplicate cell IDs; cannot align safely. "
            f"Examples: {duplicate_examples}"
        )
    
    # Read expression matrix
    adata = sc.read_10x_h5(datapath / matrix_file)
    
    # Align cell segmentations with expression data
    expr_cellids = pd.Index(adata.obs_names)
    seg_cellids = pd.Index(df['cellid'])
    matched_cellids = expr_cellids.intersection(seg_cellids)
    if len(matched_cellids) == 0:
        raise ValueError(
            "No expression barcodes match the cell segmentation IDs. "
            f"Expression examples: {expr_cellids[:5].tolist()}; "
            f"segmentation examples: {seg_cellids[:5].tolist()}"
        )
    if len(matched_cellids) < len(expr_cellids):
        missing_seg = expr_cellids.difference(seg_cellids)
        warnings.warn(
            f"{len(missing_seg)} expression cells have no segmentation geometry and will be dropped. "
            f"Examples: {missing_seg[:5].tolist()}"
        )
    extra_seg = seg_cellids.difference(expr_cellids)
    if len(extra_seg) > 0:
        warnings.warn(
            f"{len(extra_seg)} segmented cells have no expression matrix row and will be ignored. "
            f"Examples: {extra_seg[:5].tolist()}"
        )
    adata = adata[expr_cellids.isin(seg_cellids), :].copy()
    
    # Filter and reorder df to match adata.obs_names order
    df = df.set_index("cellid").reindex(adata.obs_names)
    if df["geometry"].isna().any():
        bad_cellids = df[df["geometry"].isna()].index[:5].tolist()
        raise ValueError(
            "Missing segmentation geometries after expression/segmentation alignment. "
            f"Examples: {bad_cellids}"
        )
    
    # Convert geometry strings to shapely objects if needed
    first_geometry = df["geometry"].dropna().iloc[0]
    if isinstance(first_geometry, str):
        df["geometry"] = df["geometry"].apply(wkt.loads)
    
    # Extract centroid coordinates
    df["x"] = df["geometry"].apply(lambda poly: poly.centroid.x)
    df["y"] = df["geometry"].apply(lambda poly: poly.centroid.y)
    
    # Store spatial coordinates
    adata.obsm["spatial"] = np.array(df[["x", "y"]])
    
    if 'classification' in df.columns and not df['classification'].isna().all():
        first_classification = df['classification'].dropna().iloc[0]
        if isinstance(first_classification, str):
            classifications = df['classification'].apply(
                lambda value: ast.literal_eval(value) if pd.notna(value) else None
            )
        else:
            classifications = df['classification']
        df['classification'] = classifications
        adata.obs['classification'] = [
            item.get('name') if isinstance(item, dict) else None
            for item in classifications
        ]
        adata.uns['classification_colors'] = convert_classification_to_color_dict(df, 'classification')

    # Read tissue images
    try:
        hires_img = iio.imread(datapath / hires_image_file)
        lowres_img = iio.imread(datapath / lowres_image_file)
    except FileNotFoundError as e:
        print(f"Warning: Could not load tissue images: {e}")
        hires_img = None
        lowres_img = None
    
    # Initialize spatial metadata
    adata.uns["spatial"] = {}
    adata.uns["spatial"][sample] = {}
    
    # Load scalefactors
    try:
        with open(datapath / scalefactors_file, 'r', encoding='utf-8') as file:
            scalefactor = json.load(file)
    except FileNotFoundError as e:
        print(f"Warning: Could not load scalefactors: {e}")
        scalefactor = {}
    
    # Add spot_diameter_fullres if missing (required by scanpy's sc.pl.spatial)
    # For cell segmentation data, this is not a real spot diameter but needed for plotting
    # 
    # Note: This is only used when spot_size parameter is NOT provided to sc.pl.spatial()
    # - If spot_size is provided: scanpy uses the provided value directly
    # - If spot_size is None: scanpy reads from scalefactors["spot_diameter_fullres"]
    # 
    # Setting this ensures compatibility when users don't specify spot_size parameter
    if "spot_diameter_fullres" not in scalefactor:
        # Use a reasonable default for cell segmentation data
        # Typical cell diameter in pixels (scaled to fullres)
        # If we have fiducial_diameter_fullres, we can estimate based on that
        if "fiducial_diameter_fullres" in scalefactor:
            # Estimate spot diameter as a fraction of fiducial diameter
            # Fiducial markers are typically much larger than cells/spots
            estimated_spot_diameter = scalefactor["fiducial_diameter_fullres"] / 40.0
        else:
            # Default to 20 pixels for cell segmentation (cells are smaller than spots)
            estimated_spot_diameter = 20.0
        scalefactor["spot_diameter_fullres"] = estimated_spot_diameter
    
    # Store images and scalefactors
    if hires_img is not None and lowres_img is not None:
        adata.uns["spatial"][sample]["images"] = {
            "hires": hires_img,
            "lowres": lowres_img
        }
    
    adata.uns["spatial"][sample]["scalefactors"] = scalefactor
    
    # Store geometries: GeoDataFrame for fast access and WKT strings for serialization
    # Create GeoDataFrame indexed by cellid for easy lookup
    # cellid should already be a column from the reset_index() above
    gdf_indexed = gpd.GeoDataFrame(df[["geometry"]], geometry="geometry")
    #gdf_indexed = gdf_indexed.set_index(df.index)
    adata.uns["spatial"][sample]["geometries"] = gdf_indexed
    
    # Also store WKT strings in obs for serialization compatibility
    adata.obs["geometry"] = df["geometry"].apply(lambda g: wkt.dumps(g) if g is not None else None)
    
    return adata 


def read_xenium_cellseg(
    datapath: Union[str, Path],
    sample: Optional[str] = None,
    cells_file: str = "cells.parquet",
    boundaries_file: str = "cell_boundaries.parquet",
    nucleus_boundaries_file: str = "nucleus_boundaries.parquet",
    matrix_file: str = "cell_feature_matrix.h5",
    experiment_file: str = "experiment.xenium",
    panel_file: str = "gene_panel.json",
    slice_separate: bool = False,
    slice_eps: float = 80.0,
    slice_min_samples: int = 10,
    slice_min_cells: int = 1000,
    slice_start: int = 1,
    slice_prefix: str = "S",
) -> sc.AnnData:
    """
    Read 10x Xenium cell segmentation output and create an AnnData object.
    
    Reads Xenium Analyzer output files (parquet boundaries, HDF5 expression matrix)
    and builds a squidpy-compatible AnnData with GeoDataFrame cell polygons.
    
    Parameters
    ----------
    datapath : str or Path
        Path to the Xenium output directory containing:
        - cell_feature_matrix.h5
        - cells.parquet
        - cell_boundaries.parquet
        - experiment.xenium
    sample : str, optional
        Sample name. If None, infers from directory name.
    cells_file : str
        Name of the cell metadata parquet file.
    boundaries_file : str
        Name of the cell boundary parquet file (long-table format).
    nucleus_boundaries_file : str
        Name of the nucleus boundary parquet file.
    matrix_file : str
        Name of the HDF5 expression matrix file.
    experiment_file : str
        Name of the experiment metadata JSON file.
    panel_file : str
        Name of the gene panel JSON file.
    slice_separate : bool, default False
        If True, run DBSCAN slice separation on cell centroids and add
        ``obs['slice_id']`` / ``obs['slice_cluster']`` (for Xenium TMA data).
    slice_eps : float, default 80.0
        DBSCAN neighborhood radius in µm when ``slice_separate=True``.
    slice_min_samples : int, default 10
        DBSCAN ``min_samples`` when ``slice_separate=True``.
    slice_min_cells : int, default 1000
        Minimum cells per cluster to keep as a valid slice.
    slice_start : int, default 1
        Starting index for slice IDs (S001 when ``slice_start=1``).
    slice_prefix : str, default "S"
        Prefix for slice IDs.
    
    Returns
    -------
    sc.AnnData
        AnnData object with:
        - Expression matrix in .X (CSR, cells × genes)
        - Cell metadata in .obs (centroids, area, counts, etc.)
        - Gene metadata in .var (gene_ids, feature_types)
        - Spatial coordinates in .obsm['spatial']
        - Cell polygons in .uns['spatial'][sample]['geometries'] (GeoDataFrame)
        - Cell boundary arrays in .uns['cell_boundaries'] (compact vertex arrays)
        - Nucleus boundary arrays in .uns['nucleus_boundaries'] (if available)
        - WKT geometry strings in .obs['geometry'] (for serialization)
        - Experiment metadata in .uns['experiment']
    
    Examples
    --------
    >>> import trackcell.io as tcio
    >>> adata = tcio.read_xenium_cellseg("/path/to/xenium/output")
    >>> print(adata)
    AnnData object with n_obs × n_vars = 819023 × 10029
    """
    datapath = Path(datapath).resolve()
    
    if sample is None:
        sample = datapath.name
    
    t_total = time.time()
    
    # ── 1. Read expression matrix ──
    t0 = time.time()
    h5_path = datapath / matrix_file
    if not h5_path.exists():
        raise FileNotFoundError(f"Expression matrix not found: {h5_path}")
    
    with h5py.File(h5_path, 'r') as f:
        m = f['matrix']
        shape = tuple(m['shape'][:])
        data = m['data'][:]
        indices = m['indices'][:]
        indptr = m['indptr'][:]
        barcodes = np.array([b.decode('utf-8') if isinstance(b, bytes) else str(b)
                            for b in m['barcodes'][:]])
        feat_names = np.array([n.decode('utf-8') if isinstance(n, bytes) else str(n)
                              for n in m['features']['name'][:]])
        feat_ids = np.array([i.decode('utf-8') if isinstance(i, bytes) else str(i)
                            for i in m['features']['id'][:]])
        feat_types = np.array([t.decode('utf-8') if isinstance(t, bytes) else str(t)
                              for t in m['features']['feature_type'][:]])
    
    # 10x HDF5 stores CSC (genes × cells); transpose to (cells × genes) CSR
    n_genes, n_cells = int(shape[0]), int(shape[1])
    mat_csc = sp_sparse.csc_matrix((data, indices, indptr), shape=(n_genes, n_cells))
    mat = mat_csc.tocsr().T
    del data, indices, indptr, mat_csc
    gc.collect()
    print(f'[Xenium] Expression matrix: {mat.shape[0]:,} cells × {mat.shape[1]:,} genes '
          f'({mat.nnz/1e6:.1f}M non-zero) [{time.time()-t0:.1f}s]')
    
    # ── 2. Read cell metadata ──
    t0 = time.time()
    cells_path = datapath / cells_file
    if not cells_path.exists():
        raise FileNotFoundError(f"Cell metadata not found: {cells_path}")
    cells_df = pd.read_parquet(cells_path)
    cells_df = cells_df.set_index('cell_id')
    print(f'[Xenium] Cell metadata: {len(cells_df):,} cells [{time.time()-t0:.1f}s]')
    
    # ── 3. Align expression barcodes with cell metadata ──
    t0 = time.time()
    # Xenium cell_ids already include suffix; use as-is
    obs_index = pd.Index(barcodes, name='cell_id')
    
    common_cells = obs_index.intersection(cells_df.index)
    if len(common_cells) != n_cells:
        print(f'  Warning: {len(common_cells):,}/{n_cells:,} cells matched between '
              f'expression and metadata')
    
    # Build row index map for subsetting the matrix
    row_mask = obs_index.isin(common_cells)
    row_indices = np.where(row_mask)[0]
    mat = mat[row_indices, :]
    
    # Reorder metadata to match expression order
    cells_df = cells_df.loc[common_cells]
    print(f'[Xenium] Aligned: {len(common_cells):,} cells [{time.time()-t0:.1f}s]')
    
    # ── 4. Build AnnData ──
    t0 = time.time()
    obs = pd.DataFrame(index=common_cells)
    obs_cols = ['x_centroid', 'y_centroid', 'cell_area', 'nucleus_area',
                'transcript_counts', 'total_counts', 'nucleus_count',
                'control_probe_counts', 'genomic_control_counts',
                'control_codeword_counts', 'unassigned_codeword_counts',
                'segmentation_method']
    for col in obs_cols:
        if col in cells_df.columns:
            obs[col] = cells_df[col].values
    
    var = pd.DataFrame(
        {'gene_ids': feat_ids, 'feature_types': feat_types},
        index=feat_names
    )
    var.index.name = 'gene_symbol'
    
    adata = sc.AnnData(X=mat, obs=obs, var=var)
    adata.obs.index.name = 'cell_id'
    
    if not sp_sparse.isspmatrix_csr(adata.X):
        adata.X = adata.X.tocsr()
    print(f'[Xenium] AnnData built: {adata.shape} [{time.time()-t0:.1f}s]')
    
    # ── 5. Spatial coordinates ──
    t0 = time.time()
    adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values.astype(np.float32)
    print(f'[Xenium] Spatial coords [{time.time()-t0:.1f}s]')
    
    # ── 6. Cell boundaries → GeoDataFrame polygons + compact arrays ──
    t0 = time.time()
    bd_path = datapath / boundaries_file
    if bd_path.exists():
        bd_df = pd.read_parquet(bd_path)
        print(f'  Loaded {len(bd_df):,} boundary vertices [{time.time()-t0:.1f}s]')
        
        t1 = time.time()
        cell_set = set(common_cells)
        cell_idx_map = {cid: i for i, cid in enumerate(common_cells)}
        
        # Single pass: build both compact arrays and polygons
        n_vertices_list, all_vx, all_vy = [], [], []
        cell_to_bd_idx = np.full(len(common_cells), -1, dtype=np.int32)
        polygons = {}
        cum_vertices = 0
        
        for cell_id, group in bd_df.groupby('cell_id', sort=False):
            if cell_id not in cell_set:
                continue
            vx = group['vertex_x'].values
            vy = group['vertex_y'].values
            idx = cell_idx_map[cell_id]
            
            # Compact arrays
            cell_to_bd_idx[idx] = cum_vertices
            nv = len(vx)
            n_vertices_list.append(nv)
            cum_vertices += nv
            all_vx.append(vx)
            all_vy.append(vy)
            
            # Polygon
            if nv >= 3:
                poly = geometry.Polygon(np.column_stack([vx, vy]))
                if poly.is_valid and not poly.is_empty:
                    polygons[cell_id] = poly
        
        if all_vx:
            adata.uns['cell_boundaries'] = {
                'vertex_x': np.concatenate(all_vx).astype(np.float32),
                'vertex_y': np.concatenate(all_vy).astype(np.float32),
                'n_vertices': np.array(n_vertices_list, dtype=np.int32),
                'cell_idx': cell_to_bd_idx,
            }
        del bd_df
        gc.collect()
        print(f'  Boundaries: {len(polygons):,}/{len(common_cells):,} cells '
              f'({len(polygons)/len(common_cells)*100:.1f}%) '
              f'[{time.time()-t1:.1f}s]')
    else:
        print(f'  No cell boundary file: {bd_path}')
        polygons = {}
    
    # ── 7. Nucleus boundaries (optional) ──
    nuc_path = datapath / nucleus_boundaries_file
    if nuc_path.exists():
        t0 = time.time()
        nuc_df = pd.read_parquet(nuc_path)
        nuc_grouped = nuc_df.groupby('cell_id', sort=False)
        
        nuc_nv_list, nuc_vx, nuc_vy = [], [], []
        nuc_to_idx = np.full(len(common_cells), -1, dtype=np.int32)
        cum_nuc = 0
        
        for cell_id, group in nuc_grouped:
            if cell_id in cell_set:
                vx = group['vertex_x'].values
                vy = group['vertex_y'].values
                idx = cell_idx_map[cell_id]
                nuc_to_idx[idx] = cum_nuc
                nuc_nv_list.append(len(vx))
                cum_nuc += len(vx)
                nuc_vx.append(vx)
                nuc_vy.append(vy)
        
        if nuc_vx:
            adata.uns['nucleus_boundaries'] = {
                'vertex_x': np.concatenate(nuc_vx).astype(np.float32),
                'vertex_y': np.concatenate(nuc_vy).astype(np.float32),
                'n_vertices': np.array(nuc_nv_list, dtype=np.int32),
                'cell_idx': nuc_to_idx,
            }
        del nuc_df
        gc.collect()
        print(f'[Xenium] Nucleus boundaries [{time.time()-t0:.1f}s]')
    
    # ── 8. Experiment metadata ──
    exp_path = datapath / experiment_file
    if exp_path.exists():
        with open(exp_path) as f:
            adata.uns['experiment'] = json.load(f)
    
    # ── 9. Gene panel metadata ──
    panel_path = datapath / panel_file
    if panel_path.exists():
        with open(panel_path) as f:
            adata.uns['gene_panel'] = json.load(f)
    
    # ── 10. Setup spatial metadata (squidpy-compatible) ──
    adata.uns['spatial'] = {
        sample: {
            'metadata': {
                'sample_id': sample,
                'pixel_size': adata.uns.get('experiment', {}).get('pixel_size', 0.2125),
            },
            'scalefactors': {
                'spot_diameter_fullres': 1.0,
                'tissue_hires_scalef': 1.0,
                'tissue_lowres_scalef': 1.0,
            },
        }
    }
    
    # ── 11. Store geometries (GeoDataFrame + WKT strings) ──
    if polygons:
        # Build GeoDataFrame indexed by cell_id
        gdf = gpd.GeoDataFrame(
            {'geometry': pd.Series(polygons)},
            geometry='geometry'
        )
        gdf.index.name = 'cell_id'
        adata.uns['spatial'][sample]['geometries'] = gdf
        
        # Store WKT strings in obs for serialization
        adata.obs['geometry'] = pd.Series(
            {cid: wkt.dumps(poly) for cid, poly in polygons.items()},
            name='geometry'
        )
        
        n_geo = len(polygons)
        print(f'[Xenium] Geometries: {n_geo:,}/{len(common_cells):,} cells '
              f'({n_geo/len(common_cells)*100:.1f}%)')
    
    if slice_separate:
        from ..tl.slice_separation import spatial_slice_cluster

        spatial_slice_cluster(
            adata,
            eps=slice_eps,
            min_samples=slice_min_samples,
            min_cells=slice_min_cells,
            slice_start=slice_start,
            slice_prefix=slice_prefix,
        )
        n_slices = adata.uns.get("slice_id_summary", {}).get("n_slices", 0)
        print(f'[Xenium] Slice separation: {n_slices} slices (eps={slice_eps})')

    elapsed = time.time() - t_total
    print(f'[Xenium] Total: read_xenium_cellseg done [{elapsed:.1f}s]')
    
    return adata


def convert_classification_to_color_dict(df, classification_col='classification'):
    """
    Convert classification information in DataFrame column to color dictionary.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing classification information.
    classification_col : str, default 'classification'
        Column name containing classification information.
    
    Returns
    -------
    dict
        Dictionary mapping classification names to hexadecimal color codes.
    """
    color_dict = {}
    for value in df[classification_col].dropna():
        if isinstance(value, str):
            try:
                value = ast.literal_eval(value)
            except (ValueError, SyntaxError):
                continue
        entries = value if isinstance(value, list) else [value]
        for cls in entries:
            if not isinstance(cls, dict):
                continue
            name = cls.get('name')
            rgb = cls.get('color')
            if name is None or rgb is None or len(rgb) < 3:
                continue
            # Convert RGB list to hexadecimal color code.
            hex_color = '#{:02x}{:02x}{:02x}'.format(*rgb[:3])
            color_dict[name] = hex_color
    
    return color_dict


def sync_geometries_after_subset(adata, sample: Optional[str] = None):
    """
    Synchronize geometries in adata.uns["spatial"][sample]["geometries"] 
    after subsetting the AnnData object.
    
    When you subset an AnnData object (e.g., adata[mask] or adata[adata.obs['col'] == 'value']),
    the adata.obs["geometry"] column is automatically subset, but the GeoDataFrame in
    adata.uns["spatial"][sample]["geometries"] is not automatically updated. This function
    ensures that the geometries GeoDataFrame matches the subsetted cells.
    
    Parameters
    ----------
    adata : sc.AnnData
        AnnData object that has been subsetted.
    sample : str, optional
        Sample name. If None, will use the first available sample in adata.uns["spatial"].
    
    Returns
    -------
    sc.AnnData
        The same AnnData object with synchronized geometries (modified in place).
    
    Examples
    --------
    >>> import trackcell.io as tcio
    >>> adata = tcio.read_hd_cellseg("path/to/data", sample="sample1")
    >>> # Subset the data
    >>> adata_subset = adata[adata.obs['classification'] == 'Cluster-1'].copy()
    >>> # Synchronize geometries
    >>> tcio.sync_geometries_after_subset(adata_subset, sample="sample1")
    >>> # Now adata_subset.uns["spatial"]["sample1"]["geometries"] only contains geometries for subsetted cells
    """
    if "spatial" not in adata.uns:
        warnings.warn(
            "No spatial information found in adata.uns['spatial']. "
            "This function is intended for data loaded with read_hd_cellseg()."
        )
        return adata
    
    # Determine sample name
    if sample is None:
        available_samples = list(adata.uns["spatial"].keys())
        if len(available_samples) == 0:
            warnings.warn("No samples found in adata.uns['spatial'].")
            return adata
        sample = available_samples[0]
        if len(available_samples) > 1:
            warnings.warn(
                f"Multiple samples found: {available_samples}. "
                f"Using '{sample}'. Specify `sample` explicitly to use a different one."
            )
    
    if sample not in adata.uns["spatial"]:
        warnings.warn(f"Sample '{sample}' not found in adata.uns['spatial'].")
        return adata
    
    spatial_info = adata.uns["spatial"][sample]
    
    # Check if geometries exist
    if "geometries" not in spatial_info:
        warnings.warn(
            f"No geometries found in adata.uns['spatial']['{sample}']['geometries']. "
            "This function is intended for data loaded with read_hd_cellseg()."
        )
        return adata
    
    geometries = spatial_info["geometries"]
    
    # Get the current cell IDs in the subsetted adata
    current_cell_ids = set(adata.obs_names)
    
    # Filter geometries to only include cells in the subset
    # Check if geometries is indexed by cell IDs
    if isinstance(geometries, gpd.GeoDataFrame):
        # Filter geometries based on current cell IDs
        # The geometries GeoDataFrame should be indexed by cell IDs (matching adata.obs_names)
        # This is how read_hd_cellseg stores them
        try:
            # Try to filter by index - this should work if geometries are indexed by cell IDs
            geometries_subset = geometries.loc[geometries.index.isin(current_cell_ids)]
        except (KeyError, IndexError) as e:
            # If index-based filtering fails, geometries might not be properly indexed
            # Try to match by intersecting indices
            common_indices = geometries.index.intersection(current_cell_ids)
            if len(common_indices) > 0:
                geometries_subset = geometries.loc[common_indices]
            else:
                warnings.warn(
                    f"Could not filter geometries by index. "
                    f"Geometries may not be properly indexed by cell IDs. "
                    f"Error: {e}"
                )
                return adata
        
        # Update the geometries in place
        spatial_info["geometries"] = geometries_subset
        
        # Note: adata.obs["geometry"] is automatically subset when adata is subset,
        # so we don't need to manually update it here
    else:
        warnings.warn(
            f"Unexpected geometry format in adata.uns['spatial']['{sample}']['geometries']. "
            "Expected GeoDataFrame."
        )
    
    return adata

