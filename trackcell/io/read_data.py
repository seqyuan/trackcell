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
    
    adata.obsm['spatial'] = adata.obs[coord_cols].values
    
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
    
    # Create cellid in the format expected by SpaceRanger
    df['cellid'] = df['cell_id'].apply(lambda x: f"cellid_{str(x).zfill(9)}-1")
    
    # Read expression matrix
    adata = sc.read_10x_h5(datapath / matrix_file)
    
    # Align cell segmentations with expression data
    # Filter adata to only include cells that have segmentations
    adata = adata[adata.obs_names.isin(df['cellid']),:]
    
    # Filter and reorder df to match adata.obs_names order
    # Keep cellid as a column throughout (reset_index converts index back to column)
    df = df.set_index("cellid").loc[adata.obs_names]
    
    # Handle case where reset_index creates 'index' column instead of 'cellid'
    # This can happen if the index name was lost during operations
    """
    if "cellid" not in df.columns:
        if "index" in df.columns:
            # Rename 'index' to 'cellid' if it exists
            df = df.rename(columns={"index": "cellid"})
        else:
            raise ValueError(
                f"Unexpected: cellid is not a column after reset_index(). "
                f"Index name: {df.index.name}, Columns: {list(df.columns)}"
            )
    """
    # Convert geometry strings to shapely objects if needed
    if isinstance(df["geometry"].iloc[0], str):
        df["geometry"] = df["geometry"].apply(wkt.loads)
    
    # Extract centroid coordinates
    df["x"] = df["geometry"].apply(lambda poly: poly.centroid.x)
    df["y"] = df["geometry"].apply(lambda poly: poly.centroid.y)
    
    # Store spatial coordinates
    adata.obsm["spatial"] = np.array(df[["x", "y"]])
    
    if 'classification' in df.columns:
        if isinstance(df['classification'].iloc[0], str):
            classifications = df['classification'].apply(ast.literal_eval)
        else:
            classifications = df['classification']
        adata.obs['classification'] = [i['name'] for i in classifications]
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
    # Ensure data is in dictionary format (convert to dictionary if it's a string)
    classifications = df[classification_col]
    
    # Get unique classifications
    unique_classes = classifications.explode().unique()
    
    # Create color dictionary
    color_dict = {}
    for cls in unique_classes:
        if isinstance(cls, dict):  # Ensure it's in dictionary format
            name = cls['name']
            rgb = cls['color']
            # Convert RGB list to hexadecimal color code
            hex_color = '#{:02x}{:02x}{:02x}'.format(*rgb)
            color_dict[name] = hex_color
    
    return color_dict

