"""
Spatial transcriptomics data input/output functions.

This module provides functions for reading and writing spatial transcriptomics data,
particularly from SpaceRanger output.
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
    gdf_seg = gpd.read_file(datapath / cell_segmentations_file)
    df = pd.DataFrame(gdf_seg)
    df['cellid'] = df['cell_id'].apply(lambda x: f"cellid_{str(x).zfill(9)}-1")
    
    # Read expression matrix
    adata = sc.read_10x_h5(datapath / matrix_file)
    
    # Align cell segmentations with expression data
    adata = adata[adata.obs_names.isin(df['cellid']),:]
    
    df = df.set_index("cellid").loc[adata.obs_names].reset_index()
    
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
    
    # Store images and scalefactors
    if hires_img is not None and lowres_img is not None:
        adata.uns["spatial"][sample]["images"] = {
            "hires": hires_img,
            "lowres": lowres_img
        }
    
    adata.uns["spatial"][sample]["scalefactors"] = scalefactor
    
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

