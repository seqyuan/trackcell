"""
Convert annohdcell bin2cell output to trackcell format.

This module provides functions to convert 2μm bin-level h5ad files from annohdcell
(with cell assignment labels) into trackcell-compatible h5ad files with cell-level
data and polygon geometries for spatial visualization.
"""

import scanpy as sc
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Polygon, MultiPoint
from shapely import wkt
from scipy.sparse import csr_matrix, issparse
from typing import Optional, Union
from pathlib import Path
import warnings
from multiprocessing import Pool, cpu_count
from functools import partial


def bins_to_cell_polygon(bin_coords: np.ndarray, bin_size_um: float = 2.0,
                         microns_per_pixel: Optional[float] = None) -> Polygon:
    """
    Create a polygon geometry from bin coordinates using convex hull.

    Parameters
    ----------
    bin_coords : np.ndarray
        Array of shape (n_bins, 2) containing [x, y] pixel coordinates of bins
    bin_size_um : float, default 2.0
        Size of each bin in micrometers
    microns_per_pixel : float, optional
        Conversion factor from micrometers to pixels. If provided, will expand
        the polygon slightly to account for bin size.

    Returns
    -------
    Polygon
        Shapely Polygon representing the cell boundary
    """
    if len(bin_coords) == 0:
        raise ValueError("Cannot create polygon from empty bin coordinates")

    if len(bin_coords) == 1:
        # Single bin - create a small square around it
        x, y = bin_coords[0]
        if microns_per_pixel is not None:
            half_size = (bin_size_um / microns_per_pixel) / 2
        else:
            half_size = 1.0  # Default pixel size
        return Polygon([
            (x - half_size, y - half_size),
            (x + half_size, y - half_size),
            (x + half_size, y + half_size),
            (x - half_size, y + half_size)
        ])

    # Multiple bins - use convex hull
    points = MultiPoint(bin_coords)
    polygon = points.convex_hull

    # Optionally buffer the polygon slightly to account for bin size
    if microns_per_pixel is not None and isinstance(polygon, Polygon):
        buffer_pixels = (bin_size_um / microns_per_pixel) / 2
        polygon = polygon.buffer(buffer_pixels)

    return polygon


def _create_polygon_for_cell(args):
    """
    Helper function for parallel polygon creation.

    Parameters
    ----------
    args : tuple
        (cell_id, bin_indices, spatial_coords, bin_size_um, microns_per_pixel)

    Returns
    -------
    tuple
        (cell_id, polygon)
    """
    cell_id, bin_indices, spatial_coords, bin_size_um, microns_per_pixel = args
    bin_coords = spatial_coords[bin_indices]
    polygon = bins_to_cell_polygon(bin_coords, bin_size_um, microns_per_pixel)
    return cell_id, polygon


def convert_annohdcell_to_trackcell(
    bin_h5ad_path: Union[str, Path],
    output_h5ad_path: Optional[Union[str, Path]] = None,
    sample: Optional[str] = None,
    labels_key: str = "labels_joint",
    bin_size_um: float = 2.0,
    create_polygons: bool = True,
    buffer_polygons: bool = True
) -> sc.AnnData:
    """
    Convert annohdcell 2μm bin h5ad to trackcell-compatible cell-level h5ad.

    This function reads a bin-level h5ad file from annohdcell (with cell assignment
    labels) and converts it to a cell-level h5ad file compatible with trackcell
    visualization tools. It aggregates bin counts into cells and creates polygon
    geometries for each cell.

    Parameters
    ----------
    bin_h5ad_path : str or Path
        Path to the input 2μm bin h5ad file (e.g., b2c_2um.h5ad from annohdcell)
    output_h5ad_path : str or Path, optional
        Path to save the output h5ad file. If None, will not save to disk.
    sample : str, optional
        Sample name for spatial metadata. If None, inferred from filename.
    labels_key : str, default "labels_joint"
        Column name in .obs containing cell assignment labels (0 = unassigned)
    bin_size_um : float, default 2.0
        Size of each bin in micrometers
    create_polygons : bool, default True
        Whether to create polygon geometries from bin coordinates
    buffer_polygons : bool, default True
        Whether to buffer polygons slightly to account for bin size

    Returns
    -------
    sc.AnnData
        Cell-level AnnData object with:
        - .X: Summed expression counts per cell
        - .obs: Cell metadata including "geometry" (WKT strings)
        - .obsm["spatial"]: Cell centroid coordinates
        - .uns["spatial"][sample]["geometries"]: GeoDataFrame with cell polygons
        - .uns["spatial"][sample]["images"]: Tissue images (if present in input)
        - .uns["spatial"][sample]["scalefactors"]: Scale factors

    Examples
    --------
    >>> import trackcell.io as tcio
    >>> adata = tcio.convert_annohdcell_to_trackcell(
    ...     "b2c_2um.h5ad",
    ...     output_h5ad_path="trackcell_format.h5ad",
    ...     sample="sample1"
    ... )
    >>> print(adata)
    >>> # Now can use trackcell visualization functions
    >>> import trackcell.pl as tcpl
    >>> tcpl.spatial_cell(adata, sample="sample1")

    Notes
    -----
    - Input h5ad must have .obs[labels_key] with integer cell labels (0 = unassigned)
    - Input h5ad must have .obsm["spatial"] with bin coordinates
    - Bins with label 0 are excluded (unassigned bins)
    - Gene counts are summed across all bins in each cell
    - Cell coordinates are the mean of constituent bin coordinates
    - Polygon geometries are created using convex hull of bin coordinates
    """

    # Read input h5ad
    bin_h5ad_path = Path(bin_h5ad_path)
    print(f"Reading bin-level h5ad from {bin_h5ad_path}")
    adata_bin = sc.read_h5ad(bin_h5ad_path)

    # Infer sample name if not provided
    if sample is None:
        sample = bin_h5ad_path.stem.replace("_2um", "").replace("b2c_", "")

    # Check required fields
    if labels_key not in adata_bin.obs.columns:
        raise ValueError(f"Column '{labels_key}' not found in .obs. Available columns: {adata_bin.obs.columns.tolist()}")

    if "spatial" not in adata_bin.obsm:
        raise ValueError("'spatial' not found in .obsm. Cannot extract bin coordinates.")

    # Filter out unassigned bins (label == 0)
    assigned_mask = adata_bin.obs[labels_key] != 0
    n_unassigned = (~assigned_mask).sum()
    print(f"Filtering out {n_unassigned} unassigned bins (label=0)")
    adata_bin = adata_bin[assigned_mask].copy()

    # Get unique cell labels
    cell_labels = adata_bin.obs[labels_key].unique()
    cell_labels = np.sort(cell_labels)
    n_cells = len(cell_labels)
    print(f"Found {n_cells} cells")

    # Get microns_per_pixel from scalefactors if available
    microns_per_pixel = None
    if "spatial" in adata_bin.uns:
        for key in adata_bin.uns["spatial"]:
            if "scalefactors" in adata_bin.uns["spatial"][key]:
                microns_per_pixel = adata_bin.uns["spatial"][key]["scalefactors"].get("microns_per_pixel")
                if microns_per_pixel is not None:
                    print(f"Using microns_per_pixel = {microns_per_pixel} from input h5ad")
                    break

    # Create cell-to-bin mapping matrix
    print("Creating cell-to-bin mapping...")
    bin_labels = adata_bin.obs[labels_key].values
    label_to_idx = {label: idx for idx, label in enumerate(cell_labels)}
    cell_indices = np.array([label_to_idx[label] for label in bin_labels])

    # Build sparse matrix: rows = cells, cols = bins
    n_bins = adata_bin.n_obs
    cell_to_bin = csr_matrix(
        (np.ones(n_bins), (cell_indices, np.arange(n_bins))),
        shape=(n_cells, n_bins)
    )

    # Aggregate expression matrix
    print("Aggregating expression counts...")
    if issparse(adata_bin.X):
        X_cell = cell_to_bin.dot(adata_bin.X)
    else:
        X_cell = cell_to_bin.dot(adata_bin.X)

    # Aggregate spatial coordinates (mean)
    spatial_bin = adata_bin.obsm["spatial"]
    spatial_cell = cell_to_bin.dot(spatial_bin) / cell_to_bin.sum(axis=1).A

    # Count bins per cell
    bin_counts = np.array(cell_to_bin.sum(axis=1)).flatten()

    # Create cell-level AnnData
    print("Creating cell-level AnnData...")
    adata_cell = sc.AnnData(
        X=X_cell,
        obs=pd.DataFrame({
            "cellid": [f"cellid_{int(label):09d}-1" for label in cell_labels],
            "object_id": cell_labels,
            "bin_count": bin_counts
        }),
        var=adata_bin.var.copy()
    )
    adata_cell.obs_names = adata_cell.obs["cellid"]

    # Add spatial coordinates
    adata_cell.obsm["spatial"] = spatial_cell

    # Create polygon geometries
    if create_polygons:
        print("Creating polygon geometries...")
        geometries = []
        geometry_wkt = []

        for label in cell_labels:
            # Get bins for this cell
            bin_mask = bin_labels == label
            bin_coords = spatial_bin[bin_mask]

            # Create polygon
            try:
                if buffer_polygons and microns_per_pixel is not None:
                    polygon = bins_to_cell_polygon(bin_coords, bin_size_um, microns_per_pixel)
                else:
                    polygon = bins_to_cell_polygon(bin_coords, bin_size_um, None)
                geometries.append(polygon)
                geometry_wkt.append(polygon.wkt)
            except Exception as e:
                warnings.warn(f"Failed to create polygon for cell {label}: {e}")
                # Create a point as fallback
                centroid = bin_coords.mean(axis=0)
                point = Polygon([(centroid[0], centroid[1])] * 3)  # Degenerate polygon
                geometries.append(point)
                geometry_wkt.append(point.wkt)

        # Add WKT strings to .obs
        adata_cell.obs["geometry"] = geometry_wkt

        # Create GeoDataFrame for .uns
        gdf = gpd.GeoDataFrame(
            geometry=geometries,
            index=adata_cell.obs_names
        )

        # Initialize .uns["spatial"] structure
        if "spatial" not in adata_cell.uns:
            adata_cell.uns["spatial"] = {}

        if sample not in adata_cell.uns["spatial"]:
            adata_cell.uns["spatial"][sample] = {}

        adata_cell.uns["spatial"][sample]["geometries"] = gdf
        print(f"Created {len(geometries)} polygon geometries")

    # Copy spatial metadata from input if available
    if "spatial" in adata_bin.uns:
        if "spatial" not in adata_cell.uns:
            adata_cell.uns["spatial"] = {}

        # Copy from first available sample in input
        for input_sample in adata_bin.uns["spatial"]:
            input_spatial = adata_bin.uns["spatial"][input_sample]

            if sample not in adata_cell.uns["spatial"]:
                adata_cell.uns["spatial"][sample] = {}

            # Copy images
            if "images" in input_spatial:
                adata_cell.uns["spatial"][sample]["images"] = input_spatial["images"].copy()

            # Copy scalefactors
            if "scalefactors" in input_spatial:
                adata_cell.uns["spatial"][sample]["scalefactors"] = input_spatial["scalefactors"].copy()

                # Update spot_diameter_fullres to reflect cell size
                if "spot_diameter_fullres" in adata_cell.uns["spatial"][sample]["scalefactors"]:
                    mean_bin_count = bin_counts.mean()
                    diameter_scale = np.sqrt(mean_bin_count)
                    original_diameter = adata_cell.uns["spatial"][sample]["scalefactors"]["spot_diameter_fullres"]
                    adata_cell.uns["spatial"][sample]["scalefactors"]["spot_diameter_fullres"] = original_diameter * diameter_scale
                    print(f"Scaled spot_diameter_fullres by {diameter_scale:.2f}x (mean bin_count = {mean_bin_count:.1f})")

            # Copy binsize if present
            if "binsize" in input_spatial:
                adata_cell.uns["spatial"][sample]["binsize"] = None  # Cell-level, not bin-level

            break  # Only copy from first sample

    # Calculate QC metrics
    print("Calculating QC metrics...")
    adata_cell.var_names_make_unique()
    sc.pp.calculate_qc_metrics(adata_cell, inplace=True)

    # Save if output path provided
    if output_h5ad_path is not None:
        output_h5ad_path = Path(output_h5ad_path)
        print(f"Saving to {output_h5ad_path}")
        adata_cell.write_h5ad(output_h5ad_path)

    print("Conversion complete!")
    print(f"Output: {adata_cell.n_obs} cells × {adata_cell.n_vars} genes")

    return adata_cell


def add_geometries_to_annohdcell_output(
    bin_h5ad_path: Union[str, Path],
    cell_h5ad_path: Union[str, Path],
    output_h5ad_path: Optional[Union[str, Path]] = None,
    sample: Optional[str] = None,
    labels_key: str = "labels_joint",
    bin_size_um: float = 2.0,
    buffer_polygons: bool = True,
    n_jobs: int = -1
) -> sc.AnnData:
    """
    Add polygon geometries to annohdcell's final cell h5ad using bin-level data.

    This function reads both the 2μm bin h5ad (with cell labels) and the final
    cell-level h5ad from annohdcell, then adds polygon geometries to the cell h5ad
    by creating convex hulls from the constituent bins of each cell.

    Parameters
    ----------
    bin_h5ad_path : str or Path
        Path to the 2μm bin h5ad file (e.g., b2c_2um.h5ad) with cell labels
    cell_h5ad_path : str or Path
        Path to the final cell h5ad file (e.g., b2c_cell.h5ad) from annohdcell
    output_h5ad_path : str or Path, optional
        Path to save the output h5ad file. If None, will not save to disk.
    sample : str, optional
        Sample name for spatial metadata. If None, inferred from filename.
    labels_key : str, default "labels_joint"
        Column name in bin h5ad .obs containing cell assignment labels
    bin_size_um : float, default 2.0
        Size of each bin in micrometers
    buffer_polygons : bool, default True
        Whether to buffer polygons slightly to account for bin size
    n_jobs : int, default -1
        Number of parallel workers for polygon creation.
        -1 uses all available CPU cores, 1 disables parallelization.

    Returns
    -------
    sc.AnnData
        Cell-level AnnData object (copy of input cell h5ad) with added:
        - .obs["geometry"]: WKT strings representing cell boundaries
        - .uns["spatial"][sample]["geometries"]: GeoDataFrame with cell polygons

    Examples
    --------
    >>> import trackcell.io as tcio
    >>> adata = tcio.add_geometries_to_annohdcell_output(
    ...     bin_h5ad_path="b2c_2um.h5ad",
    ...     cell_h5ad_path="b2c_cell.h5ad",
    ...     output_h5ad_path="b2c_cell_with_geom.h5ad",
    ...     sample="sample1"
    ... )
    >>> # Now can use trackcell visualization
    >>> import trackcell.pl as tcpl
    >>> tcpl.spatial_cell(adata, sample="sample1")

    Notes
    -----
    - The cell h5ad must have .obs["object_id"] matching labels in bin h5ad
    - Bin h5ad must have .obs[labels_key] with integer cell labels
    - Bin h5ad must have .obsm["spatial"] with bin coordinates
    - The function preserves all data from the input cell h5ad
    - Only adds geometry information, does not modify counts or other data
    """

    # Read both h5ad files
    bin_h5ad_path = Path(bin_h5ad_path)
    cell_h5ad_path = Path(cell_h5ad_path)

    print(f"Reading bin-level h5ad from {bin_h5ad_path}")
    adata_bin = sc.read_h5ad(bin_h5ad_path)

    print(f"Reading cell-level h5ad from {cell_h5ad_path}")
    adata_cell = sc.read_h5ad(cell_h5ad_path)

    # Infer sample name if not provided
    if sample is None:
        sample = cell_h5ad_path.stem.replace("_cell", "").replace("b2c_", "")

    # Check required fields in bin h5ad
    if labels_key not in adata_bin.obs.columns:
        raise ValueError(f"Column '{labels_key}' not found in bin h5ad .obs. Available: {adata_bin.obs.columns.tolist()}")

    if "spatial" not in adata_bin.obsm:
        raise ValueError("'spatial' not found in bin h5ad .obsm")

    # Check required fields in cell h5ad
    if "object_id" not in adata_cell.obs.columns:
        raise ValueError("Column 'object_id' not found in cell h5ad .obs. This should be the cell label ID.")

    # Get microns_per_pixel from scalefactors
    microns_per_pixel = None
    if "spatial" in adata_bin.uns:
        for key in adata_bin.uns["spatial"]:
            if "scalefactors" in adata_bin.uns["spatial"][key]:
                microns_per_pixel = adata_bin.uns["spatial"][key]["scalefactors"].get("microns_per_pixel")
                if microns_per_pixel is not None:
                    print(f"Using microns_per_pixel = {microns_per_pixel}")
                    break

    # Filter out unassigned bins (label == 0)
    assigned_mask = adata_bin.obs[labels_key] != 0
    adata_bin_filtered = adata_bin[assigned_mask].copy()

    # Get bin labels and coordinates
    bin_labels = adata_bin_filtered.obs[labels_key].values
    spatial_bin = adata_bin_filtered.obsm["spatial"]

    # Get cell labels from cell h5ad
    cell_labels = adata_cell.obs["object_id"].values
    n_cells = len(cell_labels)
    print(f"Processing {n_cells} cells")

    # Create polygon geometries for each cell
    print("Creating polygon geometries from bin coordinates...")

    # Build a mapping from cell_label to bin indices for vectorized processing
    print("Building cell-to-bin mapping...")
    cell_to_bins = {}
    for bin_idx, label in enumerate(bin_labels):
        if label not in cell_to_bins:
            cell_to_bins[label] = []
        cell_to_bins[label].append(bin_idx)

    # Convert lists to numpy arrays for faster indexing
    for label in cell_to_bins:
        cell_to_bins[label] = np.array(cell_to_bins[label], dtype=np.int32)

    # Prepare arguments for parallel processing
    args_list = []
    cell_label_to_obs_idx = {}  # Map cell label to obs index

    for idx, label in enumerate(cell_labels):
        cell_label_to_obs_idx[label] = idx
        if label not in cell_to_bins:
            warnings.warn(f"Cell {label} has no bins assigned, skipping geometry creation")
            continue

        bin_indices = cell_to_bins[label]
        mp = microns_per_pixel if buffer_polygons else None
        args_list.append((label, bin_indices, spatial_bin, bin_size_um, mp))

    print(f"Processing {len(args_list)} cells with bins...")

    # Process in parallel or serial depending on n_jobs
    if n_jobs == 1:
        # Serial processing
        results = []
        for args in args_list:
            results.append(_create_polygon_for_cell(args))
    else:
        # Parallel processing
        actual_jobs = min(n_jobs if n_jobs > 0 else cpu_count(), len(args_list))
        print(f"Using {actual_jobs} parallel workers")
        with Pool(processes=actual_jobs) as pool:
            results = pool.map(_create_polygon_for_cell, args_list)

    print(f"Created {len(results)} polygon geometries")

    # Unpack results - directly store Shapely objects, no WKT conversion
    geometries = []
    cell_ids_with_geom = []
    geometry_wkt = []

    for cell_label, polygon in results:
        obs_idx = cell_label_to_obs_idx[cell_label]
        cell_ids_with_geom.append(adata_cell.obs_names[obs_idx])
        geometries.append(polygon)
        geometry_wkt.append(polygon.wkt)  # Only convert to WKT for .obs storage

    # Add WKT strings to .obs
    # Initialize with None for all cells
    adata_cell.obs["geometry"] = None
    # Fill in geometries for cells that have them
    for cell_id, wkt_str in zip(cell_ids_with_geom, geometry_wkt):
        adata_cell.obs.loc[cell_id, "geometry"] = wkt_str

    # Create GeoDataFrame directly from Shapely objects (no WKT round-trip)
    gdf = gpd.GeoDataFrame(
        geometry=geometries,
        index=cell_ids_with_geom
    )

    # Initialize .uns["spatial"] structure if needed
    if "spatial" not in adata_cell.uns:
        adata_cell.uns["spatial"] = {}

    if sample not in adata_cell.uns["spatial"]:
        adata_cell.uns["spatial"][sample] = {}

    # Add geometries to .uns
    adata_cell.uns["spatial"][sample]["geometries"] = gdf

    # Copy spatial metadata from bin h5ad if not already present in cell h5ad
    if "spatial" in adata_bin.uns:
        for input_sample in adata_bin.uns["spatial"]:
            input_spatial = adata_bin.uns["spatial"][input_sample]

            # Copy images if not present
            if "images" in input_spatial and "images" not in adata_cell.uns["spatial"][sample]:
                adata_cell.uns["spatial"][sample]["images"] = input_spatial["images"].copy()
                print("Copied tissue images from bin h5ad")

            # Copy scalefactors if not present
            if "scalefactors" in input_spatial and "scalefactors" not in adata_cell.uns["spatial"][sample]:
                adata_cell.uns["spatial"][sample]["scalefactors"] = input_spatial["scalefactors"].copy()
                print("Copied scalefactors from bin h5ad")

            break

    # Save if output path provided
    if output_h5ad_path is not None:
        output_h5ad_path = Path(output_h5ad_path)
        print(f"Saving to {output_h5ad_path}")
        adata_cell.write_h5ad(output_h5ad_path)

    print("Geometry addition complete!")
    print(f"Output: {adata_cell.n_obs} cells × {adata_cell.n_vars} genes with geometries")

    return adata_cell
