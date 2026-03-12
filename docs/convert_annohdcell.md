# Converting annohdcell Output to TrackCell Format

This guide explains how to convert annohdcell output files into trackcell-compatible format with polygon geometries for spatial visualization.

## Two Conversion Methods

### Method 1: Convert from 2μm Bin H5AD Only
Use `convert_annohdcell_to_trackcell()` to create a new cell-level h5ad from scratch.

### Method 2: Add Geometries to Existing Cell H5AD
Use `add_geometries_to_annohdcell_output()` to add geometries to annohdcell's final cell h5ad output.

---

## Method 1: Convert from 2μm Bin H5AD Only

### Background

**annohdcell Format**
- **Input**: 2μm bin-level h5ad (e.g., `b2c_2um.h5ad`)
- **Structure**: Each observation is a 2μm bin with cell assignment labels in `.obs["labels_joint"]`
- **Missing**: No polygon geometries for cells, preventing spatial visualization with trackcell

**trackcell Format**
- **Required**: Cell-level h5ad with polygon geometries
- **Structure**:
  - `.obs["geometry"]`: WKT strings representing cell boundaries
  - `.uns["spatial"][sample]["geometries"]`: GeoDataFrame with Shapely Polygon objects
  - `.obsm["spatial"]`: Cell centroid coordinates

### Usage

#### Basic Example

```python
import trackcell.io as tcio

# Convert annohdcell format to trackcell format
adata = tcio.convert_annohdcell_to_trackcell(
    bin_h5ad_path="sample_2um.h5ad",
    output_h5ad_path="sample_trackcell.h5ad",
    sample="sample1"
)

# Now visualize with trackcell
import trackcell.pl as tcpl
tcpl.spatial_cell(adata, sample="sample1")
```

#### Advanced Options

```python
adata = tcio.convert_annohdcell_to_trackcell(
    bin_h5ad_path="b2c_2um.h5ad",
    output_h5ad_path="trackcell_format.h5ad",
    sample="my_sample",
    labels_key="labels_joint",      # Column with cell labels (default)
    bin_size_um=2.0,                # Bin size in micrometers
    create_polygons=True,           # Create polygon geometries
    buffer_polygons=True            # Buffer polygons by bin size
)
```

---

## Method 2: Add Geometries to Existing Cell H5AD

### When to Use This Method

Use this method when you want to:
- **Preserve annohdcell's exact count aggregation** (which may differ from simple summation)
- **Keep all metadata** from annohdcell's final output
- **Add only geometries** without recalculating counts

### Usage

#### Basic Example

```python
import trackcell.io as tcio

# Add geometries to annohdcell's final cell h5ad
adata = tcio.add_geometries_to_annohdcell_output(
    bin_h5ad_path="b2c_2um.h5ad",      # 2μm bin h5ad with cell labels
    cell_h5ad_path="b2c_cell.h5ad",    # Final cell h5ad from annohdcell
    output_h5ad_path="b2c_cell_with_geom.h5ad",
    sample="sample1"
)

# Now visualize with trackcell
import trackcell.pl as tcpl
tcpl.spatial_cell(adata, sample="sample1")
```

#### Advanced Options

```python
adata = tcio.add_geometries_to_annohdcell_output(
    bin_h5ad_path="b2c_2um.h5ad",
    cell_h5ad_path="b2c_cell.h5ad",
    output_h5ad_path="b2c_cell_with_geom.h5ad",
    sample="my_sample",
    labels_key="labels_joint",      # Column in bin h5ad with cell labels
    bin_size_um=2.0,                # Bin size in micrometers
    buffer_polygons=True            # Buffer polygons by bin size
)
```

### What This Method Does

1. **Reads both h5ad files**: 2μm bin h5ad and final cell h5ad
2. **Matches cells**: Uses `object_id` in cell h5ad to match with labels in bin h5ad
3. **Creates polygons**: Generates convex hull polygons from bin coordinates for each cell
4. **Adds geometries**: Adds polygon geometries to the cell h5ad without modifying counts
5. **Preserves data**: Keeps all original data from annohdcell's cell h5ad

---

## Comparison of Methods

| Feature | Method 1: `convert_annohdcell_to_trackcell()` | Method 2: `add_geometries_to_annohdcell_output()` |
|---------|-----------------------------------------------|---------------------------------------------------|
| **Input** | 2μm bin h5ad only | 2μm bin h5ad + cell h5ad |
| **Count aggregation** | Simple sum across bins | Uses annohdcell's aggregation |
| **Metadata** | Basic metadata only | Preserves all annohdcell metadata |
| **Use case** | Quick conversion, simple workflow | Preserve exact annohdcell output |
| **Speed** | Faster (one file) | Slightly slower (two files) |

---

## Parameters

### Common Parameters (Both Methods)

- **bin_h5ad_path**: Path to input 2μm bin h5ad file from annohdcell
- **output_h5ad_path**: Path to save output h5ad (optional, can be None)
- **sample**: Sample name for spatial metadata (auto-detected if None)
- **labels_key**: Column in `.obs` with cell labels (default: "labels_joint")
- **bin_size_um**: Bin size in micrometers (default: 2.0)
- **buffer_polygons**: Whether to buffer polygons by bin size (default: True)

### Method 1 Specific Parameters

- **create_polygons**: Whether to create polygon geometries (default: True)

### Method 2 Specific Parameters

- **cell_h5ad_path**: Path to final cell h5ad file from annohdcell

---

## Output Format

The output h5ad file contains:

### `.X` - Expression Matrix
- **Method 1**: Summed counts across all bins in each cell
- **Method 2**: Original counts from annohdcell's cell h5ad
- Sparse CSR matrix format

### `.obs` - Cell Metadata
- `cellid`: Cell identifier (format: `cellid_000000001-1`)
- `object_id`: Original cell label from annohdcell
- `bin_count`: Number of bins comprising each cell
- `geometry`: WKT string representation of cell polygon
- QC metrics: `n_genes_by_counts`, `total_counts`, etc.

### `.obsm["spatial"]` - Spatial Coordinates
- Cell centroid coordinates (mean of constituent bin coordinates)
- Shape: (n_cells, 2)

### `.uns["spatial"][sample]` - Spatial Metadata
- `geometries`: GeoDataFrame with Shapely Polygon objects (indexed by cellid)
- `images`: Tissue images (hires, lowres) if present in input
- `scalefactors`: Scale factors including `microns_per_pixel`, `spot_diameter_fullres`

---

## Example Workflows

### Workflow 1: Quick Conversion (Method 1)

```python
# 1. Run annohdcell to get bin-level h5ad with cell labels
# (This step is done outside Python using annohdcell CLI)
# $ annohdcell bin2cell --input spaceranger_output --output results

# 2. Convert to trackcell format
import trackcell.io as tcio
adata = tcio.convert_annohdcell_to_trackcell(
    bin_h5ad_path="results/b2c_2um.h5ad",
    output_h5ad_path="results/trackcell_format.h5ad",
    sample="sample1"
)

# 3. Visualize with trackcell
import trackcell.pl as tcpl
tcpl.spatial_cell(adata, sample="sample1", show=True)
```

### Workflow 2: Preserve Exact annohdcell Output (Method 2)

```python
# 1. Run annohdcell to get both bin and cell h5ad files
# $ annohdcell bin2cell --input spaceranger_output --output results
# This produces: b2c_2um.h5ad, b2c_cell.h5ad, b2c4rds.h5ad

# 2. Add geometries to the final cell h5ad
import trackcell.io as tcio
adata = tcio.add_geometries_to_annohdcell_output(
    bin_h5ad_path="results/b2c_2um.h5ad",
    cell_h5ad_path="results/b2c_cell.h5ad",
    output_h5ad_path="results/b2c_cell_with_geom.h5ad",
    sample="sample1"
)

# 3. Visualize with trackcell
import trackcell.pl as tcpl
tcpl.spatial_cell(adata, sample="sample1", show=True)

# 4. Perform spatial analysis
import trackcell.tl as tctl
# ... use trackcell analysis functions
```

---

## Notes

- **Bins with label 0 are excluded** (unassigned bins)
- **Polygon creation**: Uses convex hull of bin coordinates
- **Single-bin cells**: Creates a small square polygon
- **Buffer option**: Expands polygons by half the bin size for better visualization
- **Coordinate system**: Preserves pixel coordinates from input h5ad
- **Scale factors**: Method 1 automatically adjusts `spot_diameter_fullres` based on mean bin count

---

## Troubleshooting

### Error: Column 'labels_joint' not found
- Check that your input h5ad has cell assignment labels
- Try specifying a different `labels_key` parameter
- Verify the h5ad is from annohdcell's bin2cell pipeline

### Error: 'spatial' not found in .obsm
- Ensure the input h5ad has spatial coordinates
- Check that the h5ad is from annohdcell (not a generic h5ad)

### Error: Column 'object_id' not found (Method 2 only)
- Ensure the cell h5ad is from annohdcell's bin_to_cell output
- Check that the cell h5ad has the `object_id` column in `.obs`

### Warning: Failed to create polygon for cell X
- Some cells may have problematic bin coordinates
- The function creates a fallback point geometry
- Check the bin coordinates in the input h5ad

### Warning: Cell X has no bins assigned
- Some cells in the cell h5ad may not have corresponding bins in the bin h5ad
- This can happen if the files are from different runs or versions
- These cells will not have geometries added

---

## See Also

- [trackcell documentation](../README.md)
- [annohdcell documentation](https://github.com/seqyuan/annohdcell)
- [Example script - Method 1](../examples/convert_annohdcell_example.py)
- [Example script - Method 2](../examples/add_geometries_to_annohdcell_example.py)

---

## Reading H5AD Files with Geometries

When you save an h5ad file with geometries using the methods above, the GeoDataFrame is converted to WKT strings for serialization. After reading the h5ad file back, you need to restore the GeoDataFrame for spatial visualization.

### Usage

```python
import scanpy as sc
import trackcell.io as tcio

# Read h5ad file
adata = sc.read_h5ad("b2c_cell_with_geom.h5ad")

# Restore GeoDataFrame from WKT strings
adata = tcio.restore_geometries(adata)

# Now you can visualize
import trackcell.pl as tcpl
tcpl.spatial_cell(adata, sample="sample1")
```

### What This Does

- Converts WKT strings in `uns["spatial"][sample]["geometries"]` back to a GeoDataFrame
- Processes all samples by default, or specify `sample="sample1"` for a specific sample
- Required for spatial visualization after reading h5ad files

### Example

```python
# After reading
adata = sc.read_h5ad("cell_with_geom.h5ad")
print(type(adata.uns["spatial"]["sample1"]["geometries"]))
# Output: <class 'pandas.core.frame.DataFrame'>

# After restoring
adata = tcio.restore_geometries(adata, sample="sample1")
print(type(adata.uns["spatial"]["sample1"]["geometries"]))
# Output: <class 'geopandas.geodataframe.GeoDataFrame'>
```
