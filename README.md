# TrackCell

A Python package for processing and visualizing single-cell and spatial transcriptomics data.
Supports 10x Visium HD (SpaceRanger) and 10x Xenium Analyzer outputs.

## Installation

```bash
pip install trackcell -i https://pypi.org/simple
# pip install --upgrade trackcell==0.1.9 -i https://pypi.org/simple
```

## Usage

### Reading SpaceRanger Output

#### Reading Cell Segmentation Data

```python
import trackcell as tcl

# Read SpaceRanger cell segmentation output
adata = tcl.io.read_hd_cellseg(
    datapath="SpaceRanger4.0/Cse1/outs/segmented_outputs",
    sample="Cse1"
)

# The resulting AnnData object contains:
# - Expression matrix in .X
# - Cell metadata in .obs
# - Gene metadata in .var
# - Spatial coordinates in .obsm["spatial"]
# - Tissue images in .uns["spatial"][sample]["images"]
# - Scalefactors in .uns["spatial"][sample]["scalefactors"]
# - Cell geometries in .uns["spatial"][sample]["geometries"] (GeoDataFrame)
# - Cell geometries in .obs["geometry"] (WKT strings for serialization)
```

#### Subsetting Data and Synchronizing Geometries

**Important**: When you subset data loaded with `read_hd_cellseg()`, you must call `sync_geometries_after_subset()` to synchronize the geometries:

```python
import trackcell as tcl
import numpy as np

# Read data
adata = tcl.io.read_hd_cellseg(
    datapath="SpaceRanger4.0/Cse1/outs/segmented_outputs",
    sample="Cse1"
)

# Subset by spatial region
x_min, x_max = 16000, 18000
y_min, y_max = 14000, 18000

spatial_coords = adata.obsm['spatial']
mask = ((spatial_coords[:, 0] >= x_min) & (spatial_coords[:, 0] <= x_max) &
        (spatial_coords[:, 1] >= y_min) & (spatial_coords[:, 1] <= y_max))

adata_subset = adata[mask].copy()

# IMPORTANT: Synchronize geometries after subsetting
tcl.io.sync_geometries_after_subset(adata_subset, sample="Cse1")

# Now you can safely plot the subset
tcl.pl.spatial_cell(adata_subset, color="classification")
```

**Why this is necessary**: When you subset an AnnData object, `adata.obs` and `adata.obsm` are automatically subset, but `adata.uns["spatial"][sample]["geometries"]` (GeoDataFrame) is not. Without synchronization, plotting may fail with errors like `ValueError: aspect must be finite and positive`.

#### Reading Xenium Output

```python
import trackcell as tcl

# Read Xenium Analyzer cell segmentation output
adata = tcl.io.read_xenium_cellseg(
    datapath="/path/to/xenium/output",
    sample="sample1"
)

# The resulting AnnData object contains:
# - Expression matrix in .X (CSR, cells × genes)
# - Cell metadata in .obs (centroids, area, counts, segmentation_method, etc.)
# - Gene metadata in .var (gene_ids, feature_types)
# - Spatial coordinates in .obsm["spatial"]
# - Cell polygons in .uns["spatial"][sample]["geometries"] (GeoDataFrame)
# - Cell boundary arrays in .uns["cell_boundaries"] (compact vertex arrays)
# - Nucleus boundary arrays in .uns["nucleus_boundaries"] (if available)
# - WKT geometry strings in .obs["geometry"] (for serialization)
# - Experiment metadata in .uns["experiment"]
```

The function reads Xenium Analyzer output files:
| File | Format | Content |
|------|--------|---------|
| `cell_feature_matrix.h5` | HDF5 (10x) | Expression matrix (genes × cells CSC) |
| `cells.parquet` | Parquet | Cell metadata (centroids, area, counts) |
| `cell_boundaries.parquet` | Parquet long-table | Cell boundary vertices |
| `nucleus_boundaries.parquet` | Parquet long-table | Nucleus boundary vertices (optional) |
| `experiment.xenium` | JSON | Experiment metadata |
| `gene_panel.json` | JSON | Gene panel information |

**Key differences from Visium HD**:
- Cell boundaries are stored in **long-table parquet** format (one row per vertex), automatically converted to Shapely polygons
- Expression matrix is in **CSC (genes × cells)** format, automatically transposed to CSR (cells × genes)
- Cell IDs already include the `-1` suffix (no stripping needed)

After loading, the AnnData object is fully compatible with `tcl.pl.spatial_cell()`:

```python
# Visualize Xenium data with cell polygons
tcl.pl.spatial_cell(adata, color="EPCAM", cmap="Reds", edges_width=0.3)
```

#### Reading Bin-Level Data (2um/8um/16um)

```python
import trackcell as tcl

# Read SpaceRanger bin-level output (2um/8um/16um bins)
adata = tcl.io.read_hd_bin(
    datapath="SpaceRanger4.0/Cse1/binned_outputs",
    sample="Cse1",
    binsize=16  # Bin size in micrometers (default: 16, common values: 2, 8, or 16)
)

# The function automatically handles:
# - filtered_feature_bc_matrix.h5 (preferred) or filtered_feature_bc_matrix/ directory
# - tissue_positions.parquet or tissue_positions.csv
# - tissue_hires_image.png and tissue_lowres_image.png
# - scalefactors_json.json

# The resulting AnnData object contains:
# - Expression matrix in .X
# - Bin metadata in .obs (with spatial coordinates)
# - Gene metadata in .var
# - Spatial coordinates in .obsm["spatial"]
# - Tissue images in .uns["spatial"][sample]["images"]
# - Scalefactors in .uns["spatial"][sample]["scalefactors"]
# - Bin size in .uns["spatial"][sample]["binsize"] (e.g., 2, 8, or 16)

# Access the bin size information:
print(f"Bin size: {adata.uns['spatial']['Cse1']['binsize']} um")
```

#### Visualization

##### Plotting with Cell Polygons

```python
# Plot cells as polygons (requires data loaded with read_hd_cellseg)
tcl.pl.spatial_cell(
    adata, 
    color="classification",  # Color by cell type
    groups=['Cluster-2', 'Cluster-3'],  # Optional: filter specific groups
    figsize=(10, 10),
    edges_width=0.5,
    edges_color="black",
    alpha=0.8
)
```

```python
# Plot continuous values (e.g., distance to a label)
tcl.pl.spatial_cell(
    adata,
    color="Cluster-2_dist",  # Distance to Cluster-2
    cmap="Reds",
    figsize=(10, 10)
)
```

##### Dual-Color Visualization (Fill + Edge)

Use `edge_color` to color cell boundaries by a categorical column (e.g., cell type),
while `color` controls the fill (gene expression or continuous value):

```python
# Fill = gene expression, Edge = cell type
tcl.pl.spatial_cell(
    adata,
    color='EPCAM',              # Fill: gene expression (continuous)
    cmap='Reds',
    edge_color='cell_type',     # Edge: cell type (categorical)
    edge_palette={
        'T cell': '#e41a1c',
        'B cell': '#377eb8',
        'Myeloid': '#4daf4a',
    },
    edges_width=1.2,
    alpha=0.7,
)
```

```python
# Fill = continuous obs, Edge = categorical obs
tcl.pl.spatial_cell(
    adata,
    color='total_counts',       # Fill: UMIs per cell
    cmap='YlOrRd',
    edge_color='cell_type',     # Edge: cell type
    edges_width=1.5,
    alpha=0.7,
)
```

##### Traditional Point-based Visualization

```python
# Using scanpy (point-based)
sc.pl.spatial(adata, color='classification', size=2, 
              groups=['Cluster-2', 'Cluster-3'],
              legend_fontsize=12, spot_size=10, frameon=True
             )
```

```python
# Using squidpy (point-based)
sq.pl.spatial_scatter(
    adata, shape=None, color=["classification"], 
    edges_width=0, size=0.1, 
    library_id="spatial", 
    groups=['Cluster-2', 'Cluster-3'],
    figsize=(5, 4), 
    #cmap='Blues'
    #palette = mycolor
    #img_key="0.3_mpp_150_buffer", 
    #basis="spatial_cropped_150_buffer"
)
```

### Converting annohdcell Output to TrackCell Format

TrackCell provides two methods to convert annohdcell's bin2cell output into trackcell-compatible format with polygon geometries for spatial visualization.

#### Method 1: Convert from 2μm Bin H5AD Only

Create a new cell-level h5ad from annohdcell's 2μm bin h5ad with cell labels:

```python
import trackcell as tcl

# Convert annohdcell 2μm bin h5ad to trackcell format
adata = tcl.io.convert_annohdcell_to_trackcell(
    bin_h5ad_path="b2c_2um.h5ad",
    output_h5ad_path="trackcell_format.h5ad",
    sample="sample1"
)

# Now visualize with trackcell
tcl.pl.spatial_cell(adata, sample="sample1")
```

#### Method 2: Add Geometries to Existing Cell H5AD

Add polygon geometries to annohdcell's final cell h5ad output (preserves exact count aggregation):

```python
import trackcell as tcl

# Add geometries to annohdcell's final cell h5ad
adata = tcl.io.add_geometries_to_annohdcell_output(
    bin_h5ad_path="b2c_2um.h5ad",      # 2μm bin h5ad with cell labels
    cell_h5ad_path="b2c_cell.h5ad",    # Final cell h5ad from annohdcell
    output_h5ad_path="b2c_cell_with_geom.h5ad",
    sample="sample1"
)

# Now visualize with trackcell
tcl.pl.spatial_cell(adata, sample="sample1")
```

**Key differences:**
- **Method 1**: Quick conversion, simple count summation
- **Method 2**: Preserves annohdcell's exact count aggregation and all metadata

For detailed documentation, see [docs/convert_annohdcell.md](docs/convert_annohdcell.md)

### Computing Distances to a Label (10x HD)

```python
# Compute distance to a specific annotation label stored in adata.obs["group_col"]
tcl.tl.hd_labeldist(
    adata,
    groupby="classification",    # obs column containing cell type annotations
    label="Cluster-2",       # target label to measure distances from
    inplace=True          # add "{label}_px" and "{label}_dist" to adata.obs
)

# When inplace=False the function returns a DataFrame with the two columns:
dist_df = tcl.tl.hd_labeldist(adata, groupby="group_col", label="Neuron", inplace=False)
```

```python
# Visualize distance using cell polygons
tcl.pl.spatial_cell(adata, color='Cluster-2_dist', cmap='Reds', figsize=(10, 10))

# Or using traditional point-based visualization
sc.pl.spatial(adata, color='Cluster-2_dist', size=2,
              legend_fontsize=12, spot_size=10, frameon=True
             )
```

### Multi-Gene Visualization (cell2location-style)

Visualize co-expression of multiple genes in two modes:

#### Mode 1: Blended Composite (`mode='blend'`)

Maps each gene to a base color, blends by expression level → single hex color per cell.
Use with `spatial_cell` for a composite view:

```python
# Compute blended colors
adata = tcl.tl.multigene_blend(
    adata,
    genes=['EPCAM', 'PECAM1', 'VWF'],
    mode='blend',
)

# Visualize with cell polygons
tcl.pl.spatial_cell(adata, color='multigene_blend')
```

```python
# Custom colors + gamma correction
tcl.tl.multigene_blend(
    adata,
    genes=['EPCAM', 'PECAM1', 'VWF', 'ACTA2', 'PTPRC'],
    colors=['#e41a1c', '#377eb8', '#4daf4a', '#ff7f00', '#984ea3'],
    vmax_percentile=98, gamma=0.8,
)
```

#### Mode 2: Faceted Subplots (`mode='facet'`)

Each gene gets its own panel with a single-hue colormap (white → gene color),
matching the cell2location paper style:

```python
tcl.tl.multigene_blend(
    adata,
    genes=['EPCAM', 'PECAM1', 'VWF', 'ACTA2'],
    mode='facet',
    ncols=2,           # 2 columns of subplots
    edges_width=0.3,
)
```

## Development


## License



## update
```shell
git tag v0.3.8
git push origin v0.3.8

# In GitHub, go to “Releases” → “Draft a new release”.
```





