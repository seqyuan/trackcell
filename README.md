# TrackCell

A Python package for processing and vis single-cell and spatial transcriptomics data.

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

## Development


## License



## update
```shell
git tag v0.2.0
git push origin v0.2.0

# In GitHub, go to “Releases” → “Draft a new release”.
```


