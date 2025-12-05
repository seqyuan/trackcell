# TrackCell

A Python package for processing and vis single-cell and spatial transcriptomics data.

## Installation

```bash
pip install trackcell
```

## Usage

### Reading SpaceRanger Output

```python
import trackcell as tcl

# Read SpaceRanger output
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
```

#### visualisation
```python
sc.pl.spatial(adata, color='classification', size=2, 
              groups=['Cluster-2', 'Cluster-3'],
              legend_fontsize=12, spot_size=10, frameon=True
             )
```

```python
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
    group="classification",    # obs column containing cell type annotations
    label="Cluster-2",       # target label to measure distances from
    inplace=True          # add "{label}_px" and "{label}_dist" to adata.obs
)

# When inplace=False the function returns a DataFrame with the two columns:
dist_df = tcl.tl.hd_labeldist(adata, group="group_col", label="Neuron", inplace=False)
```

```python
sc.pl.spatial(adata, color='Cluster-2_dist', size=2, 
              legend_fontsize=12, spot_size=10, frameon=True
             )
```

## Development


## License



## update
```shell
git tag v0.1.5
git push origin v0.1.5

# In GitHub, go to “Releases” → “Draft a new release”.
```


