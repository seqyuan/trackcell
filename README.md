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
    path="SpaceRanger4.0/Cse1/outs/segmented_outputs",
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


## Development


## License

