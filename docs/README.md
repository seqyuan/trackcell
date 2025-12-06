# TrackCell Documentation

This directory contains the Sphinx documentation for TrackCell.

## Building Documentation Locally

1. Install documentation dependencies:

```bash
pip install -e ".[docs]"
```

2. Build the documentation:

```bash
cd docs
make html
```

The built documentation will be in `docs/build/html/`.

## Structure

- `source/` - Source RST files and configuration
- `build/` - Generated documentation (not in git)
- `requirements.txt` - Documentation build dependencies

## Adding Jupyter Notebooks

To include Jupyter notebooks in the documentation:

1. Place `.ipynb` files in `source/` or a subdirectory
2. Reference them in RST files using:

```rst
.. nbsphinx:: path/to/notebook.ipynb
```

Or include them in a toctree:

```rst
.. toctree::
   :maxdepth: 2

   notebooks/example1
   notebooks/example2
```

