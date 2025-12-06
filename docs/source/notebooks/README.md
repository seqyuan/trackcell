# Example Notebooks

Place Jupyter notebook (`.ipynb`) files here to include them in the documentation.

## Adding a Notebook

1. Create or save your notebook as `.ipynb` in this directory
2. Reference it in a RST file, for example in `examples.rst`:

```rst
.. nbsphinx:: notebooks/your_notebook.ipynb
```

Or add it to a toctree:

```rst
.. toctree::
   :maxdepth: 2

   notebooks/notebook1
   notebooks/notebook2
```

## Notebook Best Practices

- Keep notebooks focused on a single topic
- Add markdown cells to explain the code
- Ensure all outputs are cleared before committing (or set `nbsphinx_execute = 'never'` in conf.py)
- Test notebooks locally before adding to documentation

