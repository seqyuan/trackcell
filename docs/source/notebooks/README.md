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
- Keep the executed outputs in the notebook: the docs build uses
  `nbsphinx_execute = 'never'`, so the committed outputs are what readers see
- Add a gallery thumbnail to `docs/source/_static/notebooks/` and register it in
  `nbsphinx_thumbnails` in `conf.py`
- Run the notebook end-to-end before committing so no cells contain tracebacks
- Test notebooks locally before adding to documentation

