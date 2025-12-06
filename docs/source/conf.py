# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'TrackCell'
copyright = '2024, Zan Yuan'
author = 'Zan Yuan'
release = '0.2.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Base extensions
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
]

# Optional extensions
try:
    import nbsphinx
    extensions.append('nbsphinx')  # For Jupyter notebook support
except ImportError:
    pass  # nbsphinx not available, skip notebook support

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']

# PyData theme options
html_theme_options = {
    "github_url": "https://github.com/seqyuan/trackcell",
    "logo": {
        "text": "TrackCell",
    },
    "use_edit_page_button": True,
    "show_toc_level": 2,
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/seqyuan/trackcell",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/trackcell/",
            "icon": "fa-brands fa-python",
        },
    ],
}

html_context = {
    "github_user": "seqyuan",
    "github_repo": "trackcell",
    "github_version": "main",
    "doc_path": "docs/source",
}

# -- Extension configuration -------------------------------------------------

# Napoleon settings for Google/NumPy style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False

# Autodoc settings
autodoc_member_order = 'bysource'
autosummary_generate = True

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'scanpy': ('https://scanpy.readthedocs.io/en/stable/', None),
    'anndata': ('https://anndata.readthedocs.io/en/stable/', None),
}

# nbsphinx settings
nbsphinx_execute = 'never'  # Don't execute notebooks on build
nbsphinx_allow_errors = True

