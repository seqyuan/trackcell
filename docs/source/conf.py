# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'trackcell'
copyright = '2025, Zan Yuan'
author = 'Zan Yuan'
release = '0.3.6'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Source file suffixes
source_suffix = ['.rst', '.md']

# Base extensions
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx_gallery.load_style',  # For nbgallery directive support
]

# Optional extensions - nbsphinx for Jupyter notebook support
try:
    import nbsphinx
    extensions.append('nbsphinx')  # For Jupyter notebook support
    # Add IPython console highlighting support (like trackc project)
    try:
        from IPython.sphinxext import ipython_console_highlighting
        extensions.append('IPython.sphinxext.ipython_console_highlighting')
    except ImportError:
        pass  # IPython extension not available, skip
except ImportError:
    import warnings
    warnings.warn("nbsphinx not available, notebook support will be disabled")

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

# See warnings about bad links (like trackc project)
nitpicky = False  # Set to True for stricter link checking
nitpick_ignore = [
    ("", "Pygments lexer name 'ipython' is not known"),
    ("", "Pygments lexer name 'ipython3' is not known"),
]

# Suppress warnings for notebooks (docutils formatting issues are common in notebooks)
# Note: This suppresses warnings but not errors. For errors, we rely on nbsphinx_allow_errors = True
suppress_warnings = [
    'ref.docutils',  # Suppress docutils warnings in notebooks
    'misc.highlighting_failure',  # Suppress highlighting failures
    'misc.highlighting_failure.*',  # Suppress all highlighting failures
    'docutils.warnings.misc',  # Suppress docutils misc warnings
    'docutils',  # Suppress all docutils warnings
    'docutils.warnings.misc.UnexpectedIndentation',  # Suppress unexpected indentation warnings
    'docutils.warnings.misc.DefinitionListEndsWithoutBlankLine',  # Suppress definition list warnings
    'docutils.warnings.misc.BlockQuoteEndsWithoutBlankLine',  # Suppress block quote warnings
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
html_logo = '_static/images/logo_icon.svg'
html_favicon = '_static/images/logo_icon.svg'

# PyData theme options
html_theme_options = {
    "github_url": "https://github.com/seqyuan/trackcell",
    "logo": {
        "text": "TrackCell",
        "image_light": "_static/images/logo_icon.svg",
        "image_dark": "_static/images/logo_icon.svg",
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
nbsphinx_allow_errors = True  # Allow errors in notebooks
nbsphinx_timeout = 60  # Timeout for notebook execution (not used when execute='never')
# Suppress warnings about docutils formatting issues in notebooks
nbsphinx_requirejs_path = ''  # Disable requirejs if not needed
# Process raw cells as markdown (helps with formatting)
nbsphinx_process_signals = False

# nbsphinx thumbnails for gallery
# Key: notebook path relative to source directory (without .ipynb extension)
# Value: thumbnail image path relative to source directory
nbsphinx_thumbnails = {
    "notebooks/Colon_Cancer_demo": "_static/notebooks/Colon_Cancer_demo.png",
}

# Pygments configuration for ipython3 lexer
# IPython.sphinxext.ipython_console_highlighting handles this automatically
# If IPython extension is not available, nbsphinx will use python3 as fallback

# -- Options for LaTeX/PDF output -------------------------------------------------
# LaTeX configuration for Chinese character support
# Use XeLaTeX for Unicode/Chinese support (required for Chinese characters in PDF)
latex_engine = 'xelatex'

latex_elements = {
    # Use xeCJK for Chinese character support
    # This configuration supports Chinese characters in notebooks
    'preamble': r'''
\usepackage{fontspec}
\usepackage{xeCJK}
% Configure Chinese fonts with fallbacks
% Try common CJK fonts available on Linux systems (e.g., Read the Docs)
\IfFontExistsTF{Noto Sans CJK SC}{
    \setCJKmainfont{Noto Sans CJK SC}
    \setCJKsansfont{Noto Sans CJK SC}
    \setCJKmonofont{Noto Sans Mono CJK SC}
}{
    \IfFontExistsTF{Source Han Sans SC}{
        \setCJKmainfont{Source Han Sans SC}
        \setCJKsansfont{Source Han Sans SC}
        \setCJKmonofont{Source Han Mono SC}
    }{
        \IfFontExistsTF{WenQuanYi Micro Hei}{
            \setCJKmainfont{WenQuanYi Micro Hei}
            \setCJKsansfont{WenQuanYi Micro Hei}
            \setCJKmonofont{WenQuanYi Micro Hei Mono}
        }{
            % Fallback: try to use any available CJK font
            % If no CJK fonts are available, XeLaTeX will use a default font
            \setCJKmainfont{SimSun}
            \setCJKsansfont{SimHei}
            \setCJKmonofont{FangSong}
        }
    }
}
% Increase buffer size for long lines (fixes "Unable to read an entire line" error)
\maxdeadcycles=1000
% Increase input buffer size
\makeatletter
\maxdimen=16383pt
\makeatother
''',
    'maxlistdepth': '10',
    'pointsize': '10pt',
    'extraclassoptions': 'openany,oneside',
}

