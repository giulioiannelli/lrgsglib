# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from importlib.metadata import version as get_version
from pathlib import Path

# -- Path setup --------------------------------------------------------------
# Add the project root to sys.path to allow autodoc to find modules
sys.path.insert(0, str(Path(__file__).parents[2] / 'src'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'lrgsglib'
copyright = '2025, Giulio Iannelli, Pablo Villegas'
author = 'Giulio Iannelli, Pablo Villegas'
try:
    release = get_version('lrgsglib')
except Exception:
    release = '0.1.0'
version = '.'.join(release.split('.')[:2])

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    # Core Sphinx extensions
    'sphinx.ext.autodoc',           # Auto-generate docs from docstrings
    'sphinx.ext.autosummary',       # Generate summary tables
    'sphinx.ext.doctest',           # Test code snippets in docs
    'sphinx.ext.napoleon',          # Support NumPy/Google style docstrings
    'sphinx.ext.intersphinx',       # Link to other project docs
    'sphinx.ext.viewcode',          # Add links to source code
    'sphinx.ext.mathjax',           # Render math equations
    'sphinx.ext.todo',              # Support TODO items
    'sphinx.ext.coverage',          # Check documentation coverage

    # Third-party extensions
    'sphinx_rtd_theme',             # ReadTheDocs theme
    'sphinx_copybutton',            # Copy button for code blocks
    'myst_parser',                  # Markdown support
    'sphinxcontrib.bibtex',         # Bibliography support
]

# Autosummary settings
autosummary_generate = True
autosummary_imported_members = False

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__',
    'show-inheritance': True,
}
autodoc_typehints = 'description'
autodoc_typehints_description_target = 'documented'

# Suppress warnings
suppress_warnings = [
    'autodoc.duplicate',  # Suppress duplicate object description warnings
]

# Napoleon settings (NumPy style docstrings)
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# Intersphinx mapping (links to other projects)
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'networkx': ('https://networkx.org/documentation/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
}

# Doctest configuration
doctest_default_flags = (
    0
    | __import__('doctest').ELLIPSIS
    | __import__('doctest').NORMALIZE_WHITESPACE
)
doctest_global_setup = """
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from lrgsglib import *
"""

# Bibliography (add .bib files to docs/source/ when references are needed)
bibtex_bibfiles = []

# Templates path
templates_path = ['_templates']

# List of patterns to ignore when looking for source files
exclude_patterns = []

# Source file suffix
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# The master toctree document
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Theme options
html_theme_options = {
    'navigation_depth': 4,
    'collapse_navigation': False,
    'sticky_navigation': True,
    'includehidden': True,
    'titles_only': False,
    'logo_only': False,
}

# HTML context
html_context = {
    'display_github': True,
    'github_user': 'giulioiannelli',
    'github_repo': 'lrgsglib',
    'github_version': 'main',
    'conf_py_path': '/docs/source/',
}

# Custom sidebar
html_sidebars = {
    '**': [
        'globaltoc.html',
        'relations.html',
        'sourcelink.html',
        'searchbox.html',
    ]
}

# Add any paths that contain custom static files
html_logo = None  # Can add logo later
html_favicon = None  # Can add favicon later

# Output file base name for HTML help builder
htmlhelp_basename = 'lrgsglibdoc'

# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    'papersize': 'a4paper',
    'pointsize': '10pt',
}

# Grouping the document tree into LaTeX files
latex_documents = [
    (master_doc, 'lrgsglib.tex', 'lrgsglib Documentation',
     'Giulio Iannelli, Pablo Villegas', 'manual'),
]

# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'lrgsglib', 'lrgsglib Documentation',
     [author], 1)
]

# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files
texinfo_documents = [
    (master_doc, 'lrgsglib', 'lrgsglib Documentation',
     author, 'lrgsglib', 'Laplacian Renormalization Group for Signed Graphs',
     'Miscellaneous'),
]

# -- Extension configuration -------------------------------------------------

# Todo extension configuration
todo_include_todos = True

# Copybutton configuration
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
