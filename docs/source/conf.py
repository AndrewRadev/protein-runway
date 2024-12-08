# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

project = 'protein-runway'
copyright = '2024, Andrey Radev, Maria Krunic, Emma Rousseau, Evgeniya Polezhaeva'
author = 'Andrey Radev, Maria Krunic, Emma Rousseau, Evgeniya Polezhaeva'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
        'sphinx.ext.autodoc',  # For docstring generation
        'sphinx.ext.napoleon', # Support for Google/NumPy docstrings
        'sphinx.ext.viewcode', # Source code links
        'autoapi.extension',
    ]
# autosummary_generate = True  # Turn on sphinx.ext.autosummary
autoapi_dirs = ['../../lib']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
master_doc = 'index'
