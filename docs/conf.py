# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(0, os.path.abspath("../src"))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "paftacular"
copyright = "2026, Patrick Garrett"
author = "Patrick Garrett"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",  # For Google/NumPy style docstrings
    "sphinx_autodoc_typehints",  # Uses your type hints automatically
    "sphinx.ext.viewcode",  # Adds source code links
    "sphinx.ext.doctest",
]


# Optional: configure source suffixes
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# Autodoc settings
autodoc_default_options = {
    "members": True,
    "inherited-members": False,
    "show-inheritance": True,
}

autodoc_typehints = "description"  # Or 'signature' to put types in signature
napoleon_google_docstring = True
napoleon_numpy_docstring = True
