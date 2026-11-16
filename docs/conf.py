"""Sphinx configuration for the pyEnvipat documentation.

Build with ``nox -s docs`` (what CI runs) or directly with
``sphinx-build -W --keep-going docs docs/_build/html`` after
``pip install -e ".[docs]"``.
"""

from __future__ import annotations

import envipat

project = "pyEnvipat"
author = "Johannes Boog, Martin Loos, Christian Gerber"
copyright = "Eawag; GPL-2.0-only"
version = release = envipat.__version__

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "numpydoc",
    "myst_nb",
]

source_suffix = {".rst": "restructuredtext", ".md": "myst-nb", ".ipynb": "myst-nb"}
exclude_patterns = ["_build", "**.ipynb_checkpoints"]

# --- API reference ---------------------------------------------------------
autodoc_default_options = {"members": True, "member-order": "bysource"}
# Keep signatures readable: the docstrings describe types in prose.
autodoc_typehints = "none"
# Functions document their parameters; numpydoc only needs the Parameters/Returns
# tables, not an extra members table per module page.
numpydoc_show_class_members = False
numpydoc_class_members_toctree = False
# Docstrings name types in prose ("dict[str, pandas.DataFrame | str]"), which are
# not resolvable cross-reference targets.
numpydoc_xref_param_type = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

# --- Notebooks (myst-nb) ---------------------------------------------------
# Notebooks are stored without outputs and executed on every build; any cell
# error fails the build, so the examples stay runnable.
nb_execution_mode = "force"
nb_execution_raise_on_error = True
nb_execution_timeout = 120
nb_merge_streams = True

# Fail on every unresolved cross-reference, not just on broken markup.
nitpicky = True
myst_enable_extensions = ["colon_fence"]

# --- HTML ------------------------------------------------------------------
html_theme = "furo"
html_title = f"pyEnvipat {release}"
