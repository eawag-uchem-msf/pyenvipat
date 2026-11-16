"""envipat: Python port of the R package enviPat.

Isotope pattern, profile and centroid calculation for mass spectrometry.

This package houses the bundled reference datasets (see :mod:`envipat.data`)
and the compiled core extension (``envipat._core.enviPat_py``), which
implements the low-level isotope-pattern, profile and centroid algorithms.
The pure-Python pipeline (:mod:`envipat.isopattern`, :mod:`envipat.envelope`,
:mod:`envipat.vdetect`, :mod:`envipat.check_several`, :mod:`envipat.isowrap`)
wraps that extension with the same argument validation and result shaping as
the original R functions.

The public pipeline functions are re-exported at the top level, so callers
can write::

    import envipat as ep

    checked = ep.check_chemform(isotopes, "C6H12O6")
    pattern = ep.isopattern(isotopes, checked["new_formula"])
    profiles = ep.envelope(pattern, resolution=resolution)
    centroids = ep.vdetect(profiles)

instead of importing each function from its own submodule. The bundled
datasets remain a submodule import (``import envipat.data as data``), since
they are data, not callables, and importing all of them eagerly at package
import time would defeat :mod:`envipat.data`'s own lazy-loading design (see
its module docstring).
"""

from .check_several import check_several
from .envelope import envelope
from .formula import check_chemform, check_ded, mergeform, multiform, subform
from .isopattern import isopattern
from .isowrap import isowrap
from .resolution import get_resolution
from .vdetect import vdetect

# The only place the version is defined: pyproject.toml reads it from here, and the docs
# show it.
__version__ = "0.2.0"

__all__ = [
    "__version__",
    "check_chemform",
    "check_ded",
    "check_several",
    "envelope",
    "get_resolution",
    "isopattern",
    "isowrap",
    "mergeform",
    "multiform",
    "subform",
    "vdetect",
]
