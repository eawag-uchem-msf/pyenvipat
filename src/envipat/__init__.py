"""envipat: Python port of the R package enviPat.

Isotope pattern, profile and centroid calculation for mass spectrometry.

This package houses the bundled reference datasets (see :mod:`envipat.data`)
and the compiled core extension (``envipat._core.enviPat_py``), which
implements the low-level isotope-pattern, profile and centroid algorithms.
The pure-Python pipeline (:mod:`envipat.isopattern`, :mod:`envipat.envelope`,
:mod:`envipat.vdetect`, :mod:`envipat.check_several`, :mod:`envipat.isowrap`)
wraps that extension with the same argument validation and result shaping as
the original R functions.
"""

__version__ = "0.1.0"
