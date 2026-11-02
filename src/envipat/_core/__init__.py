"""envipat._core: thin CPython/NumPy C-extension wrapper around iso_core.*.

The compiled extension module lives at ``envipat._core.enviPat_py`` (built
from ``enviPat_py.c`` in this directory, plus the vendored, unmodified
algorithm sources -- ``iso_core.c``, ``combination.c``, ``profile.c``,
``parse.c``, ``element.c``, ``isotope.c``, ``peak.c`` -- against the same
pure-C core the R package's ``.Call`` entry points delegate to). It exposes
the 5 core operations 1:1 with ``iso_core.h``, with NumPy arrays in/out and
no attempt to replicate the R wrappers' SEXP-shaping logic (e.g. the
``rica``/``return_iso_calc_amount`` branch in ``main.c``'s ``iso_pattern``
wrapper) -- that belongs to a future, higher-level Python port of the R
wrapper functions themselves.

This subpackage intentionally has no pure-Python logic of its own; import
``enviPat_py`` directly for the available functions.
"""

from . import enviPat_py  # noqa: F401
