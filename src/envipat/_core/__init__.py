"""envipat._core: thin CPython/NumPy C-extension wrapper around iso_core.*.

The compiled extension module lives at ``envipat._core.enviPat_py`` (built
from ``enviPat_py.c`` in this directory, plus the vendored, unmodified
algorithm sources -- ``iso_core.c``, ``combination.c``, ``profile.c``,
``parse.c``, ``element.c``, ``isotope.c``, ``peak.c``). It exposes the 5
core operations 1:1 with ``iso_core.h``, with NumPy arrays in/out and no
result-shaping logic of its own (e.g. the ``return_iso_calc_amount``
branch) -- that lives in the higher-level Python pipeline functions
(:mod:`envipat.isopattern`, :mod:`envipat.envelope`, ...).

This subpackage intentionally has no pure-Python logic of its own; import
``enviPat_py`` directly for the available functions.
"""

import envipat._core.enviPat_py as enviPat_py  # noqa: F401
