"""Isotope pattern envelope (profile) calculation.

Port of R's ``envelope``. Convolutes an isotope pattern (from
:func:`envipat.isopattern.isopattern`) with a peak shape function (Gaussian
or Cauchy-Lorentz) into a discretized theoretical profile, by calling
``envipat._core.enviPat_py.ppm_trace`` and ``profile_with_trace`` (Python
bindings of the original ``iso_core.c`` algorithms).

Mass trace construction
------------------------
R builds the mass trace ("stick" positions) three different ways,
depending on ``ppm``/``dmz``:

- ``ppm=True``: via the C core directly (``_core.ppm_trace``) -- no
  ``seq()`` involved, so this branch is exact (same core, same inputs) up to
  the ``rtol=1e-12`` bound this module's results are verified at overall
  (see module docstring below).
- ``ppm=False, dmz`` numeric: ``seq(min(m) - 0.5, max(m) + 0.5, dmz)``.
- ``ppm=False, dmz="get"``: derives the step from ``resolution`` --
  ``dmz2 = mean(m) / resolution``; ``seq(min(m) - 1, max(m) + 0.5,
  dmz2 * frac)``.

The two ``seq()``-based branches are reproduced via :func:`envipat._rseq.r_seq`
(see that module for why this is a semantic, not bit-exact, port of R's
``seq.default`` -- this is *the* known source of floating-point divergence
from R, which is why ``envelope``/``vdetect``
are verified at ``rtol=1e-12`` rather than exact match, unlike
``isopattern``).

Result shape
------------
Returns ``dict[str, pandas.DataFrame | str]`` keyed like the input
``pattern`` dict, each DataFrame with two columns, ``"m/z"`` and
``"abundance"``. As in R (and :func:`envipat.isopattern.isopattern`), an
empty/erroring convolution result is reported as the literal string
``"error"`` for that entry, covering both a
raised ``EnviPatCoreError`` and a legitimate zero-length result.
"""

from __future__ import annotations

import sys
from collections.abc import Mapping

import numpy as np
import pandas as pd

from ._core import enviPat_py as _core
from ._rseq import r_seq

__all__ = ["envelope"]

_ENV_TYPES = {"Gaussian": 0, "CauchyLorentz": 1}


def envelope(
    pattern: Mapping,
    ppm: bool = False,
    dmz="get",
    frac: float = 0.25,
    env: str = "Gaussian",
    resolution=5e5,
    verbose: bool = True,
):
    """Convolute isotope pattern(s) into their theoretical envelope profile(s).

    Port of R's ``envelope``.

    Parameters
    ----------
    pattern
        Mapping of formula -> isotope-pattern DataFrame, as returned by
        :func:`envipat.isopattern.isopattern` (columns ``"m/z"``,
        ``"abundance"``, ... in that order).
    ppm
        If True, discretize the mass trace in ppm (``dmz`` is then a ppm
        value, must be numeric -- not ``"get"``). If False (default),
        discretize in absolute m/z, either fixed (``dmz`` numeric) or
        derived from ``resolution`` (``dmz="get"``).
    dmz
        ``"get"`` (default) to derive the step from ``resolution``/``frac``,
        or a numeric step (ppm or absolute m/z, per ``ppm``).
    frac
        Only used when ``dmz="get"``: fraction of the peak's FWHM
        (``mean(m)/resolution``) used as the trace step.
    env
        Peak shape: ``"Gaussian"`` (default) or ``"CauchyLorentz"``.
    resolution
        A single resolution value, or a sequence with one value per
        formula in ``pattern`` (length must then equal ``len(pattern)``).
    verbose
        If True, prints ``" done."`` on completion (matching R's
        ``message(" done.")``).

    Returns
    -------
    dict[str, pandas.DataFrame | str]
        Keyed like ``pattern``; see module docstring for the
        ``"error"`` sentinel.

    Notes
    -----
    R's ``plotit`` argument is intentionally not ported -- no plotting
    side effect, matching :func:`envipat.isopattern.isopattern`.
    """
    pattern = dict(pattern)
    n = len(pattern)

    if hasattr(resolution, "__len__") and not isinstance(resolution, str):
        resolution_list = list(resolution)
        if len(resolution_list) != n and len(resolution_list) > 1:
            raise ValueError("length of resolution does not match length of pattern list!\n")
    else:
        resolution_list = [resolution]

    if ppm and dmz == "get":
        raise ValueError("\n WARNING: ppm=TRUE -> dmz must be numerical!\n")
    if ppm and not isinstance(dmz, str) and dmz < 0:
        raise ValueError("\n WARNING: ppm=TRUE -> dmz must be >0\n")
    if dmz == "get" and any(r < 1 for r in resolution_list):
        raise ValueError("WARNING: invalid resolution!\n")
    if env not in _ENV_TYPES:
        raise ValueError("WARNING: invalid env argument\n")
    if ppm and dmz == "get" and frac > 1:
        raise ValueError("WARNING: invalid frac argument\n")
    if ppm and dmz == "get" and frac < 0:
        raise ValueError("WARNING: invalid frac argument\n")

    first = next(iter(pattern.values())) if pattern else None
    if (
        first is None
        or not isinstance(first, pd.DataFrame)
        or len(first.columns) < 3
        or list(first.columns[:2]) != ["m/z", "abundance"]
    ):
        raise ValueError("WARNING: is list, but has invalid entries\n")

    if verbose:
        print("Calculate profiles ...", file=sys.stderr)

    profile_type = _ENV_TYPES[env]

    if len(resolution_list) != n:
        resolution_list = resolution_list * n

    extend = 0.5
    profiles: dict = {}
    for (formula, pat), res in zip(pattern.items(), resolution_list):
        m = pat.iloc[:, 0].to_numpy(dtype=np.float64)
        a = pat.iloc[:, 1].to_numpy(dtype=np.float64)

        if ppm:
            trace = _core.ppm_trace(start=float(m.min() - extend), end=float(m.max() + extend), ppm=float(dmz))
        elif dmz != "get":
            trace = r_seq(m.min() - extend, m.max() + extend, float(dmz))
        else:
            dmz2 = float(m.mean()) / res
            trace = r_seq(m.min() - 1.0, m.max() + extend, dmz2 * frac)

        try:
            out_mz, out_abundance = _core.profile_with_trace(
                profile_type=profile_type,
                mass=m,
                abundance=a,
                trace=trace,
                resolution=int(res),
                threshold=0.0,
                filter=1,
            )
        except _core.EnviPatCoreError:
            profiles[formula] = "error"
            continue

        if out_mz.size == 0:
            profiles[formula] = "error"
            continue

        profiles[formula] = pd.DataFrame({"m/z": out_mz, "abundance": out_abundance})

    if verbose:
        print(" done.", file=sys.stderr)

    return profiles
