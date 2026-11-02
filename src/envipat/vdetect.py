"""Valley detection and centroidization.

Port of R's ``vdetect``. Calls ``envipat._core.enviPat_py.centroid`` (the
Python binding of the original ``iso_core.c`` algorithm) on each envelope
profile from
:func:`envipat.envelope.envelope`.

Result shape
------------
Returns ``dict[str, pandas.DataFrame | str]`` keyed like the input
``profiles`` dict, each DataFrame with two columns ``"m/z"`` and
``"abundance"``. As in :func:`envipat.envelope.envelope`, an empty/erroring
detection result is reported as the literal string ``"error"``,
covering both a raised ``EnviPatCoreError`` and a legitimate zero-length
result.

Verified against the golden fixtures at ``rtol=1e-12`` (not exact), for
the same reason as ``envelope``: the profiles this operates on were
themselves built from a ``seq()``-based mass trace whose floating-point
values only approximately match R's (see ``envipat._rseq``).
"""

from __future__ import annotations

import sys
from collections.abc import Mapping

import numpy as np
import pandas as pd

from ._core import enviPat_py as _core

__all__ = ["vdetect"]

_DETECT_TYPES = {"centroid": 0, "intensoid": 1, "valley": 2}


def vdetect(
    profiles: Mapping,
    detect: str = "centroid",
    verbose: bool = True,
):
    """Detect centroids, intensoids or valleys in envelope profile(s).

    Port of R's ``vdetect``.

    Parameters
    ----------
    profiles
        Mapping of formula -> envelope-profile DataFrame, as returned by
        :func:`envipat.envelope.envelope` (columns ``"m/z"``,
        ``"abundance"``).
    detect
        ``"centroid"`` (default), ``"intensoid"`` or ``"valley"`` -- see
        R's ``?vdetect`` "Definitions" section for the precise meaning of
        each.
    verbose
        If True, prints ``" done."`` on completion (matching R's
        ``message(paste0("Detect ", detect, "s... "))`` /
        ``message(" done.")`` -- the first, detect-specific message is not
        reproduced here, only the completion one, matching this module's
        siblings).

    Returns
    -------
    dict[str, pandas.DataFrame | str]
        Keyed like ``profiles``; see module docstring for the ``"error"``
        sentinel.

    Notes
    -----
    R's ``plotit`` argument is intentionally not ported -- no plotting
    side effect, matching :func:`envipat.isopattern.isopattern`.
    """
    profiles = dict(profiles)

    first = next(iter(profiles.values())) if profiles else None
    if (
        first is None
        or not isinstance(first, pd.DataFrame)
        or len(first.columns) < 2
        or list(first.columns[:2]) != ["m/z", "abundance"]
    ):
        raise ValueError("WARNING: profile has invalid entries\n")
    if detect not in _DETECT_TYPES:
        raise ValueError("WARNING: invalid detect argument")

    return_type = _DETECT_TYPES[detect]

    getlist: dict = {}
    for formula, prof in profiles.items():
        m = prof.iloc[:, 0].to_numpy(dtype=np.float64)
        a = prof.iloc[:, 1].to_numpy(dtype=np.float64)

        try:
            out_mz, out_abundance = _core.centroid(profile_mass=m, profile_abundance=a, type=return_type)
        except _core.EnviPatCoreError:
            getlist[formula] = "error"
            continue

        if out_mz.size == 0:
            getlist[formula] = "error"
            continue

        getlist[formula] = pd.DataFrame({"m/z": out_mz, "abundance": out_abundance})

    if verbose:
        print(" done.", file=sys.stderr)

    return getlist
