"""Combined (batch) calculation of isotope pattern, envelope and centroids.

Port of R's ``isowrap``. Composes :func:`envipat.isopattern.isopattern`,
:func:`envipat.envelope.envelope` and :func:`envipat.vdetect.vdetect`.

``resolution`` is explicit
--------------------------
R's ``isowrap`` takes a calibration table (``resmass``) plus ``nknots``/
``spar``, and calls ``getR`` internally to interpolate a per-formula
resolution from it whenever ``resmass`` is given (``resolution=FALSE``
means "use ``getR``"; a numeric ``resolution`` bypasses it). This port's
``isowrap`` drops ``resmass``/``nknots``/``spar`` entirely and only accepts
``resolution`` directly (a single value or a per-formula sequence, passed
straight through to :func:`envipat.envelope.envelope`) -- callers wanting
curve-based interpolation call
:func:`envipat.resolution.get_resolution` themselves and pass its result
in as ``resolution``.

Result shape
------------
Returns whatever :func:`envipat.vdetect.vdetect` returns: a
``dict[str, pandas.DataFrame | str]`` keyed by formula.
"""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from .envelope import envelope
from .isopattern import isopattern
from .vdetect import vdetect

__all__ = ["isowrap"]

_DEFAULT_EMASS = 0.00054857990924


def isowrap(
    isotopes: pd.DataFrame,
    checked: pd.DataFrame,
    resolution: float | Sequence[float],
    threshold: float = 0.1,
    charge: float | Sequence[float] | bool = 1,
    emass: float = _DEFAULT_EMASS,
    algo: int = 2,
    ppm: bool = False,
    dmz: float | str = "get",
    frac: float = 0.25,
    env: str = "Gaussian",
    detect: str = "centroid",
    verbose: bool = True,
) -> dict[str, pd.DataFrame | str]:
    """Run the full isopattern -> envelope -> vdetect pipeline.

    Port of R's ``isowrap``; see module docstring for
    why ``resolution`` is explicit here rather than interpolated from a
    calibration table (``resmass``/``nknots``/``spar`` in R -- not
    ported).

    Parameters
    ----------
    isotopes
        Isotope reference table, e.g. ``envipat.data.isotopes``. Passed
        straight through to :func:`envipat.isopattern.isopattern`.
    checked
        Output of ``check_chemform`` (``envipat.formula.check_chemform``):
        a DataFrame with (by position) a ``warning`` column and a
        ``new_formula`` column. Every row must have ``warning`` falsy,
        matching R's ``if(any(checked[, 1])) stop(...)``.
    resolution
        A single resolution value, or a sequence with one value per row
        of ``checked`` -- see module docstring.
    threshold, charge, emass, algo
        Passed through to :func:`envipat.isopattern.isopattern`.
    ppm, dmz, frac, env
        Passed through to :func:`envipat.envelope.envelope`.
    detect
        Passed through to :func:`envipat.vdetect.vdetect`.
    verbose
        Passed through to all three composed functions.

    Returns
    -------
    dict[str, pandas.DataFrame | str]
        See :func:`envipat.vdetect.vdetect`'s return value.

    Raises
    ------
    ValueError
        If any row of ``checked`` has a truthy ``warning`` (an invalid
        chemical formula).
    """
    if bool(checked.iloc[:, 0].astype(bool).any()):
        raise ValueError("WARNING: checked with incorrect chemical detected!")

    pattern = isopattern(
        isotopes,
        checked.iloc[:, 1],
        threshold=threshold,
        charge=charge,
        emass=emass,
        algo=algo,
        verbose=verbose,
    )
    profiles = envelope(
        pattern,
        ppm=ppm,
        dmz=dmz,
        frac=frac,
        env=env,
        resolution=resolution,
        verbose=verbose,
    )
    centro = vdetect(
        profiles,
        detect=detect,
        verbose=verbose,
    )
    return centro
