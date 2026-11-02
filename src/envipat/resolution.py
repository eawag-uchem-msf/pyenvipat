"""Resolution interpolation (port of R's ``getR``).

R's ``getR`` fits a smoothing spline to the instrument's mass -> resolution
calibration table and evaluates it at each formula's mass. This port fits
with :func:`scipy.interpolate.make_smoothing_spline` (smoothing parameter
chosen by GCV) instead of reproducing R's fixed-``spar``/fixed-``nknots``
``stats::smooth.spline``.

Accuracy note: results are *not* bit-compatible with R. Against the frozen R
reference the interpolated resolution deviates by up to ~17% relative error
on the ``resolution_list`` curves. The ``nknots``/``spar`` arguments of R's
``getR`` have no scipy equivalent and are therefore not supported.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.interpolate import make_smoothing_spline

__all__ = ["get_resolution"]


def _r_iqr(x: np.ndarray) -> float:
    """Interquartile range using R's default (type-7) quantile algorithm,
    which matches numpy's default ('linear') percentile interpolation."""
    q75, q25 = np.percentile(x, [75, 25], method="linear")
    return float(q75 - q25)


def _dedup_and_aggregate(x, y, w=None, tol=None):
    """Mirrors smooth.spline's binning of near-duplicate x values: sorts by
    x, groups consecutive points whose (x - mean(x))/tol rounds to the same
    integer bin, and aggregates y (weighted mean) and w (sum) per group.

    Returns (ux, ybar, wbar): the unique/aggregated x values (ascending),
    the per-group weighted-mean y, and the per-group summed weight.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    if tol is None:
        tol = 1e-6 * _r_iqr(x)
    if w is None:
        wv = np.ones(n)
    else:
        w = np.asarray(w, dtype=float)
        scale = (w * np.sum(w > 0)) / np.sum(w)
        wv = np.full(n, scale) if np.ndim(scale) == 0 else scale

    order = np.argsort(x, kind="stable")
    xs, ys, ws = x[order], y[order], wv[order]
    bins = np.round((xs - x.mean()) / tol)

    group_id = np.zeros(n, dtype=int)
    gid = 0
    for i in range(1, n):
        if bins[i] != bins[i - 1]:
            gid += 1
        group_id[i] = gid
    nx = gid + 1

    ux = np.empty(nx)
    wbar = np.zeros(nx)
    ybar = np.zeros(nx)
    for g in range(nx):
        mask = group_id == g
        ux[g] = xs[mask][0]
        wbar[g] = ws[mask].sum()
        ybar[g] = (ws[mask] * ys[mask]).sum() / wbar[g] if wbar[g] > 0 else 0.0
    return ux, ybar, wbar


def get_resolution(checked: pd.DataFrame, resmass: pd.DataFrame) -> pd.Series:
    """Interpolate instrument resolution at each checked formula's mass.

    Port of R's ``getR``, fitted with
    :func:`scipy.interpolate.make_smoothing_spline` (see the module docstring
    for the accuracy caveat versus R).

    Parameters
    ----------
    checked
        The result of ``check_chemform``: a DataFrame with (at least) boolean
        ``warning`` and numeric ``monoisotopic_mass`` columns, one row per
        formula.
    resmass
        Two-column mass/resolution calibration table, such as one entry of
        ``resolution_list`` (columns ``mz``, ``R``); only column *position*
        is used (mass first, resolution second), matching R's ``resmass[,
        1]`` / ``resmass[, 2]``. At least 10 rows are required.

    Returns
    -------
    pandas.Series
        Interpolated resolutions, indexed like ``checked`` (name
        ``"resolution"``).

    Notes
    -----
    R's ``getR`` also accepts ``nknots``, ``spar`` and ``plotit``; none of
    them are ported (no scipy equivalent / no plotting side effect).
    """
    if len(resmass) < 10:
        raise ValueError("not enough data points in resmass")
    if bool(np.any(checked["warning"])):
        raise ValueError("invalid chemical formula in checked['warning']")

    resmass_x = resmass.iloc[:, 0].to_numpy(dtype=float)
    resmass_y = resmass.iloc[:, 1].to_numpy(dtype=float)
    masses = checked["monoisotopic_mass"].to_numpy(dtype=float)

    if np.any(masses < resmass_x.min()) or np.any(masses > resmass_x.max()):
        raise ValueError("some mean_mass out of range of resmass")

    ux, ybar, wbar = _dedup_and_aggregate(resmass_x, resmass_y)
    spline = make_smoothing_spline(ux, ybar, w=wbar)
    return pd.Series(spline(masses), index=checked.index, name="resolution")
