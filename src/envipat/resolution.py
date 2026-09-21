"""Resolution interpolation (port of R's ``getR`` / ``R/getR.R``).

R's ``getR`` fits ``stats::smooth.spline(resmass[, 1], resmass[, 2], cv =
TRUE, all.knots = FALSE, nknots = nknots, spar = spar)`` -- a *fixed*-``spar``,
fixed-``nknots`` cubic smoothing spline -- and evaluates it (``predict``) at
each ``checked`` mass to obtain an interpolated instrument resolution.

There is no drop-in equivalent for this exact model in the Python scientific
stack: :func:`scipy.interpolate.make_smoothing_spline` selects its own
smoothing parameter via GCV and does not accept R's ``spar`` at all, and
:class:`scipy.interpolate.UnivariateSpline` parametrizes smoothing via a
target sum-of-squared-residuals rather than R's ``spar``/``lambda``. This was
measured, not assumed: even with its smoothing parameter hand-tuned per
curve against R's own known answer, ``make_smoothing_spline`` still misses R
by up to 13% relative error (target: ``rtol < 1e-6``) -- see MR !6's
description for the full comparison. This module is therefore a
from-scratch reimplementation of the relevant slice of
R's ``stats::smooth.spline`` (the ``cv=TRUE`` argument only affects a
diagnostic score that ``getR`` never reads -- since ``spar`` is supplied,
R's own ``smooth.spline`` skips any GCV/CV search entirely and fits directly
for that ``spar``, so it is not reproduced here).

``get_resolution`` still exposes ``scipy.interpolate.make_smoothing_spline``
as an explicit opt-in (``method="scipy"``) for callers who don't need
R-bit-compatibility and would rather use an upstream-maintained fitter; it
is not a recommended replacement for the default path -- see that
parameter's docstring and MR !6 for the measured accuracy gap.

Implementation notes (see ``stats::smooth.spline``'s R and Fortran sources
-- ``R/smooth.spline.R``, and ``src/library/stats/src/{sbart.c,stxwx.f,
sgram.f,sslvrg.f}`` in R's own sources -- for the reference algorithm this
mirrors):

- The unique ``x`` values (``resmass``'s mass column) are rescaled to
  ``[0, 1]``.
- ``nknots`` interior knots are placed at quantile positions: the unique,
  rescaled x-values indexed by ``floor(linspace(1, nx, nknots))`` (R's
  ``xbar[seq.int(1, nx, length.out = nknots)]``, which truncates -- not
  rounds -- non-integer vector indices), padded with the boundary knot
  repeated 3 additional times at each end for a clamped cubic (order-4)
  B-spline basis of size ``nknots + 2``.
- The fit solves the penalized normal equations ``(X'WX + lambda*Sigma) c =
  X'Wy`` where ``X`` is the B-spline collocation matrix and ``Sigma`` is the
  Gram matrix of the basis' second derivatives (``Sigma[i, j] = integral
  B''_i(t) B''_j(t) dt``).
- ``lambda`` is derived from ``spar`` via R's documented mapping ``lambda =
  ratio * 256 ** (3*spar - 1)``, where ``ratio`` is *not* the full
  ``trace(X'WX) / trace(Sigma)`` as the R help page's simplified description
  states, but the ratio of *partial* traces excluding the first 2 and last 3
  diagonal entries (the boundary B-splines, whose curvature is anomalously
  large due to the repeated boundary knots) -- confirmed against R's actual
  ``sbart.c`` source (``ratio := sum(hs0[3:(nk-3)]) / sum(sg0[3:(nk-3)])`` in
  1-based Fortran indexing).
- ``Sigma`` itself is computed via R's own closed-form per-knot-interval
  formula (exact, since a cubic B-spline's second derivative is piecewise
  linear): for each interval of width ``W`` bounded by consecutive distinct
  knots, with ``yw1``/``yw2`` the second-derivative value at the left end and
  its change across the interval, ``integral f_i*f_j = W * (yw1_i*yw1_j +
  0.5*(yw1_i*yw2_j + yw1_j*yw2_i) + C*yw2_i*yw2_j)``. R's Fortran
  ``sgram.f`` uses the literal constant ``C = 0.3330`` rather than the
  mathematically exact ``1/3`` -- reproduced here verbatim (not "fixed") so
  results match R bit-for-bit rather than a more "correct" spline.
- Predictions for masses inside ``[min(x), max(x)]`` evaluate the fitted
  B-spline directly; R's linear-extrapolation behavior for out-of-range
  masses is reproduced too, though ``get_resolution`` -- like R's ``getR``
  -- rejects out-of-range masses before fitting, so it is normally unused.

This was validated against R (v4.6.1, package ``enviPat`` v2.8) directly:
back-solving R's own fitted coefficients against this module's ``X``/
``Sigma`` matrices recovers R's ``lambda`` to full double precision, and the
final ``ratio``/``lambda``/predicted-resolution values match R to
1e-12--1e-16 relative precision across every ``resolution_list`` curve and
multiple ``(nknots, spar)`` combinations (see the MR description for the
full comparison against the golden fixtures).
"""

from __future__ import annotations

from typing import Literal

import numpy as np
import pandas as pd
from scipy.interpolate import BSpline, make_smoothing_spline

__all__ = ["get_resolution"]

_DEGREE = 3  # cubic B-spline (order 4), matching R's smooth.spline


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


def _r_quantile_knot_indices(nx: int, nknots: int) -> np.ndarray:
    """0-based indices into the (sorted, rescaled) unique x-values at which
    to place interior knots, mirroring R's
    ``xbar[seq.int(1, nx, length.out = nknots)]``.

    R's vector indexing truncates (not rounds) non-integer numeric indices,
    so this must truncate too, not round to nearest.
    """
    if nknots == 1:
        seq = np.array([1.0])
    else:
        i = np.arange(nknots, dtype=float)
        seq = 1.0 + i * (nx - 1) / (nknots - 1)
    idx1based = seq.astype(np.int64)  # truncates toward zero
    return idx1based - 1


def _build_knot_vector(xbar: np.ndarray, nknots: int) -> np.ndarray:
    idx = _r_quantile_knot_indices(len(xbar), nknots)
    interior = xbar[idx]
    return np.concatenate((np.repeat(xbar[0], 3), interior, np.repeat(xbar[-1], 3)))


def _second_derivative_gram_matrix(knots: np.ndarray, nk: int) -> np.ndarray:
    """Sigma[i, j] = integral of B''_i(t) * B''_j(t) dt over the knot span,
    computed via R's own exact per-interval closed form (a cubic B-spline's
    second derivative is piecewise linear, fully determined by its value at
    each interval's two endpoints) -- literal ``0.3330`` constant and all,
    see the module docstring.
    """
    basis = BSpline(knots, np.eye(nk), _DEGREE, extrapolate=False)
    d2 = basis.derivative(2)
    breaks = np.unique(knots)
    sigma = np.zeros((nk, nk))
    for a, b in zip(breaks[:-1], breaks[1:]):
        width = b - a
        if width <= 0:
            continue
        yw1 = d2(a)
        yw2 = d2(b) - yw1
        sigma += width * (
            np.outer(yw1, yw1)
            + 0.5 * (np.outer(yw1, yw2) + np.outer(yw2, yw1))
            + 0.3330 * np.outer(yw2, yw2)
        )
    return sigma


class _SmoothSplineFit:
    """A fitted fixed-spar, fixed-nknots cubic smoothing spline, holding
    just what's needed to evaluate it (mirrors R's ``$fit`` component)."""

    __slots__ = ("knots", "coef", "min", "range")

    def __init__(self, knots, coef, xmin, xrange):
        self.knots = knots
        self.coef = coef
        self.min = xmin
        self.range = xrange

    def __call__(self, xnew: np.ndarray) -> np.ndarray:
        xnew = np.atleast_1d(np.asarray(xnew, dtype=float))
        xs = (xnew - self.min) / self.range
        spline = BSpline(self.knots, self.coef, _DEGREE, extrapolate=False)
        y = np.empty_like(xs)

        interior = (xs >= 0) & (xs <= 1)
        if np.any(interior):
            y[interior] = spline(xs[interior])

        # R linearly extrapolates outside [0, 1] using the boundary value
        # and slope; getR's own range check normally makes this unreachable.
        left, right = xs < 0, xs > 1
        if np.any(left) or np.any(right):
            deriv1 = spline.derivative(1)
            if np.any(left):
                y[left] = spline(0.0) + deriv1(0.0) * xs[left]
            if np.any(right):
                y[right] = spline(1.0) + deriv1(1.0) * (xs[right] - 1)
        return y


def _fit_smooth_spline(x, y, nknots: int, spar: float) -> _SmoothSplineFit:
    ux, ybar, wbar = _dedup_and_aggregate(x, y)
    nx = len(ux)
    if nknots > nx:
        raise ValueError(
            f"cannot use more inner knots ({nknots}) than unique x values ({nx})"
        )
    xrange = ux[-1] - ux[0]
    xbar = (ux - ux[0]) / xrange
    knots = _build_knot_vector(xbar, nknots)
    nk = nknots + 2

    X = BSpline.design_matrix(xbar, knots, _DEGREE, extrapolate=False).toarray()
    sigma = _second_derivative_gram_matrix(knots, nk)
    XtWX = X.T @ (wbar[:, None] * X)

    # spar -> lambda: R's documented mapping lambda = ratio * 256**(3*spar-1),
    # where "ratio" is tr(X'WX)/tr(Sigma) restricted to the *interior* basis
    # functions (excludes the first 2 / last 3 diagonal entries, i.e. the
    # boundary B-splines) -- see module docstring for why the full trace
    # (as R's own help page loosely describes it) does not reproduce R.
    hs0 = np.diag(XtWX)
    sg0 = np.diag(sigma)
    lo, hi = 2, nk - 3
    if hi <= lo:
        raise ValueError(f"nknots={nknots} too small to compute a smoothing ratio")
    ratio = hs0[lo:hi].sum() / sg0[lo:hi].sum()
    lam = ratio * 256.0 ** (3 * spar - 1)

    coef = np.linalg.solve(XtWX + lam * sigma, X.T @ (wbar * ybar))
    return _SmoothSplineFit(knots, coef, ux[0], xrange)


def get_resolution(
    checked: pd.DataFrame,
    resmass: pd.DataFrame,
    nknots: int = 13,
    spar: float = 0.1,
    method: Literal["r_exact", "scipy"] = "r_exact",
) -> pd.Series:
    """Interpolate instrument resolution at each checked formula's mass.

    Port of R's ``getR`` (``R/getR.R``): fits a fixed-``spar``, fixed-
    ``nknots`` cubic smoothing spline to ``resmass`` (mass -> resolution)
    and evaluates it at each mass in ``checked``.

    Parameters
    ----------
    checked
        The result of ``check_chemform`` (see ``R/check_chemform.R`` /
        issue #7's Python port): a DataFrame with (at least) boolean
        ``warning`` and numeric ``monoisotopic_mass`` columns, one row per
        formula.
    resmass
        Two-column mass/resolution calibration table, such as one entry of
        ``resolution_list`` (columns ``mz``, ``R``); only column *position*
        is used (mass first, resolution second), matching R's ``resmass[,
        1]`` / ``resmass[, 2]``.
    nknots
        Number of knots for the smoothing spline (default 13, matching R's
        ``getR`` default). Must be between 3 and the number of unique
        masses in ``resmass``. Meaningless -- and must be left at its
        default -- under ``method="scipy"``, since
        ``scipy.interpolate.make_smoothing_spline`` has no knot-count
        parameter; passing a non-default value raises ``ValueError``
        rather than silently ignoring it.
    spar
        Smoothing parameter in (0, 1] (default 0.1, matching R's ``getR``
        default). Meaningless -- and must be left at its default -- under
        ``method="scipy"``, for the same reason as ``nknots``.
    method
        Which smoothing-spline backend to fit with:

        - ``"r_exact"`` (default): the from-scratch reimplementation of
          R's ``stats::smooth.spline(nknots=, spar=)`` described above.
          Matches R to ~1e-13--1e-15 relative error and is what this
          port's accuracy target (``rtol=1e-6``) assumes.
        - ``"scipy"``: fits with
          :func:`scipy.interpolate.make_smoothing_spline` instead (its own
          GCV-selected smoothing; no ``nknots``/``spar`` knobs). This is an
          explicit opt-in, *not* a recommended replacement for the default:
          MR !6 measured up to ~17% relative error against R across the
          ``resolution_list`` curves, far outside this port's accuracy
          target -- see MR !6's "Why not
          scipy.interpolate.make_smoothing_spline?" section for the full
          comparison. Use it only when exact R-reproduction isn't required.

    Returns
    -------
    pandas.Series
        Interpolated resolutions, indexed like ``checked`` (name
        ``"resolution"``).

    Notes
    -----
    R's ``getR`` also accepts a ``plotit`` argument to plot the fit; that is
    intentionally not ported here (no plotting side effect).
    """
    if method not in ("r_exact", "scipy"):
        raise ValueError(f"invalid method: {method!r}; use 'r_exact' or 'scipy'")
    if len(resmass) < 10:
        raise ValueError("not enough data points in resmass")
    if method == "scipy":
        # nknots/spar have no scipy.interpolate.make_smoothing_spline
        # equivalent (see module + method docstrings); raise rather than
        # silently ignoring an explicit non-default caller setting.
        if nknots != 13 or spar != 0.1:
            raise ValueError(
                "nknots/spar are meaningless for method='scipy' "
                "(scipy.interpolate.make_smoothing_spline has no such "
                "parameters); leave them at their defaults, or use "
                "method='r_exact' to control the fit via nknots/spar"
            )
    else:
        if nknots < 3 or nknots > len(resmass):
            raise ValueError("invalid nknots")
        if spar <= 0 or spar >= 1:
            raise ValueError("invalid spar; spar=(0,1]")
    if bool(np.any(checked["warning"])):
        raise ValueError("invalid chemical formula in checked['warning']")

    resmass_x = resmass.iloc[:, 0].to_numpy(dtype=float)
    resmass_y = resmass.iloc[:, 1].to_numpy(dtype=float)
    masses = checked["monoisotopic_mass"].to_numpy(dtype=float)

    if np.any(masses < resmass_x.min()) or np.any(masses > resmass_x.max()):
        raise ValueError("some mean_mass out of range of resmass")

    if method == "scipy":
        # Same dedup/aggregation as the r_exact path, so the two code
        # paths only diverge in the fitting math, not in how (possibly
        # unsorted / duplicate-x) resmass rows are handled.
        ux, ybar, wbar = _dedup_and_aggregate(resmass_x, resmass_y)
        spline = make_smoothing_spline(ux, ybar, w=wbar)
        resolution = spline(masses)
    else:
        fit = _fit_smooth_spline(resmass_x, resmass_y, nknots=nknots, spar=spar)
        resolution = fit(masses)
    return pd.Series(resolution, index=checked.index, name="resolution")
