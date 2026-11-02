"""A faithful re-implementation of R's three-argument ``seq(from, to, by)``.

R's ``envelope`` builds its mass trace two different ways
depending on argument ``dmz``:

- ``dmz`` numeric, ``ppm=FALSE``: ``traceit <- seq(min(m)-0.5, max(m)+0.5, dmz)``
- ``dmz="get"``: ``traceit <- seq(min(m)-1, max(m)+0.5, dmz2*frac)``

(the third, ``ppm=TRUE`` branch calls the C entry point ``iso_ppm_Call``
directly -- ``envipat._core.enviPat_py.ppm_trace`` -- and needs no seq()
port at all).

``numpy.arange``/``linspace`` do not reproduce R's ``seq.default`` exactly:
R computes the element count as ``as.integer((to - from) / by + 1e-10)``
(a fudge-factor guard against the count coming out one-short due to
floating point error in the division) and then clamps the *last* generated
element back down to ``to`` (``pmin``/``pmax``) to correct any floating
point overshoot in the final step. Both of those behaviours matter for
matching R's own element count and last value.

This is intentionally a semantic port, not a bit-exact one: seq()-based
mass trace construction is a known source of floating-point divergence
from R's evaluator, which is exactly why the golden-fixture comparison
for ``envelope``/``vdetect`` uses ``rtol=1e-12`` rather than an exact match
(see tests/python/test_pipeline.py).
Source consulted: base R's ``seq.default`` (``print(seq.default)`` in R
4.1/4.6, unchanged in this respect across versions), for the
from/to/by-all-given, non-integer branch.
"""

from __future__ import annotations

import numpy as np

__all__ = ["r_seq"]


def r_seq(from_: float, to: float, by: float) -> np.ndarray:
    """Reimplements ``seq(from, to, by)`` for numeric (double) arguments.

    Only the ``from``/``to``/``by`` (all given, none integer-typed) branch
    of R's ``seq.default`` is implemented -- the only one enviPat's R
    sources ever hit for this call shape.
    """
    from_ = float(from_)
    to = float(to)
    by = float(by)

    delta = to - from_
    if delta == 0 and to == 0:
        return np.array([to], dtype=np.float64)

    n = delta / by
    if not np.isfinite(n):
        if by == 0 and delta == 0:
            return np.array([from_], dtype=np.float64)
        raise ValueError("invalid '(to - from)/by' in seq(.)")
    if n < 0:
        raise ValueError("wrong sign in 'by' argument")

    # R: "if (finite.del && abs(del)/max(abs(to), abs(from)) < 100 *
    # .Machine$double.eps) return(from)" -- from and to are numerically
    # indistinguishable given by's magnitude; seq collapses to one point.
    if abs(delta) / max(abs(to), abs(from_)) < 100 * np.finfo(np.float64).eps:
        return np.array([from_], dtype=np.float64)

    n_int = int(n + 1e-10)
    x = from_ + np.arange(n_int + 1, dtype=np.float64) * by
    # Correct for floating point overshoot on the last element, exactly as
    # R's seq.default does (pmin/pmax against `to`).
    if by > 0:
        x = np.minimum(x, to)
    else:
        x = np.maximum(x, to)
    return x
