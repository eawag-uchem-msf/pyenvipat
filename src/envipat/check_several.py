"""Check for overlapping molecules.

Port of R's ``check_several``. Screens isotope fine structures (from
:func:`envipat.isopattern.isopattern`) or envelope centroids/intensoids
(from :func:`envipat.vdetect.vdetect`) for m/z overlaps between different
compounds, within a tolerance window (absolute m/z or ppm).

Algorithm (mirrors R exactly): every peak of every compound is flattened
into one long (compound_index, mz, peak_number_within_compound) list,
stably sorted by mz ascending, then walked pairwise. Whenever two
*adjacent* peaks in that sorted order belong to different compounds and
the lower one's mz plus the tolerance window reaches at least the higher
one's mz, both compounds are flagged, and each accumulates a ``"/"``
separated trail of which other compound it overlapped with and which of
that other compound's peaks.

Known R bug -- deliberately fixed here, not reproduced
--------------------------------------------------------
R's ``ppm==TRUE`` branch has a stray trailing comma in one of its four
``paste0()`` calls (the second ``get3[...]`` update):
``paste0(get3[getit1[i - 1]], getit3[i], "/",)``. R's parser accepts a
trailing comma as an *empty* (missing) argument rather than a syntax
error, but ``paste0`` then raises ``"argument is missing, with no
default"`` the moment it actually tries to
evaluate that missing argument -- which happens on the very first overlap
found while ``ppm=TRUE``. In other words, ``check_several(..., ppm=TRUE)``
is not merely "a bit buggy": it hard-errors on *every* input where an
overlap exists, i.e. on the one condition the function exists to report.
Golden fixture ``tests/golden/check_several/case_004.csv`` captures this:
R's own ``check_several(pat_overlap, dmz=50, ppm=TRUE)`` errors with exactly
``"argument is missing, with no default"`` instead of producing a normal
output table.

This is a plain typo (a copy-paste of the ``ppm==FALSE`` block's four
updates with one stray comma), not a deliberate-but-odd design choice
worth preserving for behavioural parity. Reproducing it would make
``ppm=True`` unusable in Python for the one scenario it is meant to
detect, for no compatibility
benefit -- nothing in the R package or its docs relies on this branch
raising. This port removes the extra comma so the ``ppm=True`` branch is
symmetric with ``ppm=False`` and returns a normal result instead of
raising. See MR description for the same disclosure, per this repo's
established norm (MRs !4-!7) of calling out such decisions explicitly
rather than silently diverging from R.

Result shape
------------
Returns a :class:`pandas.DataFrame` with columns ``"compound"``,
``"warning"``, ``"to?"`` and ``"peak_number"``, corresponding
column-for-column to R's ``data.frame(names(pattern), get1, get2,
get3)`` -- but reshaped to native Python/pandas types rather than R's
own string-serialization idioms (per reviewer decision on MR !8; see MR
description):

- ``"compound"``: ``str`` (unchanged -- already native on both sides).
- ``"warning"``: ``bool``, not R's character column of ``""``/``"TRUE"``.
  In R's original, ``""`` isn't a missing/NA marker -- ``get1`` is
  initialized to ``""`` for every compound and only ever reassigned to
  the literal string ``"TRUE"`` when an overlap involving that compound
  is found, so ``""`` unambiguously means "no overlap found" i.e.
  ``False``, not "unknown"/``NA``. Hence a plain ``bool`` column
  (``True``/``False``), not a nullable one.
- ``"to?"`` / ``"peak_number"``: ``list[int]`` per row (empty list, not
  ``""``, when there is no overlap), not R's ``"/"``-terminated string
  of 1-based row numbers / peak indices (e.g. ``"2/2/2/"`` -> ``[2, 2,
  2]``). These are always non-negative integer indices (1-based, R-style,
  preserved as-is rather than renumbered to 0-based, since they index
  into ``pattern`` / a compound's own peak table the same way R's did),
  so a list of Python ``int`` is the natural, idiomatic shape -- not
  parsed floats, and not re-encoded as another delimited string.

A caller diffing this against R's own ``data.frame`` output (e.g. the
golden fixtures, which are frozen from R and still hold the original
``""``/``"TRUE"`` and ``"2/2/2/"`` strings) must parse those fixture
strings into the same native types before comparing -- see
``tests/python/test_pipeline.py``'s ``check_several`` test.
"""

from __future__ import annotations

import sys
from collections.abc import Mapping

import numpy as np
import pandas as pd

__all__ = ["check_several"]


def check_several(pattern: Mapping, dmz: float, ppm: bool = True) -> pd.DataFrame:
    """Check for m/z overlaps between different compounds' peaks.

    Port of R's ``check_several``.

    Parameters
    ----------
    pattern
        Mapping of compound -> DataFrame (columns ``"m/z"``,
        ``"abundance"``, ...), as returned by
        :func:`envipat.isopattern.isopattern` or
        :func:`envipat.vdetect.vdetect`. Must have at least 2 entries.
    dmz
        Overlap tolerance window: an absolute m/z value if ``ppm=False``,
        or a ppm value (must be > 0) if ``ppm=True``.
    ppm
        Whether ``dmz`` is a ppm value (default True, matching R) or an
        absolute m/z value.

    Returns
    -------
    pandas.DataFrame
        Columns ``"compound"`` (str), ``"warning"`` (bool), ``"to?"``
        (list[int]) and ``"peak_number"`` (list[int]) -- see module
        docstring for what each represents and why these are the native
        types, rather than R's own string-serialized column values.
    """
    pattern = dict(pattern)

    first = next(iter(pattern.values()))
    if len(first.columns) < 2 or list(first.columns[:2]) != ["m/z", "abundance"]:
        raise ValueError("WARNING: pattern has invalid entries\n")
    if ppm and dmz < 0:
        raise ValueError("\n WARNING: ppm=TRUE -> dmz must be >0\n")
    if len(pattern) < 2:
        raise ValueError("WARNING: nothing to compare")

    names = list(pattern.keys())
    n_compounds = len(names)

    getit1: list[int] = []  # 0-based compound index, one entry per peak
    getit2: list[float] = []  # that peak's m/z
    getit3: list[int] = []  # 1-based peak-within-compound index (R-style)
    for idx, name in enumerate(names):
        mz = pattern[name].iloc[:, 0].to_numpy(dtype=np.float64)
        getit1.extend([idx] * len(mz))
        getit2.extend(mz.tolist())
        getit3.extend(range(1, len(mz) + 1))

    getit1_arr = np.asarray(getit1, dtype=np.int64)
    getit2_arr = np.asarray(getit2, dtype=np.float64)
    getit3_arr = np.asarray(getit3, dtype=np.int64)

    # R's order() default method ("radix") is a stable sort.
    order = np.argsort(getit2_arr, kind="stable")
    getit1_arr = getit1_arr[order]
    getit2_arr = getit2_arr[order]
    getit3_arr = getit3_arr[order]

    get1 = [False] * n_compounds  # warning: True iff this compound overlaps another
    get2: list[list[int]] = [[] for _ in range(n_compounds)]  # to?: 1-based compound indices
    get3: list[list[int]] = [[] for _ in range(n_compounds)]  # peak_number: 1-based peak indices

    for k in range(1, len(getit1_arr)):
        left, right = k - 1, k
        if getit1_arr[left] == getit1_arr[right]:
            continue
        if not ppm:
            overlap = (getit2_arr[left] + dmz) >= getit2_arr[right]
        else:
            overlap = (getit2_arr[left] + (getit2_arr[left] * dmz / 1e6)) >= getit2_arr[right]
        if not overlap:
            continue

        i_left = int(getit1_arr[left])
        i_right = int(getit1_arr[right])
        get1[i_right] = True
        get1[i_left] = True
        get2[i_right].append(i_left + 1)
        get2[i_left].append(i_right + 1)
        get3[i_right].append(int(getit3_arr[left]))
        get3[i_left].append(int(getit3_arr[right]))

    if any(get1):
        print("Overlaps detected!", file=sys.stderr)

    return pd.DataFrame(
        {
            "compound": names,
            "warning": pd.array(get1, dtype="bool"),
            # dtype=object (not the default list->2D-ndarray inference)
            # matters here: if every row happens to hold a same-length
            # list (e.g. all empty, or all length 2), naive conversion
            # would collapse this into a 2D array instead of a Series of
            # per-row list objects.
            "to?": pd.Series(get2, dtype=object),
            "peak_number": pd.Series(get3, dtype=object),
        }
    )
