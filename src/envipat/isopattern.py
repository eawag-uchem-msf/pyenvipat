"""Isotope pattern calculation.

Port of R's ``isopattern`` (``R/isopattern.R``) -- issue #8 of the R ->
Python port. Calls ``envipat._core.enviPat_py.pattern_algo1``/
``pattern_algo2`` (the Python bindings, from issue #5, of the same
``iso_core.c`` algorithms R's ``iso_pattern``/``iso_pattern_4`` ``.Call``
entry points delegate to -- see ``R/isopattern.R`` and
``src/envipat/_core/enviPat_py.c``).

Result shape
------------
Returns ``dict[str, pandas.DataFrame | str]`` keyed by (stringified) input
formula, matching R's ``pattern`` list (named by ``as.character(chemforms)``).
Unlike R -- where a duplicate name in ``names()<-`` simply overwrites the
earlier list entry, silently discarding its result -- this port raises
``ValueError`` as soon as a duplicate formula is detected in the input,
naming every duplicated formula (not just the first). A dict keyed by
formula cannot represent two results for the same formula anyway, so
silently keeping only the last one would be a silent data-loss bug in
Python even though it mirrors R's own list-naming semantics; see MR
description for this decision. Each value is a DataFrame with columns
``"m/z"``, ``"abundance"``,
then one integer column per isotope actually present in that formula's
peaks (e.g. ``"12C"``, ``"13C"``, ``"1H"``, ...), in the same order
``_core`` returns them.

For an invalid formula (the same "unknown element" / malformed-formula
condition ``check_chemform`` flags), R's C entry point returns
``R_NilValue`` and ``isopattern`` sets that list entry to the literal
string ``"error"`` instead of a data.frame (``R/isopattern.R``:
``if(!length(out[[1]])) pattern[[i]] <- "error"``). This is reproduced
here exactly: on ``EnviPatCoreError`` (or an empty result -- ``_core``
never actually returns one for pattern_algo1/2, but the empty-result path
is mirrored anyway for parity with the same idiom in ``envelope``/
``vdetect``), the dict value is the string ``"error"``, not a DataFrame.
Callers must check for this (``isinstance(pattern[f], str)``) before
treating a value as a DataFrame -- exactly as R callers must check
``identical(pattern[[i]], "error")`` (see ``check_several``'s R original,
which does exactly that).
"""

from __future__ import annotations

import sys
import warnings
from collections import Counter
from collections.abc import Sequence

import numpy as np
import pandas as pd

from ._core import enviPat_py as _core

__all__ = ["isopattern"]

# R default: 0.00054857990924 (electron mass, Da). Matches R/isopattern.R.
_DEFAULT_EMASS = 0.00054857990924

_ALGOS = {1: "pattern_algo1", 2: "pattern_algo2"}


def _as_charge_list(charge, n_formulas: int) -> list:
    """Mirrors R's charge-argument handling in isopattern() (R/isopattern.R):
    a scalar (including FALSE) recycled to length n_formulas, or a
    per-formula vector that must already match n_formulas in length."""
    if charge is False:
        return [False] * n_formulas
    if isinstance(charge, Sequence) and not isinstance(charge, str):
        charges = list(charge)
    else:
        charges = [charge]

    if len(charges) != n_formulas and len(charges) > 1:
        raise ValueError("length of charge does not match number of chemforms!")
    if len(charges) == 1 and n_formulas > 1:
        charges = charges * n_formulas

    for c in charges:
        if c is False:
            continue
        if c == 0:
            raise ValueError("charge=0?")
        if not isinstance(c, (int, float, np.integer, np.floating)):
            raise ValueError("charge either numeric or FALSE!")
    return charges


def isopattern(
    isotopes: pd.DataFrame,
    chemforms,
    threshold: float = 0.001,
    charge=False,
    emass: float = _DEFAULT_EMASS,
    algo: int = 1,
    rel_to: int = 0,
    verbose: bool = True,
    return_iso_calc_amount: bool = False,
):
    """Calculate isotope fine structure pattern(s) for one or more formulas.

    Port of R's ``isopattern`` (``R/isopattern.R``).

    Parameters
    ----------
    isotopes
        Isotope reference table (columns, by position: element, isotope,
        mass, abundance, ...), e.g. ``envipat.data.isotopes``. Must have
        exactly 5 columns, matching R's ``length(isotopes) != 5`` check
        (a data.frame's ``length()`` is its column count).
    chemforms
        A single formula string or a sequence of formula strings, such as
        ``envipat.data.chemforms`` or ``check_chemform(...)["new_formula"]``.
    threshold
        Probability (percentage points, 0-100) below which isotope peaks
        are pruned, as defined by ``rel_to``.
    charge
        ``False`` to omit charge/m-to-z conversion (default), a single
        int/float charge applied to every formula, or a sequence of
        per-formula charges (length must equal ``len(chemforms)``).
    emass
        Electron mass (Da), only used when ``charge`` is not ``False``.
    algo
        ``1`` (element-wise transition trees, ``iso_core_pattern_algo1``)
        or ``2`` (two-element transition trees, ``iso_core_pattern_algo4``
        -- R's "algo_4").
    rel_to
        Probability definition/pruning basis: 0=relative to the most
        intense peak, 1=relative to the "monoisotopic" peak, 2=absolute,
        3/4=percentage-of-(0/1) but returned as absolute probabilities.
        See ``R/isopattern.R``'s details section.
    verbose
        If True, prints ``" done."`` on completion (matching R's
        ``message(" done.")``; R's own per-formula progress messages from
        the C layer are not reproduced here, matching #5's core wrapper).
    return_iso_calc_amount
        If True, returns raw per-element isotope-combination counts
        (``dict[str, int]``) instead of the pattern DataFrame for every
        formula -- mirrors R's ``rica``/``return_iso_calc_amount`` branch,
        which returns a completely different (named-vector) shape.

    Returns
    -------
    dict[str, pandas.DataFrame | str | dict[str, int]]
        Keyed by (stringified) formula; see module docstring for the
        per-formula value shapes, including the ``"error"`` sentinel.

    Raises
    ------
    ValueError
        If ``chemforms`` contains the same (stringified) formula more
        than once -- see module docstring; this diverges from R, which
        silently keeps only the last duplicate's result.

    Notes
    -----
    R's ``plotit`` argument (base-graphics stick plot per formula) is
    intentionally not ported -- no plotting side effect, matching
    ``envipat.resolution.get_resolution``'s (issue #9) same omission of
    ``getR``'s own ``plotit``.
    """
    if isotopes.shape[1] != 5:
        raise ValueError("invalid isotope list")
    if threshold > 100 or threshold < 0:
        raise ValueError("invalid threshold; 0<=threshold<100.")
    if emass != _DEFAULT_EMASS:
        warnings.warn("You are sure that is the mass of an electrone?")
    if algo not in _ALGOS:
        raise ValueError("invalid algo argument!")
    if rel_to not in (0, 1, 2, 3, 4):
        raise ValueError("invalid rel_to")

    if isinstance(chemforms, str):
        chemforms = [chemforms]
    chemforms = [str(f) for f in chemforms]

    duplicates = sorted(f for f, n in Counter(chemforms).items() if n > 1)
    if duplicates:
        raise ValueError(
            "duplicate formula(s) in chemforms (each formula must be "
            "unique, since results are keyed by formula): "
            + ", ".join(duplicates)
        )

    charges = _as_charge_list(charge, len(chemforms))

    algo_fn = getattr(_core, _ALGOS[algo])

    elem = isotopes.iloc[:, 0].astype(str).tolist()
    iso_name = isotopes.iloc[:, 1].astype(str).tolist()
    iso_mass = isotopes.iloc[:, 2].to_numpy(dtype=np.float64)
    iso_abu = isotopes.iloc[:, 3].to_numpy(dtype=np.float64)

    pattern: dict = {}
    for formula, charge_i in zip(chemforms, charges):
        try:
            result = algo_fn(
                sum=formula,
                peak_limit=int(1e6),
                threshold=float(threshold),
                rel_to_mono=int(rel_to),
                iso_list_elem=elem,
                iso_list_iso=iso_name,
                iso_list_mass=iso_mass,
                iso_list_abu=iso_abu,
            )
        except _core.EnviPatCoreError:
            pattern[formula] = "error"
            continue

        if return_iso_calc_amount:
            pattern[formula] = dict(
                zip(result["element_names"], (int(c) for c in result["element_all_iso_calc_amount"]))
            )
            continue

        mz = result["mass"]
        if mz.size == 0:
            # Mirrors R's `if(!length(out[[1]]))` -- an empty-but-not-NULL
            # result would hit the same "error" branch there.
            pattern[formula] = "error"
            continue

        # R's `order(out[, 1], decreasing = FALSE)`: base R's default
        # order() method is "radix", which is a *stable* sort -- ties keep
        # their original (transition-tree emission) order.
        order = np.argsort(mz, kind="stable")
        mz_sorted = mz[order]
        abundance_sorted = result["abundance"][order]
        compound_sorted = result["compound"][order, :]

        if charge_i is not False:
            mz_sorted = mz_sorted - charge_i * emass
            if charge_i != 1:
                mz_sorted = mz_sorted / abs(charge_i)

        df = pd.DataFrame(compound_sorted, columns=result["isotope_names"])
        df.insert(0, "abundance", abundance_sorted)
        df.insert(0, "m/z", mz_sorted)
        pattern[formula] = df

    if verbose:
        # R's message() writes to stderr, not stdout.
        print(" done.", file=sys.stderr)

    return pattern
