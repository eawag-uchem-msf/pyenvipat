"""Chemical formula string helpers.

Port of R's ``check_chemform``, ``check_ded``, ``mergeform``, ``subform``
and ``multiform``.

``check_chemform`` validates and normalizes a chemical formula string (an
optional ``[NN]`` isotope-mass-number prefix, one capital letter, zero or
more lowercase letters, and an optional integer count -- repeated, and
optionally grouped in round brackets with their own trailing count, e.g.
``"(NH4)2SO4"``) and computes its monoisotopic mass from an ``isotopes``
reference table. ``check_ded``, ``mergeform``, ``subform`` and
``multiform`` are formula-algebra helpers built on the same element/count
grammar, but -- like their R originals -- they operate on *already
normalized, fully counted* formulas (no round brackets, every element
followed by an explicit count) such as ``check_chemform``'s own
``new_formula`` output; they are typically used downstream of it (e.g. for
adduct arithmetic, see ``adducts``).

Tokenizer, not a character-by-character port
----------------------------------------------
R's originals re-scan the formula character by character (``substr(x, j,
j)``) through several sequential passes (validate characters, validate
bracket balance, insert implicit ``1``s, expand ``(...)N`` groups
textually, then disassemble into element/count pairs). This module instead
parses each formula once with a small recursive-descent tokenizer
(``_parse_flat`` for the four flat-formula helpers, ``_parse_nested`` for
``check_chemform``), but reproduces R's exact normalization rules:

- an omitted count after an element or a ``(...)`` group defaults to 1
  ("implicit 1 insertion");
- ``(...)N`` groups multiply every element count inside them by N, with
  correct nesting (``"((CH2)2O)3"`` -> ``"C6H12O3"``);
- ``D`` is substituted for ``[2]H`` (a literal text substitution) in all
  five functions, consistently: ``check_chemform`` (both at the top level
  and inside a ``(...)`` group), ``check_ded`` (both its *deduct* and its
  *formulas* argument), ``mergeform``, ``subform`` and ``multiform``.
  Every one of these substitutions uses a negative lookahead
  (``D(?![a-z])``) so as not to clobber a longer symbol starting with
  "D" -- e.g. "Dy" (dysprosium, present in ``isotopes``).

  This is a deliberate, disclosed divergence from R (issue #14): R's
  five functions disagree with each other about ``D`` --
  ``check_chemform``'s top level looks ``"D"`` up directly as its own
  ``isotopes`` row instead of substituting it (so
  ``check_chemform(isotopes, "D2O1")`` reports formula ``"D2O1"``, not
  ``"[2]H2O1"``, with a very slightly different mass, because the ``D``
  and ``[2]H`` rows of ``isotopes`` carry mass values differing in their
  last digit), only its internal bracket-multiplier helper substitutes
  ``D`` -> ``[2]H`` (so R's own ``check_chemform`` already normalizes
  ``"(D2)3"`` to ``"[2]H6"``, not ``"D6"``); and R's ``check_ded``
  substitutes ``D`` -> ``[2]H`` in ``deduct`` but ``D`` -> plain ``H`` in
  ``formula1``. Issue #14 resolved this by canonicalizing: ``D`` now
  means deuterium (``[2]H``) everywhere in this Python port, in every
  function and every argument. This is a real, permanent behavioural
  difference from R for any formula containing a literal ``D`` -- not
  just a bug fix like the "Dy" guard below.

  Separately, R's own ``gsub("D", ...)`` calls lack the "Dy" guard above
  in ``check_ded``/``mergeform``/``subform`` and silently corrupt any
  formula containing "Dy" there (only R's ``multiform`` guards against
  it); this port applies the guard everywhere a ``D`` substitution
  happens, which -- unlike the canonicalization above -- *is* a plain bug
  fix, flagged as a divergence from R like the ``adducts.csv`` fix.
- ``check_chemform`` merges repeated elements in a formula, summing their
  counts and keeping the position of first occurrence (``"C2C3"`` ->
  ``"C5"``). ``mergeform`` does *not* do this for you: a repeated element
  within a single input formula (already-normalized formulas are not
  supposed to have one) raises ``ValueError`` instead of being silently
  merged -- see :func:`mergeform`.

Reference mass lookup
----------------------
Each element symbol's reference mass for the monoisotopic-mass sum is the
mass of its *most abundant* isotope row in ``isotopes`` (this is what
makes an isotope-labelled symbol like ``"[13]C"`` -- which appears in
``isotopes`` as its own single row with abundance 1 -- resolve to that
row's own mass, while a plain ``"C"`` resolves to natural carbon's
most-abundant isotope, ``12C``). This mirrors R's ``check_chemform``
internals exactly.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Literal, overload

import pandas as pd

__all__ = ["check_chemform", "check_ded", "mergeform", "multiform", "subform"]

_ALLOWED_CHARS_RE = re.compile(r"^[A-Za-z0-9()\[\]]*$")
# Matches a literal "D" meant as the deuterium shorthand, not as the first
# letter of a longer element symbol (e.g. "Dy" dysprosium, or a
# hypothetical future "Db" dubnium row) -- R's own gsub("D", ...) calls in
# check_ded/mergeform/subform lack this guard and silently corrupt any
# formula containing "Dy" (see MR !7 discussion); multiform's R original
# already guards against this and this port follows that pattern
# everywhere D is substituted.
_D_NOT_BEFORE_LOWER_RE = re.compile(r"D(?![a-z])")

_INVALID_MASS = -9999


def _format_r_number(x: float) -> str:
    """Format a number the way R's ``paste0``/``as.character`` would.

    Only used by :func:`multiform`, whose ``fact`` may be non-integral
    (R's own docs only show integer factors, but the function does not
    require one -- ``count * fact`` is pasted as-is). Whole-valued floats
    print without a trailing ``.0`` (R: ``as.character(5)`` -> ``"5"``,
    not ``"5.0"``); otherwise Python's ``repr`` already produces the same
    shortest round-tripping decimal representation R's default formatter
    does for values in this range (atom counts), verified against R for a
    range of fractional factors.
    """
    if x == int(x):
        return str(int(x))
    return repr(x)


class _FormulaSyntaxError(ValueError):
    """Raised internally when a formula does not conform to the grammar.

    Callers that need R's "invalid formula" *row* semantics (check_chemform)
    catch this; callers that mirror R's original hard-error-on-malformed-
    input behaviour (check_ded, mergeform, subform, multiform -- whose R
    originals loop forever on a formula missing an explicit trailing count,
    see module docstring / MR description) let it propagate.
    """


# ---------------------------------------------------------------------------
# Element/count tokenizer shared by check_ded / mergeform / subform /
# multiform: a *flat* sequence of (optional "[NN]" prefix + element symbol)
# tokens, each followed by a REQUIRED integer count. No round brackets.
# ---------------------------------------------------------------------------

_FLAT_TOKEN_RE = re.compile(r"(\[\d+\][A-Z][a-z]*|[A-Z][a-z]*)(\d+)")


def _parse_flat(formula: str) -> list[tuple[str, int]]:
    """Parse an already-normalized, fully-counted formula (no round
    brackets) into an ordered list of (element_symbol, count) pairs.

    Mirrors the element/number disassembly loop shared verbatim by
    check_ded, mergeform, subform and multiform in R -- which requires
    every element to be followed by an explicit digit count (an element
    with an omitted count makes R's while loop spin forever, see this
    module's docstring); this raises ``_FormulaSyntaxError`` promptly
    instead.
    """
    pos = 0
    n = len(formula)
    tokens: list[tuple[str, int]] = []
    while pos < n:
        m = _FLAT_TOKEN_RE.match(formula, pos)
        if not m:
            raise _FormulaSyntaxError(
                f"malformed formula {formula!r} at position {pos}: "
                "expected an element symbol followed by an explicit count "
                "(R's tokenizer never terminates on formulas with an "
                "implicit/omitted count here -- normalize with "
                "check_chemform first)"
            )
        symbol, count = m.group(1), int(m.group(2))
        tokens.append((symbol, count))
        pos = m.end()
    return tokens


def _merge_counts(tokens: list[tuple[str, int]]) -> dict[str, int]:
    """Sum counts per element symbol, in order of first occurrence."""
    merged: dict[str, int] = {}
    for symbol, count in tokens:
        merged[symbol] = merged.get(symbol, 0) + count
    return merged


def _reject_repeated_symbol(tokens: list[tuple[str, int]], formula: str) -> None:
    """Raise if ``formula`` (already tokenized into ``tokens``) contains
    the same element symbol more than once.

    Used by :func:`mergeform` on each of its two input formulas *before*
    combining them: a formula's own repeated, un-merged element (e.g.
    R's ``adducts`` ``"M+H+NH4"`` -> ``"H1N1H4"``) is what trips
    R's ``mergeform`` vector-recycling bug (see module docstring / MR
    description). An earlier version of this port silently merged such
    duplicates before combining, "fixing" the corruption rather than
    surfacing it; per reviewer feedback that is a real error the caller
    should see, not something mergeform should paper over, so this now
    fails loud instead.
    """
    seen: set[str] = set()
    for symbol, _ in tokens:
        if symbol in seen:
            raise _FormulaSyntaxError(
                f"formula {formula!r} contains element {symbol!r} more "
                "than once as separate tokens (e.g. 'H1...H4'); mergeform "
                "does not merge a formula's own repeated tokens for you -- "
                "normalize with check_chemform first"
            )
        seen.add(symbol)


# ---------------------------------------------------------------------------
# check_chemform's nested tokenizer: element/count tokens, optionally
# grouped in round brackets with their own trailing count, arbitrarily
# nested (e.g. "((CH2)2O)3"). Counts (element or group) default to 1 when
# omitted -- "implicit 1 insertion".
# ---------------------------------------------------------------------------


def _parse_nested(formula: str) -> list[tuple[str, int]]:
    """Parse a (possibly bracketed/nested) formula into a flat, ordered
    list of (element_symbol, count) pairs with all group multipliers
    already applied. Raises _FormulaSyntaxError for anything R's
    check_chemform would flag with warnit <- TRUE."""
    if formula == "":
        raise _FormulaSyntaxError("empty formula")
    if not _ALLOWED_CHARS_RE.match(formula):
        raise _FormulaSyntaxError(f"formula {formula!r} contains illegal characters")

    n = len(formula)
    pos = 0

    def parse_count(default: int = 1) -> int:
        nonlocal pos
        start = pos
        while pos < n and formula[pos].isdigit():
            pos += 1
        if pos == start:
            return default
        return int(formula[start:pos])

    def parse_element() -> tuple[str, int]:
        nonlocal pos
        prefix = ""
        if formula[pos] == "[":
            start = pos
            pos += 1
            digits_start = pos
            while pos < n and formula[pos].isdigit():
                pos += 1
            if pos == digits_start or pos >= n or formula[pos] != "]":
                raise _FormulaSyntaxError(f"malformed isotope bracket in {formula!r}")
            pos += 1
            prefix = formula[start:pos]
        if pos >= n or not formula[pos].isupper():
            raise _FormulaSyntaxError(
                f"expected an element symbol in {formula!r} at position {pos}"
            )
        start = pos
        pos += 1
        while pos < n and formula[pos].islower():
            pos += 1
        symbol = prefix + formula[start:pos]
        # D -> [2]H, unconditionally (issue #14): unlike R, which only
        # substitutes inside a round-bracket group at the top level (see
        # module docstring), this port canonicalizes D to deuterium
        # everywhere, so top-level and in-group formulas are treated the
        # same way.
        if symbol == "D":
            symbol = "[2]H"
        count = parse_count()
        return symbol, count

    def parse_group(top_level: bool) -> list[tuple[str, int]]:
        nonlocal pos
        terms: list[tuple[str, int]] = []
        started = False
        while pos < n and formula[pos] != ")":
            started = True
            ch = formula[pos]
            if ch == "(":
                pos += 1
                inner = parse_group(top_level=False)
                if pos >= n or formula[pos] != ")":
                    raise _FormulaSyntaxError(f"unbalanced '(' in {formula!r}")
                pos += 1
                if not inner:
                    raise _FormulaSyntaxError(f"empty group '()' in {formula!r}")
                count = parse_count()
                terms.extend((sym, cnt * count) for sym, cnt in inner)
            elif ch == "[" or ch.isupper():
                terms.append(parse_element())
            else:
                raise _FormulaSyntaxError(f"unexpected character {ch!r} in {formula!r}")
        if top_level and not started:
            raise _FormulaSyntaxError(f"empty formula body in {formula!r}")
        return terms

    terms = parse_group(top_level=True)
    if pos != n:
        raise _FormulaSyntaxError(f"unbalanced ')' in {formula!r}")
    return terms


@dataclass
class _IsotopeRef:
    order: dict[str, int] = field(default_factory=dict)
    mass: dict[str, float] = field(default_factory=dict)


def _build_isotope_ref(isotopes: pd.DataFrame) -> _IsotopeRef:
    """For each unique element symbol (isotopes' 1st column), in order of
    first appearance, record its reference (most-abundant-isotope) mass:
    the mass (3rd column) of the row with the highest abundance (4th
    column) among all rows sharing that symbol. Mirrors R's check_chemform
    (positional, not name-based, matching R's isotopes[,1]/
    [,3]/[,4] column access -- so a caller-supplied isotopes table need
    only match column *order*, not names)."""
    elem_col = isotopes.iloc[:, 0]
    mass_col = isotopes.iloc[:, 2].astype(float)
    abund_col = isotopes.iloc[:, 3].astype(float)

    ref = _IsotopeRef()
    best_abund: dict[str, float] = {}
    for symbol, mass, abund in zip(elem_col, mass_col, abund_col, strict=True):
        if symbol not in ref.order:
            ref.order[symbol] = len(ref.order)
            ref.mass[symbol] = mass
            best_abund[symbol] = abund
        elif abund > best_abund[symbol]:
            ref.mass[symbol] = mass
            best_abund[symbol] = abund
    return ref


@overload
def check_chemform(
    isotopes: pd.DataFrame,
    chemforms: str | Iterable[str],
    get_sorted: bool = ...,
    get_list: Literal[False] = ...,
) -> pd.DataFrame: ...
@overload
def check_chemform(
    isotopes: pd.DataFrame,
    chemforms: str | Iterable[str],
    get_sorted: bool = ...,
    *,
    get_list: Literal[True],
) -> list[dict[str, int]]: ...
@overload
def check_chemform(
    isotopes: pd.DataFrame,
    chemforms: str | Iterable[str],
    get_sorted: bool = ...,
    get_list: bool = ...,
) -> pd.DataFrame | list[dict[str, int]]: ...
def check_chemform(
    isotopes: pd.DataFrame,
    chemforms: str | Iterable[str],
    get_sorted: bool = False,
    get_list: bool = False,
) -> pd.DataFrame | list[dict[str, int]]:
    """Check and normalize chemical formulas; compute monoisotopic mass.

    Port of R's ``check_chemform``.

    Parameters
    ----------
    isotopes
        Isotope reference table (columns, by position: element, isotope,
        mass, abundance, ...), e.g. ``envipat.data.isotopes``.
    chemforms
        A single formula string or an iterable of formula strings.
    get_sorted
        If True, reorder each formula's elements into the order they
        first appear in ``isotopes`` (both in the returned formula string
        and, for ``get_list=True``, in the returned element/count mapping)
        instead of the order they first appear in the input formula. This
        also changes floating-point summation order and therefore, in
        rare cases, the last bit(s) of ``monoisotopic_mass``.
    get_list
        If True, return a ``list[dict[str, int]]`` (one dict of
        element -> count per input formula, in input order; an empty dict
        for an invalid formula) instead of a DataFrame.

    Returns
    -------
    pandas.DataFrame
        Columns ``warning`` (bool), ``new_formula`` (str) and
        ``monoisotopic_mass`` (float), one row per input formula, when
        ``get_list=False``. For an invalid formula, ``warning`` is True,
        ``new_formula`` echoes the original (unmodified) input string and
        ``monoisotopic_mass`` is -9999 -- matching R exactly.
    list[dict[str, int]]
        When ``get_list=True``.

    Raises
    ------
    ValueError
        If any formula contains a space (matches R's hard ``stop()``).
    """
    if isinstance(chemforms, str):
        chemforms = [chemforms]
    chemforms = list(chemforms)

    if any(" " in f for f in chemforms):
        raise ValueError("Some or all of strings in chemforms contain empty spaces - please revise")

    ref = _build_isotope_ref(isotopes)

    warn: list[bool] = []
    new_formula: list[str] = []
    mass: list[float] = []
    listed: list[dict[str, int]] = []

    for original in chemforms:
        try:
            terms = _parse_nested(original)
            counts = _merge_counts(terms)
            for symbol in counts:
                if symbol not in ref.mass:
                    raise _FormulaSyntaxError(f"unknown element {symbol!r} in {original!r}")
            order = list(counts.items())
            if get_sorted:
                order.sort(key=lambda item: ref.order[item[0]])
        except _FormulaSyntaxError:
            warn.append(True)
            new_formula.append(original)
            mass.append(float(_INVALID_MASS))
            listed.append({})
            continue

        formula_str = ""
        total_mass = 0.0
        result_counts: dict[str, int] = {}
        for symbol, count in order:
            formula_str += f"{symbol}{count}"
            total_mass += ref.mass[symbol] * count
            result_counts[symbol] = count

        warn.append(False)
        new_formula.append(formula_str)
        mass.append(total_mass)
        listed.append(result_counts)

    if get_list:
        return listed

    return pd.DataFrame({"warning": warn, "new_formula": new_formula, "monoisotopic_mass": mass})


def check_ded(formulas: str | Iterable[str], deduct: str) -> list[bool]:
    """Check whether ``deduct`` is fully contained in each of ``formulas``.

    Port of R's ``check_ded``. ``formulas`` and
    ``deduct`` must already be fully counted, unbracketed formulas (as
    produced by :func:`check_chemform`).

    ``D`` is substituted for ``[2]H`` before tokenizing, in both
    ``deduct`` and each of ``formulas``. R's original substitutes
    ``D`` -> ``[2]H`` in ``deduct`` but ``D`` -> plain ``H`` in
    ``formula1``; this port
    canonicalizes both sides to deuterium instead (issue #14), so e.g.
    ``check_ded(["C4D10"], "C4H10")`` now reports ``True`` (deduct's
    ``H10`` is *not* contained in the formula's ``[2]H10``), unlike R.

    Parameters
    ----------
    formulas
        A single formula string or an iterable of formula strings (the
        "containing" side).
    deduct
        A single formula string (the "deduct").

    Returns
    -------
    list[bool]
        One entry per formula in ``formulas``: True if ``deduct`` is
        *not* fully contained in that formula, False if it is.
    """
    if isinstance(formulas, str):
        formulas = [formulas]

    deduct_str = _D_NOT_BEFORE_LOWER_RE.sub("[2]H", str(deduct))
    deduct_counts = _merge_counts(_parse_flat(deduct_str))

    out: list[bool] = []
    for formula in formulas:
        formula_str = _D_NOT_BEFORE_LOWER_RE.sub("[2]H", str(formula))
        formula_counts = _merge_counts(_parse_flat(formula_str))
        contained = all(
            symbol in formula_counts and formula_counts[symbol] >= count
            for symbol, count in deduct_counts.items()
        )
        out.append(not contained)
    return out


def mergeform(formula1: str | Iterable[str], formula2: str | Iterable[str]) -> list[str]:
    """Combine chemical formulas by summing element counts.

    Port of R's ``mergeform``. Both arguments must
    already be fully counted, unbracketed formulas.

    Unlike R's original -- which, given a formula with an internally
    repeated, un-merged element (e.g. R's ``adducts``
    ``"M+H+NH4"`` -> ``"H1N1H4"``), silently corrupts every element
    processed afterwards via vector recycling (see MR !7 for the full
    trace) -- this raises ``ValueError`` promptly instead of "fixing" the
    corruption by merging the duplicate for the caller. Normalize with
    :func:`check_chemform` first if a formula may contain repeats.

    Parameters
    ----------
    formula1
        A single formula string or an iterable of formula strings.
    formula2
        A single formula string, merged into each of ``formula1``.

    Returns
    -------
    list[str]
        One merged formula string per entry in ``formula1``.

    Raises
    ------
    ValueError
        If a formula omits an explicit count after some element, or if
        either ``formula1`` or ``formula2`` contains the same element
        symbol more than once as separate tokens.
    """
    if isinstance(formula1, str):
        formula1 = [formula1]

    formula2_str = _D_NOT_BEFORE_LOWER_RE.sub("[2]H", str(formula2))
    tokens2 = _parse_flat(formula2_str)
    _reject_repeated_symbol(tokens2, formula2_str)
    counts2 = _merge_counts(tokens2)

    out: list[str] = []
    for f1 in formula1:
        f1_str = _D_NOT_BEFORE_LOWER_RE.sub("[2]H", str(f1))
        tokens1 = _parse_flat(f1_str)
        _reject_repeated_symbol(tokens1, f1_str)
        counts1 = _merge_counts(tokens1)

        merged: dict[str, int] = dict(counts1)
        for symbol, count in counts2.items():
            merged[symbol] = merged.get(symbol, 0) + count
        # R iterates `both <- unique(c(element1, element2))`: formula1's
        # elements first (in formula1's order), then any *new* elements
        # from formula2 (in formula2's order) that aren't already in
        # formula1.
        ordered = list(counts1.keys()) + [s for s in counts2 if s not in counts1]
        out.append("".join(f"{s}{merged[s]}" for s in ordered))
    return out


def subform(formula1: str | Iterable[str], formula2: str | Iterable[str]) -> list[str]:
    """Subtract ``formula2`` from each of ``formula1``.

    Port of R's ``subform``. Both arguments must already
    be fully counted, unbracketed formulas. Elements whose count drops to
    (or was already) zero are dropped from the result, matching R.

    Parameters
    ----------
    formula1
        A single formula string or an iterable of formula strings (the
        "minuend" side).
    formula2
        A single formula string, subtracted from each of ``formula1``.

    Returns
    -------
    list[str]
        One entry per formula in ``formula1``: the subtracted formula
        string if ``formula2`` was fully contained in it, otherwise a
        human-readable message describing why not (matching R, which
        returns a message string instead of raising) -- e.g.
        ``"Cl1 from formula 2 not part of formula1"`` or ``"Atom number of
        C from formula 2 not fully subset of formula1 atom number"``.
    """
    if isinstance(formula1, str):
        formula1 = [formula1]

    formula2_str = _D_NOT_BEFORE_LOWER_RE.sub("[2]H", str(formula2))
    counts2 = _merge_counts(_parse_flat(formula2_str))

    out: list[str] = []
    for f1 in formula1:
        f1_str = _D_NOT_BEFORE_LOWER_RE.sub("[2]H", str(f1))
        counts1 = _merge_counts(_parse_flat(f1_str))

        failure = None
        for symbol, count in counts2.items():
            if symbol not in counts1:
                failure = f"{symbol} from formula 2 not part of formula1"
                break
            if count > counts1[symbol]:
                failure = (
                    f"Atom number of {symbol} from formula 2 not fully "
                    "subset of formula1 atom number"
                )
                break
        if failure is not None:
            out.append(failure)
            continue

        result = ""
        for symbol, count in counts1.items():
            remaining = count - counts2.get(symbol, 0)
            if remaining > 0:
                result += f"{symbol}{remaining}"
        out.append(result)
    return out


def multiform(formula_in: str | Iterable[str], fact: float) -> list[str]:
    """Multiply every element count in a chemical formula by a factor.

    Port of R's ``multiform``. Must already be a
    fully counted, unbracketed formula.

    Parameters
    ----------
    formula_in
        A single formula string or an iterable of formula strings.
    fact
        Scalar factor to multiply every element count by (need not be an
        integer; R's own examples pass ``fact`` through unchanged into
        ``count * fact``, so e.g. ``fact=0.5`` halves every count without
        any rounding).

    Returns
    -------
    list[str]
        One multiplied formula string per entry in ``formula_in``.
    """
    if isinstance(formula_in, str):
        formula_in = [formula_in]

    out: list[str] = []
    for f in formula_in:
        f_str = _D_NOT_BEFORE_LOWER_RE.sub("[2]H", str(f))
        tokens = _parse_flat(f_str)
        out.append("".join(f"{symbol}{_format_r_number(count * fact)}" for symbol, count in tokens))
    return out
