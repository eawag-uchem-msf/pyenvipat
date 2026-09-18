"""Chemical formula string helpers.

Port of R's ``check_chemform`` (``R/check_chemform.R``), ``check_ded``
(``R/check_ded.R``), ``mergeform`` (``R/mergeform.R``), ``subform``
(``R/subform.R``) and ``multiform`` (``R/multiform.R``) -- issue #7 of the
R -> Python port.

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
- ``D`` is substituted for ``[2]H`` (a literal text substitution, matching
  R's ``gsub("D", "[2]H", ...)``) in ``check_ded``'s *deduct* formula,
  ``mergeform``, ``subform`` and (via a negative lookahead so as not to
  clobber a hypothetical "Db" symbol) ``multiform``; ``check_ded``'s
  *formula1* substitutes ``D`` for plain ``H`` instead -- a deliberate
  asymmetry present in R (see ``R/check_ded.R`` lines 29 vs 87) that this
  port preserves rather than "fixing". ``check_chemform`` itself performs
  *no* top-level D substitution: ``"D"`` is looked up directly as its own
  element symbol in ``isotopes`` (it is present there, with abundance 1),
  so e.g. ``check_chemform(isotopes, "D2O1")`` reports formula ``"D2O1"``
  (not ``"[2]H2O1"``) and a very slightly different mass than the
  chemically-equivalent ``"[2]H2O1"``, because the ``D`` and ``[2]H``
  rows of ``isotopes`` carry mass values that differ in their last
  digit. This is reproduced exactly, not smoothed over. There is one
  exception: R's ``check_chemform`` *does* apply the ``D`` -> ``[2]H``
  substitution inside its internal bracket-multiplier helper, so a
  formula like ``"(D2)3"`` normalizes to ``"[2]H6"``, not ``"D6"`` -- an
  inconsistency in the R original that is reproduced here too (see the
  ``_in_group`` flag inside ``_parse_nested`` below).
- merging repeated elements in a formula sums their counts, keeping the
  position of first occurrence (``"C2C3"`` -> ``"C5"``).

Reference mass lookup
----------------------
Each element symbol's reference mass for the monoisotopic-mass sum is the
mass of its *most abundant* isotope row in ``isotopes`` (this is what
makes an isotope-labelled symbol like ``"[13]C"`` -- which appears in
``isotopes`` as its own single row with abundance 1 -- resolve to that
row's own mass, while a plain ``"C"`` resolves to natural carbon's
most-abundant isotope, ``12C``). This mirrors R's ``check_chemform``
internals exactly (see ``R/check_chemform.R`` lines 133-147).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import pandas as pd

__all__ = ["check_chemform", "check_ded", "mergeform", "subform", "multiform"]

_ALLOWED_CHARS_RE = re.compile(r"^[A-Za-z0-9()\[\]]*$")
_D_RE = re.compile(r"D")
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
            raise _FormulaSyntaxError(f"expected an element symbol in {formula!r} at position {pos}")
        start = pos
        pos += 1
        while pos < n and formula[pos].islower():
            pos += 1
        symbol = prefix + formula[start:pos]
        # R's internal bracket-multiplier helper (multif) substitutes
        # D -> [2]H before tokenizing; the top-level disassembly does not.
        # Reproduced verbatim: only elements parsed *inside* a round-bracket
        # group go through this substitution.
        if _in_group[0] and symbol == "D":
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
                was_in_group = _in_group[0]
                _in_group[0] = True
                inner = parse_group(top_level=False)
                _in_group[0] = was_in_group
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

    _in_group = [False]
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
    column) among all rows sharing that symbol. Mirrors R/check_chemform.R
    lines 133-147 (positional, not name-based, matching R's isotopes[,1]/
    [,3]/[,4] column access -- so a caller-supplied isotopes table need
    only match column *order*, not names)."""
    elem_col = isotopes.iloc[:, 0]
    mass_col = isotopes.iloc[:, 2].astype(float)
    abund_col = isotopes.iloc[:, 3].astype(float)

    ref = _IsotopeRef()
    best_abund: dict[str, float] = {}
    for symbol, mass, abund in zip(elem_col, mass_col, abund_col):
        if symbol not in ref.order:
            ref.order[symbol] = len(ref.order)
            ref.mass[symbol] = mass
            best_abund[symbol] = abund
        elif abund > best_abund[symbol]:
            ref.mass[symbol] = mass
            best_abund[symbol] = abund
    return ref


def check_chemform(
    isotopes: pd.DataFrame,
    chemforms,
    get_sorted: bool = False,
    get_list: bool = False,
):
    """Check and normalize chemical formulas; compute monoisotopic mass.

    Port of R's ``check_chemform`` (``R/check_chemform.R``).

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
        raise ValueError(
            "Some or all of strings in chemforms contain empty spaces - please revise"
        )

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

    return pd.DataFrame(
        {"warning": warn, "new_formula": new_formula, "monoisotopic_mass": mass}
    )


def check_ded(formulas, deduct) -> list:
    """Check whether ``deduct`` is fully contained in each of ``formulas``.

    Port of R's ``check_ded`` (``R/check_ded.R``). ``formulas`` and
    ``deduct`` must already be fully counted, unbracketed formulas (as
    produced by :func:`check_chemform`).

    Note the intentional asymmetry preserved from R: within ``deduct``,
    ``D`` is substituted for ``[2]H`` before tokenizing; within each of
    ``formulas``, ``D`` is instead substituted for plain ``H`` (see
    ``R/check_ded.R`` lines 29 and 87) -- so e.g.
    ``check_ded(["C4D10"], "C4H10")`` reports ``False`` (not contained is
    False, i.e. it *is* considered contained), because the formula's "D"
    is read as ordinary H here, not as deuterium.

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

    deduct_str = _D_RE.sub("[2]H", str(deduct))
    deduct_counts = _merge_counts(_parse_flat(deduct_str))

    out: list[bool] = []
    for formula in formulas:
        formula_str = _D_RE.sub("H", str(formula))
        formula_counts = _merge_counts(_parse_flat(formula_str))
        contained = all(
            symbol in formula_counts and formula_counts[symbol] >= count
            for symbol, count in deduct_counts.items()
        )
        out.append(not contained)
    return out


def mergeform(formula1, formula2) -> list:
    """Combine chemical formulas by summing element counts.

    Port of R's ``mergeform`` (``R/mergeform.R``). Both arguments must
    already be fully counted, unbracketed formulas.

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
    """
    if isinstance(formula1, str):
        formula1 = [formula1]

    formula2_str = _D_RE.sub("[2]H", str(formula2))
    counts2 = _merge_counts(_parse_flat(formula2_str))

    out: list[str] = []
    for f1 in formula1:
        f1_str = _D_RE.sub("[2]H", str(f1))
        counts1 = _merge_counts(_parse_flat(f1_str))

        merged: dict[str, int] = dict(counts1)
        for symbol, count in counts2.items():
            merged[symbol] = merged.get(symbol, 0) + count
        # R iterates `both <- unique(c(element1, element2))`: formula1's
        # elements first (in formula1's order), then any *new* elements
        # from formula2 (in formula2's order) that aren't already in
        # formula1.
        ordered = list(counts1.keys()) + [
            s for s in counts2.keys() if s not in counts1
        ]
        out.append("".join(f"{s}{merged[s]}" for s in ordered))
    return out


def subform(formula1, formula2) -> list:
    """Subtract ``formula2`` from each of ``formula1``.

    Port of R's ``subform`` (``R/subform.R``). Both arguments must already
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

    formula2_str = _D_RE.sub("[2]H", str(formula2))
    counts2 = _merge_counts(_parse_flat(formula2_str))

    out: list[str] = []
    for f1 in formula1:
        f1_str = _D_RE.sub("[2]H", str(f1))
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


def multiform(formula_in, fact) -> list:
    """Multiply every element count in a chemical formula by a factor.

    Port of R's ``multiform`` (``R/multiform.R``). Must already be a
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
        out.append(
            "".join(f"{symbol}{_format_r_number(count * fact)}" for symbol, count in tokens)
        )
    return out
