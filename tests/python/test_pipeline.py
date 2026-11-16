"""Golden-fixture regression tests for the pipeline functions:
isopattern, envelope, vdetect, check_several, isowrap.

Golden fixtures were frozen from the real R package (v2.8) in issue #3
(``tests/golden/<fn>/cases.csv`` + per-case output CSVs). Each fixture
file was written by stacking every input formula's own (possibly
differently-shaped) result table into one CSV -- a
union of every formula's columns, NA-filled for formulas that don't have
a given column. :func:`_split_by_compound` below reverses that: it
regroups the stacked table back into one frame per compound (columns
aligned by *name*, not position, since the union/NA-fill process can
interleave columns in a different relative order than any single
formula's own call would produce -- see that function's docstring).

Tolerance policy (per issue #8's own scope):

- ``isopattern``: ``rtol=1e-12`` on m/z and abundance. Same compiled
  ``iso_core`` algorithm as R (issue #5's Python bindings), same inputs, so
  the only slack needed is for pandas 2.x parsing some 17-digit CSV values
  ~1 ULP off (pandas 3 is exact). Isotope counts and strings match exactly.
- ``envelope``/``vdetect``: ``rtol=1e-12`` (see below for vdetect). The R idiom this ports
  builds mass traces with base R's ``seq(from, to, by)``
  (:func:`envipat._rseq.r_seq`); that is a semantic, not bit-exact, port
  (see that module's docstring) -- issue #8 explicitly calls this out as
  "a known source of floating-point divergence from R's evaluator".
  Measured deviations across the full envelope/vdetect case matrix are
  reported in the MR description (all comfortably under 1e-12, generally
  1e-14 to 1e-15 -- i.e. a handful of ULPs, not a real tolerance need).
  On pandas 2.x the golden/data CSV values themselves are parsed ~1 ULP
  off, which vdetect's valley arithmetic amplifies to ~2e-12 on a few
  small abundances, so ``vdetect`` uses ``rtol=1e-11``.
- ``check_several``: exact match, except the one case that exercises R's
  known ``ppm=TRUE`` trailing-comma bug (see
  ``envipat.check_several``'s module docstring for the fix-vs-reproduce
  decision) -- that case has no numeric golden target (R itself errored
  instead of producing one) and is asserted not to raise instead.
- ``isowrap``: structure and coarse agreement only (same compounds and
  columns, finite values, m/z span within 0.5 Da).
  Per-formula resolutions come from ``envipat.get_resolution``, which fits
  scipy's smoothing spline and deviates from R's by up to ~17% (issue #17),
  so values cannot be compared to the R golden outputs.

Run after `pip install -e .`:
    pytest tests/python/test_pipeline.py
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import envipat.data as data
from envipat.check_several import check_several
from envipat.envelope import envelope
from envipat.formula import check_chemform
from envipat.isopattern import isopattern
from envipat.isowrap import isowrap
from envipat.resolution import get_resolution
from envipat.vdetect import vdetect

REPO_ROOT = Path(__file__).resolve().parents[2]
GOLDEN = REPO_ROOT / "tests" / "golden"

ISOTOPES = data.isotopes
CHEMFORMS = list(data.chemforms)


# --------------------------------------------------------------------- #
# golden fixture I/O helpers
# --------------------------------------------------------------------- #


def _parse_params(params: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in params.split(";"):
        if not part:
            continue
        k, v = part.split("=", 1)
        out[k] = v
    return out


def _load_cases(section: str) -> dict[str, dict[str, str]]:
    with (GOLDEN / section / "cases.csv").open(newline="") as f:
        return {row["case_id"]: row for row in csv.DictReader(f)}


def _read_stacked_csv(path: Path) -> pd.DataFrame:
    # A formula that mixes a labelled isotope with that same element's
    # natural form makes rbind_fill() emit two genuinely duplicate-named
    # columns (e.g. two "13C" columns -- see issue #15 and
    # envipat.isopattern's own module docstring, which does the same:
    # `pd.DataFrame(..., columns=result["isotope_names"])` keeps literal
    # duplicate labels, not "13C"/"13C.1"). pandas.read_csv's header
    # parser always mangles duplicate header names into "13C"/"13C.1" --
    # that would desync the golden columns from what isopattern() itself
    # returns, so the header is read separately and reapplied as-is.
    header = next(csv.reader(path.open(newline="")))
    df = pd.read_csv(path, header=0)
    df.columns = header
    return df


def _split_by_compound(df: pd.DataFrame) -> dict[str, object]:
    """Reverses the fixture generator's row stacking (a column union,
    NA-filled per formula): regroups a stacked golden fixture table back
    into one value per compound, in first-appearance order.

    Columns are aligned by *name*: the column union can end up
    interleaving different formulas' isotope columns in a different
    relative order than any single formula's own isopattern() call would
    produce for itself (e.g. a later formula whose own element-parse
    order disagrees with the order an earlier formula already fixed in
    the union) -- callers must compare by column name, not position, and
    this returns each compound's own (non-NA) columns already reindexed
    that way (in the union's left-to-right order, which is a superset of,
    not necessarily identical to, that compound's own natural order --
    hence "align by name" on both sides of every comparison below).
    """
    out: dict[str, object] = {}
    for compound, group in df.groupby("compound", sort=False):
        if bool(group["error"].iloc[0]):
            out[compound] = "error"
            continue
        # dict.fromkeys: dedupe by label while preserving first-seen order.
        # A duplicate-named label (e.g. two "13C" columns) must appear only
        # once in the selector -- selecting a duplicate label N times from
        # a frame that already has N matching physical columns would index
        # by the cross product (N*N columns), not return them once each.
        cols = list(dict.fromkeys(c for c in df.columns if c not in ("compound", "error")))
        sub = group[cols].dropna(axis=1, how="all").reset_index(drop=True)
        out[compound] = sub
    return out


def _assert_profile_dict_matches(got: dict, expected: dict, *, rtol: float) -> None:
    assert set(got) == set(expected)
    for name, exp in expected.items():
        g = got[name]
        if isinstance(exp, str):
            assert g == exp, f"{name}: expected {exp!r}, got {g!r}"
            continue
        assert not isinstance(g, str), f"{name}: expected a DataFrame, got {g!r}"
        np.testing.assert_allclose(g["m/z"].to_numpy(), exp["m/z"].to_numpy(), rtol=rtol)
        np.testing.assert_allclose(
            g["abundance"].to_numpy(), exp["abundance"].to_numpy(), rtol=rtol
        )


# --------------------------------------------------------------------- #
# Shared formula pools -- replicates the fixture generator's pool10/pool3/
# pool5 construction using this package's own (golden-verified exact)
# check_chemform -- cross-checked once to confirm this reproduces the same
# 10 formulas in the same order the golden fixtures were built from.
# --------------------------------------------------------------------- #


def _build_pool10() -> list[str]:
    ck_all = check_chemform(ISOTOPES, CHEMFORMS)
    mass = ck_all["monoisotopic_mass"].to_numpy(dtype=np.float64)
    order = np.argsort(mass, kind="stable")  # R's order(): stable radix sort
    n = len(CHEMFORMS)

    idx_min = order[0]
    idx_max = order[-1]
    idx_med = order[round(n / 2) - 1]
    idx_q1 = order[round(n / 4) - 1]
    cl_heavy = next(f for f in CHEMFORMS if re.search("Cl[2-9]", f))
    br_heavy = next(f for f in CHEMFORMS if re.search("Br[1-9]", f))

    candidates = [
        CHEMFORMS[idx_min],
        CHEMFORMS[idx_max],
        CHEMFORMS[idx_med],
        CHEMFORMS[idx_q1],
        cl_heavy,
        br_heavy,
        "C1H4",
        "C50H80N10O15S2",
        "[13]C2C8D2H8[15]N2N2[18]O2O3",
        "C20Cl10",
    ]
    return list(dict.fromkeys(candidates))  # R's unique(): first-occurrence order


POOL10 = _build_pool10()
POOL3 = [POOL10[0], POOL10[6], POOL10[8]]
POOL5 = POOL10[:5]


# --------------------------------------------------------------------- #
# isopattern
# --------------------------------------------------------------------- #

ISOPATTERN_CASES = _load_cases("isopattern")


def _formulas_for_isopattern_case(row: dict[str, str]) -> list[str]:
    p = _parse_params(row["params"])
    n = int(p["n_formulas"])
    if row["notes"] == "threshold_0_full_fine_structure":
        return ["C1H4", POOL10[0]]
    if n == len(POOL10):
        return POOL10
    if n == len(CHEMFORMS):
        return CHEMFORMS
    raise AssertionError(f"{row['case_id']}: unexpected n_formulas={n}")


def _charge_for_isopattern_case(p: dict[str, str]):
    c = p["charge"]
    if c == "FALSE":
        return False
    if "/" in c:
        return [int(x) for x in c.split("/")]
    return int(c)


@pytest.mark.parametrize("case_id", sorted(ISOPATTERN_CASES))
def test_isopattern_matches_golden(case_id):
    row = ISOPATTERN_CASES[case_id]
    p = _parse_params(row["params"])
    formulas = _formulas_for_isopattern_case(row)
    charge = _charge_for_isopattern_case(p)

    got = isopattern(
        ISOTOPES,
        formulas,
        threshold=float(p["threshold"]),
        charge=charge,
        algo=int(p["algo"]),
        rel_to=int(p["rel_to"]),
        verbose=False,
    )

    expected = _split_by_compound(_read_stacked_csv(GOLDEN / "isopattern" / row["output_file"]))
    assert set(got) == set(expected)

    for name, exp in expected.items():
        g = got[name]
        if isinstance(exp, str):
            assert g == exp
            continue
        assert not isinstance(g, str), f"{name}: expected a DataFrame, got {g!r}"
        np.testing.assert_allclose(g["m/z"].to_numpy(), exp["m/z"].to_numpy(), rtol=1e-12, atol=0)
        np.testing.assert_allclose(
            g["abundance"].to_numpy(), exp["abundance"].to_numpy(), rtol=1e-12, atol=0
        )

        # dict.fromkeys: a formula mixing a labelled isotope with that same
        # element's natural form gives rbind_fill() (and isopattern() itself)
        # two genuinely duplicate-named columns (e.g. two "13C" columns --
        # see issue #15); dedupe by label so each name is only selected
        # once below (`g[col]`/`exp[col]` already returns *all* same-named
        # columns for a duplicate label, in that formula's own tokenizing
        # order on both sides -- see _read_stacked_csv/_split_by_compound).
        isotope_cols = list(dict.fromkeys(c for c in exp.columns if c not in ("m/z", "abundance")))
        assert set(isotope_cols) == {c for c in g.columns if c not in ("m/z", "abundance")}
        for col in isotope_cols:
            np.testing.assert_array_equal(
                g[col].to_numpy(dtype=np.int64), exp[col].to_numpy(dtype=np.int64)
            )


def test_isopattern_duplicate_formula_raises():
    # Regression test for reviewer decision on MR !8: a duplicate input
    # formula must raise, not silently overwrite the earlier dict entry
    # (see envipat.isopattern's module/Raises docstring).
    with pytest.raises(ValueError, match="C1H4"):
        isopattern(ISOTOPES, ["C1H4", "C2H6", "C1H4"], threshold=0.1, verbose=False)


def test_isopattern_multiple_duplicate_formulas_all_named():
    # All duplicated formulas must be named in the error, not just the
    # first one encountered.
    with pytest.raises(ValueError) as excinfo:
        isopattern(ISOTOPES, ["C1H4", "C2H6", "C1H4", "C2H6", "C3H8"], threshold=0.1, verbose=False)
    assert "C1H4" in str(excinfo.value)
    assert "C2H6" in str(excinfo.value)
    assert "C3H8" not in str(excinfo.value)


# --------------------------------------------------------------------- #
# envelope
# --------------------------------------------------------------------- #

ENVELOPE_CASES = _load_cases("envelope")
# Hardcoded to match the fixture generator's dmz_modes list exactly (its
# "abs" and "ppm" numeric dmz values -- the manifest's params only record
# the *mode name*, not the value that generated it).
_ENVELOPE_DMZ_ABS = 0.001
_ENVELOPE_DMZ_PPM = 5.0


@pytest.fixture(scope="module")
def pattern3():
    return isopattern(ISOTOPES, POOL3, threshold=0.1, charge=False, algo=1, verbose=False)


def _envelope_kwargs_for_case(p: dict[str, str]) -> dict:
    dmz_name = p["dmz"]
    if dmz_name == "get":
        dmz = "get"
    elif dmz_name == "abs":
        dmz = _ENVELOPE_DMZ_ABS
    elif dmz_name == "ppm":
        dmz = _ENVELOPE_DMZ_PPM
    else:
        raise AssertionError(f"unexpected dmz mode {dmz_name!r}")
    return {
        "ppm": p["ppm"] == "TRUE",
        "dmz": dmz,
        "env": p["env"],
        "resolution": float(p["resolution"]),
    }


@pytest.mark.parametrize("case_id", sorted(ENVELOPE_CASES))
def test_envelope_matches_golden(case_id, pattern3):
    row = ENVELOPE_CASES[case_id]
    p = _parse_params(row["params"])
    kwargs = _envelope_kwargs_for_case(p)

    got = envelope(pattern3, frac=0.25, verbose=False, **kwargs)
    expected = _split_by_compound(_read_stacked_csv(GOLDEN / "envelope" / row["output_file"]))

    _assert_profile_dict_matches(got, expected, rtol=1e-12)


# --------------------------------------------------------------------- #
# vdetect -- operates directly on the golden *envelope* fixture (per its
# own cases.csv "input_file" column), isolating vdetect's own fidelity
# from envelope's mass-trace construction.
# --------------------------------------------------------------------- #

VDETECT_CASES = _load_cases("vdetect")


@pytest.mark.parametrize("case_id", sorted(VDETECT_CASES))
def test_vdetect_matches_golden(case_id):
    row = VDETECT_CASES[case_id]
    p = _parse_params(row["params"])

    input_path = (GOLDEN / "vdetect" / row["input_file"]).resolve()
    profiles = _split_by_compound(_read_stacked_csv(input_path))

    got = vdetect(profiles, detect=p["detect"], verbose=False)
    expected = _split_by_compound(_read_stacked_csv(GOLDEN / "vdetect" / row["output_file"]))

    _assert_profile_dict_matches(got, expected, rtol=1e-11)


# --------------------------------------------------------------------- #
# check_several
# --------------------------------------------------------------------- #

CHECK_SEVERAL_CASES = _load_cases("check_several")


@pytest.fixture(scope="module")
def pattern_overlap():
    return isopattern(ISOTOPES, ["C10H10", "C9H10[13]C1"], threshold=0.1, verbose=False)


@pytest.fixture(scope="module")
def pattern_no_overlap():
    return isopattern(ISOTOPES, POOL3, threshold=0.1, verbose=False)


def _parse_golden_check_several_warning(s: str) -> bool:
    """R's golden ``"warning"`` column holds ``""``/``"TRUE"`` strings;
    envipat.check_several reshapes this to native bool (see its module
    docstring) -- ``""`` means "no overlap" i.e. False, not NA."""
    assert s in ("", "TRUE")
    return s == "TRUE"


def _parse_golden_check_several_index_list(s: str) -> list[int]:
    """R's golden ``"to?"``/``"peak_number"`` columns hold ``"/"``-
    terminated index-list strings (e.g. ``"2/2/2/"``, or ``""`` for no
    overlap); envipat.check_several reshapes these to native
    ``list[int]`` (``[]`` for no overlap) -- see its module docstring."""
    if s == "":
        return []
    return [int(x) for x in s.split("/") if x != ""]


@pytest.mark.parametrize("case_id", sorted(CHECK_SEVERAL_CASES))
def test_check_several_matches_golden(case_id, pattern_overlap, pattern_no_overlap):
    row = CHECK_SEVERAL_CASES[case_id]
    p = _parse_params(row["params"])
    scenario = p["scenario"]
    dmz = float(p["dmz"])
    ppm = p["ppm"] == "TRUE"
    pattern = (
        pattern_overlap
        if "overlap" in scenario and "no_overlap" not in scenario
        else pattern_no_overlap
    )

    if scenario == "ppmTRUE_overlap":
        # Known R bug, deliberately fixed rather than reproduced -- see
        # envipat.check_several's module docstring. R itself errors here
        # ("argument is missing, with no default";
        # tests/golden/check_several/case_004.csv captures that error
        # message, not a numeric result), so there is no golden output
        # to compare against: assert this port does *not* raise and
        # instead returns a normal, symmetric result.
        result = check_several(pattern, dmz=dmz, ppm=ppm)
        assert list(result.columns) == ["compound", "warning", "to?", "peak_number"]
        assert pd.api.types.is_bool_dtype(result["warning"])
        assert result["warning"].any()
        return

    result = check_several(pattern, dmz=dmz, ppm=ppm)
    assert pd.api.types.is_bool_dtype(result["warning"])
    assert result["to?"].map(lambda v: isinstance(v, list)).all()
    assert result["peak_number"].map(lambda v: isinstance(v, list)).all()

    raw = pd.read_csv(
        GOLDEN / "check_several" / row["output_file"], keep_default_na=False, dtype=str
    )
    expected = pd.DataFrame(
        {
            "compound": raw["compound"],
            "warning": raw["warning"].map(_parse_golden_check_several_warning),
            "to?": raw["to?"].map(_parse_golden_check_several_index_list),
            "peak_number": raw["peak_number"].map(_parse_golden_check_several_index_list),
        }
    )
    pd.testing.assert_frame_equal(
        result.reset_index(drop=True), expected.reset_index(drop=True), check_dtype=False
    )


# --------------------------------------------------------------------- #
# isowrap
# --------------------------------------------------------------------- #

# ``get_resolution`` fits scipy's make_smoothing_spline, which deviates from
# R's fixed-spar smooth.spline by up to ~17% (issue #17). Profile sampling and
# peak shapes depend on the resolution, so the isowrap golden outputs can no
# longer be compared value by value. We only check structure and coarse
# agreement: same compounds, same string-vs-table outcome and columns, finite
# values, and the same overall m/z span. (The strongest peak is not compared: near-equal
# isotope peaks can swap rank under a different resolution.)
_ISOWRAP_MZ_ATOL = 0.5  # Da; loose on purpose


def _assert_profile_structure_matches(got: dict, expected: dict) -> None:
    assert set(got) == set(expected)
    for name, exp in expected.items():
        g = got[name]
        if isinstance(exp, str):
            assert g == exp, f"{name}: expected {exp!r}, got {g!r}"
            continue
        assert not isinstance(g, str), f"{name}: expected a DataFrame, got {g!r}"
        assert list(g.columns) == list(exp.columns)
        assert len(g) > 0
        assert np.all(np.isfinite(g.to_numpy(dtype=float)))
        gm, em = g["m/z"].to_numpy(), exp["m/z"].to_numpy()
        assert abs(gm.min() - em.min()) < _ISOWRAP_MZ_ATOL
        assert abs(gm.max() - em.max()) < _ISOWRAP_MZ_ATOL


ISOWRAP_CASES = _load_cases("isowrap")


@pytest.mark.parametrize("case_id", sorted(ISOWRAP_CASES))
def test_isowrap_matches_golden(case_id):
    row = ISOWRAP_CASES[case_id]
    p = _parse_params(row["params"])

    checked5 = check_chemform(ISOTOPES, POOL5)
    resolution = get_resolution(checked5, data.resolution_list[p["curve"]]).to_numpy(
        dtype=np.float64
    )

    got = isowrap(
        ISOTOPES,
        checked5,
        resolution,
        threshold=float(p["threshold"]),
        charge=int(p["charge"]),
        algo=int(p["algo"]),
        ppm=False,
        dmz="get",
        frac=0.25,
        env="Gaussian",
        detect=p["detect"],
        verbose=False,
    )
    expected = _split_by_compound(_read_stacked_csv(GOLDEN / "isowrap" / row["output_file"]))

    _assert_profile_structure_matches(got, expected)
