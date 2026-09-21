"""Golden-fixture regression tests for issue #8's pipeline port:
isopattern, envelope, vdetect, check_several, isowrap.

Golden fixtures were frozen from the real R package (v2.8) in issue #3
(``tools/make_golden.R``, ``tests/golden/<fn>/cases.csv`` + per-case
output CSVs). Each fixture file was written by stacking every input
formula's own (possibly differently-shaped) result table into one CSV via
that script's ``stack_pattern_list()``/``rbind_fill()`` helpers -- a
union of every formula's columns, NA-filled for formulas that don't have
a given column. :func:`_split_by_compound` below reverses that: it
regroups the stacked table back into one frame per compound (columns
aligned by *name*, not position, since the union/NA-fill process can
interleave columns in a different relative order than any single
formula's own call would produce -- see that function's docstring).

Tolerance policy (per issue #8's own scope):

- ``isopattern``: exact match. Same compiled ``iso_core`` algorithm as R
  (issue #5's Python bindings), same inputs -- no floating-point
  operation in between other than R-identical charge/emass arithmetic
  and a stable sort, so there is no principled reason to allow any
  slack, and none is given (``np.testing.assert_array_equal``/exact
  string equality).
- ``envelope``/``vdetect``: ``rtol=1e-12``. The R idiom this porrts
  builds mass traces with base R's ``seq(from, to, by)``
  (:func:`envipat._rseq.r_seq`); that is a semantic, not bit-exact, port
  (see that module's docstring) -- issue #8 explicitly calls this out as
  "a known source of floating-point divergence from R's evaluator".
  Measured deviations across the full envelope/vdetect case matrix are
  reported in the MR description (all comfortably under 1e-12, generally
  1e-14 to 1e-15 -- i.e. a handful of ULPs, not a real tolerance need).
- ``check_several``: exact match, except the one case that exercises R's
  known ``ppm=TRUE`` trailing-comma bug (see
  ``envipat.check_several``'s module docstring for the fix-vs-reproduce
  decision) -- that case has no numeric golden target (R itself errored
  instead of producing one) and is asserted not to raise instead.
- ``isowrap``: same ``rtol=1e-12`` as envelope/vdetect (it composes both).
  Since issue #8's ``isowrap`` takes ``resolution`` explicitly rather than
  interpolating it from a calibration table (``getR`` integration is
  issue #9's job, not merged into this branch -- see
  ``envipat.isowrap``'s module docstring), the golden isowrap fixtures'
  per-formula resolutions (originally computed by R's ``isowrap`` calling
  ``getR`` internally) are recomputed once, out of band, by asking the
  *installed R package itself* for them (``tools/dump_isowrap_resolution.R``,
  same pattern as issue #5's ``tools/dump_isotopes_precise.R``) and fed
  in as the explicit ``resolution`` argument -- this does not import or
  depend on issue #9's Python port at all. Skips (not fails) when
  Rscript / the installed R package aren't available, exactly like
  ``tests/test_core.py``'s ``full_isotopes`` fixture.

Run inside the dev container, after `pip install -e .`:
    pytest tests/python/test_pipeline.py
"""

from __future__ import annotations

import csv
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import envipat.data as data
from envipat._core import enviPat_py as _core
from envipat.check_several import check_several
from envipat.envelope import envelope
from envipat.formula import check_chemform
from envipat.isopattern import isopattern
from envipat.isowrap import isowrap
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
    # float_precision="round_trip": these fixtures were written with
    # sprintf("%.17g", .) to round-trip every double bit-exactly (see
    # tools/make_golden.R); pandas' default C float parser silently drops
    # the last 1-2 bits on some rows for 17-significant-digit strings --
    # exactly the pitfall already documented and fixed for envipat.data's
    # isotopes.csv/adducts.csv loaders (see that module). isopattern's
    # golden comparison is an *exact* match, so this matters here more
    # than anywhere else in this file.
    return pd.read_csv(path, float_precision="round_trip")


def _split_by_compound(df: pd.DataFrame) -> dict[str, object]:
    """Reverses tools/make_golden.R's stack_pattern_list()/rbind_fill():
    regroups a stacked golden fixture table back into one value per
    compound, in first-appearance order.

    Columns are aligned by *name*: rbind_fill's column union can end up
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
        cols = [c for c in df.columns if c not in ("compound", "error")]
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
        np.testing.assert_allclose(g["abundance"].to_numpy(), exp["abundance"].to_numpy(), rtol=rtol)


# --------------------------------------------------------------------- #
# Shared formula pools -- replicates tools/make_golden.R's pool10/pool3/
# pool5 construction (see that script's "Shared formula pools" section)
# using this port's own (issue #7, already golden-verified exact) check_
# chemform, rather than asking R for them -- cross-checked once against
# a direct Rscript dump to confirm this reproduces the same 10 formulas
# in the same order (see MR description).
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
# Known issue #3 golden-fixture bug: tools/make_golden.R's rbind_fill()
# silently drops data for any formula that mixes a labelled isotope (a
# "[NN]Elem" or "D" token) with that *same* element's natural, multi-
# isotope form -- e.g. "C14H10N2O1[13]C1D2" (natural C14/H10 alongside
# labelled [13]C1/D2). Both the labelled pseudo-element and the natural
# element contribute a column with the *same* isotope-name string (here,
# "13C" and "2H" each twice) -- confirmed directly against a raw R 4.6.1
# .Call (`colnames(as.data.frame(out2))` lists "13C" twice for this
# formula) -- so this is not an envipat._core or Python-port bug: our
# bindings return the exact same duplicate-named list (see
# _duplicate_isotope_columns below and envipat.isopattern's
# compound_sorted columns for any of these formulas).
#
# The bug is in the *fixture generator*: rbind_fill() builds
# `all_cols <- unique(unlist(lapply(df_list, names)))` and then selects
# `df[all_cols]` for every formula. `unique()` collapses each duplicate
# name to one entry in `all_cols`, and R's `df[, name]`-by-name selection
# then silently returns only the *first* matching column -- discarding
# the natural element's (varying) substitution-count column and
# duplicating the labelled column's (constant) value in its place, for
# every peak. Confirmed by hand against
# tests/golden/isopattern/case_001.csv: this compound's "13C"/"2H"/"15N"/
# "18O" columns are all constant (2) across every peak, while the
# algorithm's own (correct, non-stacked) output for the same formula has
# the natural-element columns genuinely varying (0/1/2...) -- see MR
# description for the full trace, including the direct-.Call comparison
# used to confirm this is a fixture bug, not a core-algorithm or
# Python-port bug.
#
# Rather than hand-list every affected formula (there turn out to be many
# -- any deuterium-labelled formula in the base 107 data(chemforms)
# collides on "2H" the same way), this is detected programmatically below
# for every formula this test file actually uses (POOL10 and the full
# data(chemforms) list), by asking envipat._core itself which isotope
# names it returns more than once for a given formula. m/z and abundance
# are untouched (they aren't name-collided) and are still checked
# normally for every affected compound; only the specific colliding
# columns' *values* are skipped for that compound -- everything else,
# including that the columns exist at all, is still checked.
# --------------------------------------------------------------------- #


def _duplicate_isotope_columns(formula: str) -> set[str]:
    """Isotope-name columns envipat._core.pattern_algo1 returns more than
    once for `formula` -- see the block comment above."""
    elem = ISOTOPES.iloc[:, 0].astype(str).tolist()
    iso_name = ISOTOPES.iloc[:, 1].astype(str).tolist()
    iso_mass = ISOTOPES.iloc[:, 2].to_numpy(dtype=np.float64)
    iso_abu = ISOTOPES.iloc[:, 3].to_numpy(dtype=np.float64)
    result = _core.pattern_algo1(
        sum=formula,
        peak_limit=int(1e6),
        threshold=0.1,
        rel_to_mono=0,
        iso_list_elem=elem,
        iso_list_iso=iso_name,
        iso_list_mass=iso_mass,
        iso_list_abu=iso_abu,
    )
    counts = Counter(result["isotope_names"])
    return {name for name, n in counts.items() if n > 1}


_KNOWN_DUPLICATE_NAME_BUG = {
    formula: dup
    for formula in dict.fromkeys([*POOL10, *CHEMFORMS])
    if (dup := _duplicate_isotope_columns(formula))
}


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
        np.testing.assert_array_equal(g["m/z"].to_numpy(), exp["m/z"].to_numpy())
        np.testing.assert_array_equal(g["abundance"].to_numpy(), exp["abundance"].to_numpy())

        isotope_cols = [c for c in exp.columns if c not in ("m/z", "abundance")]
        assert set(isotope_cols) == set(c for c in g.columns if c not in ("m/z", "abundance"))
        skip_cols = _KNOWN_DUPLICATE_NAME_BUG.get(name, set())
        for col in isotope_cols:
            if col in skip_cols:
                continue
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
        isopattern(
            ISOTOPES, ["C1H4", "C2H6", "C1H4", "C2H6", "C3H8"], threshold=0.1, verbose=False
        )
    assert "C1H4" in str(excinfo.value)
    assert "C2H6" in str(excinfo.value)
    assert "C3H8" not in str(excinfo.value)


# --------------------------------------------------------------------- #
# envelope
# --------------------------------------------------------------------- #

ENVELOPE_CASES = _load_cases("envelope")
# Hardcoded to match tools/make_golden.R's dmz_modes list exactly (its
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
    return dict(ppm=p["ppm"] == "TRUE", dmz=dmz, env=p["env"], resolution=float(p["resolution"]))


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

    _assert_profile_dict_matches(got, expected, rtol=1e-12)


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


@pytest.mark.parametrize("case_id", sorted(CHECK_SEVERAL_CASES))
def test_check_several_matches_golden(case_id, pattern_overlap, pattern_no_overlap):
    row = CHECK_SEVERAL_CASES[case_id]
    p = _parse_params(row["params"])
    scenario = p["scenario"]
    dmz = float(p["dmz"])
    ppm = p["ppm"] == "TRUE"
    pattern = pattern_overlap if "overlap" in scenario and "no_overlap" not in scenario else pattern_no_overlap

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
        assert (result["warning"] == "TRUE").any()
        return

    result = check_several(pattern, dmz=dmz, ppm=ppm)
    expected = pd.read_csv(
        GOLDEN / "check_several" / row["output_file"], keep_default_na=False, dtype=str
    )
    pd.testing.assert_frame_equal(
        result.reset_index(drop=True), expected.reset_index(drop=True), check_dtype=False
    )


# --------------------------------------------------------------------- #
# isowrap
# --------------------------------------------------------------------- #

ISOWRAP_CASES = _load_cases("isowrap")
DUMP_ISOWRAP_RESOLUTION_SCRIPT = REPO_ROOT / "tools" / "dump_isowrap_resolution.R"


@pytest.fixture(scope="session")
def isowrap_resolution_table(tmp_path_factory):
    """getR()-interpolated resolutions for POOL5, per curve -- computed by
    the installed R package itself (see this module's docstring and
    tools/dump_isowrap_resolution.R for why: issue #9's getR port is not
    merged into this branch). Skips if R / the installed package aren't
    available, mirroring tests/test_core.py's full_isotopes fixture."""
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")
    out = tmp_path_factory.mktemp("isowrap") / "resolution.csv"
    proc = subprocess.run(
        ["Rscript", str(DUMP_ISOWRAP_RESOLUTION_SCRIPT), str(out)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        pytest.skip(f"could not dump isowrap resolution via R: {proc.stderr.strip()}")
    return pd.read_csv(out, float_precision="round_trip")


@pytest.mark.parametrize("case_id", sorted(ISOWRAP_CASES))
def test_isowrap_matches_golden(case_id, isowrap_resolution_table):
    row = ISOWRAP_CASES[case_id]
    p = _parse_params(row["params"])
    curve = p["curve"]

    curve_rows = isowrap_resolution_table[isowrap_resolution_table["curve"] == curve]
    resolution = curve_rows.set_index("formula").loc[POOL5, "resolution"].to_numpy(dtype=np.float64)

    checked5 = check_chemform(ISOTOPES, POOL5)

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

    _assert_profile_dict_matches(got, expected, rtol=1e-12)
