"""Golden-fixture regression tests for envipat._core.enviPat_py.

Compares the compiled C extension's output against the .Call-level golden
fixtures frozen in tests/golden/c_direct/ (issue #3): same inputs, same
compiled algorithm code (iso_core.c and the vendored algorithm files under
src/envipat/_core/, unmodified except for one disclosed preferences.h edit
-- see that file), only the calling language differs. So every comparison
here is an exact equality check (np.testing.assert_array_equal, no
tolerance), not an allclose().

Five of the iso_pattern/iso_pattern_4 cases (case_006..case_010, formula
C10H15Cl2) call the C core with the *full* isotopes table rather than a
committed, reduced one; that table is the packaged
src/envipat/data/isotopes.csv (a full-precision export of R's
data(isotopes)).
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from envipat._core import enviPat_py as core

REPO_ROOT = Path(__file__).resolve().parents[1]
GOLDEN_DIR = REPO_ROOT / "tests" / "golden" / "c_direct"
ISOTOPES_CSV = REPO_ROOT / "src" / "envipat" / "data" / "isotopes.csv"


# --------------------------------------------------------------------- #
# golden fixture I/O helpers
# --------------------------------------------------------------------- #


def _read_csv_columns(path: Path) -> dict[str, list[str]]:
    with path.open(newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        cols: dict[str, list[str]] = {h: [] for h in header}
        for row in reader:
            for h, v in zip(header, row):
                cols[h].append(v)
    return cols


def _float_col(cols: dict[str, list[str]], name: str) -> np.ndarray:
    return np.array([float(v) for v in cols[name]], dtype=np.float64)


def _int_col(cols: dict[str, list[str]], name: str) -> np.ndarray:
    return np.array([int(round(float(v))) for v in cols[name]], dtype=np.int32)


def _parse_params(params: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in params.split(";"):
        if not part:
            continue
        k, v = part.split("=", 1)
        out[k] = v
    return out


def _load_cases() -> dict[str, dict[str, str]]:
    with (GOLDEN_DIR / "cases.csv").open(newline="") as f:
        return {row["case_id"]: row for row in csv.DictReader(f)}


CASES = _load_cases()


def _case_ids(fn: str) -> list[str]:
    return [case_id for case_id, row in CASES.items() if row["fn"] == fn]


# --------------------------------------------------------------------- #
# isotopes table fixtures
# --------------------------------------------------------------------- #


def _load_isotopes_csv(path: Path):
    cols = _read_csv_columns(path)
    return cols["element"], cols["isotope"], _float_col(cols, "mass"), _float_col(cols, "abundance")


@pytest.fixture(scope="session")
def reduced_isotopes():
    """H/C/O only -- committed as case_003_input_isotopes.csv, no R needed."""
    return _load_isotopes_csv(GOLDEN_DIR / "case_003_input_isotopes.csv")


@pytest.fixture(scope="session")
def full_isotopes():
    """The full isotopes table, as packaged with envipat."""
    return _load_isotopes_csv(ISOTOPES_CSV)


# --------------------------------------------------------------------- #
# iso_ppm_Call -> ppm_trace
# --------------------------------------------------------------------- #


@pytest.mark.parametrize("case_id", _case_ids("iso_ppm_Call"))
def test_ppm_trace_matches_golden(case_id):
    row = CASES[case_id]
    p = _parse_params(row["params"])
    expected = _float_col(_read_csv_columns(GOLDEN_DIR / row["output_file"]), "mz")

    got = core.ppm_trace(start=float(p["start"]), end=float(p["end"]), ppm=float(p["ppm"]))

    np.testing.assert_array_equal(got, expected)


# --------------------------------------------------------------------- #
# iso_profile_with_trace_Call -> profile_with_trace
# --------------------------------------------------------------------- #


@pytest.mark.parametrize("case_id", _case_ids("iso_profile_with_trace_Call"))
def test_profile_with_trace_matches_golden(case_id):
    row = CASES[case_id]
    p = _parse_params(row["params"])
    out_cols = _read_csv_columns(GOLDEN_DIR / row["output_file"])
    expected_mz = _float_col(out_cols, "mz")
    expected_abundance = _float_col(out_cols, "abundance")
    trace = _float_col(_read_csv_columns(GOLDEN_DIR / row["input_file"]), "mz")

    got_mz, got_abundance = core.profile_with_trace(
        profile_type=int(p["type1"]),
        mass=np.array([float(p["f1"])]),
        abundance=np.array([float(p["a1"])]),
        trace=trace,
        resolution=int(float(p["r1"])),
        threshold=float(p["t1"]),
        filter=int(p["filter"]),
    )

    np.testing.assert_array_equal(got_mz, expected_mz)
    np.testing.assert_array_equal(got_abundance, expected_abundance)


# --------------------------------------------------------------------- #
# iso_centroid_Call -> centroid
# --------------------------------------------------------------------- #

# case_023..case_025's manifest "input_file" (case_011.csv, a bare ppm mz
# trace with no abundance column at all) does not match what
# tools/make_golden.R actually fed iso_centroid_Call for those cases -- the
# in-memory profile from its own case_017 run (type1=0/filter=1/t1=0). This
# is a provenance-metadata mislabel in the #3 fixture generator script, not
# a defect in the golden output values themselves (those were written from
# the correct in-memory object). Routed around here using case_017.csv
# directly, whose own bit-exactness is independently covered by
# test_profile_with_trace_matches_golden. case_027..case_029's
# "case_026.csv" is correct as recorded.
_CENTROID_PROFILE_INPUT = {
    "gaussian_filter1_t0": "case_017.csv",
    "two_equal_peaks_valley": "case_026.csv",
}


@pytest.mark.parametrize("case_id", _case_ids("iso_centroid_Call"))
def test_centroid_matches_golden(case_id):
    row = CASES[case_id]
    p = _parse_params(row["params"])
    expected_cols = _read_csv_columns(GOLDEN_DIR / row["output_file"])
    expected_mz = _float_col(expected_cols, "mz")
    expected_abundance = _float_col(expected_cols, "abundance")

    input_file = _CENTROID_PROFILE_INPUT[p["profile"]]
    in_cols = _read_csv_columns(GOLDEN_DIR / input_file)
    profile_mass = _float_col(in_cols, "mz")
    profile_abundance = _float_col(in_cols, "abundance")

    got_mz, got_abundance = core.centroid(
        profile_mass=profile_mass, profile_abundance=profile_abundance, type=int(p["t1"])
    )

    np.testing.assert_array_equal(got_mz, expected_mz)
    np.testing.assert_array_equal(got_abundance, expected_abundance)


# --------------------------------------------------------------------- #
# iso_pattern / iso_pattern_4 -> pattern_algo1 / pattern_algo2
# --------------------------------------------------------------------- #

_ALGO = {"iso_pattern": core.pattern_algo1, "iso_pattern_4": core.pattern_algo2}
_PATTERN_CASE_IDS = [case_id for case_id, row in CASES.items() if row["fn"] in _ALGO]


@pytest.mark.parametrize("case_id", _PATTERN_CASE_IDS)
def test_pattern_matches_golden(case_id, reduced_isotopes, full_isotopes):
    row = CASES[case_id]
    p = _parse_params(row["params"])
    algo = _ALGO[row["fn"]]

    if row["input_file"].endswith("_input_isotopes.csv"):
        elem, iso, mass, abu = reduced_isotopes
    elif p["s1"] in ("C2H4", "ZZZ99"):
        # C2H4 only needs H/C (covered by the committed H/C/O reduced
        # table); ZZZ99 is an invalid element symbol and fails regardless
        # of table content. Neither needs the full table.
        elem, iso, mass, abu = reduced_isotopes
    else:
        elem, iso, mass, abu = full_isotopes

    call_kwargs = dict(
        sum=p["s1"],
        peak_limit=int(float(p["pl"])),
        threshold=float(p["t1"]),
        rel_to_mono=int(p["rtm"]),
        iso_list_elem=elem,
        iso_list_iso=iso,
        iso_list_mass=mass,
        iso_list_abu=abu,
    )

    if p["s1"] == "ZZZ99":
        with pytest.raises(core.EnviPatCoreError):
            algo(**call_kwargs)
        expected_cols = _read_csv_columns(GOLDEN_DIR / row["output_file"])
        assert expected_cols["is_null"] == ["TRUE"]
        return

    result = algo(**call_kwargs)

    if p["rica"] == "1":
        expected_cols = _read_csv_columns(GOLDEN_DIR / row["output_file"])
        expected_counts = _int_col(expected_cols, "count")
        assert result["element_names"] == expected_cols["element"]
        np.testing.assert_array_equal(result["element_all_iso_calc_amount"], expected_counts)
        return

    out_cols = _read_csv_columns(GOLDEN_DIR / row["output_file"])
    isotope_names = list(out_cols.keys())[2:]

    expected_mass = _float_col(out_cols, "mass")
    expected_abundance = _float_col(out_cols, "abundance")
    expected_compound = np.array(
        [_int_col(out_cols, name) for name in isotope_names], dtype=np.int32
    ).T

    np.testing.assert_array_equal(result["mass"], expected_mass)
    np.testing.assert_array_equal(result["abundance"], expected_abundance)
    assert result["isotope_names"] == isotope_names
    np.testing.assert_array_equal(result["compound"], expected_compound)
