"""Golden-fixture regression tests for envipat.resolution (issue #9).

get_resolution is a from-scratch Python reimplementation of R's
stats::smooth.spline (fixed spar, fixed nknots) as used by R's getR
(R/getR.R) -- see src/envipat/resolution.py's module docstring for why no
existing scipy function is a drop-in substitute, and for the exact R
internals (Fortran sbart.c/stxwx.f/sgram.f/sslvrg.f) this mirrors.

These tests compare against tests/golden/getR/*.csv, frozen from the real R
package (v2.8) in issue #3 (see tools/make_golden.R): for each of the 31
resolution_list curves and four (nknots, spar) combinations, R's getR was
evaluated at a 15-point mass grid spanning that curve's own mass range.

Run inside the dev container, after `pip install -e .`:
    pytest tests/python
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import envipat.data as data
from envipat.resolution import get_resolution

GOLDEN_DIR = Path(__file__).resolve().parents[2] / "tests" / "golden" / "getR"
RTOL = 1e-6

_COMBO_RE = re.compile(r"resolution_nknots(\d+)_spar([\d.]+)")
_CURVE_RE = re.compile(r"curve=([^;]+);")


def _load_golden_cases():
    """Yields (curve, nknots, spar, mass_grid, expected) for every
    (curve, combo) pair in the frozen golden fixtures."""
    cases = pd.read_csv(GOLDEN_DIR / "cases.csv")
    out = []
    for _, row in cases.iterrows():
        curve = _CURVE_RE.search(row["params"]).group(1)
        df = pd.read_csv(GOLDEN_DIR / row["output_file"])
        mass_grid = df["mass"].to_numpy(dtype=float)
        for col in df.columns:
            if col == "mass":
                continue
            m = _COMBO_RE.match(col)
            nknots, spar = int(m.group(1)), float(m.group(2))
            expected = df[col].to_numpy(dtype=float)
            out.append(
                pytest.param(
                    curve, nknots, spar, mass_grid, expected,
                    id=f"{curve}/nknots={nknots}/spar={spar}",
                )
            )
    return out


_GOLDEN_CASES = _load_golden_cases()


def _make_checked(masses: np.ndarray) -> pd.DataFrame:
    """A minimal stand-in for check_chemform's output (issue #7), with just
    the columns get_resolution reads: boolean 'warning' and numeric
    'monoisotopic_mass'."""
    return pd.DataFrame({
        "warning": False,
        "new_formula": pd.NA,
        "monoisotopic_mass": masses,
    })


@pytest.mark.parametrize(
    "curve,nknots,spar,mass_grid,expected", _GOLDEN_CASES
)
def test_matches_r_golden_fixture(curve, nknots, spar, mass_grid, expected):
    resmass = data.resolution_list[curve]
    checked = _make_checked(mass_grid)
    result = get_resolution(checked, resmass, nknots=nknots, spar=spar)

    assert isinstance(result, pd.Series)
    assert result.name == "resolution"
    assert list(result.index) == list(checked.index)

    got = result.to_numpy()
    rel = np.abs(got - expected) / np.abs(expected)
    assert rel.max() < RTOL, (
        f"{curve} nknots={nknots} spar={spar}: max relative error "
        f"{rel.max():.3e} >= {RTOL:.0e}"
    )


def test_default_nknots_and_spar_match_getR_defaults():
    # R/getR.R: function(checked, resmass, nknots = 13, spar = 0.1, ...)
    resmass = data.resolution_list["Elite_R240000@400"]
    checked = _make_checked(np.array([100.0, 500.0, 1000.0]))
    explicit = get_resolution(checked, resmass, nknots=13, spar=0.1)
    default = get_resolution(checked, resmass)
    pd.testing.assert_series_equal(explicit, default)


def test_rejects_too_few_resmass_points():
    resmass = pd.DataFrame({"mz": range(5), "R": range(5)})
    checked = _make_checked(np.array([2.0]))
    with pytest.raises(ValueError, match="not enough data points"):
        get_resolution(checked, resmass)


def test_rejects_out_of_range_mass():
    resmass = data.resolution_list["Elite_R240000@400"]
    checked = _make_checked(np.array([resmass["mz"].max() + 1.0]))
    with pytest.raises(ValueError, match="out of range"):
        get_resolution(checked, resmass)


def test_rejects_invalid_nknots():
    resmass = data.resolution_list["Elite_R240000@400"]
    checked = _make_checked(np.array([100.0]))
    with pytest.raises(ValueError, match="invalid nknots"):
        get_resolution(checked, resmass, nknots=2)
    with pytest.raises(ValueError, match="invalid nknots"):
        get_resolution(checked, resmass, nknots=len(resmass) + 1)


def test_rejects_invalid_spar():
    resmass = data.resolution_list["Elite_R240000@400"]
    checked = _make_checked(np.array([100.0]))
    with pytest.raises(ValueError, match="invalid spar"):
        get_resolution(checked, resmass, spar=0.0)
    with pytest.raises(ValueError, match="invalid spar"):
        get_resolution(checked, resmass, spar=1.0)


def test_rejects_checked_with_warning_formulas():
    resmass = data.resolution_list["Elite_R240000@400"]
    checked = _make_checked(np.array([100.0, 200.0]))
    checked.loc[1, "warning"] = True
    with pytest.raises(ValueError, match="warning"):
        get_resolution(checked, resmass)
