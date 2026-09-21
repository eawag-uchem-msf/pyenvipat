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
import envipat.resolution as resolution_module
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
        # float_precision="round_trip": pandas' default C float parser does
        # not always round-trip the 17-significant-digit strings these
        # fixtures were written with bit-exactly (same class of bug fixed
        # in envipat.data (#7) for isotopes.csv/adducts.csv). Without this,
        # two curves' (Sciex_QTOFX500R_R25000@200,
        # Agilent_QTOF6550_low_extended_highSens_R25000@200) own boundary
        # mass -- parsed here from this fixture file -- comes out 1 ULP
        # higher than the identical value parsed from resolution_list.csv
        # via envipat.data (which already uses round_trip), so
        # get_resolution's own in-range check spuriously rejects the exact
        # boundary mass. Discovered while building issue #10's end-to-end
        # comparison.
        df = pd.read_csv(GOLDEN_DIR / row["output_file"], float_precision="round_trip")
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


# --- method="scipy" opt-in backend --------------------------------------
#
# These deliberately do NOT compare against tests/golden/getR/*.csv at the
# rtol=1e-6 (or tighter) bar used above: MR !6 measured
# scipy.interpolate.make_smoothing_spline missing R's ground truth by up to
# ~17% relative error (mean ~6.7%), so it cannot meet that bar by design,
# not by a bug in this port. The regression check below uses a loose,
# explicitly-labeled tolerance (sanity-check magnitude only), not the
# port's real accuracy target.
_SCIPY_LOOSE_RTOL = 0.25  # loose on purpose -- see comment above


def test_scipy_method_returns_series_with_right_shape_and_index():
    resmass = data.resolution_list["Elite_R240000@400"]
    checked = _make_checked(np.array([100.0, 500.0, 1000.0, 1500.0]))
    result = get_resolution(checked, resmass, method="scipy")

    assert isinstance(result, pd.Series)
    assert result.name == "resolution"
    assert list(result.index) == list(checked.index)
    assert len(result) == len(checked)
    assert np.all(np.isfinite(result.to_numpy()))


def test_scipy_method_uses_scipy_make_smoothing_spline(monkeypatch):
    calls = []
    real = resolution_module.make_smoothing_spline

    def spy(*args, **kwargs):
        calls.append((args, kwargs))
        return real(*args, **kwargs)

    monkeypatch.setattr(resolution_module, "make_smoothing_spline", spy)

    resmass = data.resolution_list["Elite_R240000@400"]
    checked = _make_checked(np.array([100.0, 500.0]))
    get_resolution(checked, resmass, method="scipy")

    assert len(calls) == 1


def test_scipy_method_default_path_is_unaffected():
    # method="r_exact" (the default) must never touch make_smoothing_spline.
    resmass = data.resolution_list["Elite_R240000@400"]
    checked = _make_checked(np.array([100.0, 500.0]))

    def boom(*args, **kwargs):
        raise AssertionError("make_smoothing_spline must not be called for r_exact")

    orig = resolution_module.make_smoothing_spline
    resolution_module.make_smoothing_spline = boom
    try:
        get_resolution(checked, resmass)  # default method="r_exact"
    finally:
        resolution_module.make_smoothing_spline = orig


@pytest.mark.parametrize(
    "curve,nknots,spar,mass_grid,expected",
    # Reuse a slice of the golden cases as a source of (curve, mass_grid,
    # expected-from-R) triples; nknots/spar are ignored by method="scipy".
    [c for c in _GOLDEN_CASES if c.values[1] == 13 and c.values[2] == 0.1][:5],
)
def test_scipy_method_is_in_the_right_ballpark_vs_r(
    curve, nknots, spar, mass_grid, expected
):
    """Not a tight-tolerance regression test (see module comment above): just
    confirms the scipy backend produces plausible, same-order-of-magnitude
    resolutions, since MR !6 already established it cannot hit rtol=1e-6."""
    resmass = data.resolution_list[curve]
    checked = _make_checked(mass_grid)
    result = get_resolution(checked, resmass, method="scipy")

    got = result.to_numpy()
    rel = np.abs(got - expected) / np.abs(expected)
    assert rel.max() < _SCIPY_LOOSE_RTOL, (
        f"{curve}: scipy backend max relative error {rel.max():.3e} "
        f"exceeds even the loose {_SCIPY_LOOSE_RTOL:.0%} sanity bound"
    )


def test_scipy_method_rejects_non_default_nknots():
    resmass = data.resolution_list["Elite_R240000@400"]
    checked = _make_checked(np.array([100.0]))
    with pytest.raises(ValueError, match="meaningless for method='scipy'"):
        get_resolution(checked, resmass, method="scipy", nknots=6)


def test_scipy_method_rejects_non_default_spar():
    resmass = data.resolution_list["Elite_R240000@400"]
    checked = _make_checked(np.array([100.0]))
    with pytest.raises(ValueError, match="meaningless for method='scipy'"):
        get_resolution(checked, resmass, method="scipy", spar=0.5)


def test_scipy_method_still_rejects_out_of_range_mass():
    resmass = data.resolution_list["Elite_R240000@400"]
    checked = _make_checked(np.array([resmass["mz"].max() + 1.0]))
    with pytest.raises(ValueError, match="out of range"):
        get_resolution(checked, resmass, method="scipy")


def test_scipy_method_still_rejects_checked_with_warning_formulas():
    resmass = data.resolution_list["Elite_R240000@400"]
    checked = _make_checked(np.array([100.0, 200.0]))
    checked.loc[1, "warning"] = True
    with pytest.raises(ValueError, match="warning"):
        get_resolution(checked, resmass, method="scipy")


def test_invalid_method_raises():
    resmass = data.resolution_list["Elite_R240000@400"]
    checked = _make_checked(np.array([100.0]))
    with pytest.raises(ValueError, match="invalid method"):
        get_resolution(checked, resmass, method="not_a_real_method")
