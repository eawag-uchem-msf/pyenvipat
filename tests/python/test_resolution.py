"""Tests for envipat.resolution.get_resolution (issue #17).

get_resolution fits scipy.interpolate.make_smoothing_spline (GCV-chosen
smoothing) to a resolution_list curve. It is *not* bit-compatible with R's
getR (stats::smooth.spline with fixed nknots/spar, removed in #17): the
comparison against tests/golden/getR/*.csv, frozen from R package v2.8, uses
a loose, explicitly-labeled tolerance (sanity-check magnitude only), not a
tight regression bar.

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
# Loose on purpose: scipy's GCV smoothing deviates from R's fixed-spar fit by
# up to ~17% relative error on the resolution_list curves.
LOOSE_RTOL = 0.25

_COMBO_RE = re.compile(r"resolution_nknots(\d+)_spar([\d.]+)")
_CURVE_RE = re.compile(r"curve=([^;]+);")


def _load_golden_cases():
    """Yields (curve, nknots, spar, mass_grid, expected) for every
    (curve, combo) pair in the frozen golden fixtures."""
    cases = pd.read_csv(GOLDEN_DIR / "cases.csv")
    out = []
    for _, row in cases.iterrows():
        curve = _CURVE_RE.search(row["params"]).group(1)
        # Parse with the same (default) float parser as envipat.data so a
        # curve's boundary mass is identical on both sides; otherwise
        # get_resolution's in-range check can spuriously reject it by 1 ULP.
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
                    curve,
                    nknots,
                    spar,
                    mass_grid,
                    expected,
                    id=f"{curve}/nknots={nknots}/spar={spar}",
                )
            )
    return out


_GOLDEN_CASES = _load_golden_cases()


def _make_checked(masses: np.ndarray) -> pd.DataFrame:
    """A minimal stand-in for check_chemform's output (issue #7), with just
    the columns get_resolution reads: boolean 'warning' and numeric
    'monoisotopic_mass'."""
    return pd.DataFrame(
        {
            "warning": False,
            "new_formula": pd.NA,
            "monoisotopic_mass": masses,
        }
    )


@pytest.mark.parametrize(
    "curve,nknots,spar,mass_grid,expected",
    # One (curve, mass grid, R result) triple per curve; the (nknots, spar)
    # combination only selects which R column to compare to.
    [c for c in _GOLDEN_CASES if c.values[1] == 13 and c.values[2] == 0.1],
)
def test_is_in_the_right_ballpark_vs_r(curve, nknots, spar, mass_grid, expected):
    resmass = data.resolution_list[curve]
    checked = _make_checked(mass_grid)
    result = get_resolution(checked, resmass)

    assert isinstance(result, pd.Series)
    assert result.name == "resolution"
    assert list(result.index) == list(checked.index)

    got = result.to_numpy()
    assert np.all(np.isfinite(got))
    rel = np.abs(got - expected) / np.abs(expected)
    assert rel.max() < LOOSE_RTOL, (
        f"{curve}: max relative error {rel.max():.3e} exceeds the loose "
        f"{LOOSE_RTOL:.0%} sanity bound"
    )


def test_uses_scipy_make_smoothing_spline(monkeypatch):
    calls = []
    real = resolution_module.make_smoothing_spline

    def spy(*args, **kwargs):
        calls.append((args, kwargs))
        return real(*args, **kwargs)

    monkeypatch.setattr(resolution_module, "make_smoothing_spline", spy)

    resmass = data.resolution_list["Elite_R240000@400"]
    get_resolution(_make_checked(np.array([100.0, 500.0])), resmass)
    assert len(calls) == 1


def test_rejects_removed_keyword_arguments():
    resmass = data.resolution_list["Elite_R240000@400"]
    checked = _make_checked(np.array([100.0]))
    for kwargs in ({"nknots": 13}, {"spar": 0.1}, {"method": "scipy"}):
        with pytest.raises(TypeError):
            get_resolution(checked, resmass, **kwargs)


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


def test_rejects_checked_with_warning_formulas():
    resmass = data.resolution_list["Elite_R240000@400"]
    checked = _make_checked(np.array([100.0, 200.0]))
    checked.loc[1, "warning"] = True
    with pytest.raises(ValueError, match="warning"):
        get_resolution(checked, resmass)
