"""End-to-end ``isowrap`` pipeline parity tests (issue #10).

This is the first comparison in this port that exercises the *whole*
``check_chemform -> get_resolution -> isopattern -> envelope -> vdetect``
chain with each side deriving its own resolution independently:

- R: ``tools/dump_full_pipeline.R`` runs the real R package's
  ``isowrap()``, which calls ``getR()`` *internally* from a
  ``resolution_list`` calibration table.
- Python: this file calls :func:`envipat.resolution.get_resolution`
  itself (issue #9's from-scratch ``getR`` port) and feeds its result
  into :func:`envipat.envelope.envelope` via ``isowrap``'s own documented
  extension point (``envipat.isowrap``'s module docstring: "callers
  wanting curve-based interpolation call ``get_resolution`` themselves
  ... and pass its result in as ``resolution``").

Earlier golden-fixture tests (``tests/python/test_pipeline.py``'s
``isowrap`` cases, issue #8) fed the *same*, R-dumped resolution number
into both languages' ``envelope()`` call -- so they only verified
``isowrap``'s own composition. This file is the first to compare fully
independent resolution derivations end to end, closing the loop issue
#8's own scope note deferred to issue #9 ("isowrap and getR" in MR !8).

Cross-product covered (per issue #10's own scope): 8 formulas (spanning
the mass range common to *every* ``resolution_list`` curve, see
``tools/dump_full_pipeline.R``'s own comment for why these were picked)
x all 31 ``resolution_list`` curves x ``algo`` in {1, 2} x ``detect`` in
{"centroid", "intensoid", "valley"} = 1488 (formula, curve, algo, detect)
combinations.

Tolerances
----------
- ``resolution``: ``rtol=1e-6``, the bound already established for
  :func:`envipat.resolution.get_resolution` against R's ``getR`` in
  issue #9. Measured max relative deviation across all 248
  (curve, formula) pairs: **~2.6e-14** (see MR description for the exact
  figure) -- comfortably inside the bound, consistent with issue #9's
  own measurements.
- ``m/z``: ``rtol=1e-12``, the bound issue #8 established for
  ``envelope``/``vdetect`` given a *shared* resolution value. Measured
  max relative deviation across all 1488 cases: **~8.9e-16** -- still
  comfortably inside the tighter bound even with independently-derived
  resolutions on each side, since m/z (peak *position*) is far less
  sensitive to the mass-trace grid than abundance is.
- ``abundance``: **not** held to issue #8's ``rtol=1e-12`` here, and
  that is disclosed rather than silently loosened. Unlike issue #8's
  isowrap tests, the two sides' resolution values differ by a real
  (tiny, ~1e-14 relative, well inside the established ``rtol=1e-6``)
  amount; ``envelope``'s ``dmz="get"`` branch derives its mass-sampling
  step from that resolution (``dmz2 = mean(m) / resolution``), and
  Gaussian profile abundance is locally steep enough near peak
  maxima/minima that this residual perturbs the *sampled* abundance
  more than it perturbs peak position. Measured worst case among
  non-degenerate (not near-double-underflow) abundances: **~2.9e-11**
  relative deviation -- see MR description for the full distribution.
  ``detect="valley"`` additionally produces some genuinely
  near-zero local-minimum abundances (down to ~1e-93 in this matrix),
  where *relative* error is meaningless (a ~1e-102 absolute difference
  between two ~1e-93 numbers is not a real disagreement). Comparison
  therefore uses a combined tolerance (``numpy.isclose`` semantics:
  ``|diff| <= atol + rtol * |expected|``) with ``rtol=1e-9,
  atol=1e-9`` -- comfortably tighter than "loose" (>30x margin over the
  measured non-degenerate worst case) while correctly handling the
  near-zero ``valley`` tail.

Adducts.csv/adducts.rda divergence (issue #10's own caveat)
-------------------------------------------------------------
This comparison's formula set is 8 hand-picked ``chemforms`` entries
(chosen for mass-range coverage, not adduct arithmetic) -- it never
touches ``envipat.data.adducts`` or ``mergeform``/``subform``/etc., so
the known ``adducts.csv`` vs. ``adducts.rda`` divergence for the 3
corrected ``Formula_add`` rows (``"F"``, ``"H2CO3"``, ``"H1N1H4"`` --
see issue #10's own comment thread and MR !7) cannot appear in *this*
file's comparison. That divergence is already covered by
``tests/python/test_formula.py``'s
``test_flat_helpers_reject_missing_explicit_count`` and
``test_mergeform_rejects_repeated_element_in_one_formula`` (issue #7);
``test_adducts_csv_carries_the_three_documented_patches`` below is a
light confirmation, against the real bundled ``adducts.csv`` (not
hardcoded strings), that those 3 rows still carry the patched values and
that Python's ``mergeform`` still succeeds on them -- without invoking R
on them here (R hangs on the unpatched ``"F"`` row, see the caveat).

Skips (not fails) when R is unavailable
-----------------------------------------
Same pattern as ``tests/test_core.py``'s ``full_isotopes`` fixture and
``tests/python/test_pipeline.py``'s ``isowrap_resolution_table``: the
``r_reference`` fixture below calls ``shutil.which("Rscript")`` and
``pytest.skip``s if it's ``None`` or the dump script errors, rather than
failing. Every other test in this repository's suite (``pytest
tests/test_core.py tests/python``, 303 cases as of this issue) needs no R
at all -- see issue #10's MR description for how that "no R installed"
claim was verified concretely (a separate, R-less container).

Run inside the dev container, after `pip install -e .`:
    pytest tests/python/test_end_to_end.py
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import envipat.data as data
from envipat.envelope import envelope
from envipat.formula import check_chemform, mergeform
from envipat.isopattern import isopattern
from envipat.resolution import get_resolution
from envipat.vdetect import vdetect

REPO_ROOT = Path(__file__).resolve().parents[2]
DUMP_SCRIPT = REPO_ROOT / "tools" / "dump_full_pipeline.R"

# Must stay in sync with tools/dump_full_pipeline.R's own `formulas` --
# see that script's comment for why these 8 were picked.
FORMULAS = [
    "C4H5N5D6",
    "C6H10Cl1N2[15]N3",
    "C7H7D5Cl1N5",
    "C10H8Cl1N3O1",
    "C14H10N2O1[13]C1D2",
    "C14H16Cl1N3O1",
    "C14H13N1O7S1",
    "C38H72N2O12",
]
ALGOS = (1, 2)
DETECT_MODES = ("centroid", "intensoid", "valley")

RESOLUTION_RTOL = 1e-6  # established in issue #9
MZ_RTOL = 1e-12  # established in issue #8, for a *shared* resolution
ABUNDANCE_RTOL = 1e-9  # see module docstring for why this differs from #8
ABUNDANCE_ATOL = 1e-9  # handles detect="valley"'s near-zero local minima


@pytest.fixture(scope="session")
def r_reference(tmp_path_factory):
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not available")
    outdir = tmp_path_factory.mktemp("full_pipeline")
    proc = subprocess.run(
        ["Rscript", str(DUMP_SCRIPT), str(outdir)],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        pytest.skip(f"could not dump full pipeline reference via R: {proc.stderr.strip()}")
    resolution = pd.read_csv(outdir / "resolution.csv", float_precision="round_trip")
    results = pd.read_csv(outdir / "results.csv", float_precision="round_trip")
    return resolution, results


@pytest.fixture(scope="session")
def checked_formulas():
    checked = check_chemform(data.isotopes, FORMULAS)
    assert not checked["warning"].any(), "one or more FORMULAS failed check_chemform"
    return checked


@pytest.fixture(scope="session")
def patterns(checked_formulas):
    return {
        algo: isopattern(
            data.isotopes, checked_formulas["new_formula"],
            threshold=0.1, charge=1, algo=algo, verbose=False,
        )
        for algo in ALGOS
    }


@pytest.fixture(scope="session")
def comparison(r_reference, checked_formulas, patterns):
    """Runs the real Python pipeline (get_resolution -> isopattern ->
    envelope -> vdetect) once for the full cross-product and returns a
    long DataFrame of per-case results/deviations against the R
    reference. See module docstring for the tolerance rationale.
    """
    resolution_r, results_r = r_reference
    new_formulas = checked_formulas["new_formula"].tolist()
    curves = resolution_r["curve"].unique().tolist()
    assert len(curves) == 31, f"expected all 31 resolution_list curves, got {len(curves)}"

    rows = []
    for curve in curves:
        resmass = data.resolution_list[curve]
        # getR (both languages) rejects nknots > nrow(resmass); 2 of the
        # 31 curves have only 12 calibration rows, below getR's own
        # nknots=13 default -- see tools/dump_full_pipeline.R's comment.
        nknots = min(13, len(resmass))
        resolution_py = get_resolution(checked_formulas, resmass, nknots=nknots, spar=0.1)

        rsub = resolution_r[resolution_r["curve"] == curve].set_index("formula")["resolution"]
        for idx, fm in zip(checked_formulas.index, new_formulas):
            r_val = float(rsub.loc[fm])
            p_val = float(resolution_py.loc[idx])
            rows.append({
                "kind": "resolution", "curve": curve, "algo": None, "detect": None,
                "formula": fm, "rel_dev": abs(p_val - r_val) / abs(r_val),
                "abs_dev": abs(p_val - r_val), "error_mismatch": False, "len_mismatch": False,
            })

        for algo in ALGOS:
            profiles = envelope(
                patterns[algo], ppm=False, dmz="get", frac=0.25, env="Gaussian",
                resolution=resolution_py, verbose=False,
            )
            for detect in DETECT_MODES:
                centro = vdetect(profiles, detect=detect, verbose=False)
                for fm in new_formulas:
                    sub_r = results_r[
                        (results_r["curve"] == curve)
                        & (results_r["algo"] == algo)
                        & (results_r["detect"] == detect)
                        & (results_r["formula"] == fm)
                    ]
                    is_error_r = bool(sub_r["is_error"].iloc[0]) if len(sub_r) else True
                    py_val = centro[fm]
                    is_error_py = isinstance(py_val, str)
                    base = {"curve": curve, "algo": algo, "detect": detect, "formula": fm}
                    if is_error_r or is_error_py:
                        rows.append({
                            **base, "kind": "peaks", "rel_dev": np.nan, "abs_dev": np.nan,
                            "error_mismatch": is_error_r != is_error_py, "len_mismatch": False,
                        })
                        continue

                    exp_mz = sub_r.sort_values("row")["mz"].to_numpy(dtype=float)
                    exp_ab = sub_r.sort_values("row")["abundance"].to_numpy(dtype=float)
                    got_mz = py_val.iloc[:, 0].to_numpy(dtype=float)
                    got_ab = py_val.iloc[:, 1].to_numpy(dtype=float)

                    if len(exp_mz) != len(got_mz):
                        rows.append({
                            **base, "kind": "peaks", "rel_dev": np.nan, "abs_dev": np.nan,
                            "error_mismatch": False, "len_mismatch": True,
                        })
                        continue

                    mz_rel = np.abs(got_mz - exp_mz) / np.abs(exp_mz)
                    ab_abs = np.abs(got_ab - exp_ab)
                    ab_rel = ab_abs / np.where(exp_ab == 0, 1.0, np.abs(exp_ab))
                    ab_close = np.isclose(got_ab, exp_ab, rtol=ABUNDANCE_RTOL, atol=ABUNDANCE_ATOL)

                    rows.append({
                        **base, "kind": "mz", "rel_dev": float(mz_rel.max()),
                        "abs_dev": float(np.abs(got_mz - exp_mz).max()),
                        "error_mismatch": False, "len_mismatch": False,
                    })
                    rows.append({
                        **base, "kind": "abundance", "rel_dev": float(ab_rel.max()),
                        "abs_dev": float(ab_abs.max()), "error_mismatch": False,
                        "len_mismatch": False, "abundance_close": bool(ab_close.all()),
                    })
    return pd.DataFrame(rows)


def test_full_cross_product_covers_every_combination(comparison):
    n_formulas = len(FORMULAS)
    assert len(comparison[comparison["kind"] == "resolution"]) == 31 * n_formulas
    assert len(comparison[comparison["kind"].isin(["mz", "peaks"])]) == 31 * len(ALGOS) * len(DETECT_MODES) * n_formulas


def test_no_error_or_peak_count_mismatches(comparison):
    mism = comparison[comparison["error_mismatch"] | comparison["len_mismatch"]]
    assert mism.empty, (
        f"{len(mism)} (curve, algo, detect, formula) case(s) disagree on whether R/Python "
        f"raised an error, or on peak count:\n{mism.to_string()}"
    )


def test_resolution_matches_r_getR(comparison):
    res = comparison[comparison["kind"] == "resolution"]
    worst = res["rel_dev"].max()
    assert worst < RESOLUTION_RTOL, (
        f"max resolution relative deviation {worst:.3e} >= established bound {RESOLUTION_RTOL:.0e}"
    )


def test_mz_matches_r(comparison):
    mz = comparison[comparison["kind"] == "mz"]
    worst = mz["rel_dev"].max()
    assert worst < MZ_RTOL, (
        f"max m/z relative deviation {worst:.3e} >= established bound {MZ_RTOL:.0e}"
    )


def test_abundance_matches_r_within_justified_tolerance(comparison):
    ab = comparison[comparison["kind"] == "abundance"]
    fails = ab[ab["abundance_close"] == False]  # noqa: E712 (mixed dtype column; avoid `~` on object dtype)
    assert fails.empty, (
        f"{len(fails)} case(s) exceed rtol={ABUNDANCE_RTOL:.0e}/atol={ABUNDANCE_ATOL:.0e} "
        f"(numpy.isclose semantics):\n"
        f"{fails.sort_values('rel_dev', ascending=False).head(10).to_string()}"
    )


def test_adducts_csv_carries_the_three_documented_patches():
    """Confirms the real bundled adducts.csv still carries #7's 3 patched
    Formula_add values, and that mergeform succeeds on them -- see module
    docstring for why this doesn't need (and deliberately doesn't run) a
    live R comparison here."""
    adducts = data.adducts
    patched = {
        "M+F": "F1",
        "M+H2CO3-H": "H2C1O3",
        "M+H+NH4": "H5N1",
    }
    for name, expected_formula_add in patched.items():
        row = adducts[adducts["Name"] == name]
        assert len(row) == 1, f"expected exactly one {name!r} row in adducts.csv"
        actual = row["Formula_add"].iloc[0]
        assert actual == expected_formula_add, (
            f"{name!r} Formula_add: expected patched {expected_formula_add!r}, got {actual!r}"
        )
        # Doesn't hang/corrupt like R's unpatched data/adducts.rda would
        # (see module docstring) -- just needs to not raise.
        mergeform("C1H4", actual)
