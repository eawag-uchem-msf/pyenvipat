"""Bundled-data sanity checks.

The live R-vs-Python end-to-end comparison (issue #10) was removed together
with the R implementation (issue #11); it stays available in git history
(``tools/dump_full_pipeline.R``). The frozen R reference lives on in
``tests/golden/``.

The bundled ``adducts.csv`` deliberately differs from the original
``adducts.rda`` for 3 rows (``"F"``, ``"H2CO3"``, ``"H1N1H4"``, which R's
``mergeform`` mishandles); this checks the patched values are intact.

Run after `pip install -e .`:
    pytest tests/python/test_end_to_end.py
"""

from __future__ import annotations

import envipat.data as data
from envipat.formula import mergeform


def test_adducts_csv_carries_the_three_documented_patches():
    """Confirms the real bundled adducts.csv still carries #7's 3 patched
    Formula_add values, and that mergeform succeeds on them."""
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
        # Just needs to not raise.
        mergeform("C1H4", actual)
