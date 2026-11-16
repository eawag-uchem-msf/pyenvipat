"""Spot-check tests for envipat.data (issue #6).

These are ad hoc verification tests for the bundled CSV/TXT datasets and
the Python loaders (src/envipat/data.py). They
check dataset shapes and a handful of exact values that were independently
read off the original R data() objects (see the commit message / PR
description for the R-side sprintf("%.17g", .) values these were compared
against).

Run inside the dev container, after `pip install -e .`:
    pytest tests/python
"""

import math

import pandas as pd

import envipat.data as data


def test_isotopes_shape():
    assert data.isotopes.shape == (309, 5)
    assert list(data.isotopes.columns) == [
        "element",
        "isotope",
        "mass",
        "abundance",
        "ratioC",
    ]


def test_isotopes_exact_values():
    # R: isotopes[1, ] == H / 1H / mass=1.007825032 / abundance=0.99988500000000002 / ratioC=6
    row0 = data.isotopes.iloc[0]
    assert row0["element"] == "H"
    assert row0["isotope"] == "1H"
    assert row0["mass"] == 1.007825032
    assert row0["abundance"] == 0.99988500000000002
    assert int(row0["ratioC"]) == 6

    # R: isotopes[309, "mass"] == 34.96885271 (last row, Cl; [81]Br was inserted before it)
    assert data.isotopes.iloc[308]["mass"] == 34.96885271


def test_chemforms_shape_and_values():
    assert data.chemforms.shape == (107,)
    assert isinstance(data.chemforms, pd.Series)
    # R: chemforms[1] == "C4H4Cl1N3O1"
    assert data.chemforms.iloc[0] == "C4H4Cl1N3O1"


def test_adducts_shape_and_values():
    assert data.adducts.shape == (51, 9)
    assert list(data.adducts.columns) == [
        "Name",
        "calc",
        "Charge",
        "Mult",
        "Mass",
        "Ion_mode",
        "Formula_add",
        "Formula_ded",
        "Multi",
    ]
    # R: adducts[1, "Name"] == "2M-H", adducts[1, "Mass"] == -1.0072760000000001
    row0 = data.adducts.iloc[0]
    assert row0["Name"] == "2M-H"
    assert row0["Mass"] == -1.0072760000000001


def test_resolution_list_shape_and_values():
    assert isinstance(data.resolution_list, dict)
    assert len(data.resolution_list) == 31
    for df in data.resolution_list.values():
        assert list(df.columns) == ["mz", "R"]
        assert 12 <= len(df) <= 196

    # R: names(resolution_list)[19] == "OTFusion,QExactiveHF_480000@200"
    # (name contains a literal comma - exercises CSV quoting/round-trip)
    curve = data.resolution_list["OTFusion,QExactiveHF_480000@200"]
    assert len(curve) == 34
    # R: resolution_list[[19]][1, ] == mz=50, R=957196
    assert curve.iloc[0]["mz"] == 50
    assert curve.iloc[0]["R"] == 957196


def test_float_precision_survives_export():
    # abundance is a value that would be silently truncated by write.csv's
    # default 7-significant-digit formatting (0.9998850 vs the true
    # 0.99988500000000002) - assert full round-trip precision.
    abundance = data.isotopes.iloc[0]["abundance"]
    assert math.isclose(abundance, 0.99988500000000002, rel_tol=0, abs_tol=0)


def test_lazy_loading_returns_same_cached_object():
    # __getattr__ + lru_cache: repeated access should not re-parse the file.
    assert data.isotopes is data.isotopes
    assert data.adducts is data.adducts


def test_isotopes_br81_label():
    # Labelled single-isotope [81]Br (issue #29), inserted before the final Cl row.
    row = data.isotopes[data.isotopes["element"] == "[81]Br"]
    assert len(row) == 1
    row = row.iloc[0]
    assert row["isotope"] == "81Br"
    assert row["mass"] == 80.9162897
    assert row["abundance"] == 1
    assert int(row["ratioC"]) == 2


def test_isopattern_br81_label():
    from envipat import isopattern

    res = isopattern(data.isotopes, ["C6H5[81]Br1"], threshold=0.0001, verbose=False)
    pat = res["C6H5[81]Br1"]
    # Single isotope per element position => the monoisotopic peak is the 81Br one.
    assert abs(pat["m/z"].iloc[0] - (6 * 12 + 5 * 1.007825032 + 80.9162897)) < 1e-3
