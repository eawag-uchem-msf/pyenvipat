"""The version is defined once, as envipat.__version__ (see CONTRIBUTING.md)."""

import re
from importlib.metadata import version
from pathlib import Path

import pytest

import envipat

CITATION = Path(__file__).resolve().parents[2] / "CITATION.cff"


def test_version_matches_distribution_metadata():
    # Distribution name (pyEnvipat) and import name (envipat) differ.
    assert envipat.__version__ == version("pyEnvipat")


def test_citation_version_matches():
    if not CITATION.is_file():
        pytest.skip("CITATION.cff not available (tests run outside the source tree)")
    match = re.search(r"^version:\s*(\S+)\s*$", CITATION.read_text(encoding="utf-8"), re.M)
    assert match is not None, "no version field in CITATION.cff"
    assert match.group(1) == envipat.__version__
