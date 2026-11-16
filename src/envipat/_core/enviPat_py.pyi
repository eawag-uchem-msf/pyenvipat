"""Type stub for the compiled extension ``envipat._core.enviPat_py``.

Signatures mirror the ``PyArg_ParseTupleAndKeywords`` format strings and the
return values built in ``enviPat_py.c``; keep both in sync.
"""

from collections.abc import Sequence
from typing import TypedDict

import numpy as np
import numpy.typing as npt

_F64 = npt.NDArray[np.float64]
_I32 = npt.NDArray[np.int32]
_ArrayLikeF64 = npt.ArrayLike

class EnviPatCoreError(Exception):
    """Raised when an underlying ``iso_core_*`` call reports a failure."""

class PatternResult(TypedDict):
    """Dict returned by :func:`pattern_algo1` / :func:`pattern_algo2`."""

    mass: _F64
    abundance: _F64
    compound: _I32  # 2-D: one row per peak, one column per isotope
    isotope_names: list[str]
    element_names: list[str]
    element_all_iso_calc_amount: _I32

def ppm_trace(start: float, end: float, ppm: float) -> _F64:
    """ppm-spaced mass trace between ``start`` and ``end``."""

def centroid(
    profile_mass: _ArrayLikeF64, profile_abundance: _ArrayLikeF64, type: int
) -> tuple[_F64, _F64]:
    """Centroid (0), intensoid (1) or valley (2) detection on a profile."""

def profile_with_trace(
    profile_type: int,
    mass: _ArrayLikeF64,
    abundance: _ArrayLikeF64,
    trace: _ArrayLikeF64,
    resolution: int,
    threshold: float,
    filter: int,
) -> tuple[_F64, _F64]:
    """Convolve a pattern onto a mass trace (0=Gaussian, 1=Cauchy-Lorentz)."""

def pattern_algo1(
    sum: str,
    peak_limit: int,
    threshold: float,
    rel_to_mono: int,
    iso_list_elem: Sequence[str],
    iso_list_iso: Sequence[str],
    iso_list_mass: _ArrayLikeF64,
    iso_list_abu: _ArrayLikeF64,
) -> PatternResult:
    """Isotope pattern calculation, algorithm 1."""

def pattern_algo2(
    sum: str,
    peak_limit: int,
    threshold: float,
    rel_to_mono: int,
    iso_list_elem: Sequence[str],
    iso_list_iso: Sequence[str],
    iso_list_mass: _ArrayLikeF64,
    iso_list_abu: _ArrayLikeF64,
) -> PatternResult:
    """Isotope pattern calculation, algorithm 2."""
