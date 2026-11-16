"""Lazy-loaded access to enviPat's bundled reference datasets.

Mirrors the R package's ``data(isotopes)``, ``data(chemforms)``,
``data(adducts)`` and ``data(resolution_list)`` objects. The plain-text
CSV/TXT files backing this module live under ``envipat/data/`` (a plain
resource directory, not a sub-package) and were exported once from the R
package's datasets, with floating point columns written at full precision
("%.17g") so every double round-trips bit-exactly.

Each dataset is exposed as a module-level attribute (``isotopes``,
``chemforms``, ``adducts``, ``resolution_list``) but is only read and parsed
from disk on first access, via module ``__getattr__`` (PEP 562) - not at
import time.

Data files and shapes
----------------------
- ``isotopes.csv``: columns ``element, isotope, mass, abundance, ratioC``,
  loaded as a :class:`pandas.DataFrame` (309 rows).
- ``chemforms.txt``: one chemical formula per line, loaded as a
  :class:`pandas.Series` of strings (107 elements).
- ``adducts.csv``: 9 columns matching R's ``adducts`` -
  ``Name, calc, Charge, Mult, Mass, Ion_mode, Formula_add, Formula_ded,
  Multi`` - loaded as a :class:`pandas.DataFrame` (51 rows).
- ``resolution_list.csv``: long format, columns ``name, mz, R`` (one row per
  curve point; the R side's literal ``m/z`` column name is renamed to ``mz``
  for CSV/pandas friendliness). Loaded as a ``dict[str, pandas.DataFrame]``
  keyed by instrument/resolution curve name, each value a two-column
  (``mz``, ``R``) DataFrame, in the same curve order as the R
  ``resolution_list`` (31 curves).

Examples
--------
>>> import envipat.data as data
>>> data.isotopes.shape
(309, 5)
>>> data.chemforms.iloc[0]
'C4H4Cl1N3O1'
>>> list(data.adducts.columns)
['Name', 'calc', 'Charge', 'Mult', 'Mass', 'Ion_mode', 'Formula_add', 'Formula_ded', 'Multi']
>>> data.resolution_list["Elite_R240000@400"].columns.tolist()
['mz', 'R']
"""

from __future__ import annotations

from functools import lru_cache
from importlib import resources

import pandas as pd

_ANCHOR = "envipat"
_DATA_SUBDIR = "data"


def _data_path(filename: str):
    """Return an importlib.resources Traversable for a file under
    ``envipat/data/``. Works both from a repo checkout and an installed
    (including zipped) wheel."""
    return resources.files(_ANCHOR) / _DATA_SUBDIR / filename


@lru_cache(maxsize=1)
def _load_isotopes() -> pd.DataFrame:
    # float_precision="round_trip": these CSVs were exported with sprintf
    # ("%.17g", .) specifically to round-trip every double bit-exactly;
    # pandas' default C float parser ("low") does not
    # guarantee that for 17-significant-digit strings and silently drops the
    # last 1-2 bits on some rows (discovered while verifying issue #7's
    # check_chemform against R: isotopes.csv had 6/308 mass and 81/308
    # abundance values off by 1 ULP, adducts.csv 1/51 Mass value).
    with resources.as_file(_data_path("isotopes.csv")) as path:
        return pd.read_csv(path, float_precision="round_trip")


@lru_cache(maxsize=1)
def _load_chemforms() -> pd.Series:
    with (
        resources.as_file(_data_path("chemforms.txt")) as path,
        open(path, encoding="utf-8") as fh,
    ):
        lines = [line.rstrip("\n") for line in fh]
    return pd.Series(lines, name="chemform")


@lru_cache(maxsize=1)
def _load_adducts() -> pd.DataFrame:
    with resources.as_file(_data_path("adducts.csv")) as path:
        return pd.read_csv(path, float_precision="round_trip")


@lru_cache(maxsize=1)
def _load_resolution_list() -> dict:
    with resources.as_file(_data_path("resolution_list.csv")) as path:
        long_df = pd.read_csv(path, float_precision="round_trip")
    curves: dict = {}
    for name, group in long_df.groupby("name", sort=False):
        curves[name] = group[["mz", "R"]].reset_index(drop=True)
    return curves


_LOADERS = {
    "isotopes": _load_isotopes,
    "chemforms": _load_chemforms,
    "adducts": _load_adducts,
    "resolution_list": _load_resolution_list,
}


def __getattr__(name: str):
    try:
        loader = _LOADERS[name]
    except KeyError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
    return loader()


def __dir__():
    return sorted(set(globals()) | set(_LOADERS))
