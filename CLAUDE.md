# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

**pyEnvipat** (import name `envipat`) is a Python port of the R package enviPat (Loos, Gerber et al., EAWAG; GPL-2.0-only, see `LICENSE`, `AUTHORS`, `NOTICE`) for isotope pattern and mass spectrometry envelope calculations. The original algorithms are C code, compiled into a CPython extension; the Python layer handles formula parsing, pandas I/O and orchestration. No R code is part of this repository; the R package survives only as frozen reference outputs in `tests/golden/`.

## Common Commands

```bash
# Create an environment and install in editable mode with dev tooling (compiles the C extension)
python -m venv .venv && . .venv/bin/activate
pip install -e ".[dev]"
pre-commit install

# Quick loop on the editable install (no R needed)
pytest
pytest tests/python/test_pipeline.py -k isowrap

# What CI runs (.github/workflows/ calls these nox sessions)
nox -s lint format typecheck    # ruff check, ruff format (--check in CI), mypy
nox -s tests                    # builds + installs a wheel, pytest + coverage on 3.10-3.13
nox -s cwarnings                # own C files (iso_core.c, enviPat_py.c) with -Wall -Wextra -Werror
nox -s sanitize                 # pytest against an ASan/UBSan build of the extension (GCC, Linux)
nox -s build                    # sdist + wheel, twine check, sdist must ship all .c/.h and data files
nox -s version_bump             # __version__ must be higher than on origin/main (CI: MR pipelines only, not renovate/* branches)
nox -s audit                    # pip-audit the resolved runtime dependencies (not dev tooling)
nox -s release_check -- vX.Y.Z  # tag == v + __version__ == CITATION.cff version (CI: v* tag pipelines)
nox -s renovate_check           # validate renovate.json (needs npx); pinning policy in CONTRIBUTING.md
nox -s docs                     # Sphinx site (-W, nitpicky), executes docs/tutorial.ipynb; not yet published by CI

# Wheels as CI builds them (config in [tool.cibuildwheel]; Linux needs docker).
# CI (wheels.yml) builds Linux, Windows (win_amd64) and macOS (x86_64, arm64) wheels;
# the `tests` job runs on Linux, Windows and macOS.
cibuildwheel --platform linux
```

Tool configuration (ruff, mypy, pytest, coverage threshold) lives in `pyproject.toml`; `.pre-commit-config.yaml` mirrors the quality checks and must never touch `tests/golden/` or the vendored C files. On a host filesystem without symlink support (e.g. a Windows mount), pass `nox --envdir <dir outside the mount>`.

After changing any `.c`/`.h` file under `src/envipat/_core/`, re-run `pip install -e .` to rebuild the extension.

## Architecture

### Layout

- `src/envipat/` — Python package
  - `formula.py` — `check_chemform`, `check_ded`, `mergeform`, `multiform`, `subform`
  - `isopattern.py`, `envelope.py`, `vdetect.py`, `check_several.py`, `isowrap.py` — pipeline functions
  - `resolution.py` — `get_resolution` (port of R's `getR`, fitted with scipy's `make_smoothing_spline`; not bit-compatible with R); `_rseq.py` holds R's `seq` semantics
  - `data.py` + `data/` — lazily loaded `isotopes`, `chemforms`, `adducts`, `resolution_list` (CSV/TXT)
- `src/envipat/_core/` — C code: `iso_core.c/h` (pure-C core), `enviPat_py.c` (CPython bindings), and the vendored original algorithm files (`combination`, `profile`, `parse`, `element`, `isotope`, `peak`, `data.h`, `preferences.h`). `setup.py` builds them into `envipat._core.enviPat_py`, with flags per compiler (`-O2 -ffp-contract=off -Wall -Wextra` for GCC/Clang, `/O2` for MSVC); the vendored files' remaining warnings are triaged there.
- `docs/` — Sphinx site (furo, numpydoc, myst-nb): `api/` pages autodoc one module each; `data.rst` documents the `envipat.data` datasets by hand (they are lazy module attributes autodoc cannot see); `tutorial.ipynb` is committed without outputs and executed on every build. Cross-references must resolve (`nitpicky`), so private names go in double backticks, not `:func:`.

### Typical workflow

```
envipat.data.isotopes / chemforms / adducts
      ↓
check_chemform()   ← validates & normalizes formula strings
      ↓
isopattern()       ← isotope fine structure (algo=1 or algo=2)
      ↓
envelope()         ← profile convolution (Gaussian / Cauchy-Lorentz)
      ↓
vdetect()          ← centroid / intensoid / valley detection
```

`isowrap()` runs the pipeline end to end; `get_resolution()` interpolates instrument resolution from a `resolution_list` curve.

### Tests

- `tests/test_core.py` — C extension against `tests/golden/c_direct/` (exact equality).
- `tests/python/` — data loaders, formula helpers, pipeline functions, `get_resolution`.
- `tests/golden/` — fixtures frozen from the original R package v2.8. **Never regenerate or edit them**: the R scripts that generated them are not part of this repository (only in git history). Tolerances per module are documented in each test file's docstring.

## Conventions

- Keep the vendored C algorithm files unmodified so floating-point behaviour matches the original; changes belong in `iso_core.c`/`enviPat_py.c`.
- Commit messages follow Conventional Commits.
- Every MR bumps `__version__` in `src/envipat/__init__.py`, the single source of the version (pyproject reads it via `[tool.setuptools.dynamic]`, the docs via `envipat.__version__`). CI job `version-bump` (`nox -s version_bump`) enforces it on MR pipelines. `CITATION.cff` repeats the version (a test checks they match). User-facing changes go into `CHANGELOG.md` under `[Unreleased]`.
- Releases: a protected `vX.Y.Z` tag runs `release-check`, builds the sdist and cibuildwheel wheels (GitHub Release assets), and publishes to PyPI via trusted publishing, gated by the `pypi` environment's required reviewer.
