# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

**pyEnvipat** (import name `envipat`) is a Python port of the R package enviPat (Loos, Gerber et al., EAWAG; GPL-2.0-only, see `LICENSE`, `AUTHORS`, `NOTICE`) for isotope pattern and mass spectrometry envelope calculations. The original algorithms are C code, compiled into a CPython extension; the Python layer handles formula parsing, pandas I/O and orchestration. The R implementation has been removed; it lives on only in git history and as frozen reference outputs in `tests/golden/`.

## Common Commands

```bash
# Create an environment and install in editable mode (compiles the C extension)
python -m venv .venv && . .venv/bin/activate
pip install -e . pytest

# Run the whole test suite (no R needed)
pytest tests

# Run one module / test
pytest tests/python/test_pipeline.py -k isowrap
```

After changing any `.c`/`.h` file under `src/envipat/_core/`, re-run `pip install -e .` to rebuild the extension.

## Architecture

### Layout

- `src/envipat/` — Python package
  - `formula.py` — `check_chemform`, `check_ded`, `mergeform`, `multiform`, `subform`
  - `isopattern.py`, `envelope.py`, `vdetect.py`, `check_several.py`, `isowrap.py` — pipeline functions
  - `resolution.py` — `get_resolution` (port of R's `getR`, fitted with scipy's `make_smoothing_spline`; not bit-compatible with R); `_rseq.py` holds R's `seq` semantics
  - `data.py` + `data/` — lazily loaded `isotopes`, `chemforms`, `adducts`, `resolution_list` (CSV/TXT)
- `src/envipat/_core/` — C code: `iso_core.c/h` (pure-C core), `enviPat_py.c` (CPython bindings), and the vendored original algorithm files (`combination`, `profile`, `parse`, `element`, `isotope`, `peak`, `data.h`, `preferences.h`). `setup.py` builds them into `envipat._core.enviPat_py` with `-O2`.

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
- `tests/golden/` — fixtures frozen from the original R package v2.8. **Never regenerate or edit them**: the generating R scripts are no longer in the tree (`git log -- tools/make_golden.R`). Tolerances per module are documented in each test file's docstring.

## Conventions

- Keep the vendored C algorithm files unmodified so floating-point behaviour matches the original; changes belong in `iso_core.c`/`enviPat_py.c`.
- Commit messages follow Conventional Commits.
