# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

**enviPat** is an R package (v2.8) for isotope pattern and mass spectrometry envelope calculations. It uses a hybrid R/C architecture: user-facing R functions call performance-critical C algorithms via `.Call()`.

## Common Commands

```bash
# Install the package locally (compiles C code)
R CMD INSTALL .

# Build a distributable tarball
R CMD build .

# Run full CRAN package checks
R CMD check enviPat_*.tar.gz
# or interactively:
Rscript -e 'rcmdcheck::rcmdcheck()'

# Regenerate man/ documentation from roxygen2 comments
Rscript -e 'roxygen2::roxygenise()'

# Run a single function's examples
Rscript -e 'library(enviPat); example(isopattern)'
```

There is no automated test suite. All function documentation includes runnable `@examples` blocks that serve as integration tests.

## Architecture

### R/C Split

The package has two computational layers:

- **R layer** (`R/`): User-facing functions handle argument validation, looping over batches, and result formatting.
- **C layer** (`src/`): Core algorithms called via `.Call()`. The five registered C entry points are declared in [src/enviPat_init.c](src/enviPat_init.c):
  - `iso_pattern` / `iso_pattern_4` — two isotope fine structure algorithms (algo=1 and algo=2)
  - `iso_ppm_Call` — PPM-based m/z discretization
  - `iso_profile_with_trace_Call` — peak profile convolution (Gaussian / Cauchy-Lorentz)
  - `iso_centroid_Call` — centroid / intensoid / valley detection

The bulk of the C logic lives in [src/main.c](src/main.c) (~2700 lines) and [src/combination.c](src/combination.c) (~3500 lines). Static isotope and element data are embedded in [src/data.h](src/data.h) (~4600 lines).

### Typical User Workflow

```
data(isotopes) / data(chemforms) / data(adducts)
      ↓
check_chemform()   ← validates & normalizes formula strings
      ↓
isopattern()       ← calls iso_pattern or iso_pattern_4 in C
      ↓
envelope()         ← calls iso_profile_with_trace_Call in C
      ↓
vdetect()          ← calls iso_centroid_Call in C
```

`isowrap()` is a convenience wrapper that runs this full pipeline including resolution interpolation via `getR()`.

### Bundled Data (`data/`)

| Object | Description |
|---|---|
| `isotopes` | Reference table of stable isotopes (mass, abundance) |
| `chemforms` | Example chemical formulas for demos |
| `adducts` | 47 ESI ionization adduct definitions |
| `resolution_list` | Empirical resolution curves for Orbitrap/QExactive instruments |

These `.rda` files are documented in `R/data-*.R` and must be loaded with `data(...)` before use.

### Documentation

All exported functions use roxygen2 (`@param`, `@return`, `@examples`, etc.). Edit the roxygen block in the corresponding `R/*.R` file and regenerate with `roxygen2::roxygenise()` — do **not** edit `man/*.Rd` files directly.

## Development Container

A `.devcontainer/` setup is included (rocker/verse base, RStudio Server on port 8787). The container includes devtools, roxygen2, testthat, and a C/C++ compiler. VS Code extensions for R and C++ are preconfigured.

**Always build, install, and run tests inside the dev container, not on the host.** The host's glibc version may not match what the container's toolchain produces, causing `.so` load failures (`GLIBC_x.xx not found`) if compiled objects cross between the two. If the container isn't already running, start it (e.g. `docker ps` to check, or use the VS Code "Reopen in Container" flow) before running `R CMD INSTALL`, `R CMD build`, `R CMD check`, or any test suite. Use `docker exec <container> bash -lc '...'` to run commands inside it from the CLI. Run `R CMD INSTALL --install-tests .` so `testthat.R` can find `tests/testthat` in the installed package, and invoke it as `cd tests && Rscript testthat.R` (mirrors how `R CMD check` sets the working directory).
