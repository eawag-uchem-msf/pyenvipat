# Changelog

All notable changes to pyEnvipat are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/). Releases are tagged `vX.Y.Z`; see
"Releasing" in [CONTRIBUTING.md](CONTRIBUTING.md) for how this file is updated.

## [Unreleased]

### Fixed

- Renovate tracks the Renovate version that validates `renovate.json` in CI
  (`renovate-config.yml`), which the move to GitHub Actions had left
  unmanaged; the custom manager for `pip install x==y` pins in workflows,
  which no longer matched anything, is removed.
- CI no longer fails on `main` for lack of a macOS runner: the `tests` job
  runs macOS only on Python 3.10 and 3.13 (the `wheels` workflow still tests
  every Python version on macOS), and runs on `main` and on tags are no longer
  cancelled by a later push.

## [1.0.5] - 2026-11-24

### Added

- Release process: a pipeline for a protected `vX.Y.Z` tag checks that the tag
  matches `envipat.__version__` and `CITATION.cff`, and keeps the sdist and the
  Linux wheels as release artifacts (no upload to PyPI yet).
- `pip-audit` CI job and `nox -s audit` for known vulnerabilities in the runtime
  dependencies.
- PyPI metadata: classifiers, keywords, license files and well-known project
  URLs. The `Private :: Do Not Upload` classifier blocks uploads until
  publishing is set up.
- This changelog.
- Renovate opens weekly MRs for dev tools, pre-commit hooks and CI images;
  the pinning policy is documented in `CONTRIBUTING.md`.
- GitHub Actions CI (`.github/workflows/`): lint, type check, tests, C warning
  check, sanitizer run, build and Linux wheels on every push/PR, plus a
  `release` workflow on `v*` tags that publishes a GitHub Release and uploads
  to PyPI via trusted publishing (gated by the `pypi` environment's required
  reviewer). Run manually, the `release` workflow is a dry run that uploads
  the same sdist and wheels to TestPyPI instead.
- `.github/ISSUE_TEMPLATE/` and `.github/PULL_REQUEST_TEMPLATE.md`.
- Documentation site on GitHub Pages
  (<https://eawag-uchem-msf.github.io/pyenvipat/>): the `docs` workflow builds
  it on every PR and deploys it from `main`; `Documentation` added to the
  project URLs.
- Windows and macOS support (#22, #23): the test suite runs on Windows (MSVC)
  and macOS (Intel and Apple Silicon) for Python 3.10-3.13, and binary wheels
  are built for `win_amd64`, macOS `x86_64` and macOS `arm64`; releases publish
  them next to the Linux wheels. Building from source on Windows uses
  `clang-cl` (the vendored C sources use variable-length arrays, which MSVC's
  `cl.exe` does not support); on Windows and arm64, golden-fixture abundances
  are compared with a 4-ULP relative tolerance instead of exactly.

### Changed

- The build requires `setuptools>=77.0.3` (PEP 639 license metadata) and builds
  against `numpy>=2.0`; the built extension still runs on numpy 1.24 and newer.
- Repository moved from `gitlab.eawag.ch` to `github.com/eawag-uchem-msf/pyenvipat`:
  CI, issue/PR templates and all documentation references now point at GitHub;
  Renovate runs via the GitHub App instead of a self-hosted CI job.

### Removed

- The `Private :: Do Not Upload` classifier, now that PyPI publishing is wired
  up (uploads still only happen from a `v*` tag, gated by a required reviewer).

### Fixed

- The C extension is compiled with `-ffp-contract=off` on GCC/Clang/clang-cl, so
  compilers that fuse multiply-adds by default (Apple clang on arm64, GCC on
  aarch64) produce the same results as the x86_64 reference builds.
- Release uploads to PyPI now contain the sdist and all cibuildwheel wheels
  instead of the sdist plus an unpublishable plain `linux_x86_64` wheel.
- A failing test run in CI no longer passes because its output is piped
  through `tee`.
- `isopattern` warns about a non-default `emass` with an explicit `UserWarning`
  that names the given and the default value (it read "electrone" before).

## [0.2.2] - 2026-11-17

### Added

- Contributing guide, security policy, `CITATION.cff`, and GitLab issue and
  merge request templates.

## [0.2.1] - 2026-11-16

### Added

- Type hints for the public API, `py.typed` and a stub for the C extension.
- `[81]Br` labelled isotope in the isotope list.
- Sphinx documentation site with API reference and tutorial, published to
  GitLab Pages.
- CI on Linux: lint, type check, tests on Python 3.10-3.13, C warning check,
  AddressSanitizer/UBSan run, sdist and manylinux/musllinux wheels.

### Changed

- The version is defined only in `envipat.__version__`; every merge request
  must raise it.

### Fixed

- Elements are sorted stably, so m/z no longer depends on the C library's
  `qsort`.
- The sdist ships the C headers, so it can be built into a wheel.
- CSV loaders no longer pass pandas' deprecated `float_precision` argument.

## [0.2.0] - 2026-11-02

### Removed

- The R implementation; the Python port is complete.
- **Breaking:** the `r_exact` method of `get_resolution`.

## [0.1.0] - 2026-09-18

### Added

- First Python port of enviPat: C extension (`envipat._core.enviPat_py`),
  bundled data, formula helpers, `isopattern`, `envelope`, `vdetect`,
  `check_several`, `isowrap` and `get_resolution`, tested against golden
  fixtures frozen from R enviPat v2.8.

[Unreleased]: https://github.com/eawag-uchem-msf/pyenvipat/compare/v1.0.5...HEAD
[1.0.5]: https://github.com/eawag-uchem-msf/pyenvipat/compare/a485e1f...v1.0.5
[0.2.2]: https://gitlab.eawag.ch/uchem/msf-server/pyenvipat/-/commit/a485e1f
[0.2.1]: https://gitlab.eawag.ch/uchem/msf-server/pyenvipat/-/commit/977d999
[0.2.0]: https://gitlab.eawag.ch/uchem/msf-server/pyenvipat/-/tags/0.2.0
