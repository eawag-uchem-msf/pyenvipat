# Installation

pyEnvipat is not on PyPI yet. Install it from a source checkout:

```bash
git clone https://github.com/eawag-uchem-msf/pyenvipat.git
cd pyenvipat
pip install .
```

Requirements:

- Python >= 3.10
- a C compiler, to build the bundled `envipat._core.enviPat_py` extension at
  install time; on Windows that is `clang-cl` (Visual Studio component "C++
  Clang tools for Windows", or LLVM), because the vendored C sources use
  variable-length arrays, which `cl.exe` cannot compile
- `numpy` at build time (installed automatically by pip's build isolation)

CI builds and tests the package on Linux, Windows (MSVC) and macOS (Intel and
Apple Silicon), and releases ship binary wheels for these platforms, which need
no C compiler.

The runtime dependencies (`pandas`, `numpy`, `scipy`) are installed
automatically. R is not needed: the reference data in {mod}`envipat.data` ships
as plain CSV/TXT files.

## Development install

```bash
pip install -e ".[dev]"   # editable install with test, lint and type-check tooling
pre-commit install
```

To build this documentation locally:

```bash
pip install -e ".[docs]"
sphinx-build -W --keep-going docs docs/_build/html
# or, as CI does:
nox -s docs
```
