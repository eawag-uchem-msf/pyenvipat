# Installation

Install pyEnvipat from [PyPI](https://pypi.org/project/pyEnvipat/):

```bash
pip install pyEnvipat
```

It needs Python >= 3.10. The runtime dependencies (`pandas`, `numpy`, `scipy`)
are installed automatically. R is not needed: the reference data in
{mod}`envipat.data` ships as plain CSV/TXT files.

## Binary wheels

Releases ship binary wheels, which need no C compiler, for CPython 3.10-3.13 on:

- Linux x86_64 (manylinux_2_28 and musllinux_1_2)
- Windows x86_64
- macOS x86_64 (Intel) and arm64 (Apple Silicon)

Every wheel passes the full test suite in CI before it is published.

## Building from source

On other platforms or Python versions, pip builds the package from the sdist.
That needs:

- a C compiler, to build the bundled `envipat._core.enviPat_py` extension; on
  Windows that is `clang-cl` (Visual Studio component "C++ Clang tools for
  Windows", or LLVM), because the vendored C sources use variable-length
  arrays, which `cl.exe` cannot compile
- `numpy` at build time (installed automatically by pip's build isolation)

The same applies to installing from a source checkout:

```bash
git clone https://github.com/eawag-uchem-msf/pyenvipat.git
cd pyenvipat
pip install .
```

## Development install

From a source checkout:

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
