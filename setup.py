"""Build script for the envipat._core.enviPat_py C extension.

pyproject.toml carries build-system requirements and project metadata, but
setuptools has no declarative pyproject.toml syntax for compiled Extension
objects -- those still need to live here.
"""

import numpy
from setuptools import Extension, setup

CORE_DIR = "src/envipat/_core"

# Matches the optimization level the original enviPat C sources were built
# with (-O2), so floating-point behavior does not drift from the original
# results the tests compare against (see iso_core.c/combination.c/etc.,
# vendored here unmodified). GCC and Clang both take the *last* -O flag
# seen on the command line, and setuptools appends extra_compile_args after
# its own sysconfig-derived flags (commonly -O3 on this platform's Python
# build), so appending "-O2" here is sufficient to win.
EXTRA_COMPILE_ARGS = ["-O2"]

enviPat_py = Extension(
    "envipat._core.enviPat_py",
    sources=[
        f"{CORE_DIR}/enviPat_py.c",
        f"{CORE_DIR}/iso_core.c",
        f"{CORE_DIR}/combination.c",
        f"{CORE_DIR}/profile.c",
        f"{CORE_DIR}/parse.c",
        f"{CORE_DIR}/element.c",
        f"{CORE_DIR}/isotope.c",
        f"{CORE_DIR}/peak.c",
    ],
    include_dirs=[CORE_DIR, numpy.get_include()],
    extra_compile_args=EXTRA_COMPILE_ARGS,
)

setup(ext_modules=[enviPat_py])
