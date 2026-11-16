"""Build script for the envipat._core.enviPat_py C extension.

pyproject.toml carries build-system requirements and project metadata, but
setuptools has no declarative pyproject.toml syntax for compiled Extension
objects -- those still need to live here.
"""

import numpy
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

CORE_DIR = "src/envipat/_core"

# Compiler flags per compiler family (setuptools' compiler_type), chosen when
# the compiler is known, i.e. at build time rather than import time.
#
# -O2 / /O2 matches the optimization level the original enviPat C sources were
# built with, so floating-point behavior does not drift from the original
# results the tests compare against. GCC, Clang and MSVC all honor the *last*
# optimization flag on the command line, and setuptools appends
# extra_compile_args after its own defaults (e.g. -O3 from sysconfig, /Ox on
# MSVC), so appending it here is sufficient to win.
#
# -Wall -Wextra: the vendored algorithm files (combination.c, profile.c,
# parse.c) emit warnings that were triaged and left as is, because those files
# must stay unmodified (see CLAUDE.md):
#   -Wsign-compare           ptrdiff_t/int vs. size_t loop bounds; all values
#                            are non-negative counts far below 2**63.
#   -Wunused-parameter       parameters kept for the original call signatures.
#   -Wunterminated-string-initialization
#                            parse.c: char test[10] = "123456789\0" -- only the
#                            implicit second NUL is dropped; still terminated.
# Our own C files (iso_core.c, enviPat_py.c) must stay warning-free; the
# `cwarnings` nox session / CI job enforces that with -Werror.
COMPILE_ARGS = {
    "msvc": ["/O2"],
    "unix": ["-O2", "-Wall", "-Wextra"],  # GCC and Clang, incl. macOS and MinGW
    "mingw32": ["-O2", "-Wall", "-Wextra"],
}


class BuildExt(build_ext):
    """Add the compiler-specific flags from COMPILE_ARGS to every extension."""

    def build_extensions(self):
        args = COMPILE_ARGS.get(self.compiler.compiler_type, [])
        for ext in self.extensions:
            ext.extra_compile_args = [*args, *ext.extra_compile_args]
        super().build_extensions()


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
)

setup(ext_modules=[enviPat_py], cmdclass={"build_ext": BuildExt})
