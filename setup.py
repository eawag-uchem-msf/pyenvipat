"""Build script for the envipat._core.enviPat_py C extension.

pyproject.toml carries build-system requirements and project metadata, but
setuptools has no declarative pyproject.toml syntax for compiled Extension
objects -- those still need to live here.
"""

import os
import shutil

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
# -ffp-contract=off: never fuse a*b+c into one FMA instruction. Clang (Apple
# clang on arm64 in particular) and GCC on aarch64 contract by default, which
# rounds once instead of twice and breaks the exact-equality golden fixtures
# (tests/test_core.py), frozen from x86_64 builds without FMA. On x86_64 without
# -mfma there is no FMA to emit, so Linux/Windows x86_64 builds are unchanged.
# On Windows clang-cl gets the same option via /clang:-ffp-contract=off.
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
    # clang-cl, not cl.exe (see BuildExt); /clang: passes a plain clang option.
    "msvc": ["/O2", "/clang:-ffp-contract=off"],
    "unix": ["-O2", "-ffp-contract=off", "-Wall", "-Wextra"],  # GCC and Clang, incl. macOS
    "mingw32": ["-O2", "-ffp-contract=off", "-Wall", "-Wextra"],
}


def _use_clang_cl(compiler) -> None:
    """Compile with clang-cl instead of cl.exe in the MSVC toolchain.

    The vendored combination.c uses C99 variable-length arrays, which cl.exe does
    not support, and the vendored files stay unmodified (see CLAUDE.md). clang-cl
    accepts them and is a drop-in for cl.exe: same options, same MSVC runtime and
    link.exe, so the extension is ABI-compatible with the python.org CPython.
    It ships with Visual Studio ("C++ Clang tools for Windows") and LLVM; set
    CLANG_CL to its path if it is not on PATH.
    """
    clang_cl = os.environ.get("CLANG_CL") or shutil.which("clang-cl")
    if not clang_cl:
        raise RuntimeError(
            "Building pyEnvipat on Windows needs clang-cl (the vendored C sources use "
            "variable-length arrays, which cl.exe cannot compile). Install the Visual "
            "Studio component 'C++ Clang tools for Windows' or LLVM, and put clang-cl "
            "on PATH or set CLANG_CL to its full path."
        )
    if not compiler.initialized:
        compiler.initialize()
    compiler.cc = clang_cl
    # /GL emits cl.exe-specific whole-program objects; clang-cl's are plain COFF.
    for attr in ("compile_options", "compile_options_debug"):
        opts = getattr(compiler, attr, None)
        if opts is not None:
            setattr(compiler, attr, [o for o in opts if o != "/GL"])


class BuildExt(build_ext):
    """Add the compiler-specific flags from COMPILE_ARGS to every extension."""

    def build_extensions(self):
        if self.compiler.compiler_type == "msvc":
            _use_clang_cl(self.compiler)
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
