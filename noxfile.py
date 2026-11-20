"""Task runner. CI (.gitlab-ci.yml) calls these sessions, so local and CI runs match.

    nox                       # lint, format, typecheck, tests on every Python
    nox -s tests-3.12         # tests on one Python version
    nox -s format -- --check  # what CI runs (check only, no rewrite)
    nox -s build              # sdist + wheel into dist/, checks sdist contents
    nox -s cwarnings          # our own C files must compile warning-free (-Werror)
    nox -s sanitize           # test suite against an ASan/UBSan build (GCC, Linux)
    nox -s docs               # Sphinx site into docs/_build/html, warnings are errors
    nox -s version_bump       # version must be higher than on origin/main (CI: every MR)
    nox -s renovate_check     # renovate.json against CI's Renovate version (needs npx)

`tests` builds a wheel and installs it non-editably, so the suite runs against
the installed package (compiled extension, package data), not the source tree.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tarfile
from pathlib import Path

import nox
from packaging.version import Version  # a dependency of nox itself

PYTHONS = ["3.10", "3.11", "3.12", "3.13"]
LINT_PATHS = ["src", "tests", "docs", "noxfile.py", "setup.py"]
ROOT = Path(__file__).parent
CORE_DIR = ROOT / "src" / "envipat" / "_core"
DOCS_DIR = ROOT / "docs"
VERSION_FILE = "src/envipat/__init__.py"

# C files we maintain. The vendored algorithm files stay unmodified, so their
# triaged warnings (see setup.py) are not held to -Werror.
OWN_C_FILES = ["iso_core.c", "enviPat_py.c"]

# GCC: -fno-sanitize-recover makes any UBSan finding abort instead of logging and
# continuing, so the test run fails on it. setup.py still appends -O2.
SANITIZE_FLAGS = "-fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer -g"

# The test import must resolve to site-packages, never to ./src.
ASSERT_INSTALLED = """
import pathlib, sys, envipat
p = pathlib.Path(envipat.__file__).resolve()
print("envipat imported from", p)
if "site-packages" not in p.parts:
    sys.exit("envipat was imported from the source tree, not the installed wheel")
"""

nox.needs_version = ">=2025.2.9"
nox.options.default_venv_backend = "uv|virtualenv"
nox.options.sessions = ["lint", "format", "typecheck", "tests"]
nox.options.error_on_missing_interpreters = bool(os.environ.get("CI"))


@nox.session
def lint(session: nox.Session) -> None:
    """ruff check."""
    session.install("ruff==0.17.5")
    session.run("ruff", "check", *LINT_PATHS, *session.posargs)


@nox.session(name="format")
def format_(session: nox.Session) -> None:
    """ruff format (pass `-- --check` or set CI to only check)."""
    session.install("ruff==0.17.5")
    check = ["--check", "--diff"] if "--check" in session.posargs or os.environ.get("CI") else []
    session.run("ruff", "format", *check, *LINT_PATHS)


@nox.session
def typecheck(session: nox.Session) -> None:
    """mypy (config in pyproject.toml)."""
    session.install("-e", ".[dev]")
    session.run("mypy", *session.posargs)


def _install_fresh_wheel(session: nox.Session, env: dict[str, str] | None = None) -> None:
    """Build a wheel from a clean build/ dir and install it with the test extra.

    setuptools skips compiling when build/ already holds an up-to-date extension,
    which would mix e.g. a sanitizer build into a normal one; so start clean.
    """
    wheelhouse = Path(session.create_tmp()) / "wheelhouse"
    shutil.rmtree(wheelhouse, ignore_errors=True)
    shutil.rmtree(ROOT / "build", ignore_errors=True)
    session.install("build")
    session.run("python", "-m", "build", "--wheel", "--outdir", str(wheelhouse), str(ROOT), env=env)
    (wheel,) = wheelhouse.glob("*.whl")
    session.install(f"{wheel}[test]")
    # Same version as a previous run's wheel would count as already installed.
    session.install("--force-reinstall", "--no-deps", str(wheel))


@nox.session(python=PYTHONS)
def tests(session: nox.Session) -> None:
    """pytest + coverage against the installed wheel."""
    _install_fresh_wheel(session)
    session.run("python", "-c", ASSERT_INSTALLED)

    session.run("python", "-c", ASSERT_INSTALLED)

    session.run(
        "pytest",
        "--cov=envipat",
        "--cov-report=term-missing",
        "--cov-report=xml:coverage.xml",
        "--junitxml=report.xml",
        *session.posargs,
    )


@nox.session
def cwarnings(session: nox.Session) -> None:
    """Compile our own C files with -Wall -Wextra -Werror (GCC/Clang)."""
    session.install("numpy")
    cc = os.environ.get("CC", "cc")
    includes = session.run(
        "python",
        "-c",
        "import numpy, sysconfig; print(sysconfig.get_paths()['include']); "
        "print(numpy.get_include())",
        silent=True,
    )
    assert isinstance(includes, str)
    inc_flags = [f"-I{d}" for d in [str(CORE_DIR), *includes.split()]]
    for name in OWN_C_FILES:
        session.run(
            cc, "-c", "-O2", "-Wall", "-Wextra", "-Werror", *inc_flags,
            str(CORE_DIR / name), "-o", os.devnull, external=True,
        )  # fmt: skip


@nox.session
def sanitize(session: nox.Session) -> None:
    """Run the test suite against an ASan/UBSan build of the extension (GCC, Linux).

    The interpreter itself is not instrumented, so the ASan runtime is preloaded.
    Leak detection is off: CPython and NumPy keep allocations alive at exit.
    """
    cc = os.environ.get("CC", "gcc")
    libasan = session.run(cc, "-print-file-name=libasan.so", silent=True, external=True)
    assert isinstance(libasan, str)
    build_env = {"CC": cc, "CFLAGS": SANITIZE_FLAGS, "LDFLAGS": SANITIZE_FLAGS}
    try:
        _install_fresh_wheel(session, env=build_env)
        run_env = {
            "LD_PRELOAD": libasan.strip(),
            "ASAN_OPTIONS": "detect_leaks=0:abort_on_error=1",
            "UBSAN_OPTIONS": "print_stacktrace=1:halt_on_error=1",
        }
        session.run("python", "-c", ASSERT_INSTALLED, env=run_env)
        session.run("pytest", "-p", "no:cacheprovider", *session.posargs, env=run_env)
    finally:
        # Never leave sanitizer objects in build/ for a later normal build.
        shutil.rmtree(ROOT / "build", ignore_errors=True)


@nox.session
def docs(session: nox.Session) -> None:
    """Build the Sphinx site into docs/_build/html; any warning fails the build.

    The tutorial notebook is executed during the build, so it also checks that the
    examples still run. Pass e.g. `-- -b linkcheck` to add sphinx-build options.
    """
    session.install("-e", ".[docs]")
    out = DOCS_DIR / "_build" / "html"
    shutil.rmtree(out, ignore_errors=True)
    session.run("sphinx-build", "-W", "--keep-going", *session.posargs, str(DOCS_DIR), str(out))


def _read_version(source: str) -> str:
    match = re.search(r'^__version__ = "([^"]+)"$', source, re.MULTILINE)
    if match is None:
        raise RuntimeError(f"no __version__ assignment found in {VERSION_FILE}")
    return match.group(1)


@nox.session(python=False)
def version_bump(session: nox.Session) -> None:
    """Fail unless envipat.__version__ is higher than on the base branch.

    The base ref defaults to origin/main; CI passes the MR's target branch:
    `nox -s version_bump -- origin/<target>`.
    """
    base_ref = session.posargs[0] if session.posargs else "origin/main"
    base_source = subprocess.run(
        ["git", "show", f"{base_ref}:{VERSION_FILE}"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    base = Version(_read_version(base_source))
    head = Version(_read_version((ROOT / VERSION_FILE).read_text()))
    if head <= base:
        session.error(
            f"version {head} is not higher than {base} on {base_ref}; "
            f"bump __version__ in {VERSION_FILE}"
        )
    print(f"version bumped: {base} ({base_ref}) -> {head}")


@nox.session(python=False)
def renovate_check(session: nox.Session) -> None:
    """Validate renovate.json with the Renovate version CI runs (needs Node.js/npx)."""
    ci = (ROOT / ".gitlab-ci.yml").read_text()
    match = re.search(r"ghcr\.io/renovatebot/renovate:(\S+)", ci)
    if match is None:
        session.error("no ghcr.io/renovatebot/renovate image in .gitlab-ci.yml")
    session.run(
        "npx",
        "--yes",
        "--package",
        f"renovate@{match.group(1)}",
        "--",
        "renovate-config-validator",
        "--strict",
        external=True,
    )


def _check_sdist(sdist: Path) -> None:
    """Fail unless the sdist ships every C source/header and data file."""
    required = [
        *CORE_DIR.glob("*.c"),
        *CORE_DIR.glob("*.h"),
        *(ROOT / "src" / "envipat" / "data").glob("*"),
    ]
    with tarfile.open(sdist) as tar:
        # Members are "<name>-<version>/<path>"; strip the top-level directory.
        members = {m.split("/", 1)[1] for m in tar.getnames() if "/" in m}
    missing = sorted(
        p.relative_to(ROOT).as_posix()
        for p in required
        if p.relative_to(ROOT).as_posix() not in members
    )
    if missing:
        raise RuntimeError(f"sdist {sdist.name} is missing: {', '.join(missing)}")
    print(f"sdist {sdist.name}: all {len(required)} C/header/data files present")


@nox.session
def build(session: nox.Session) -> None:
    """Build sdist + wheel into dist/, validate metadata and sdist contents."""
    session.install("build", "twine")
    for f in (ROOT / "dist").glob("*"):
        f.unlink()
    shutil.rmtree(ROOT / "build", ignore_errors=True)
    session.run("python", "-m", "build", "--outdir", "dist", str(ROOT))
    session.run("twine", "check", "--strict", *map(str, (ROOT / "dist").glob("*")))
    (sdist,) = (ROOT / "dist").glob("*.tar.gz")
    _check_sdist(sdist)
