"""Task runner. CI (.gitlab-ci.yml) calls these sessions, so local and CI runs match.

    nox                       # lint, format, typecheck, tests on every Python
    nox -s tests-3.12         # tests on one Python version
    nox -s format -- --check  # what CI runs (check only, no rewrite)
    nox -s build              # sdist + wheel into dist/

`tests` builds a wheel and installs it non-editably, so the suite runs against
the installed package (compiled extension, package data), not the source tree.
"""

from __future__ import annotations

import os
from pathlib import Path

import nox

PYTHONS = ["3.10", "3.11", "3.12", "3.13"]
LINT_PATHS = ["src", "tests", "noxfile.py", "setup.py"]
ROOT = Path(__file__).parent

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


@nox.session(python=PYTHONS)
def tests(session: nox.Session) -> None:
    """pytest + coverage against the installed wheel."""
    wheelhouse = Path(session.create_tmp()) / "wheelhouse"
    session.install("build")
    session.run("python", "-m", "build", "--wheel", "--outdir", str(wheelhouse), str(ROOT))
    wheels = sorted(wheelhouse.glob("*.whl"), key=lambda p: p.stat().st_mtime)
    session.install(f"{wheels[-1]}[test]")

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
def build(session: nox.Session) -> None:
    """Build sdist + wheel into dist/ and validate their metadata."""
    session.install("build", "twine")
    for f in (ROOT / "dist").glob("*"):
        f.unlink()
    session.run("python", "-m", "build", "--outdir", "dist", str(ROOT))
    session.run("twine", "check", "--strict", *map(str, (ROOT / "dist").glob("*")))
