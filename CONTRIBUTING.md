# Contributing to pyEnvipat

Thanks for helping. pyEnvipat computes isotope patterns and envelopes for mass spectrometry.
Please read this page before opening a merge request (MR).

## Development setup

Requirements: Python >= 3.10 and a C compiler (the `envipat._core.enviPat_py`
extension is built at install time).

```bash
git clone git@gitlab.eawag.ch:joboog/pyenvipat.git
cd pyenvipat
python -m venv .venv && . .venv/bin/activate
pip install -e ".[dev]"
pre-commit install
```

After changing any `.c`/`.h` file under `src/envipat/_core/`, run
`pip install -e .` again to rebuild the extension.

## Checks and tests

CI (`.gitlab-ci.yml`) runs the same [nox](https://nox.thea.codes/) sessions
you can run locally:

```bash
pytest                    # quick loop on the editable install
nox -s lint format typecheck
nox -s tests              # builds + installs a wheel, runs the suite on Python 3.10-3.13
nox -s build              # sdist + wheel, twine check
pre-commit run --all-files
```

Tool configuration (ruff, mypy, pytest, coverage threshold) lives in
`pyproject.toml`. An MR must pass lint, format, type check and tests, and new
behavior needs new tests.

## Golden-fixture rule

`tests/golden/` holds reference outputs frozen from enviPat v2.8. They are the
regression baseline for the calculations.

- Never edit, reformat or regenerate files under `tests/golden/`. Pre-commit
  hooks are configured to leave them alone; do not work around that.
- A failing golden test means a result changed. Treat it as a regression and
  fix the code unless the change is intended.
- If you believe a fixture is wrong, or a result should deliberately change,
  open an issue first and explain the change in the MR description.
  Tolerances and known divergences are documented in `tests/`.
- The vendored original C sources (`combination`, `profile`, `parse`,
  `element`, `isotope`, `peak`, `data.h`, `preferences.h`) are likewise not
  reformatted or altered algorithmically.

## Commit conventions

Commits follow [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>(<optional scope>)<!>: <imperative summary, capitalized>

<body explaining why, wrapped at ~72 columns>
```

Types in use: `feat`, `fix`, `docs`, `test`, `refactor`, `style`, `build`, `ci`,
`chore`. Mark breaking changes with `!`. Keep commits small and logical, and
name the issue in the MR (`Closes #<n>`).

## Workflow

1. Pick or open an issue (use the templates in `.gitlab/issue_templates`).
2. Branch from `main` as `<issue-number>-<short-slug>`.
3. Commit, push, and open an MR using the MR template; one concern per MR.
4. Increase the version number (`version` in `pyproject.toml`) in every MR,
   following [semantic versioning](https://semver.org/). This is required
   for every new MR, including docs-only ones.
5. Make sure CI is green and respond to review comments.

## Licensing

pyEnvipat is GPL-2.0-only, the same as upstream enviPat. By contributing you
agree that your contribution is released under that license. See `LICENSE`,
`AUTHORS` and `NOTICE`.

## Security

Report vulnerabilities as described in [SECURITY.md](SECURITY.md), not in
public issues.
