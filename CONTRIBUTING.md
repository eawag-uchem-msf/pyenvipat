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
4. Increase the version number in every MR, following
   [semantic versioning](https://semver.org/), including docs-only ones.
   Change `__version__` in `src/envipat/__init__.py` (`pyproject.toml` reads
   it from there) and `version` in `CITATION.cff`. The CI job `version-bump`
   (`nox -s version_bump -- <target-ref>`) fails if the version was not raised.
   Renovate MRs (branches `renovate/*`) are exempt: they only change tooling
   and CI images, never the package.
5. Make sure CI is green and respond to review comments.

## Dependencies and pinning

pyEnvipat is a library, so its requirements must not block the environments it
is installed into. Which dependencies are pinned, and how:

| What | Where | Policy | Updated by |
|------|-------|--------|------------|
| Runtime dependencies | `[project.dependencies]` | lower bound only, no upper cap | hand, only when the code needs a newer version |
| Build requirements | `[build-system].requires` | lower bound only | hand |
| ruff, mypy | `dev` extra **and** `.pre-commit-config.yaml` | exact pin, same version in both files | Renovate, one MR per tool |
| Other test/dev/docs extras | `[project.optional-dependencies]` | lower bound or unpinned | hand |
| pre-commit hooks | `.pre-commit-config.yaml` (`rev:`) | exact tag | Renovate |
| CI images | `.gitlab-ci.yml` (`image:`) | major or minor tag (`python:3.12`, `alpine:3.22`, `docker:29`) | Renovate |
| CI tools | `.gitlab-ci.yml` (`pip install x==y`) | exact pin | Renovate |

The test matrix (`python:$PYTHON_VERSION`) and the Python versions in
`noxfile.py` and `[tool.cibuildwheel]` change only by hand.

[Renovate](https://docs.renovatebot.com/) reads `renovate.json`. Once a week it
opens MRs: one each for ruff and for mypy, one for a new Python in the CI default
image, and one for all other minor and patch updates. Major updates get their own MRs.
Nothing is merged automatically. The Dependency Dashboard issue lists pending and
rejected updates. Check a change to `renovate.json` with `nox -s renovate_check`
(needs Node.js); CI runs the same validator in the `renovate-config` job.

Renovate runs as the `renovate` CI job. Maintainers set it up once:

1. Project access token (Settings > Access tokens): role Developer, scopes `api`
   and `write_repository`. Store it as the masked, protected CI/CD variable
   `RENOVATE_TOKEN`.
2. A github.com personal access token without any scopes (pre-commit hooks are looked up
   as GitHub tags, and release notes come from GitHub). Store it as the masked,
   protected CI/CD variable `GITHUB_COM_TOKEN`.
3. Pipeline schedule (Build > Pipeline schedules) on `main`, weekly, with the
   variable `RENOVATE` = `true`. That pipeline runs only the `renovate` job.

## Licensing

pyEnvipat is GPL-2.0-only, the same as upstream enviPat. By contributing you
agree that your contribution is released under that license. See `LICENSE`,
`AUTHORS` and `NOTICE`.

## Security

Report vulnerabilities as described in [SECURITY.md](SECURITY.md), not in
public issues.
