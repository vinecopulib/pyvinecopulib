# Contributing

This guide covers the dev-side workflow: setting up an environment, the
build pipeline, the Makefile + pre-commit conventions, and the release
flow. End-user install instructions live in [README.md](README.md).

## Quick start

The Python environment is managed by [`uv`](https://docs.astral.sh/uv/) on
top of a Python interpreter you provide. Native build deps (Boost, Eigen,
libclang) are easiest to install via conda/mamba, but anything that
supplies them works.

1. Clone with submodules:

   ```bash
   git clone --recursive https://github.com/vinecopulib/pyvinecopulib.git
   cd pyvinecopulib
   ```

2. Provide native build deps and a Python ≥ 3.9 interpreter. The recommended
   recipe (mirrors the maintainers' setup):

   ```bash
   mamba create -n pyvinecopulib python=3.11 boost eigen 'python-clang=18.*'
   mamba activate pyvinecopulib
   ```

3. Install `uv` (one-off; if it's not already on your `PATH`):

   ```bash
   pip install uv
   ```

4. Sync all dependency groups + extras, then perform the editable install
   (the editable install must run after the sync because scikit-build-core
   needs the build-system deps in place):

   ```bash
   make dev-setup                # = make sync + pre-commit install
   ```

   Under the hood `make sync` runs:

   ```bash
   uv sync --all-extras --group dev --group test --group notebooks
   uv pip install -e . --no-build-isolation
   ```

5. Iterate:

   ```bash
   make quick-check              # ruff + ty + fast tests
   make check-all                # lint + type-check + tests with coverage
   ```

   Or call uv directly when you don't want the Make wrapper:

   ```bash
   uv run pytest
   uv run ruff format --check
   uv run ruff check
   uv run ty check
   uv build                      # sdist + wheel
   ```

## Repository layout

| Path | What lives there |
|---|---|
| `src/pyvinecopulib/` | Pure-Python package source (`__init__.py`, helpers). |
| `src/pyvinecopulib_ext.cpp`, `src/include/` | nanobind bindings + binding-side headers. |
| `lib/` | Vendored C++ submodules: `vinecopulib`, `wdm`, `kde1d`. |
| `tests/` | pytest test suite. |
| `examples/` | Jupyter notebooks rendered into the docs by nbsphinx. |
| `docs/` | Sphinx source. `conf.py` stages README/CHANGELOG/examples and generates `features.rst`/`examples.rst` at build time. |
| `scripts/` | Maintainer-facing utilities (see below). |
| `scripts/build/` | Build-internal tooling invoked by CMake. **Don't run these directly.** |

## How the build works

`pip install -e . --no-build-isolation` (the editable install) drives the
full pipeline:

1. CMake configure: `scripts/build/find_libclang.py` locates a libclang
   shared library (PyPI package, conda env, or `LIBCLANG_PATH` override).
2. CMake build, pre-compile: `scripts/build/generate_docstring.py` parses
   the C++ headers in `lib/` and writes `src/include/docstr.hpp`.
3. CMake build, compile: `nanobind_add_module` links
   `src/pyvinecopulib_ext.cpp` against the freshly generated `docstr.hpp`.
4. CMake POST_BUILD: `scripts/build/generate_stubs.py` stages the
   assembled package in a tempdir and writes
   `src/pyvinecopulib/__init__.pyi` from the live extension.

Both `docstr.hpp` and `__init__.pyi` are gitignored — the build is the
single source of truth. With `editable.rebuild = true` in
`pyproject.toml`, importing `pyvinecopulib` after a C++ edit triggers a
rebuild and refresh on the spot.

## Dependency layout

User-facing extras live under `[project.optional-dependencies]` in
`pyproject.toml` (`doc`, `examples`, `sklearn`) — these ship with the wheel
and end users install them via `pip install pyvinecopulib[<extra>]`.

Developer-only deps live under PEP 735 `[dependency-groups]`:

| Group | What it pulls |
|---|---|
| `test` | `pytest`, `pytest-cov`, `pytest-rerunfailures`, `pytest-xdist`, `coverage[toml]`. |
| `notebooks` | `nbmake`, `jupyter`, `nbconvert`, `nbformat` (for the notebook-test + regen workflows). |
| `dev` | `ruff`, `ty`, `pre-commit`, `twine`. |

Activate them with `uv sync --group <name>` (combine with `--all-extras` /
`--extra <name>` as needed). The generated `requirements.txt` /
`environment.yml` workflow was retired — `uv sync` is the single source of
truth.

## Makefile cheatsheet

`make help` lists everything. Most-used:

| Command | Purpose |
|---|---|
| `make sync` | `uv sync` all groups/extras + editable install. |
| `make dev-setup` | `make sync` + install pre-commit hooks. |
| `make quick-check` | Fast feedback: lint + type-check + tests without coverage. |
| `make check-all` | Pre-push: lint + type-check + tests with coverage. |
| `make test` / `make test-fast` / `make test-examples` | Test variants. |
| `make lint` / `make format` | Ruff. |
| `make type-check` | ty. |
| `make docs` / `make docs-serve` | Build / live-serve HTML documentation. |
| `make build` / `make sdist` / `make wheel` | `uv build` artifacts. |
| `make examples` | Re-execute notebooks (used by the docs/release flow; runs `scripts/regenerate_notebooks.py`). |
| `make clean` | Wipe Python caches and build dirs. |

## Pre-commit hooks

Hooks run on every commit. Manage with `make pre-commit-install` /
`make pre-commit`. The configured set:

- **ruff** — Python lint + format.
- **ty** — Astral's type checker (local hook; needs `ty` on PATH, which
  the dev env install provides).
- **clang-format** — C++ in `src/` (Google style; `docstr.hpp` excluded).
- **cmake-format** — CMake formatting.
- General whitespace / YAML / TOML / JSON checks.

## CI overview

`.github/workflows/pypi.yml` is the single workflow. Jobs:

- **`build`** — cibuildwheel matrix (16 wheels: Linux glibc/musl, macOS arm64,
  Windows × cp39/cp310/cp311/cp312-ABI3).
- **`check_wheels`** — counts and twine-checks wheel artifacts.
- **`verify_docs_build`** — RTD-equivalent doc build with `-W` (warnings as
  errors). Uploads `docs-html` artifact for PR review.
- **`install_and_unit_test`** — installs each wheel and runs pytest +
  notebook tests.
- **`regenerate_notebooks`** — fires on PRs to `main` (or any PR labelled
  `regenerate-notebooks`). Re-executes notebooks via the wheel and
  auto-commits the refreshed outputs back to the PR branch with `[skip ci]`.
- **`build_sdist`** — lints, type-checks, builds the sdist, installs from
  it as a sanity check.
- **`upload_to_pypi`** — publishes on tag push.

## Release process

Day-to-day flow:

1. Feature branches → PR → `dev`. CI runs everything except notebook
   regen.
2. When `dev` is ready: open a PR `dev → main`. The
   `regenerate_notebooks` job re-executes notebooks and commits the
   refreshed outputs to the PR branch automatically.
3. Merge `dev → main`, then tag the merge commit on `main`. The
   `upload_to_pypi` job publishes wheels + sdist to PyPI.
4. Read the Docs picks up the new tag automatically and rebuilds the
   `stable` version against the published wheel. No manual gh-pages
   deploy step.

`make release-check` runs the relevant local subset before opening the
release PR.

## Code style

- **Python**: PEP 8 (ruff-enforced).
- **C++**: Google style (clang-format-enforced). `src/include/docstr.hpp`
  is excluded — it's a build artifact.
- **Type hints**: required on Python source; ty checks them.
- **Docstrings**: required on public functions.

## Troubleshooting

- **`docstr.hpp` looks stale or import fails after C++ changes**: re-run
  `uv pip install -e . --no-build-isolation` — `editable.rebuild = true`
  handles most cases, but a fresh manual install fixes anything weird.
- **`uv` can't find the native build deps**: confirm `boost`, `eigen`, and
  a libclang ≤ 18 are reachable from your active env (`mamba list | grep
  -E 'boost|eigen|clang'`).
- **Nuclear option**: `git clean -fdx` (⚠️ removes everything not tracked
  by git).

## Tips

- Run `make quick-check` early and often.
- Pre-commit hooks fix most formatting issues automatically; just
  re-stage and commit again.
- Notebook output noise: don't commit re-executed notebooks from a local
  editable build (they capture absolute paths from the editable rebuild
  print). Let the CI `regenerate_notebooks` job own that.
