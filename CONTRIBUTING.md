# Contributing to pyvinecopulib

This guide covers the dev-side workflow: setting up an environment, the
build pipeline, the Makefile + pre-commit conventions, and the release
flow. End-user install instructions live in [README.md](README.md).

## Quick start

1. Clone with submodules and create a conda env:

   ```bash
   git clone --recursive https://github.com/vinecopulib/pyvinecopulib.git
   cd pyvinecopulib
   make env-conda                # creates the `pyvinecopulib` conda env
   conda activate pyvinecopulib
   ```

2. Install development dependencies and the package in editable mode, plus
   pre-commit hooks:

   ```bash
   make dev-setup
   ```

3. Iterate:

   ```bash
   make quick-check              # ruff + ty + fast tests
   make check-all                # lint + type-check + tests with coverage
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

## Makefile cheatsheet

`make help` lists everything. Most-used:

| Command | Purpose |
|---|---|
| `make dev-setup` | One-shot: install dev/doc/examples extras + pre-commit hooks. |
| `make quick-check` | Fast feedback: lint + type-check + tests without coverage. |
| `make check-all` | Pre-push: lint + type-check + tests with coverage. |
| `make test` / `make test-fast` / `make test-examples` | Test variants. |
| `make lint` / `make format` | Ruff. |
| `make type-check` | ty. |
| `make docs` | Build HTML documentation. |
| `make docs-serve` | Live-reload server (`sphinx-autobuild`). |
| `make clear-cache` | Wipe Python caches and build dirs. |
| `make examples` | Re-execute notebooks (used by the docs/release flow; runs `scripts/regenerate_notebooks.py`). |

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
  `pip install -e . --no-build-isolation` — `editable.rebuild = true`
  handles most cases, but a fresh manual install fixes anything weird.
- **Build issues**: `make debug-build` for verbose output.
- **Project status**: `make status` shows git + Python + installed deps.
- **Nuclear option**: `make git-clean` (⚠️ removes everything not
  tracked by git).

## Tips

- Run `make quick-check` early and often.
- Pre-commit hooks fix most formatting issues automatically; just
  re-stage and commit again.
- Notebook output noise: don't commit re-executed notebooks from a local
  editable build (they capture absolute paths from the editable rebuild
  print). Let the CI `regenerate_notebooks` job own that.
