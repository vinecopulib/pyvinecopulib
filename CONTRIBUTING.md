# Contributing

This guide covers the dev-side workflow: setting up an environment, the
build pipeline, the Makefile, pre-commit conventions, and the release flow.
End-user install instructions live in [README.md](README.md).

## Quick start

```bash
git clone --recursive https://github.com/vinecopulib/pyvinecopulib.git
cd pyvinecopulib
mamba create -n pyvinecopulib python=3.11 boost eigen 'python-clang=18.*' uv
mamba activate pyvinecopulib
make sync
```

`make sync` runs `uv sync --all-extras --group dev --group test --group notebooks`,
performs the editable install, and installs pre-commit hooks.

Iterate:

```bash
make check          # ruff + format-check + ty
make test           # pytest
make docs           # sphinx
```

Or drop the make wrapper and call `uv` directly:

```bash
uv run pytest
uv run ruff check
uv run ty check
uv build            # sdist + wheel
```

## Repository layout

| Path | What lives there |
|---|---|
| `src/pyvinecopulib/` | Pure-Python package source. |
| `src/pyvinecopulib_ext.cpp`, `src/include/` | nanobind bindings. |
| `lib/` | Vendored C++ submodules: `vinecopulib`, `wdm`, `kde1d`. |
| `tests/` | pytest suite. |
| `examples/` | Notebooks rendered into the docs by nbsphinx. |
| `docs/` | Sphinx source. |
| `scripts/` | Build helpers (`find_libclang.py`, `generate_docstring.py`, `generate_stubs.py`) + `regenerate_notebooks.py`. |

## How the build works

`uv pip install -e . --no-build-isolation` drives the full pipeline:

1. CMake configure: `scripts/find_libclang.py` locates a libclang shared library.
2. Pre-compile: `scripts/generate_docstring.py` writes `src/include/docstr.hpp`.
3. Compile: `nanobind_add_module` links the extension.
4. POST_BUILD: `scripts/generate_stubs.py` writes one `__init__.pyi` per
   subpackage from the live extension.

`docstr.hpp` and the `.pyi` files are gitignored — the build is the single
source of truth. `editable.rebuild = true` re-runs the build on import after
C++ edits.

Because `ty` reads those `.pyi` files for type info, run an editable build
before `ty check`:

```bash
uv pip install -ve . --no-build-isolation
uv run ty check src tests
```

`make sync` does both in one shot; `make check` then calls `ty` directly.

## Dependency layout

User-facing extras (`doc`, `examples`, `sklearn`) live under
`[project.optional-dependencies]` — installable via
`pip install pyvinecopulib[<extra>]`.

Developer-only deps live under PEP 735 `[dependency-groups]`:

| Group | Pulls |
|---|---|
| `build` | `[build-system].requires` + `ninja` + `cmake`. |
| `test` | pytest stack. |
| `notebooks` | jupyter, nbmake, nbconvert, nbformat. |
| `dev` | ruff, ty, pre-commit, twine (`include-group = "build"`). |

Combine groups and extras as needed: `uv sync --group test --extra examples`.

## Makefile targets

`make help` lists everything. The full set:

| Command | Purpose |
|---|---|
| `make sync` | Install all deps + editable build + pre-commit hooks. |
| `make check` | Read-only lint + format-check + ty. |
| `make format` | Apply ruff autofixes + format. |
| `make test` / `make test-examples` | Pytest suite / notebook tests. |
| `make docs` | Sphinx build (-W). |
| `make build` / `make sdist` | Artifacts via `uv build`. |
| `make notebooks` | Re-execute example notebooks. |
| `make clean` | Wipe build artifacts and Python caches. |

## Pre-commit hooks

Configured set: ruff (lint + format), ty, clang-format, cmake-format, and
general whitespace/TOML/JSON checks. `make sync` installs them.

## CI overview

`.github/workflows/pypi.yml` is the single workflow. Jobs:

- **`build`** — cibuildwheel matrix (12 wheels: Linux glibc/musl, macOS arm64,
  Windows × cp310/cp311/cp312-ABI3).
- **`check_wheels`** — counts and twine-checks artifacts.
- **`verify_docs_build`** — RTD-mirror doc build with `-W`.
- **`install_and_unit_test`** — installs each wheel and runs pytest + notebook
  tests.
- **`regenerate_notebooks`** — fires on PRs to `main` (or any PR labelled
  `regenerate-notebooks`); auto-commits refreshed outputs.
- **`build_sdist`** — runs `make check && make sdist` and tests the installed
  sdist.
- **`upload_to_pypi`** — publishes on tag push.

## Release process

1. Feature branches → PR → `dev`.
2. PR `dev → main` — `regenerate_notebooks` refreshes outputs automatically.
3. Merge and tag on `main` — `upload_to_pypi` ships to PyPI and Read the Docs
   picks up the tag.

## Code style

- Python: PEP 8 (ruff-enforced).
- C++: Google style (clang-format).
- Type hints required on Python source; `ty` checks them.

## Docstring convention

All public-API docstrings follow the **numpydoc** convention
(<https://numpydoc.readthedocs.io/en/latest/format.html>). Section
headings (`Parameters`, `Returns`, `Raises`, `Notes`, `Warnings`,
`See Also`, `References`, `Examples`), parameter and return
formatting, and reference numbering use numpydoc's strict grammar.
Every parameter and return value carries a type (use the
`name : type` form, e.g. `u : ndarray, shape (n, 2), dtype float`).

C++-derived docstrings inherit the convention via
`scripts/generate_docstring.py`, which translates Doxygen tags
(`@param`, `@return`, `@note`,
`@warning`, `@exception`, `@see`) to the matching numpydoc sections
and emits Python type annotations on each parameter / return.
Edit the C++ source comment to change the rendered Python
docstring; do not hand-edit the generated `src/include/docstr.hpp`.

Enforcement runs through `numpydoc.validation` at pre-commit time
(and ad-hoc via `uv run python -m numpydoc lint <path>`); the
active rule set is configured in `[tool.numpydoc_validation]` in
`pyproject.toml`. Internal modules (`_python_helpers`, `_base`,
`_batched`, `_fit_tll`, `_util`, `_interp`, `_deprecations`) are
excluded by path —
they're off-contract per AGENTS.md §"Module boundaries".

## Troubleshooting

- **Build fails after C++ changes**: re-run `make sync` (or `uv pip install
  -e . --no-build-isolation` directly).
- **uv can't find native build deps**: confirm `boost`, `eigen`, and a
  libclang ≤ 18 are reachable in the active env.
- **Nuclear option**: `git clean -fdx` (⚠️ wipes everything untracked).
