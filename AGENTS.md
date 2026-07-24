# AGENTS.md

Normative engineering spec for contributors and coding agents working on
this repository: scope, stability tiers, module boundaries, conventions,
and where to look for what.

For the **dev-side workflow** (environment setup, Makefile, pre-commit,
CI, release flow) see [CONTRIBUTING.md](CONTRIBUTING.md); for the
**user-facing pitch** and install instructions see [README.md](README.md).
This file does not duplicate those — it concentrates on engineering
invariants that survive across PRs.

## Project overview

`pyvinecopulib` is the Python interface to
[vinecopulib](https://github.com/vinecopulib/vinecopulib) — a
header-only C++ library for vine-copula and bivariate-copula
inference, built on Eigen. Two sister C++ libraries ship in the same
wheel:

- [`wdm`](https://github.com/tnagler/wdm) — weighted Kendall's τ /
  Spearman's ρ / Pearson etc.
- [`kde1d`](https://github.com/vinecopulib/kde1d-cpp) — 1-d kernel
  density estimation with boundary correction and discrete support.

All three are vendored as **git submodules under `lib/`** (see
`.gitmodules`). The Python package wraps them through a single
nanobind extension (`pyvinecopulib_ext.cpp`) and adds Python-only
extensions on top:

1. `pyvinecopulib.core`, `pyvinecopulib.families`, `pyvinecopulib.utils`
   — re-exports of the bound C++ surface, organised by topic; `core`
   additionally ships a backend-neutral pair-copula / vine abstraction
   layer (`BicopLike` / `VinecopLike` protocols, `BicopBase` /
   `VinecopBase` canonical bases, `ConditioningContext` policies) that
   custom NumPy / PyTorch backends subclass.
2. `pyvinecopulib.sklearn` — scikit-learn-compatible estimators
   (`VineDensity`, `VineRegressor`, `VineForestDensity`,
   `VineForestRegressor`) on top of the core, with a pluggable
   backend layer.
3. `pyvinecopulib.torch` — pure-PyTorch port of the evaluation
   cascade for GPU and autograd workflows.

Three design principles inform the rest of this file:

- **The C++ libraries are upstream.** `lib/vinecopulib`, `lib/wdm`, and
  `lib/kde1d` are git submodules. Behaviour changes belong upstream;
  this repo bumps the submodule pin and adjusts the bindings.
- **Generated files are build artefacts.** `src/include/docstr.hpp`
  (libclang-extracted C++ docstrings) and every `.pyi` stub under
  `src/pyvinecopulib/**/__init__.pyi` are gitignored. The build is the
  single source of truth — do not hand-edit, do not commit.
- **Code is quantitatively sensitive.** Pseudo-observation transforms,
  h-functions, Rosenblatt cascades, family parameterisations,
  pickling round-trips, and TLL grids all encode mathematically
  precise behaviour. Small "obvious-looking" changes can silently
  break copula identities. Treat numerical paths as
  correctness-critical and prefer round-trip / parity tests over
  structural ones.

### Stability tiers

Different subpackages have different change policies. Honour the tier
when proposing API changes:

| Surface | Tier | Policy |
|---|---|---|
| `pyvinecopulib.core`, `pyvinecopulib.families`, `pyvinecopulib.utils`, top-level `pyvinecopulib` (core class re-exports) | **Stable-ish** | Solid user base. Prefer deprecation aliases over breaks; document migrations in `CHANGELOG.md`. PR #207 is the model: the reorg kept old import paths working via `_deprecations.py` + `DeprecationWarning`. Breaks are allowed (e.g. the pybind11→nanobind migration; the #207 cleanup) but must be intentional, documented, and worth the churn. |
| `pyvinecopulib.sklearn` | **Active development** | API may change in breaking ways between minor releases. The latest break is the `#218` public backend system (estimators now take a single `backend=` instead of loose `controls=`/`structure=`/`seed=` kwargs); forest estimators (`#213`) are also new. |
| `pyvinecopulib.torch` | **Active development** | Same status. Defaults are still being tuned (cf. `990f997` device-aware `batched`, `cache_integrals=True`); the torch↔C++ cascade parity is a hard guarantee, but the `FitControlsTorchVinecop` surface and `TorchVinecop` method signatures may still shift. |
| `pyvinecopulib._python_helpers`, `pyvinecopulib._deprecations` | **Internal** | Underscore-prefixed. Not part of any contract; rename / restructure freely. `_deprecations.py` itself is slated for removal on the next major release. |

The "Solid user base" claim refers to the released `main` branch (see
the [GitHub project](https://github.com/vinecopulib/pyvinecopulib) for
the latest tag). Pre-release work on `dev` (and feature branches
ahead of it) is allowed to break sklearn/torch APIs as needed.

## Scope

### Included

- **Bivariate copula modelling** — every family bound from
  `lib/vinecopulib`: `indep`, `gaussian`, `student`, `clayton`,
  `gumbel`, `frank`, `joe`, `bb1/6/7/8`, `tawn`, `tll`; with their
  rotations, mixed-discrete handling, and family-set constraints via
  `FitControlsBicop`.
- **Vine copula modelling** — `Vinecop` with Dissmann selection
  (`mst_prim`), Wilson-weighted random structures
  (`random_weighted`), and user-supplied `RVineStructure` /
  `CVineStructure` / `DVineStructure`. Truncation, threading,
  bootstrap, pre-fit selection criteria, and family sets are exposed
  through `FitControlsVinecop`.
- **Univariate marginals** — `Kde1d` (`lib/kde1d`) with continuous,
  ordered-discrete, and unordered-categorical support.
- **Dependence measures** — `wdm` (`lib/wdm`).
- **Quasi-random sampling** — `sobol`, `ghalton`, `simulate_uniform`.
- **Pseudo-observations** — `to_pseudo_obs`.
- **Scikit-learn-compatible estimators** — `VineDensity`,
  `VineRegressor`, `VineForestDensity`, `VineForestRegressor` with a
  pluggable backend (`VinecopBackend` / `TorchVinecopBackend`).
- **PyTorch evaluator** — `TorchBicop`, `TorchVinecop` (pure-torch
  cascade with GPU placement, autograd, and an optional `batched`
  evaluation fast path; byte-for-byte parity with the C++ cascade).
- **Backend-neutral extension layer** — the `BicopLike` / `VinecopLike`
  contracts and canonical `BicopBase` / `VinecopBase` bases (NumPy or
  PyTorch) for hosting custom pair copulas in a vine, including
  **non-simplified / conditional** vines via a `ConditioningContext`
  (walk-through: `examples/11_extending_pyvinecopulib.ipynb`).

### Excluded (explicit)

- **Copula families outside the bound set.** New parametric families
  belong upstream in `lib/vinecopulib`; bindings then follow.
- **Custom C++ forks.** The repo always tracks the upstream
  `lib/vinecopulib` submodule pin; local C++ patches under
  `lib/` are not accepted.
- **Discrete margins on the torch backend.** `TorchVinecopBackend`
  is continuous-only (raises `NotImplementedError` when any
  `var_types[i] != "c"`). Use `VinecopBackend` for discrete /
  mixed-type problems.
- **Density estimators outside the vine framework.** General-purpose
  multivariate density models (normalising flows, Gaussian mixtures,
  …) are not in scope; `pyvinecopulib` is a vine-copula library.
- **Pinned legacy alias for every old import path forever.**
  Deprecation aliases live in `_deprecations.py` and warn on access;
  they are removed on the next major release.

## Package structure

```text
pyvinecopulib/
  AGENTS.md, CLAUDE.md           # this file + thin pointer (`@AGENTS.md`)
  README.md, CONTRIBUTING.md, CHANGELOG.md, LICENSE
  pyproject.toml                 # uv / scikit-build-core / ruff / ty / pytest / coverage config
  CMakeLists.txt                 # nanobind build, libclang docstring + .pyi stub generation
  Makefile                       # thin wrapper over `uv` (see CONTRIBUTING.md)
  .pre-commit-config.yaml        # ruff + ty + clang-format + cmake-format hooks

  lib/                           # upstream C++ — git submodules; do not patch locally
    vinecopulib/                 # core copula library
    wdm/                         # weighted dependence measures
    kde1d/                       # 1-d KDE with boundary correction

  src/
    pyvinecopulib_ext.cpp        # nanobind binding entry point (single .so)
    include/                     # binding-side C++ headers
      pyvinecopulib.hpp          # init_* declarations
      bicop/, vinecop/, kde1d/, misc/   # per-topic binding headers
      docstr.hpp                 # AUTO-GENERATED via scripts/generate_docstring.py (gitignored)

    pyvinecopulib/
      __init__.py                # top-level: core re-exports + lazy sklearn import + deprecation shim
      _deprecations.py           # warn-on-access aliases for pre-#207 top-level names
      py.typed                   # PEP 561 marker (built by scripts/generate_stubs.py)

      core/__init__.py           # Bicop, Vinecop, *VineStructure, FitControls* (re-exports from ext)
        protocols.py             # BicopLike / VinecopLike backend-neutral contracts
        bicop_base.py            # BicopBase (canonical BicopLike partial impl)
        vinecop_base.py          # VinecopBase (array-agnostic cascades + fit/select)
        context.py               # ConditioningContext / Simplified / NonSimplified
        _rootfind.py             # solve_increasing (monotone bisection; internal)
      families/__init__.py       # BicopFamily enum + 13 family constants + 15 group constants
      utils/__init__.py          # Kde1d, to_pseudo_obs, wdm, sobol, ghalton, simulate_uniform, benchmark
        _pair_plots.py           # pairs_copula_data plotting helper (pure Python)

      sklearn/__init__.py        # VineDensity, VineRegressor, VineForest{Density,Regressor}, backends
        backends.py              # VinecopBackend / TorchVinecopBackend + resolve_backend
        _base.py                 # VineBase (parameter-constraints, schema, 3-step pipeline)
        _mcs.py                  # Model Confidence Set survivor selection (da_mcs_marg / da_mcs_unif)
        density.py               # VineDensity
        regressor.py             # VineRegressor
        _forest_base.py          # VineForestBase ensemble machinery
        forest_density.py        # VineForestDensity
        forest_regressor.py      # VineForestRegressor

      torch/__init__.py          # TorchBicop, TorchVinecop, FitControlsTorch*
        bicop.py, vinecop.py     # nn.Module evaluators
        controls.py              # FitControlsTorchBicop / FitControlsTorchVinecop dataclasses
        _interp.py               # InterpolationGrid2D (bilinear; Sinkhorn margin renormalisation) — internal
        _fit_tll.py              # pure-torch TLL kernel
        _batched.py              # batched evaluation variants

      _python_helpers/           # internal; pure-Python wrappers used by the binding
        bicop.py, vinecop.py, kde1d.py, stats.py
      pyvinecopulib_ext.*.so     # compiled extension (gitignored build artefact)
      **/__init__.pyi            # type stubs AUTO-GENERATED via scripts/generate_stubs.py (gitignored)

  tests/                         # flat layout; one file per topic; shared fixtures in conftest.py
  docs/                          # Sphinx; conf.py drives features.rst via autosummary
  examples/                      # Jupyter notebooks (10), executed in CI and embedded via nbsphinx
  scripts/                       # build helpers + benchmarks
```

## Tooling

Full dev-workflow detail lives in [CONTRIBUTING.md](CONTRIBUTING.md).
The minimum an agent needs:

```bash
make sync           # editable build + all extras + pre-commit hooks
make check          # ruff lint + ruff format-check + ty
make test           # pytest tests/                 (serial — see note)
make docs           # sphinx -W
```

Conventions baked in by the toolchain:

- **`uv` is canonical.** Every `make` target dispatches to `uv run …`.
  Use `uv run pytest …` directly when iterating on a single test.
- **Editable installs need `--no-build-isolation`** so the active env's
  libclang ≤ 18 is reused. `scikit-build` reads `editable.rebuild =
  true` from `pyproject.toml` and re-runs the C++ build on import.
- **`ty check` requires a fresh editable build first**, because it
  reads the auto-generated `__init__.pyi` files for type info from the
  compiled extension. `make sync` does both; bare `make check` assumes
  the stubs are current.
- **Do not add `-n auto` to pytest.** xdist workers intermittently
  crash on GHA when re-importing native extensions ("node down: Not
  properly terminated"). The note is in
  [pyproject.toml](pyproject.toml) next to `addopts`.
- **`ruff` is pinned at `0.11.6`** (formatter output stability); only
  bump deliberately. Line length 80; indent width 2.
- **Tests for own deprecations are loud:** `filterwarnings` in
  `[tool.pytest.ini_options]` promotes `pyvinecopulib.*` deprecation
  warnings to errors, so internal code that still calls a deprecated
  path will fail CI.

For performance work: profile first, optimise demonstrated hotspots
only, and preserve every quantitative invariant (round-trip identities,
parity with the C++ cascade, pickling stability).

### Running heavy commands when the agent shell shares the host GPU

When the agent shell is on the same host whose GPU also runs the
maintainer's X session (or any interactive desktop session), a CUDA
lock-up in a pyvinecopulib subprocess can take the desktop driver
down with it — observed twice in prior sessions, both times the box
disappeared from ssh and required a hard reboot. The full pytest
suite is allowed and reliable on its own; what is **not** allowed is
queueing GPU-touching work alongside still-alive background tasks
from earlier turns. Specifically:

1. **After every `Bash(run_in_background=true)` call** (including
   monitoring shells like a backgrounded `nvidia-smi -l` or a
   `bench_torch_*` sweep), wait for the completion notification and
   then `TaskStop` the task explicitly if it did not exit. Never let
   a background task survive across turns.
2. **Before any GPU-touching command** (the full `pytest tests/`,
   any `scripts/bench_torch_*.py`, any test that drives
   `pyvinecopulib.torch`), run a clean-state check:
   `pgrep -af "uv run|pytest|bench_torch" | grep -v vscode` and
   `nvidia-smi`. Proceed only if no leftover Python/uv processes are
   running and VRAM is back at the X-session baseline.
3. **Never queue two GPU subprocesses concurrently.** No parallel
   `--devices cuda` benches via overlapping background bashes; one
   in flight at a time, strictly serial.
4. **Defer to the maintainer's own terminal for long bench sweeps.**
   When uncertain about cumulative resource use, ask the maintainer
   to run the sweep themselves where the process lifecycle is theirs
   to manage.

`make docs` (Sphinx) is unaffected — no torch involved.

## Working on this repo

### Inspection order

Before changing code, read in this order:

1. `AGENTS.md` (this file) — invariants and boundaries.
2. `docs/` — high-level intent, including the Sphinx `concepts.rst`
   primer on Sklar's theorem, pair copulas, R-vines, and TLL.
3. `src/pyvinecopulib/<subpackage>/__init__.py` — the module docstring
   is the canonical short description.
4. The implementation file you're about to touch, then the matching
   `tests/test_<topic>.py` for expected behaviour.

Match existing local patterns rather than introducing new ones.

### Definition of done

For any behaviour change:

- Diffs are scoped to the task; no opportunistic refactors that span
  unrelated files.
- Honour the [stability tier](#stability-tiers): for `core` /
  `families` / `utils`, prefer a deprecation alias over a hard break;
  for `sklearn` / `torch`, breaks are allowed but must be flagged in
  `CHANGELOG.md` (the `## Unreleased` section).
- Tests added or extended. Prefer extending an existing parametrised
  test over duplicating logic; share fixtures via `conftest.py`.
- Public-API changes update the module docstring (re-rendered in
  `docs/features.rst` via autosummary) and the matching example
  notebook when one exists.
- Run the [validation sequence](#tooling).

### Coding conventions

- **Indentation: 2 spaces. Line length: 80.** Enforced by ruff.
- **Type hints are required** on public Python source. `ty` checks
  them; the only allowed unresolved import is
  `pyvinecopulib.pyvinecopulib_ext` (the compiled `.so`).
- **`__init__.py` files use explicit `__all__`** to define the public
  surface; ruff's per-file ignore (`F403`/`F405`) covers the
  re-export pattern. No wildcard re-exports elsewhere.
- **Tests import from public namespaces** (`from
  pyvinecopulib.sklearn import VineDensity`), not deep internals.
  `_python_helpers` and other underscore-prefixed modules are off
  limits to tests.
- **Generated files stay generated.** `docstr.hpp` and every
  `__init__.pyi` are produced by `scripts/generate_docstring.py` and
  `scripts/generate_stubs.py` respectively. Do not hand-edit; do not
  commit. If a docstring or stub is wrong, fix the C++ source or the
  binding code, then rebuild.
- **Underscore-prefixed modules are internal.** Move helpers into
  `_python_helpers/` or a leading-underscore file inside the subpackage
  rather than exposing them.
- **Numpydoc docstring convention.** Public-API docstrings follow the
  [numpydoc spec](https://numpydoc.readthedocs.io/en/latest/format.html):
  short summary as the first line, `Parameters` / `Returns` /
  `Raises` / `Notes` / `Warnings` / `See Also` / `References` /
  `Examples` sections in that order, with every parameter and return
  value typed (`name : type` form, e.g.
  `u : ndarray, shape (n, 2), dtype float`). C++-derived docstrings
  inherit the convention through
  [scripts/generate_docstring.py](scripts/generate_docstring.py), which
  translates Doxygen tags to numpydoc sections and emits Python type
  annotations. `numpydoc.validation` is enabled as a pre-commit check;
  rule set + path exclusions live in `[tool.numpydoc_validation]` in
  `pyproject.toml`.
- **American English in prose.** Docstrings, comments, and docs use
  American spelling (`behavior`, `finalize`, `serialize`, `normalize`,
  `color` — not `behaviour` / `finalise` / `serialise` / …). Some legacy
  text still uses British spelling; new and edited prose should be
  American.
- **Reserve "backend" for the sklearn layer.** In `core` (and other
  user-facing) prose the word "backend" means the sklearn
  `VinecopBackend` / `TorchVinecopBackend` context, which most core users
  never touch — so don't use it for the NumPy-vs-PyTorch array
  distinction there. Say "array namespace", "array library", or just
  "NumPy or PyTorch" / "array-agnostic" instead.
- **Docs cross-references must resolve (nitpicky).** The Sphinx build runs
  with `nitpicky = True` (`docs/conf.py`), so `make docs` — which passes `-W`
  and is enforced by the `verify_docs_build` CI job — fails on *any* unresolved
  cross-reference. To keep it green: reference documented classes with **bare
  double backticks** (`` ``Vinecop`` ``), which `process_cross_references`
  (`docs/conf.py`) rewrites to a fully-qualified
  `:class:`~pyvinecopulib.core.Vinecop`` that resolves from any page — *not*
  `` :class:`pyvinecopulib.Vinecop` `` (wrong path; the documented target is
  `pyvinecopulib.core.Vinecop`) and *not* single backticks (which render
  italic, not a link). Refs to private methods / internal helpers must be
  plain ``literals`` (no `:meth:` / `:func:` role — they have no doc page).
  External types (`numpy.*`, `torch.*`, builtins) resolve via
  `intersphinx_mapping`; the autosummary class template
  (`docs/_templates/autosummary/class.rst`) gives *Attributes* a `:toctree:`
  so property refs get pages; nanobind's `numpy.ndarray[dtype=…]` signatures
  are collapsed to `numpy.ndarray` by an `autodoc-process-signature` hook. The
  only sanctioned suppression is the short `nitpick_ignore_regex` list
  (upstream-C++ getter-name mismatches, `BicopFamily` value aliases,
  scikit-learn-generated methods); prefer fixing a reference over extending it.

### Maintaining this file

If a coding agent or reviewer keeps repeating the same correction —
about a convention this repo enforces — update `AGENTS.md` rather
than relying on tribal knowledge. Do not add ephemeral, user-specific,
or machine-local preferences here. The `CHANGELOG.md` is the place for
release-by-release context; this file is for invariants.

## Module boundaries

### Upstream C++ (`lib/`)

`lib/vinecopulib`, `lib/wdm`, and `lib/kde1d` are header-only C++
libraries pinned as **shallow git submodules** (`.gitmodules`).
Behaviour and API changes belong upstream. The Python repo only:

1. Bumps the submodule SHA (clearly motivated in the PR).
2. Adjusts `src/include/**` and `src/pyvinecopulib_ext.cpp` to track
   binding-relevant upstream changes.
3. Updates the Python re-exports / docstrings if the surface changed.

Cloning requires `--recursive` (see README); CI does this
automatically.

### `pyvinecopulib_ext` (the nanobind extension)

- Single binding module compiled from `src/pyvinecopulib_ext.cpp` +
  `src/include/**`.
- `nanobind` (>= 2.7) is the binding system, **not** pybind11. The
  pybind11 → nanobind switch was a deliberate API-affecting change;
  do not partially revert it. Use `nb::` types and `nanobind_add_module`
  conventions.
- `src/include/docstr.hpp` is generated by `scripts/generate_docstring.py`
  via libclang at CMake configure time. To change a docstring, edit the
  matching C++ comment in `lib/<library>` (upstream) and rebuild.
- The module's pickled identifier is `pyvinecopulib.core.<Class>`
  etc., not `pyvinecopulib.pyvinecopulib_ext.<Class>`. Scoped module
  overrides in the binding ensure this; pickling round-trip is
  guaranteed across the canonical paths.

### `pyvinecopulib.core`

- Pure re-exports of `Bicop`, `Vinecop`, `RVineStructure`,
  `CVineStructure`, `DVineStructure`, `FitControlsBicop`,
  `FitControlsVinecop` from `pyvinecopulib_ext`.
- Use `Bicop.from_family(...)` / `Bicop.from_data(...)` and
  `Vinecop.from_data(...)` (or `RVineStructure.simulate(...)`) — these
  factories are the documented entry points; the raw constructor
  signatures are kept for nanobind-level access only.
- `tree_algorithm` on `FitControlsVinecop`: `"mst_prim"` (default,
  Dissmann) and `"random_weighted"` (Wilson-weighted random tree;
  used by the sklearn forests).
- **Backend-neutral abstraction layer** (pure Python; `core` imports
  without PyTorch). The extension point for custom (e.g. neural,
  conditional) pair copulas and vines:
  - `BicopLike[ArrayT]` / `VinecopLike[ArrayT]` (`protocols.py`) —
    generic, `runtime_checkable` protocols mirroring the `Bicop` /
    `Vinecop` evaluation surface on any array backend (NumPy or
    PyTorch); `Bicop` / `Vinecop` satisfy them structurally.
  - `BicopBase` (`bicop_base.py`) / `VinecopBase` (`vinecop_base.py`) —
    canonical partial implementations to subclass. A `BicopBase`
    subclass defines `pdf` / `hfunc1` / `hfunc2` and inherits `hinv1` /
    `hinv2` (bisection), `simulate`, `loglik`, `plot` (`flip` — needed
    only to host the pair in structure *selection* — defaults to
    raising); a `VinecopBase` subclass defines the one hook
    `_get_pair_copula` and inherits the whole tree-by-tree cascade plus
    the public `fit` and `select` engines. `select` is an
    exact port of `Vinecop`'s Dissmann / Wilson structure selection
    (same matrix encoding, selection-time pairs reused via `flip`, no
    re-fit; parity is a hard guarantee). `TorchBicop` / `TorchVinecop`
    are the torch subclasses.
  - `ConditioningContext` / `SimplifiedContext` (default) /
    `NonSimplifiedContext` (`context.py`) — the per-edge policy that
    turns the simplified cascade into a **non-simplified / conditional**
    vine (each pair also sees its conditioning-set values `u_D` and any
    external covariates `x`). Walk-through:
    `examples/11_extending_pyvinecopulib.ipynb`.
  - `solve_increasing` (`_rootfind.py`) — vectorized monotone bisection
    behind the default `hinv1` / `hinv2` (internal; not re-exported).

### `pyvinecopulib.families`

- Re-exports the `BicopFamily` enum, 13 family-tag constants, and
  15 family-group lists from `pyvinecopulib_ext`. See the module
  docstring for the full group definitions (`itau`, `lt`, `ut`,
  `rotationless`, `analytic_derivs`, …).
- The family lists are the canonical way to constrain the fit
  search space: pass them to `FitControlsBicop(family_set=...)` or
  `FitControlsVinecop(family_set=...)`. Pre-PR-#207 top-level
  aliases (`pyvinecopulib.gaussian`) still resolve but emit
  `DeprecationWarning` via `_deprecations.py`.

### `pyvinecopulib.utils`

- Re-exports `Kde1d`, `to_pseudo_obs`, `wdm`, `sobol`, `ghalton`,
  `simulate_uniform`, `benchmark` (all C++) plus the pure-Python
  `pairs_copula_data` helper from `_pair_plots.py`.
- `Kde1d` is used internally by the sklearn estimators as the
  marginal estimator; it also stands alone for any 1-d KDE problem.
- `to_pseudo_obs(data)` is the canonical input transform for
  copula fitting (rank-normalise to the unit hypercube).

### `pyvinecopulib.sklearn`

User-facing estimators on top of the core, organised so that the
3-step pipeline — fit 1-d KDE marginals → transform to
pseudo-observations → fit a vine on the copula data — happens once,
in `VineBase` (`_base.py`).

Class hierarchy:

```text
sklearn.base.BaseEstimator
├── VineBase (_base.py)            # shared pipeline + DataFrame schema
│   ├── VineDensity   (+ DensityMixin)
│   └── VineRegressor (+ RegressorMixin)
└── VineForestBase (_forest_base.py)  # ensemble + MCS survivor selection
    ├── VineForestDensity   (+ DensityMixin)
    └── VineForestRegressor (+ RegressorMixin)
```

#### `pyvinecopulib.sklearn.backends`

Public extension point introduced in `#218`. Estimators do not call
`pyvinecopulib.Vinecop.from_data(...)` directly — they go through a
backend object. Two concrete backends ship, both subclassing the
private `_VinecopBackendBase` (which owns `structure_of` and the
copy-on-write forest plumbing; concrete backends override only the
divergent members + hooks):

- `VinecopBackend(controls=None, structure=None)` — default. Wraps
  `pyvinecopulib.Vinecop`; no PyTorch dependency.
- `TorchVinecopBackend(controls=None, structure=None)` — wraps
  `pyvinecopulib.torch.TorchVinecop`. Constructing this class imports
  `torch` — it is the explicit opt-in signal that PyTorch is
  required. The sklearn subpackage itself does not import torch. Since
  `TorchVinecop.from_data` now auto-selects a structure (mirroring
  `pv.Vinecop.from_data`), this backend delegates structure selection
  to the vine rather than branching itself.

Both backends expose:

```python
backend.fit_vine(U, var_types=...) -> VinecopLike
backend.pdf(vine, U)        -> np.ndarray
backend.cdf(vine, U, N=..., seeds=...) -> np.ndarray
backend.simulate(vine, n, seeds=...)   -> np.ndarray
backend.structure_of(vine)  -> RVineStructure
backend.with_random_structure(d, seeds)  -> Backend  # copy-on-write
backend.with_local_random(seeds)         -> Backend  # ditto
backend.with_num_threads(n)              -> Backend  # ditto (torch: no-op)
```

`pyvinecopulib.core.VinecopLike` is the canonical `runtime_checkable`
Protocol describing the post-fit vine surface (`pdf` / `cdf` /
`rosenblatt` / `inverse_rosenblatt` / `simulate`, plus a `structure`
attribute); both `pv.Vinecop` and `pv.torch.TorchVinecop` satisfy it
structurally (no inheritance). `fit_vine` returns conforming vines, so
downstream code that only needs evaluation can type against
`pyvinecopulib.core.VinecopLike` instead of either concrete class. The
backend layer no longer defines its own copy.

`resolve_backend(backend)` is the dispatch helper: `None` →
default-constructed `VinecopBackend`; any other value is returned
as-is. Estimators call this once at `fit()` time and pin the resolved
backend as `self.backend_`.

#### Estimator conventions (scikit-learn developer guide)

`#218` aligned every estimator with the
[scikit-learn third-party-estimator developer guide](https://scikit-learn.org/stable/developers/develop.html).
The rules that bind:

- `__init__` performs **no** validation — it stores arguments
  verbatim. All checks live in `_parameter_constraints` +
  `_validate_params()`, called at the top of `fit()`.
- The single backend parameter is `backend=` (no
  `backend="cpp"` / `"torch"` shortcuts; pass an instance).
- Fitted attributes follow sklearn naming with a trailing underscore:
  `feature_names_in_`, `schema_`, `random_state_` (resolved RNG),
  `backend_` (resolved backend pinned at fit time).
- Every post-fit method calls
  `sklearn.utils.validation.check_is_fitted` before doing anything.
- `base_params` on the forests is defensively-copied at fit time so
  callers can't mutate the forest's view of it; this also makes
  `sklearn.base.clone()` round-trip cleanly.
- `random_state` is the canonical RNG kwarg name throughout (no
  legacy `seed=` / `seeds=` on `VineDensity.sample`,
  `VineDensity.cdf`, `VineForestDensity.cdf`, or the forest
  constructors).

#### Forest survivor selection (`_mcs.py`)

Forests sample `n_vines` candidate vine structures (Joe's uniform
algorithm by default; `"local"` adds Kendall-τ weighting via
Wilson's walk), fit them, then prune by **Model Confidence Set**
on a held-out fraction (`val_fraction=0.25`):

- `da_mcs_marg(losses, alpha)` — dual-argmin survivor test with
  per-model coverage. Cheaper; default.
- `da_mcs_unif(losses, alpha)` — uniform / familywise coverage via
  adaptive pre-screening. Use when the candidate count is large.

### `pyvinecopulib.torch`

Pure-PyTorch port of the evaluation cascade. Every public class is a
`torch.nn.Module`, so `.to("cuda")`, autograd, and composition with
other torch models are first-class. The submodule **hard-requires**
torch at import time (raises `ImportError` with an install hint).

Key surface:

- `TorchBicop` / `TorchVinecop` — evaluators.
  - `TorchBicop` is a density on a grid; constructors:
    `TorchBicop(grid_points, values, cache_integrals=True, ...)`,
    `TorchBicop.from_bicop(cop, ...)` (lift a C++ `Bicop`),
    `TorchBicop.from_data(u, controls=None, ...)` (fit; dispatches on
    `controls.method`).
  - `TorchVinecop` mirrors `pv.Vinecop`'s `pdf` / `cdf` /
    `rosenblatt` / `inverse_rosenblatt` / `simulate` signatures.
- `FitControlsTorchBicop` / `FitControlsTorchVinecop` — fit-time
  dataclasses. Notable knobs:
  - `method` — `"tll"` (the only fitter; kept as the dispatch seam
    for future torch fitters).
  - `cache_integrals` — default `True` (set in `990f997`); precomputes
    integral grids for ~80–300× evaluation speed-up with mean IAE
    `< 1e-3`.
  - `batched` — fires a single batched bicop call per tree level
    (available on `pdf` / `rosenblatt`, not on `inverse_rosenblatt`).
    The non-batched cascade is a byte-for-byte port of the C++
    evaluator; the batched path agrees with it to floating-point
    tolerance (`tests/test_torch_vinecop.py`).
  - `device`, `dtype` — propagate to every tensor on construction;
    fitted modules respect `.to(device)` afterwards.
  - `structure_controls` — `FitControlsTorchVinecop` only; the
    structure-selection controls used when `TorchVinecop.from_data`
    is called with `structure=None` (the R-vine is selected on the
    compiled `pv.Vinecop`, then lifted). `None` defaults to TLL with
    `trunc_lvl=20`.
- `InterpolationGrid2D` (`torch/_interp.py`) — the 2-d bilinear grid
  backing `TorchBicop`; **internal** (not re-exported). Margin
  normalisation uses Sinkhorn iterations to drive marginals to uniform.

### Top-level `pyvinecopulib`

`src/pyvinecopulib/__init__.py` is *not* empty: it re-exports the core
class surface so the long-standing
`from pyvinecopulib import Bicop, Vinecop, to_pseudo_obs` pattern keeps
working. Specifically:

- `__all__` covers the eight core classes plus `to_pseudo_obs` and the
  three always-loaded subpackages (`core`, `families`, `utils`); plus
  the lazily-loaded `sklearn`.
- `__getattr__` provides two things: lazy import of `sklearn` (the
  extra is only triggered on `import pyvinecopulib.sklearn` or
  attribute access) and a deprecation shim for the 35 pre-#207
  top-level names (every family constant, every utility function from
  `utils`). Each access emits a `DeprecationWarning` pointing at the
  new canonical path.
- For new code, **always import from the canonical subpackage**
  (`from pyvinecopulib.families import gaussian`,
  `from pyvinecopulib.utils import Kde1d`), not the top-level alias.

`torch` is **not** re-exported at the top level — `import
pyvinecopulib.torch` is the only entry. Same for `sklearn.backends`.

### Internal: `_python_helpers`, `_deprecations`

- `_python_helpers/{bicop,vinecop,kde1d,stats}.py` — pure-Python
  helpers (DataFrame conversions, plotting glue) called by the
  nanobind extension via `nb::module_::def(...)` lookup. Not part of
  any public contract. Move new internal helpers here rather than
  exposing them.
- `_deprecations.py` — `_DEPRECATED_TOP_LEVEL` dict + `_resolve_deprecated`
  helper for the top-level `__getattr__` shim. Slated for removal on
  the next major release; new deprecation aliases can be added here
  in the meantime, but each entry is a debt to be paid down.

## Public APIs

The canonical, exhaustive API listing is `docs/features.rst`, which is
generated at Sphinx-build time via `autosummary` from the
`DOCSTRING_SUBPACKAGES` table in `docs/conf.py`. The Sphinx HTML build
is the source of truth for signatures and parameter docs; the listings
below are a quick orientation.

- **`pyvinecopulib.core`** — `Bicop`, `Vinecop`, `RVineStructure`,
  `CVineStructure`, `DVineStructure`, `FitControlsBicop`,
  `FitControlsVinecop`; plus the backend-neutral abstraction layer
  `BicopLike`, `VinecopLike`, `BicopBase`, `VinecopBase`,
  `ConditioningContext`, `SimplifiedContext`, `NonSimplifiedContext`.
- **`pyvinecopulib.families`** — `BicopFamily` enum; per-family
  constants (`indep`, `gaussian`, `student`, `clayton`, `gumbel`,
  `frank`, `joe`, `bb1`, `bb6`, `bb7`, `bb8`, `tawn`, `tll`); group
  lists (`all`, `parametric`, `nonparametric`, `one_par`, `two_par`,
  `three_par`, `elliptical`, `archimedean`, `extreme_value`, `bb`,
  `rotationless`, `lt`, `ut`, `itau`, `analytic_derivs`).
- **`pyvinecopulib.utils`** — `Kde1d`, `to_pseudo_obs`, `wdm`,
  `sobol`, `ghalton`, `simulate_uniform`, `benchmark`,
  `pairs_copula_data`.
- **`pyvinecopulib.sklearn`** — `VineDensity`, `VineRegressor`,
  `VineForestDensity`, `VineForestRegressor`; plus the `backends`
  submodule (`VinecopBackend`, `TorchVinecopBackend`,
  `resolve_backend`).
- **`pyvinecopulib.torch`** — `TorchBicop`, `TorchVinecop`,
  `FitControlsTorchBicop`, `FitControlsTorchVinecop`.

Top-level `pyvinecopulib` re-exports the eight core classes and
`to_pseudo_obs`; everything else — including the `core` abstraction
layer (`BicopBase` / `VinecopBase`, the protocols, the contexts) — is
reachable only through the subpackages.

## Tests

The suite lives flat under `tests/` (one file per topic; no per-module
subdirectories). Shared fixtures in `tests/conftest.py`:

- `random_state` — fixed `RandomState(42)` for reproducibility.
- `sample_array_data` — 300-sample 2-d multivariate-normal NumPy array.
- `sample_dataframe_data` — mixed numeric / ordered-categorical /
  unordered-categorical pandas DataFrame.
- `regression_data` — linear-regression fixture
  `(X, y, true_coef, noise_std)` for `VineRegressor` tests.
- `unique_json_path` — per-test path in `tmp_path` for serialisation
  round-trips.

Conventions:

- **Tests import from public namespaces.** `from
  pyvinecopulib.sklearn import VineDensity` — not `from
  pyvinecopulib.sklearn._base import VineBase`. Internal helpers are
  validated indirectly through the public surface.
- **Coverage** is collected via `--cov=src/pyvinecopulib` with
  `term-missing`, `html`, and `xml` reports
  (configured in `pyproject.toml`). No threshold is enforced as a CI
  blocker, but new code should keep parity with the existing roughly-
  full coverage on every `src/pyvinecopulib/**` file.
- **Serial pytest only.** Do not pass `-n auto` to `pytest`; xdist
  workers crash on native-extension re-import. Local iterative runs
  can use `pytest tests/test_<topic>.py::<name>` directly.
- **Test extras are required to run their tests.** `make sync`
  installs `[sklearn]` and `[torch]` so the suite runs in full;
  CI installs the same. `test_sklearn_*` requires `scikit-learn`;
  `test_torch_*` requires `torch`.
- **DeprecationWarnings from this package are errors.** The
  `filterwarnings` entry in `[tool.pytest.ini_options]` enforces
  this — internal code paths must already be on the post-`#207`
  imports.
- **Example notebooks are test targets.** `make test-examples`
  re-executes `examples/*.ipynb` via `pytest --nbmake`. CI runs
  this on every PR; notebook outputs are auto-regenerated on
  PRs targeting `main` (the `regenerate_notebooks` workflow job).

Round-trip / parity properties to preserve when touching numerics:

- Pickling round-trip: `Bicop`, `Vinecop`, `Kde1d`, the sklearn
  estimators, and the torch evaluators all round-trip through
  `pickle.dumps` / `pickle.loads`.
- C++ ↔ torch parity: `TorchVinecop.from_vinecop(cpp_vine)` yields a
  module whose `pdf` / `cdf` / `rosenblatt` outputs match the C++
  vine to within floating-point tolerance, on the operational grid.
- `batched=True` ↔ `batched=False`: numerically equivalent `pdf` /
  `rosenblatt` outputs on the same fitted vine.
- `sklearn.base.clone()` round-trip: every estimator clones cleanly
  with all `__init__` parameters preserved verbatim.

## Extension points

- **New copula families.** Add the C++ implementation upstream in
  `lib/vinecopulib`, bump the submodule, then extend
  `src/include/bicop/family.hpp` + `src/pyvinecopulib/families/__init__.py`
  to bind the new tag and group it appropriately.
- **Custom pair copulas / vines (`pyvinecopulib.core`).** Subclass
  `BicopBase` (define `pdf` / `hfunc1` / `hfunc2`) for a custom pair
  copula, and host it by subclassing `VinecopBase` (define the one hook
  `_get_pair_copula`); both run on NumPy or PyTorch and inherit the full
  evaluation surface. Implement `BicopLike` / `VinecopLike` directly for
  an immutable / functional backend. For a **non-simplified /
  conditional** vine, pass a `NonSimplifiedContext` and drive
  `VinecopBase.fit` with a `fit_edge` callback; see
  `examples/11_extending_pyvinecopulib.ipynb`. `TorchBicop` /
  `TorchVinecop` are the reference torch subclasses.
- **New `pyvinecopulib.sklearn` backends.** Subclass the private
  `_VinecopBackendBase`, which already provides `structure_of` and the
  copy-on-write forest plumbing (`with_random_structure` /
  `with_local_random` / `with_num_threads`); override the divergent
  members (`fit_vine`, `pdf`, `cdf`, `simulate`, `_default_controls`,
  and the `_tree_selection_controls` / `_install_tree_selection_controls`
  hooks the local-random forest walks through). `resolve_backend`
  accepts any such object (it only defaults `None`). Consider whether
  the underlying vine satisfies the `pyvinecopulib.core.VinecopLike`
  Protocol so downstream code can stay type-stable.
- **New torch fit methods.** Add the implementation under
  `pyvinecopulib/torch/_fit_<name>.py`, add `"<name>"` to `METHODS`
  in `controls.py`, and add the dispatch branch in
  `TorchBicop.from_data` (the lone `"tll"` path is the template).
- **New sklearn-style ensembles.** Subclass `VineForestBase` (or
  `VineBase` for single-vine variants) and override the loss /
  prediction methods. Stay inside the
  `_parameter_constraints` / `_validate_params()` pattern and pin
  fitted attributes with trailing underscores.
- **Wider quasi-random integration.** `simulate_uniform` is the
  high-level driver behind `Vinecop.cdf` Monte-Carlo evaluation and
  random structure sampling; new low-discrepancy methods slot in
  alongside `sobol` / `ghalton` in `lib/vinecopulib`'s `misc/stats.hpp`.
