# AGENTS.md

Normative engineering spec for contributors and coding agents working on
this repository: scope, stability tiers, module boundaries, conventions,
and where to look for what.

For the **mechanics** of the dev workflow (environment setup, Makefile
targets, pre-commit, the CI job graph) see
[CONTRIBUTING.md](CONTRIBUTING.md); for the **user-facing pitch** and
install instructions see [README.md](README.md); for the
release-by-release history see [CHANGELOG.md](CHANGELOG.md). This file does
not duplicate those — it holds the engineering invariants that survive
across pull requests, including the policy half of the branching and
release model.

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
   — re-exports of the bound C++ surface, organized by topic; `core`
   additionally ships a backend-neutral pair-copula / vine abstraction
   layer (`BicopLike` / `VinecopLike` protocols, `BicopBase` /
   `VinecopBase` canonical bases, `ConditioningContext` policies) that
   custom NumPy / PyTorch backends subclass, and `Vinedist` — a vine
   copula combined with univariate margins, i.e. a full multivariate
   distribution on the data scale.
2. `pyvinecopulib.margins` — the univariate marginal layer `Vinedist`
   composes: the built-in margins, family selection, and an adapter
   registry that presents a SciPy / PyTorch / other-ecosystem
   distribution object as a margin.
3. `pyvinecopulib.sklearn` — scikit-learn-compatible estimators
   (`VineDensity`, `VineRegressor`) on top of the core, with a
   pluggable
   backend layer.
4. `pyvinecopulib.torch` — pure-PyTorch port of the evaluation
   cascade for GPU and autograd workflows.

Three design principles inform the rest of this file:

- **The C++ libraries are upstream.** `lib/vinecopulib`, `lib/wdm`, and
  `lib/kde1d` are git submodules. Behavior changes belong upstream;
  this repo bumps the submodule pin and adjusts the bindings.
- **Generated files are build artifacts.** `src/include/docstr.hpp`
  (libclang-extracted C++ docstrings) and every `.pyi` stub under
  `src/pyvinecopulib/**/__init__.pyi` are gitignored. The build is the
  single source of truth — do not hand-edit, do not commit.
- **Code is quantitatively sensitive.** Pseudo-observation transforms,
  h-functions, Rosenblatt cascades, family parameterizations,
  pickling round-trips, and TLL grids all encode mathematically
  precise behavior. Small "obvious-looking" changes can silently
  break copula identities. Treat numerical paths as
  correctness-critical and prefer round-trip / parity tests over
  structural ones.

### Stability tiers

Different subpackages have different change policies. Honor the tier
when proposing API changes:

| Surface | Tier | Policy |
|---|---|---|
| `pyvinecopulib.core`, `pyvinecopulib.families`, `pyvinecopulib.utils`, top-level `pyvinecopulib` (core class re-exports) | **Stable-ish** | Solid user base. Prefer deprecation aliases over breaks; document migrations in `CHANGELOG.md`. PR #207 is the model: the reorg kept old import paths working via `_deprecations.py` + `DeprecationWarning`. Breaks are allowed (e.g. the pybind11→nanobind migration; the #207 cleanup) but must be intentional, documented, and worth the churn. |
| `pyvinecopulib.margins` | **Active development** | New in the vine-distribution work. The margin contract (`MarginLike` / `MarginBase` in `core`) is the part to treat as load-bearing; the curated parametric family registry, the selection criteria, and the report schema are all expected to move as they meet real data. |
| `pyvinecopulib.sklearn` | **Active development** | API may change in breaking ways between minor releases. The latest break is the `#218` public backend system (estimators now take a single `backend=` instead of loose `controls=`/`structure=`/`seed=` kwargs). |
| `pyvinecopulib.torch` | **Active development** | Same status. Defaults are still being tuned (cf. `990f997` device-aware `batched`, `cache_integrals=True`); the torch↔C++ cascade parity is a hard guarantee, but the `FitControlsTorchVinecop` surface and `TorchVinecop` method signatures may still shift. |
| `pyvinecopulib._python_helpers`, `pyvinecopulib._deprecations` | **Internal** | Underscore-prefixed. Not part of any contract; rename / restructure freely. `_deprecations.py` itself is slated for removal in 2.0. |

The "Solid user base" claim refers to the newest tag (see the
[GitHub project](https://github.com/vinecopulib/pyvinecopulib)).
Unreleased work on `main` is allowed to break sklearn/torch APIs as needed.

### Branching and releases

Pull requests go to `main`. Releases are tags on `main`; there is no
long-lived development branch. Read the Docs' `latest` follows `main` and
`stable` follows the newest tag.

- **Squash on merge: one pull request, one commit on `main`.** The
  iterations taken to reach a working state are review history, not project
  history, so the squashed message — not the intermediate commits — has to
  explain the change: what changed, why, and anything a future reader needs.
  The `(#NNN)` suffix keeps the pull request discoverable from `git log`,
  which is where `CHANGELOG.md` entries are sourced.
- **A pull request is a feature, not a commit.** Because squashing collapses
  the branch to a single commit, the unit of a pull request is the unit a future
  reader wants to find in `git log`: a change that stands on its own and earns
  roughly one changelog bullet. Several commits on the branch are expected and
  cost nothing. The failure mode to avoid is the opposite one — opening a pull
  request per commit, which turns one feature into a dozen entries on `main` and
  asks a reviewer to hold the whole chain in their head to judge any part of it.
  In particular, **the implementation sequence in a design document is not the
  pull-request list**: it is an order for writing the code, and several of its
  steps usually belong in the same pull request.
- **Commit subjects are `type(scope): subject`**, with `!` marking a
  breaking change. Scopes are the subpackages and areas: `core`, `families`,
  `utils`, `sklearn`, `torch`, `bicop`, `vinecop`, `build`, `ci`, `docs`,
  `deps`, `examples`.
- **Stack dependent work** rather than merging to unblock yourself: each
  pull request branches off the previous one and targets it. `gh stack`
  (`init` / `add` / `submit` / `sync` / `rebase`) manages the chain. Two
  consequences: CI's `pull_request` trigger must stay unfiltered, because a
  stacked child targets the branch below it; and nothing may push to a
  branch in a stack outside `gh stack` — including bots, which is why
  `regenerate_notebooks` is label-gated.
- **Never merge to `main` without express consent.** Open the pull request,
  get it green, and stop. A green matrix and an approved plan are not
  authorization. This applies equally to pushing tags and to changing
  repository or Read the Docs settings — and it matters more here than
  upstream, because **a `v*` tag push publishes to PyPI, and PyPI never
  re-accepts a version number**. A mistagged release is permanent.

Note for anyone reasoning about branch topology: a shallow clone makes
`git merge-base` report a divergence that does not exist. Run
`git fetch --unshallow` first.

### Changelog

`CHANGELOG.md` is newest-first. The top heading carries `(unreleased)`
while a cycle is open, is dated when the release ships, and a fresh
`(unreleased)` heading is opened immediately after tagging — so a released
version is never indistinguishable from an unreleased one.

Each change is **one bullet, one to three lines, four at the very most**:
imperative present, identifiers in backticks, no bold, and a trailing
`(#NNN)` naming the pull request (upstream work is cited as
`([vinecopulib#NNN](…))`). Anything that needs more room than that belongs
in the migration guide, not in a bullet.

Entries are sourced from `git log` — which is why the squashed commit
message has to stand on its own. Source them **by commit**, not by pull
request number: upstream's numbering has gaps where a pull request was
closed unmerged, and some commits carry no number at all.

## Scope

### Included

- **Bivariate copula modeling** — every family bound from
  `lib/vinecopulib`: `indep`, `gaussian`, `student`, `clayton`,
  `gumbel`, `frank`, `joe`, `bb1/6/7/8`, `tawn`, `tll`; with their
  rotations, mixed-discrete handling, analytic parameter/argument
  derivatives, tail-dependence / Blomqvist-β summaries, log-likelihood
  scores / gradient / Hessian / score-covariance, and family-set
  constraints via `FitControlsBicop`.
- **Vine copula modeling** — `Vinecop` with Dissmann selection
  (`mst_prim`, `mst_kruskal`), random spanning trees
  (`random_weighted`, `random_unweighted`), and user-supplied `RVineStructure` /
  `CVineStructure` / `DVineStructure`. Truncation, threading,
  bootstrap, pre-fit selection criteria, and family sets are exposed
  through `FitControlsVinecop`. Also: conditional sampling
  (`sample_conditional`), conditioning-aware selection
  (`FitControlsVinecop.conditioning_set`) and `reorient`, the
  list-of-trees round-trip (`get_trees` / `RVineStructure.from_trees`),
  the gradient/diagnostics surface (`scores` / `gradient` / `hessian` /
  `scores_cov`, with per-observation-parameter overloads).
- **Univariate marginals** — `Kde1d` (`lib/kde1d`) with continuous,
  ordered-discrete, and unordered-categorical support; plus the
  `pyvinecopulib.margins` layer on top of it: `Kde1d` *is* the
  nonparametric margin, `ParametricMargin` the parametric one,
  with AIC / BIC / AICc family selection over a curated candidate
  set (`MarginSelector`), and an adapter registry (`as_margin` /
  `register_margin_adapter`) that accepts a SciPy or PyTorch
  distribution object as a margin.
- **Vine distributions** — `Vinedist`: any `VinecopLike` combined with
  one margin per variable, giving `pdf` / `logpdf` / `cdf` / `loglik` /
  `sample` / `sample_conditional` / `rosenblatt` /
  `inverse_rosenblatt` on the **data** scale rather than the copula
  scale, for continuous, discrete and mixed margins alike.
- **Dependence measures** — `wdm` (`lib/wdm`).
- **Quasi-random sampling** — `sobol`, `ghalton`, `sample_uniform`.
- **Pseudo-observations** — `to_pseudo_obs`.
- **Estimator ensembling / model averaging.** Combining several
  fitted vines — bagging, averaging over candidate structures,
  post-hoc selection among them — is left to downstream packages.
  The library ships single-vine estimators plus the seams such a
  package needs: the copy-on-write `with_*` backend derivations, a
  pre-settable `schema_`, `VineRegressor.normalize_weights`, and the
  `_weights_for_batch` / `_predict_from_iter` split. Do not remove
  those as "unused" — they are deliberate, tested extension points.

- **Scikit-learn-compatible estimators** — `VineDensity`,
  `VineRegressor` with a pluggable backend (`VinecopBackend` /
  `TorchVinecopBackend`).
- **PyTorch evaluator** — `TorchBicop`, `TorchVinecop` (pure-torch
  cascade with GPU placement, autograd, and an optional `batched`
  evaluation fast path; byte-for-byte parity with the C++ cascade).
- **Backend-neutral extension layer** — the `BicopLike` / `VinecopLike`
  contracts and canonical `BicopBase` / `VinecopBase` bases (NumPy or
  PyTorch) for hosting custom pair copulas in a vine, including
  **non-simplified / conditional** vines via a `ConditioningContext`
  (walk-through: `examples/10_extending_pyvinecopulib.ipynb`).

### Excluded (explicit)

- **Copula families outside the bound set.** New parametric families
  belong upstream in `lib/vinecopulib`; bindings then follow.
- **Custom C++ forks.** The repo always tracks the upstream
  `lib/vinecopulib` submodule pin; local C++ patches under
  `lib/` are not accepted.
- **Discrete margins on the torch *marginal* layer.** `TorchVinecop` and
  `TorchVinecopBackend` handle discrete variables (the copula half), but
  `TorchMargin` rejects a discrete family: a margin with atoms needs a
  left-limit `cdf`, which `torch.distributions` does not expose. Use
  `Vinedist` with `pyvinecopulib.margins` for the marginal half.
- **Density estimators outside the vine framework.** General-purpose
  multivariate density models (normalizing flows, Gaussian mixtures,
  …) are not in scope; `pyvinecopulib` is a vine-copula library.
- **Pinned legacy alias for every old import path forever.**
  Deprecation aliases live in `_deprecations.py` and warn on access;
  they are removed in 2.0. They survive 1.0.0 deliberately — that release
  already breaks enough — but the reprieve is one cycle, not indefinite.

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

      core/__init__.py           # Bicop, Vinecop, *VineStructure, FitControls*, Kde1d (re-exports from ext)
        protocols.py             # BicopLike / VinecopLike backend-neutral contracts
        bicop_base.py            # BicopBase (canonical BicopLike partial impl)
        vinecop_base.py          # VinecopBase (array-agnostic cascades + fit/select)
        context.py               # ConditioningContext / Simplified / NonSimplified
        margin_base.py           # MarginBase (canonical MarginLike partial impl)
        vinedist.py              # Vinedist (copula + margins = a distribution)
        _discrete.py             # DiscretePair + the discrete layouts / per-edge types
        _reorient.py             # relabel a structure onto a chosen order tail (internal)
        _rootfind.py             # solve_increasing (monotone bisection; internal)
      families/__init__.py       # BicopFamily enum + 13 family constants + 15 group constants
      utils/__init__.py          # to_pseudo_obs, wdm, sobol, ghalton, sample_uniform, benchmark
        _pair_plots.py           # pairs_copula_data plotting helper (pure Python)

      margins/__init__.py        # as_margin, register_margin_adapter, resolve_margins, the margins
        parametric.py            # ParametricMargin (one SciPy family) — needs the [scipy] extra
        selection.py             # MarginSelector (fit every admissible candidate, keep the best)
        _openturns.py            # OpenTURNSMargin / OpenTURNSSelector — needs the [openturns] extra
        _adapters.py             # the coercion registry + per-ecosystem adapters — internal
        _resolve.py              # resolve_margins / fit_margin — internal

      sklearn/__init__.py        # VineDensity, VineRegressor, backends
        backends.py              # VinecopBackend / TorchVinecopBackend + resolve_backend
        _base.py                 # VineBase (parameter-constraints, schema, 3-step pipeline)
        density.py               # VineDensity
        regressor.py             # VineRegressor

      torch/__init__.py          # TorchBicop, TorchVinecop, TorchKde1d, TorchMargin, TorchVinedist, FitControlsTorch*
        bicop.py, vinecop.py     # nn.Module evaluators
        margin.py, vinedist.py   # TorchMargin / TorchVinedist (nn.Module margins + distribution)
        kde1d.py                 # TorchKde1d (the torch marginal estimator)
        _kde1d_interp.py         # kde1d's InterpolationGrid, ported — internal
        controls.py              # FitControlsTorchBicop / FitControlsTorchVinecop dataclasses
        _interp.py               # InterpolationGrid2D (bilinear; Sinkhorn margin renormalization) — internal
        _fit_tll.py              # pure-torch TLL kernel
        _batched.py              # batched evaluation variants

      _python_helpers/           # internal; pure-Python wrappers used by the binding
        bicop.py, vinecop.py, kde1d.py, stats.py
      pyvinecopulib_ext.*.so     # compiled extension (gitignored build artifact)
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

For performance work: profile first, optimize demonstrated hotspots
only, and preserve every quantitative invariant (round-trip identities,
parity with the C++ cascade, pickling stability).

### Which CI leg covers what

Local runs cannot substitute for the matrix, so know what each leg does
before deciding a change is verified. The 10 `build` legs run
`pytest tests/` under cibuildwheel with only `pytest-cov` and
`pytest-rerunfailures` installed, so every test guarded by
`importorskip` **skips** there. The 12 `install_and_unit_test` legs sync
`--extra examples --extra sklearn --extra torch`, install the wheel, and
run both `pytest tests/` and the notebooks — that is the matrix covering
the optional extras. `gh pr checks` and
`gh run view --job <id> --log-failed` read the result, though a run
triggered by `gh workflow run` attaches its checks to the commit without
`gh pr checks` listing them.

To check that a test survives without an optional extra, block the import
in-process (`sys.modules["scipy"] = None`) rather than reaching for a
matrix run.

## Working on this repo

### Inspection order

Before changing code, read in this order:

1. `AGENTS.md` (this file) — invariants and boundaries.
2. `docs/` — high-level intent, including the Sphinx `concepts.rst`
   primer on Sklar's theorem, pair copulas, R-vines, and TLL.
3. `src/pyvinecopulib/<subpackage>/__init__.py` — the module docstring
   is the canonical short description.
4. The implementation file you're about to touch, then the matching
   `tests/test_<topic>.py` for expected behavior.

Match existing local patterns rather than introducing new ones.

### Definition of done

For any behavior change:

- Diffs are scoped to the task; no opportunistic refactors that span
  unrelated files.
- Honor the [stability tier](#stability-tiers): for `core` /
  `families` / `utils`, prefer a deprecation alias over a hard break;
  for `sklearn` / `torch`, breaks are allowed but must be flagged in
  `CHANGELOG.md` (the top `(unreleased)` section).
- Tests added or extended. Prefer extending an existing parametrized
  test over duplicating logic; share fixtures via `conftest.py`.
- Public-API changes update the module docstring (re-rendered in
  `docs/features.rst` via autosummary) and the matching example
  notebook when one exists.
- Run the [validation sequence](#tooling), and for anything touching the
  packaging path, an sdist build and install — `build_sdist` is the only
  CI leg that runs `make check`.
- **A submodule bump additionally runs the numerics gate**, which differs by
  submodule. For `lib/vinecopulib` and `lib/wdm`: `tests/test_torch_bicop.py`,
  `tests/test_torch_vinecop.py` and `tests/test_structure_selection.py`. For
  `lib/kde1d`: `tests/test_kde1d.py`, `tests/test_torch_kde1d.py`,
  `tests/test_margins.py`, `tests/test_sklearn_margins.py` and
  `tests/test_plots.py`. They hold the torch↔C++ and NumPy↔C++ parity
  tolerances. When a number moves, **regenerate the expected value rather than
  widening the tolerance**, and check the direction of the change against what
  upstream says it fixed.
- **Fix what you find.** A defect uncovered along the way is fixed, not
  annotated, worked around, or left behind an explanatory comment. When the
  real fix belongs elsewhere — upstream, or a separate change — say so in
  the pull request description and open the issue; a comment is not a
  substitute.

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
- **American English** in code, comments, documentation, commit messages,
  and changelog entries: *behavior*, *normalize*, *serialize*, *finalize*,
  *center*, *modeling*, *honored*, *color*. There is no legacy exemption.
  `codespell` enforces it in `make lint` and in pre-commit, using its
  `en-GB_to_en-US` dictionary, so the rule covers every British spelling
  rather than a list this repository happened to drift on. It catches
  ordinary typos in the same pass. Configuration -- the skip list for
  generated and vendored files, and the domain vocabulary it would otherwise
  flag -- lives in `[tool.codespell]` in `pyproject.toml`; add a word there
  only when it is genuinely a term of art, never to silence a real
  misspelling.

- **Write for the caller, not the implementer.** A docstring says what a
  method does, what its arguments mean, and what it returns — never how it
  is computed. The algorithm it delegates to, which helper it calls, why a
  branch exists, and what it allocates are implementation details; they
  belong in the code, the commit message, or nowhere.

  This binds harder here than upstream:
  [scripts/generate_docstring.py](scripts/generate_docstring.py) lifts the
  C++ `//!` text **verbatim** into the Python docstring, so an
  implementation detail written upstream becomes user-facing Python API
  text on the rendered site. **When an upstream docstring reads as
  implementation detail, fix it upstream and bump the pin** — do not patch
  it in `src/include/**` and do not hand-write a replacement, except at
  the sites where libclang genuinely cannot disambiguate (they are
  enumerated in `src/include/vinecop/class.hpp`).
- **Doxygen upstream, numpydoc downstream.** The generator translates one
  into the other. Never "fix" an upstream `//!` comment into numpydoc.
- **Bind alternative constructors as named factories, not overloads.** C++
  overloads a constructor; Python names it. Every alternative way to build an
  object is a `def_static` — `Bicop.from_family` / `from_data` / `from_file` /
  `from_json`, `Vinecop.from_structure`, `RVineStructure.from_trees`,
  `FitControlsVinecop.from_bicop_controls` — and the binding surface carries no
  overloaded `__init__` at all. This is the Pythonic shape (each entry point is
  discoverable by name, with its own signature and docstring), and it is also
  the only shape the toolchain renders: nanobind concatenates the docstrings of
  overloaded bindings, so two `Parameters` sections collide and numpydoc fails
  the docs build, while the `.pyi` generator emits only the first signature.
  Where two forms genuinely are one operation, prefer a single method that
  dispatches on an optional argument — as `Bicop.pdf` does for per-row
  `parameters` — over two bound overloads.
- **Do not restate what a sibling documents.** A method that differs from a
  near-twin in one argument gets a short summary and only the text specific to
  it. Copying the twin's details, references, or edge cases is duplication that
  will drift.
- **Comments are documentation, not history.** Keep them aimed at whoever
  reads the code next: the constraint, the invariant, or the reason a
  non-obvious choice is required. Previous bugs, benchmark numbers, and
  review discussion belong in the commit message. The test: *would this
  comment still make sense in a file that had never had the bug?* If it
  only reads as a contrast with what the code used to do, it is history.
  This applies equally to CMake, CI workflows, and `pyproject.toml`.
- **Do not suppress diagnostics.** `# noqa`, `# type: ignore`,
  `nitpick_ignore_regex`, pytest `filterwarnings`, and bandit skips are
  narrow and justified in a comment, or the underlying code gets fixed.
  Every `nitpick_ignore_regex` entry names the upstream cause that will
  retire it.
- **Every code example is executed.** An example in prose is either a cell
  in an `examples/*.ipynb` notebook (run by `pytest --nbmake` in CI) or a
  doctest — never a paste into `docs/*.rst` or a docstring. doctest is not
  currently wired up, so prose pages link a notebook cell rather than
  inlining code.
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
Behavior and API changes belong upstream. The Python repo only:

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
  `Vinecop.from_data(...)` (or `RVineStructure.sample(...)`) — these
  factories are the documented entry points; the raw constructor
  signatures are kept for nanobind-level access only.
- `tree_algorithm` on `FitControlsVinecop`: `"mst_prim"` (default,
  Dissmann), `"mst_kruskal"`, `"random_weighted"` (Wilson-weighted random
  tree; reachable from the sklearn layer via
  `VinecopBackend.with_local_random`) and `"random_unweighted"`.
- `FitControlsVinecop.conditioning_set` (property + pickled, not a
  positional ctor arg) drives conditioning-aware selection — the fitted
  order ends with the given 1-based variables so
  `Vinecop.sample_conditional` / `reorient` can condition on them.
- `RVineStructure.from_trees(d, trees)` is the **faithful** inverse of
  `RVineStructure.get_trees()` (identity diagonal policy — each edge's
  `conditioned[0]` on the diagonal — so `from_trees(s.dim, s.get_trees()) == s`).
  Upstream `Vinecop.select` finalizes with the *same* (flip-free)
  convention, so `VinecopBase.select` assembles its selected trees through this
  same `from_trees` and matches the compiled selector's matrix byte-for-byte —
  one diagonal convention throughout.
- **Backend-neutral abstraction layer** (pure Python; `core` imports
  without PyTorch). The extension point for custom (e.g. neural,
  conditional) pair copulas and vines:
  - `BicopLike[ArrayT]` / `VinecopLike[ArrayT]` (`protocols.py`) —
    generic, `runtime_checkable` protocols mirroring the `Bicop` /
    `Vinecop` evaluation surface on any array backend (NumPy or
    PyTorch). `Bicop` / `Vinecop` satisfy them *nominally* — `isinstance` is
    `True`, because a `runtime_checkable` Protocol compares method names and
    nothing else. It is not a statement that the signatures agree: the
    compiled classes take per-row `parameters` where the protocols take a
    conditioning matrix `x`, which is why `x` is keyword-only on both.
  - `BicopBase` (`bicop_base.py`) / `VinecopBase` (`vinecop_base.py`) —
    canonical partial implementations to subclass. A `BicopBase`
    subclass defines `pdf` / `hfunc1` / `hfunc2` and inherits `hinv1` /
    `hinv2` (bisection), `sample`, `loglik`, `plot` (`flip` — needed
    only to host the pair in structure *selection* — defaults to
    raising); a `VinecopBase` subclass defines the one hook
    `_get_pair_copula` and inherits the whole tree-by-tree cascade plus
    the public `fit` and `select` engines. `select` is an
    exact port of `Vinecop`'s Dissmann / Wilson structure selection
    (same matrix encoding, selection-time pairs reused via `flip`, no
    re-fit; parity is a hard guarantee). `TorchBicop` / `TorchVinecop`
    are the torch subclasses.
  - `DiscretePair` (`_discrete.py`) — a *continuous* pair copula evaluated on a
    discrete or mixed edge. **The vine owns the discrete layouts, the pair
    copulas stay continuous**: `_bind_vine(..., var_types=)` declares which
    variables have atoms, `pair_var_types(tree, edge)` derives the types each
    slot sees from the structure alone, and the cascades hand a four-column
    `[u1, u2, u1^-, u2^-]` argument to any pair whose types include `"d"`.
    `BicopLike` is therefore unchanged — it stays a two-column continuous
    contract — and a custom pair copula opts in by implementing `cdf` and
    wrapping itself in `DiscretePair`. `fit` / `select` take `var_types` too and
    forward each edge's types to `fit_edge` as a keyword, only on the edges that
    have one (the rule `_pair_eval` applies to `x`). The parity test that binds
    is the **normalization identity** `Σ_atoms c(u₁,u₂)·(u₁ − u₁⁻) = 1`: the
    quotients telescope, so it holds exactly and needs no reference
    implementation and no tolerance argument. It is what established that
    `DiscretePair` was right and the compiled `tll` pair was wrong
    (fixed upstream in vinecopulib#739 and pinned since). Parametrize
    pair-level parity over **every** family, not a representative couple —
    covering only `gaussian` and `clayton` is why that class of defect stayed
    invisible on both sides for as long as it did. Note the identity **cannot**
    catch a cache regression: it telescopes to the four corners, so it reads
    `1 − 2e-10` for a correct density and for a 38%-wrong one alike.
    A rectangle's probability is read by differencing four `cdf` values, which
    is what the compiled pair does, so a `DiscretePair` is bit-identical to it.
    `TorchBicop.rect_mass` would be more accurate — 1.2e-15 against 9.2e-15 at
    a `1/8`-wide atom, measured against exact rational truth, and far more at
    narrower ones — but it is **deliberately not used**: the density divides by
    the atom's area, and the discrete cascade then turns a 1e-15 pair-level
    difference into `8.5e-8` at the vine, a visible divergence from
    `Vinecop`. The torch↔C++ cascade parity is a documented guarantee, and it
    outranks the accuracy here; revisit only together.
  - `sample_conditional` / `reorient` (`_reorient.py`) — conditional sampling and
    the value-preserving relabeling it rests on. A **truncated** model relabels
    like any other: the trees above the truncation are independence, so the peel
    has nothing to move there and the slot map covers only `trunc_lvl` trees. At
    `trunc_lvl == 0` every set is admissible and the relabeled order is a stable
    partition. `reorient` **returns** the
    relabeled `(structure, pair_copulas)` rather than mutating, since the base
    class leaves pair storage to the subclass; `conditioning_set` on
    `rosenblatt` / `inverse_rosenblatt` / `sample_conditional` evaluates through
    an internal reoriented view instead. The peel that steers a chosen set to
    the order tail is borrowed from the compiled `Vinecop.reorient` (run on a
    throwaway independence vine) and the slot map is then matched up in Python,
    so admissibility and the error messages are exactly `Vinecop`'s — the same
    trade `select` makes with `_select_spanning_tree`. A relabeling is refused
    on a non-simplified vine: it can permute the columns of each edge's `x_e`,
    which makes the result a different model. `select` takes `conditioning_set`
    too (the `+d` MST penalty, then the relabeling).
  - `ConditioningContext` / `SimplifiedContext` (default) /
    `NonSimplifiedContext` (`context.py`) — the per-edge policy that
    turns the simplified cascade into a **non-simplified / conditional**
    vine (each pair also sees its conditioning-set values `u_D` and any
    external covariates `x`). Walk-through:
    `examples/10_extending_pyvinecopulib.ipynb`.
  - `solve_increasing` (`_rootfind.py`) — vectorized monotone bisection
    behind the default `hinv1` / `hinv2` (internal; not re-exported).
    Brackets may be array-valued and unbounded, since it also backs
    `MarginBase.icdf` on an infinite support.
- **The marginal layer.** `MarginLike[ArrayT]` (`protocols.py`) is
  `{pdf, cdf, icdf}` and declares no attributes, the same discipline as
  `BicopLike`. `pdf` means *the density with respect to the margin's own
  reference measure* — a Lebesgue density for a continuous margin, a
  probability mass at an atom — which is what makes
  `log f(x) = log c(u) + Σ_j log pdf_j(x_j)` hold verbatim for
  continuous, discrete and mixed margins with no branch in the
  likelihood path. `MarginBase` (`margin_base.py`) needs only `pdf` /
  `cdf` and supplies `icdf` (bisection), `logpdf`, `cdf_left`, `loglik`,
  `sample`, `var_type`, `support`, `is_fitted` and a raising `fit`.
  Everything beyond `{pdf, cdf, icdf}` is an **optional capability**
  read with `getattr` (`var_type` ∈ `{"c","d","zi"}`, `cdf_left`,
  `logpdf`, `sample`, `support`, `supports_covariates`), per the house
  precedent — each
  member added to the protocol is one a foreign object must happen to
  have.
- **`Vinedist`** (`vinedist.py`) composes a `VinecopLike` with one
  margin per variable. It owns the copula-scale layout: `copula_data`
  builds the compact `(n, d + k)` matrix from `cdf` and `cdf_left`,
  clamps once, and checks `cdf_left <= cdf`, so callers never `hstack`
  a left-limit block by hand. `logpdf` sums logs rather than
  accumulating a product, because the marginal term carries the scale
  and a `d = 50` product underflows. `sample` is
  `marginal_icdf(copula.sample(n, ...))`, so it inherits the copula's
  quasi-random and seeding options and never calls a margin's own
  sampler. `sample_conditional` is the same sandwich around the
  copula's, with one difference the caller feels: a discrete
  conditioner needs no left-limit column, since it is derived from that
  variable's own margin. Both scales resolve an omitted
  `conditioning_set` through one `infer_conditioning_set`, so the
  column-to-variable rule cannot drift between them. Every method also takes optional exogenous covariates `x`,
  forwarded to each margin that declares `supports_covariates` and to a
  copula that declares it too. Two seams keep the array namespace
  coherent: `_prep` (identity here, `torch.as_tensor` on
  `TorchVinedist`) coerces one input array onto the parts' namespace, so a
  caller may hand the type they have; and `copula_data` /
  `_conditioning_data` / `marginal_cdf` / `marginal_icdf` take `xp` from
  the *columns the margins returned*, never from the input — a margin may
  legitimately answer in another array type, and stacking that through the
  input's namespace either raises or silently detaches.

### `pyvinecopulib.margins`

The univariate half of `Vinedist`, kept out of `core` so that `core`
imports without SciPy. Three groups:

- **Built-in margins** — `Kde1d` *is* the default margin, needing no
  wrapper; it takes `xmin` / `xmax` so a bounded variable is not fitted
  past its support — the sklearn estimators fill those in from a
  categorical's declared levels, since otherwise the density grid is
  padded past the data. What a bound *means* differs by variable type, and
  `docs/concepts.rst` (`concepts-kde-margins`) is where that is written down:
  for a discrete variable it is the **integer support**, so the bound and the
  data must both be integers and the fitted grid runs half a unit wider at each
  end. A categorical whose levels are not integers is therefore refused, by
  name, at fit time rather than by a bare `invalid_argument` from C++. Then `ParametricMargin` (one SciPy family) and
  `MarginSelector` (fit every admissible candidate, keep the best on
  `selected_`, report the rest on `report_` — the `GridSearchCV` shape).
- **Coercion** — `as_margin(obj)` is idempotent and routes **every**
  margin `Vinedist` receives, so a discrete SciPy object cannot slip
  past on a bare `pdf` (in SciPy's new API `pdf` is `+∞` at an atom;
  the mass is `pmf`). `register_margin_adapter(predicate, adapter)` is
  how another ecosystem is added without touching this package.
- **Resolution** — `resolve_margins(spec, ...)` mirrors
  `resolve_backend`: a string alias, one instance broadcast per column,
  a length-`d` sequence, or a dict keyed by column. Margins follow the
  library's own **construct-then-`fit`** pattern (`fit` returns `self`),
  so one class is both the specification and the fitted object, and a
  spec may freely mix already-fitted margins with unfitted ones —
  `from_data` fits only the latter.

Conventions that bind: the fit is **two-step (IFM)** — margins first,
then the copula on the resulting pseudo-observations — never fit all of
SciPy (a blind sweep ranks `vonmises` above the true `gamma` because its
reported support lies), and never silently skip a failed candidate: every
rejection gets a row in `report_` with a reason, and a column where
everything fails **raises**, naming each family and its cause.
`on_failure="fallback"` substitutes `Kde1d` with one warning instead --
available, but not the default, because answering a parametric request
nonparametrically is the same class of silent downgrade the weights
contract already refuses.

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

- Re-exports `Kde1d`, `to_pseudo_obs`, `wdm`, `find_latent_sample`,
  `sobol`, `ghalton`, `sample_uniform`, `benchmark` (all C++) plus the
  pure-Python `pairs_copula_data` helper from `_pair_plots.py`.
- `wdm`'s `method` includes Chatterjee's ξ (`"chatterjee"` / `"cxi"` /
  `"xi"`), the one **asymmetric** measure in the list — it measures how far
  `y` is a function of `x`. `FitControlsVinecop.tree_criterion` accepts it
  too, spelled **`"cxi"` only** (the full accepted set is `tau`, `rho`,
  `hoeffd`, `mcor`, `cxi`, `joe`, `custom`), so any Python-side selector that
  computes the criterion itself must accept it *and* symmetrize it the way
  `pairwise_cxi` does — `max(ξ₁₂, ξ₂₁)` — or silently diverge from
  `Vinecop.select`.
- ξ breaks **predictor ties** at random, since ordering them by the response
  would manufacture dependence. The seeds default to a constant, so ξ is a
  function of its arguments; `wdm(..., seeds=)` varies that ordering for a
  caller who wants to average over it. Untied predictors never construct the
  generator, so continuous data is unaffected.
- `wdm` **raises** on weights whose sum is not finite and positive, rather
  than returning `NaN`.
- The per-entry structure accessors (`struct_array`, `min_array`,
  `needed_hfunc1` / `needed_hfunc2`) are wrapped in a bounds check: upstream
  indexes the triangular array without one, so reading a tree above
  `trunc_lvl` **segfaulted**. The real fix belongs upstream; the guard is here
  because a crash takes the interpreter down.
- `find_latent_sample(u, b, niter=3)` recovers a continuous sample from
  interval-censored copula data — the transform a nonparametric fit on
  discrete margins runs on. The draw is deterministic and invariant to
  argument order, so a pair reused with its arguments flipped recovers the
  same latent sample.
- `Kde1d` is used internally by the sklearn estimators as the
  marginal estimator; it also stands alone for any 1-d KDE problem.
- `to_pseudo_obs(data)` is the canonical input transform for
  copula fitting (rank-normalize to the unit hypercube).

### `pyvinecopulib.sklearn`

User-facing estimators on top of the core, organized so that the
3-step pipeline — fit 1-d KDE marginals → transform to
pseudo-observations → fit a vine on the copula data — happens once,
in `VineBase` (`_base.py`).

Class hierarchy:

```text
sklearn.base.BaseEstimator
├── VineBase (_base.py)            # shared pipeline + DataFrame schema
│   ├── VineDensity   (+ DensityMixin)
│   └── VineRegressor (+ RegressorMixin)
```

#### `pyvinecopulib.sklearn.backends`

Public extension point introduced in `#218`. Estimators do not call
`pyvinecopulib.Vinecop.from_data(...)` directly — they go through a
backend object. Two concrete backends ship, both subclassing the
private `_VinecopBackendBase` (which owns `structure_of` and the
copy-on-write `with_*` derivations; concrete backends override only the
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
backend.sample(vine, n, seeds=...) -> np.ndarray
backend.structure_of(vine)  -> RVineStructure
backend.default_margin(var_type, bounds) -> MarginLike
backend.bind_distribution(vine, margins) -> Vinedist
backend.with_random_structure(d, seeds)  -> Backend  # copy-on-write
backend.with_local_random(seeds)         -> Backend  # ditto
backend.with_num_threads(n)              -> Backend  # ditto (torch: no-op)
```

`default_margin` and `bind_distribution` are what keep the *whole*
distribution on one array namespace: the backend names the margin class an
estimator fits when the caller named none, and it assembles the fitted parts.
The default backend wraps its copula in `_BackendVinecop` so `distribution_`
evaluates with the same threading and batching arguments the estimator uses;
`TorchVinecopBackend` fits `TorchKde1d` margins (carrying `device` / `dtype`
from its controls) and publishes a **`TorchVinedist` holding the raw
`TorchVinecop`**, so `.to(device)`, `state_dict` and autograd reach the whole
object. A compiled `Kde1d` the caller supplied is lifted with
`TorchKde1d.from_kde1d` rather than refused, which is what makes
`margins="kde"` behave like `margins=None` there. Two consequences: on the
torch backend nothing passes `batched` at all -- every call lets the vine
resolve it per device, which is what makes the sklearn path take the stacked
cascade on CUDA like every other caller -- and every public
estimator method converts back to NumPy through `_base._as_ndarray` —
`np.asarray` alone raises on a tensor that requires grad or lives on an
accelerator.

`pyvinecopulib.core.VinecopLike` is the canonical `runtime_checkable`
Protocol describing the post-fit vine surface (`pdf` / `cdf` /
`rosenblatt` / `inverse_rosenblatt` / `sample`, plus a `structure`
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
- `random_state` is the canonical RNG kwarg name throughout (no
  legacy `seed=` / `seeds=` on `VineDensity.sample`,
  `VineDensity.cdf`, or the estimator constructors).

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
    `rosenblatt` / `inverse_rosenblatt` / `sample` signatures.
- `TorchMargin` / `TorchVinedist` — the marginal and joint halves.
  - `TorchMargin` is a `MarginBase[Tensor]` that is *also* an
    `nn.Module`: `torch.distributions.Distribution` has no
    `.to(device)` and contributes nothing to `state_dict` as a plain
    attribute, so the parameters are registered and the distribution is
    **rebuilt per call** — the same shape `TorchBicop` uses for its
    grid. `TorchMargin.from_distribution(factory, parameters=...)` is
    the general entry point; `icdf` bisects `cdf` over `support` for the
    families that implement one but not the other (`Gamma`, `Chi2`).
  - `TorchVinedist` is `Vinedist[Tensor]` plus `nn.Module`, with margins
    in a `ModuleList` and `log_prob` as an alias for `logpdf`. Every
    margin **must** be an `nn.Module`: SciPy raises on gradient-carrying
    tensors and returns a plain `ndarray` without them, so accepting a
    SciPy margin here would detach the graph silently. `from_data` fits
    end to end in torch — a `TorchKde1d` per column, then
    `TorchVinecop.from_data` on the copula data they produce — so the whole
    joint distribution lands on one device in one dtype. It takes no
    covariates: no torch margin reads them, and an unconditional fit
    behind a conditional-looking call is worse than a refusal.
  - `TorchKde1d` is the torch marginal estimator, and the only one that
    handles **discrete** and **zero-inflated** variables. Fitting delegates
    to the compiled `Kde1d`; every evaluation is pure torch, because
    `grid_points` / `values` / `type` / `prob0` are the only state the
    compiled `pdf` / `cdf` / `icdf` read. The grid is a **buffer**, not a
    parameter — the density is fitted, not learned — so optimizing it is the
    opt-in `values.requires_grad_(True)`, the `TorchBicop` precedent.
    `icdf` reproduces the C++ inversion exactly — a bracketed Newton within
    the cell holding the requested mass, bisecting where the density is flat,
    with the C++ early exit reproduced as a frozen-once-converged mask — and
    then reattaches an exact gradient by one Newton step, so the correction
    does not move the value while `dq/dp` and `dq/d values` are right. The
    residual is written in units of mass, so the total mass carries its share
    of `dq/d values`. This is the one parity claim in the port that is a
    tolerance rather than an equality, and deliberately so: the compiled
    quantile is not portable to a few ULPs -- rebuilding kde1d with
    `-march=native` alone moves it 19 -- so no port can equal every build of
    it. The correction is gated on
    grad being enabled, not on the grid being learned — a fitted fixed grid
    still has to differentiate the quantile in `p`. Two of `Kde1d`'s attribute names could not
    be reused: `type` is `nn.Module`'s dtype cast (read `kde_type`) and
    `loglik` is the contract's method.
  - `torch/_kde1d_interp.py` ports `kde1d`'s `InterpolationGrid`, and its
    contract is *fidelity*, not improvement. One C++ behavior looks like a bug
    and must be reproduced: `integrate` adds no tail contribution, so the
    unnormalized integral saturates at the grid's mass even though the density
    beyond it is nonzero. Two others used to be on this list — the interpolant
    dropping by `exp(-0.5)` at `grid_points[-1]`, and `pdf_discrete` dividing
    by the raw interpolation — and kde1d fixed both, at which point
    reproducing them *was* the defect. That is the failure mode this contract
    invites: check the list against upstream on every bump, because a quirk
    that has been fixed reads exactly like a quirk that has not. Nothing is cached: coefficients and cell
    integrals are recomputed in the graph, which is the opposite call from
    `TorchBicop.cache_integrals` and for a stated reason — here the cached
    quantity would be an `O(m)` vector shared by the batch, not an `O(m^2)`
    integral per query.
- Discrete variables are declared with `var_types` on `TorchVinecop`'s three
  constructors. The stored pair copulas stay continuous interpolation grids and
  `_get_pair_copula` wraps a discrete edge in `DiscretePair`, so `state_dict` /
  `.to()` / pickling see only real `nn.Module` parameters. `TorchBicop.from_data`
  takes the four-column layout and reuses the compiled `find_latent_sample`,
  which is what `TllBicop::fit` now consumes for a discrete edge; the jittered
  ranks only seed the bandwidth. Two things a discrete torch vine refuses: the
  **integral cache** (`cache_integrals=None` resolves to `False`, an explicit
  `True` raises), because differencing a bilinearly interpolated `cdf` gives 38%
  error on a `("d","d")` density; and the **batched fast path**, whose stacked
  per-level grids carry no distribution function at all.
- `FitControlsTorchBicop` / `FitControlsTorchVinecop` — fit-time
  dataclasses. Notable knobs:
  - `method` — `"tll"` (the only fitter; kept as the dispatch seam
    for future torch fitters).
  - `cache_integrals` — default `True`; precomputes three `(m, m)`
    cumulative-trapezoid **prefix** tables (`sy`, `sx`, `p`), from which
    `cdf` / `hfunc*` read their value in closed form. The reconstruction is
    **exact**, not an approximation: `chat` is bilinear, so along a grid line
    it is piecewise linear and its integral is piecewise linear across cells.
    So the cache costs nothing in accuracy — it agrees with the on-the-fly
    path to summation-order noise — and it carries an exact gradient in
    `values` as well as in `u`. `hinv*` get less from them than `hfunc*` do but
    not nothing: there is no `O(1)` lookup, because locating the bracketing
    cell needs the conditional cumulative along the whole free axis, but that
    cumulative is exactly what a prefix table holds. Integration is linear, so
    blending two of its lines is the same quantity as integrating the blended
    knots -- a gather instead of a trapezoid and a scan, agreeing to 2e-16. So
    the two cache modes run the same closed-form inversion on a cumulative
    they reach differently, and agree to rounding rather than exactly.
    The tables are buffers, so `_tables` rebuilds them in-graph when `values`
    starts tracking grad after construction.
  - `rect_mass(a1, b1, a2, b2)` is available in **both** cache modes: the exact
    probability of a rectangle — the value a four-corner `cdf` difference
    defines, arranged so that almost none of it cancels. That difference turns
    an absolute error `ε` into `≈4ε/(w₁w₂)` in the atom widths; `rect_mass`
    amplifies by `1/w₂` alone, since only its `λ(b₂) − λ(a₂)` term cancels and
    that multiplies a term of order `w₁`. Measured on a `1.2e-4`-wide
    rectangle: `2.9e-12` against `8.7e-9`. `values >= 0` is a constructor
    precondition precisely because the nonnegative-weight bound depends on it.
    Note it is the **probability**, not the density's mass: `cdf` renormalizes
    each grid line by its own total, so the two differ, and a discrete edge is
    defined against the distribution function.
  - `compile` — runs the batched cascades through `torch.compile`, on CUDA
    replayed as a CUDA graph. Off by default: the first call at each input
    shape pays tens of seconds of Inductor, so it is worth it for a cascade
    called repeatedly and not for a single evaluation. Torch caps how many
    variants of one code object it will hold (`cache_size_limit`, 8), and a
    vine is a variant — so a process compiling more vines than that falls
    back to eager. The `batched` flag is **not** a control: it is resolved
    per device on each call, and overridable per call.
  - `device`, `dtype` — propagate to every tensor on construction;
    fitted modules respect `.to(device)` afterwards.
  - `trunc_lvl`, `tree_criterion`, `threshold`, `tree_algorithm`, `seeds`
    — `FitControlsTorchVinecop` only; the structure-selection knobs used
    when `TorchVinecop.from_data` is called with `structure=None`.
    Selection runs through `VinecopBase.select`, so it stays on the array
    namespace rather than round-tripping a compiled `pv.Vinecop`.
- `InterpolationGrid2D` (`torch/_interp.py`) — the 2-d bilinear grid
  backing `TorchBicop`; **internal** (not re-exported). Margin
  normalization uses Sinkhorn iterations to drive marginals to uniform.

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
  `from pyvinecopulib.core import Kde1d`), not the top-level alias.

`torch` is **not** re-exported at the top level — `import
pyvinecopulib.torch` is the only entry. Same for `sklearn.backends`.

### Internal: `_python_helpers`, `_deprecations`

- `_python_helpers/{bicop,vinecop,kde1d,stats}.py` — pure-Python
  helpers (DataFrame conversions, plotting glue) called by the
  nanobind extension via `nb::module_::def(...)` lookup. Not part of
  any public contract. Move new internal helpers here rather than
  exposing them.
- `_deprecations.py` — `_DEPRECATED_TOP_LEVEL` dict + `_resolve_deprecated`
  helper for the top-level `__getattr__` shim. Slated for removal in 2.0;
  new deprecation aliases can be added here in the meantime, but each entry
  is a debt to be paid down.

## Public APIs

The canonical, exhaustive API listing is `docs/features.rst`, which is
generated at Sphinx-build time via `autosummary` from the
`DOCSTRING_SUBPACKAGES` table in `docs/conf.py`. The Sphinx HTML build
is the source of truth for signatures and parameter docs; the listings
below are a quick orientation.

- **`pyvinecopulib.core`** — `Bicop`, `Vinecop`, `RVineStructure`,
  `CVineStructure`, `DVineStructure`, `FitControlsBicop`,
  `FitControlsVinecop`; plus the backend-neutral abstraction layer
  `BicopLike`, `VinecopLike`, `BicopBase`, `VinecopBase`, `DiscretePair`,
  `ConditioningContext`, `SimplifiedContext`, `NonSimplifiedContext`;
  plus the marginal layer `MarginLike`, `MarginBase` and the joint
  object `Vinedist`.
- **`pyvinecopulib.families`** — `BicopFamily` enum; per-family
  constants (`indep`, `gaussian`, `student`, `clayton`, `gumbel`,
  `frank`, `joe`, `bb1`, `bb6`, `bb7`, `bb8`, `tawn`, `tll`); group
  lists (`all`, `parametric`, `nonparametric`, `one_par`, `two_par`,
  `three_par`, `elliptical`, `archimedean`, `extreme_value`, `bb`,
  `rotationless`, `lt`, `ut`, `itau`, `analytic_derivs`).
- **`pyvinecopulib.utils`** — `to_pseudo_obs`, `wdm`,
  `find_latent_sample`, `sobol`, `ghalton`, `sample_uniform`,
  `benchmark`, `pairs_copula_data`.
- **`pyvinecopulib.margins`** — `ParametricMargin`, `MarginSelector`,
  `OpenTURNSMargin`, `OpenTURNSSelector`, `as_margin`,
  `register_margin_adapter`, `resolve_margins`.
- **`pyvinecopulib.sklearn`** — `VineDensity`, `VineRegressor`,
 plus the `backends`
  submodule (`VinecopBackend`, `TorchVinecopBackend`,
  `resolve_backend`).
- **`pyvinecopulib.torch`** — `TorchBicop`, `TorchVinecop`, `TorchKde1d`,
  `TorchMargin`, `TorchVinedist`, `FitControlsTorchBicop`,
  `FitControlsTorchVinecop`.

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
- `unique_json_path` — per-test path in `tmp_path` for serialization
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
  re-executes `examples/*.ipynb` via `pytest --nbmake`. CI runs this on
  every pull request. The `regenerate_notebooks` job refreshes the stored
  outputs, and fires only on the `regenerate-notebooks` label: it pushes a
  commit onto the head branch, which would rewrite the base of every pull
  request stacked above it.

Round-trip / parity properties to preserve when touching numerics:

- Pickling round-trip: `Bicop`, `Vinecop`, `Kde1d`, the sklearn
  estimators, and the torch evaluators all round-trip through
  `pickle.dumps` / `pickle.loads`.
- C++ ↔ torch parity: `TorchVinecop.from_vinecop(cpp_vine)` yields a
  module whose `pdf` / `cdf` / `rosenblatt` outputs match the C++
  vine to within floating-point tolerance, on the operational grid.
- `batched=True` ↔ `batched=False`: numerically equivalent `pdf` /
  `rosenblatt` on the same fitted vine, and **bit-identical**
  `inverse_rosenblatt` -- its waves reorder the cells without changing what
  any one of them computes, so that one is pinned at `atol=rtol=0`.
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
  an immutable / functional backend. To put that pair on a **discrete**
  edge, add a `cdf` and return `DiscretePair(pair, self.pair_var_types(t, e))`
  from `_get_pair_copula` (and from `fit_edge`, which receives the edge's
  `var_types`); the vine supplies the left-limit columns. For a
  **non-simplified / conditional** vine, pass a `NonSimplifiedContext` and drive
  `VinecopBase.fit` with a `fit_edge` callback. To condition on a subset of
  variables, implement `flip` as well and use `sample_conditional` /
  `select(conditioning_set=)`; see
  `examples/10_extending_pyvinecopulib.ipynb`. `TorchBicop` /
  `TorchVinecop` are the reference torch subclasses.
- **Custom margins (`pyvinecopulib.core`).** Subclass `MarginBase` and
  define `pdf` / `cdf`; `icdf`, `logpdf`, `cdf_left`, `loglik`,
  `sample` and `support` come with it. Add `fit(y, weights=None) ->
  Self` to make it an estimator, or leave it out for a fixed margin —
  `is_fitted` is what `resolve_margins` dispatches on. Override
  `cdf_left` whenever the family has an exact left limit (Poisson's
  `gammaincc(k, μ)`, a categorical's `cumsum(probs)[k-1]`): the derived
  `cdf(x) - pdf(x)` cancels in the right tail and `cdf(x - 1)` is
  meaningless off an integer lattice.
- **Another ecosystem's distributions (`pyvinecopulib.margins`).** Call
  `register_margin_adapter(predicate, adapter)` — from a package or a
  notebook cell — rather than editing `_adapters.py`. That is what keeps
  OpenTURNS, TFP and NumPyro out of `core` while remaining usable, and
  what lets `as_margin` stay the single funnel every margin passes
  through.
- **New `pyvinecopulib.sklearn` backends.** Subclass the private
  `_VinecopBackendBase`, which already provides `structure_of`,
  `default_margin`, `bind_distribution` and the
  copy-on-write `with_*` derivations (`with_random_structure` /
  `with_local_random` / `with_num_threads`); override the divergent
  members (`fit_vine`, `pdf`, `cdf`, `sample`, `_default_controls`,
  and `_default_controls`, which `_effective_controls` resolves
  lazily). `resolve_backend`
  accepts any such object (it only defaults `None`). Consider whether
  the underlying vine satisfies the `pyvinecopulib.core.VinecopLike`
  Protocol so downstream code can stay type-stable.
- **New torch fit methods.** Add the implementation under
  `pyvinecopulib/torch/_fit_<name>.py`, add `"<name>"` to `METHODS`
  in `controls.py`, and add the dispatch branch in
  `TorchBicop.from_data` (the lone `"tll"` path is the template).
- **New sklearn-style estimators.** Subclass `VineBase` and add the
  mixin that matches the task (`DensityMixin` / `RegressorMixin`),
  reusing the 3-step pipeline (`_validate_input` / `_fit_marginals` /
  `_to_u_scale` / `_fit_vine`). Stay inside the
  `_parameter_constraints` / `_validate_params()` pattern and pin
  fitted attributes with trailing underscores.
- **Wider quasi-random integration.** `sample_uniform` is the
  high-level driver behind `Vinecop.cdf` Monte-Carlo evaluation and
  random structure sampling; new low-discrepancy methods slot in
  alongside `sobol` / `ghalton` in `lib/vinecopulib`'s `misc/stats.hpp`.
