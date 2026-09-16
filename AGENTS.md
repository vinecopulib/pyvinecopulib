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
   additionally ships an array-agnostic abstraction layer
   (`BicopLike` / `VinecopLike` / `MarginLike` / `VinedistLike` protocols,
   `BicopBase` / `VinecopBase` / `MarginBase` / `VinedistBase` canonical
   bases, `ControlsLike` for fit configuration, `ConditioningContext`
   policies) that custom NumPy / PyTorch parts subclass, and `Vinedist`
   — a vine copula combined with univariate margins, i.e. a full
   multivariate distribution on the data scale.
2. `pyvinecopulib.margins` — the univariate marginal layer `Vinedist`
   composes: the built-in margins, family selection, and an adapter
   registry that presents a SciPy / PyTorch / other-ecosystem
   distribution object as a margin.
3. `pyvinecopulib.sklearn` — scikit-learn-compatible estimators
   (`VineDensity`, `VineRegressor`) on top of the core, fitting whichever
   `VinedistBase` subclass they are handed.
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
| The four contracts and their bases in `core` (`BicopLike` / `BicopBase`, `VinecopLike` / `VinecopBase`, `MarginLike` / `MarginBase`, `VinedistLike` / `VinedistBase`) | **Stable from 1.0.0** | `README.md` tells users to subclass these "with the same confidence as on `Vinecop`", so they are covered by the same policy as the rest of `core` from the 1.0.0 tag onward. New in 1.0.0, which is why the argument order and the optional-capability split were settled *before* it shipped rather than after. A protocol may still gain an optional capability -- that widens it -- but not a required member. |
| `pyvinecopulib.margins` | **Active development** | New in the vine-distribution work. The margin *contract* is stable (see the row above, where it belongs); the curated parametric family registry, the selection criteria, and the report schema are all expected to move as they meet real data. |
| `pyvinecopulib.sklearn` | **Active development** | API may change in breaking ways between minor releases. The latest break is `#339`, which deleted the `#218` backend layer: estimators now take `distribution=` / `controls=` / `structure=` directly. |
| `pyvinecopulib.torch` | **Active development** | Same status. Defaults are still being tuned (cf. `990f997` device-aware `batched`, `cache_integrals=True`); the torch↔C++ cascade parity is a hard guarantee, but the `FitControlsTorchVinecop` surface and `TorchVinecop` method signatures may still shift. |
| Every underscore-prefixed module, and `pyvinecopulib._deprecations` | **Internal** | Not part of any contract; rename / restructure freely. `_deprecations.py` itself is slated for removal in 2.0. |

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
    - **`!` is measured from the newest tag**, the same baseline
      `CHANGELOG.md` uses, so inside an unreleased cycle a signature no
      release shipped is a surface nothing can break. The test is the
      changelog's own: if no bullet belongs under *Breaking API changes*, the
      `!` does not belong either.
- **Stack dependent work** rather than merging to unblock yourself: each
  pull request branches off the previous one and targets it. `gh stack`
  (`init` / `add` / `submit` / `sync` / `rebase`) manages the chain. Two
  consequences: CI's `pull_request` trigger must stay unfiltered, because a
  stacked child targets the branch below it; and nothing may push to a
  branch in a stack outside `gh stack` — including bots, which is why
  `regenerate_notebooks` runs only on its label.
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

**The changelog is the diff from one release to the next, not a log of the
cycle.** So the baseline for every entry is the newest tag, and a change that
only moves code added since that tag is not an entry: it **folds into the
bullet that introduced the feature**, extending it and adding its `(#NNN)` to
the citation. That covers more than it first looks like. A signature settled
twice before shipping, a rename of something never released, a helper made
public and then narrowed, and a bug fixed in code the same open cycle added
are all invisible from the last release — a reader upgrading sees one feature,
arrived at once. Writing them up separately produces a section that reads as
the pull-request order rather than the release, and inflates one feature into
a dozen bullets. The test is the same one the `!` marker uses: state the
entry as a sentence about the newest tag, and if it is not a difference from
it, there is no entry to write.

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
  nonparametric margin, `SciPyMargin` the parametric one — named, it fits
  that family; unnamed, `select` chooses one by AIC / BIC / AICc over a
  curated candidate set — plus `FitControlsMargin` to configure either and
  `FitControlsKde1d` for the kernel knobs, and an adapter registry
  (`as_margin` / `register_margin_adapter`) that accepts a SciPy or PyTorch
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
  The library ships single-vine estimators, and what such a package needs is
  the *public* surface: `distribution=` / `controls=` / `structure=` name the
  model to fit, and `sklearn.base.clone` deep-copies a non-estimator
  parameter, so `clone(est).set_params(structure=...)` derives a member
  without touching the parent. `#339` deleted the copy-on-write `with_*`
  builders that predated that realization.

  Two behaviors here exist for ensembling and are not incidental: a
  pre-settable `schema_` is honored across a *refit*, not only the first
  `fit` (a wrapper refits its survivors, and re-inferring the types there
  silently changes the model), and `VineRegressor.normalize_weights` is a
  real `__init__` parameter so it survives `clone`.

  **No underscore-prefixed name is protected.** An earlier version of this
  file listed private hooks that must not be removed "as unused" because a
  downstream package reached for them. That rule is withdrawn: it froze
  internals by decree and let a private hook substitute for API worth
  designing. An internal with no in-library caller is a candidate for
  deletion, and a downstream need is an argument for making something
  *public*, not for pinning a private name.

- **Scikit-learn-compatible estimators** — `VineDensity`,
  `VineRegressor`, fitting whichever `VinedistBase` subclass they are
  handed (`Vinedist`, or `TorchVinedist` for the PyTorch lane).
- **PyTorch evaluator** — `TorchTllBicop`, `TorchVinecop` (pure-torch
  cascade with GPU placement, autograd, and an optional `batched`
  evaluation fast path; parity with the `Vinecop` cascade to
  floating-point tolerance).
- **Array-agnostic extension layer** — the `BicopLike` / `VinecopLike`
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
- **A copula-family registry to adapt.** The three parametric margin classes
  (`SciPyMargin`, `OpenTURNSMargin`, `TorchDistributionMargin`) each adapt one
  ecosystem's family registry. There is no pair-copula counterpart, and that is
  settled rather than pending: scipy and `torch.distributions` ship no
  copulas at all, and OpenTURNS' 19 are **not** adapted either — pair copulas
  are `lib/vinecopulib`'s own domain (rotations, h-functions, discrete
  handling, tau maps). A learnable pair copula is written by subclassing
  `BicopBase`, which `TorchVinecop` hosts like any other; see
  `examples/10_extending_pyvinecopulib.ipynb`.
- **Density estimators outside the vine framework.** General-purpose
  multivariate density models (normalizing flows, Gaussian mixtures,
  …) are not in scope; `pyvinecopulib` is a vine-copula library.
- **Pinned legacy alias for every old import path forever.**
  Deprecation aliases live in `_deprecations.py` and warn on access;
  they are removed in 2.0. They survive 1.0.0 because that release
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
      kde1d/controls.hpp         # FitControlsKde1d — the one type defined here, not upstream
      docstr.hpp                 # AUTO-GENERATED via scripts/generate_docstring.py (gitignored)

    pyvinecopulib/
      __init__.py                # top-level: core re-exports + lazy sklearn import + deprecation shim
      _deprecations.py           # warn-on-access aliases for pre-#207 top-level names
      py.typed                   # PEP 561 marker (built by scripts/generate_stubs.py)

      core/__init__.py           # Bicop, Vinecop, *VineStructure, FitControls*, Kde1d (re-exports from ext)
        protocols.py             # every contract: the array and its namespace, then Bicop/Vinecop/Margin/Vinedist/Controls
        bicop_base.py            # BicopBase (canonical BicopLike partial impl)
        vinecop_base.py          # VinecopBase (array-agnostic cascades + fit/select)
        vinecop_context.py       # ConditioningContext / Simplified / NonSimplified
        margin_base.py           # MarginBase (canonical MarginLike partial impl)
        vinedist_base.py         # VinedistBase (array-agnostic cascade + IFM fit)
        vinedist.py              # Vinedist (NumPy + compiled Vinecop)
        margin_controls.py       # FitControlsMargin (the family search a margin runs)
        _covariates.py           # the two `x`-forwarding rules + `prepare_covariates`
        _vinecop_discrete.py     # the discrete layouts / per-edge types
        _vinecop_fit_engines.py  # fit_parts / select_parts — the two fit engines (internal)
        bicop_independence.py    # IndependenceBicop
        _placement.py            # place / reference_array / to_numpy + the `_prep` and `_sample_uniform` hooks
        _vinecop_reorient.py     # relabel a structure onto a chosen order tail (internal)
        _rootfind.py             # solve_increasing (monotone bisection; internal)
        _loglik.py               # safe_log / sum_loglik — the two likelihood steps (internal)
        _json.py                 # how a model payload is encoded and written (internal)
        _margins.py              # everything about a margin but its contract: coercion, resolution, JSON (internal)
        extend.py                # the extension surface: placement, covariates, and the batching sentinel
        _trim.py                 # trim — the domain step of the input pipeline
        _validation.py           # the layout / weights / covariate validators (internal)
        _bicop_plot.py           # what `Bicop.plot` / `BicopBase.plot` draw (internal)
        _vinecop_plot.py         # what `Vinecop.plot` / `VinecopBase.plot` draw (internal)
        _margin_plot.py          # what `Kde1d.plot` / `MarginBase.plot` draw (internal)
        _normal.py               # SciPy-free normal / exponential scales for those plots (internal)
      families/__init__.py       # BicopFamily enum + 13 family constants + 15 group constants
      utils/__init__.py          # to_pseudo_obs, wdm, sobol, ghalton, sample_uniform
        _pair_plots.py           # pairs_copula_data plotting helper (pure Python)

      margins/__init__.py        # the two ecosystem adapters + re-exports of core's margin internals
        scipy.py                 # SciPyMargin (one SciPy family, or select one) — needs the [scipy] extra
        openturns.py             # OpenTURNSMargin — needs the [openturns] extra

      sklearn/__init__.py        # VineDensity, VineRegressor
        _base.py                 # VineBase (parameter-constraints, schema, 3-step pipeline)
        _sklearn_private.py      # the scikit-learn internals `_parameter_constraints` needs
        density.py               # VineDensity
        regressor.py             # VineRegressor

      torch/__init__.py          # TorchTllBicop, TorchVinecop, TorchKde1d, TorchDistributionMargin, TorchVinedist, FitControlsTorch*
        tll_bicop.py, vinecop.py # nn.Module evaluators
        distribution_margin.py   # TorchDistributionMargin (torch.distributions adapter)
        vinedist.py              # TorchVinedist (nn.Module margins + distribution)
        kde1d.py                 # TorchKde1d (the torch marginal estimator)
        _margin_kde1d_interp.py  # kde1d's InterpolationGrid, ported — internal
        controls.py              # FitControlsTorchBicop / FitControlsTorchVinecop dataclasses
        _bicop_interp.py         # InterpolationGrid2D (bilinear; Sinkhorn margin renormalization) — internal
        _bicop_fit_tll.py        # pure-torch TLL kernel
        _vinecop_batched.py      # batched evaluation variants
        _placement.py            # the torch lane's `_prep` hook / reference tensor — internal

      _build_info.py             # build provenance, read by `__version__` reporting
      _cpu.py                    # the AVX2 / FMA check the x86-64 wheels need
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

Conventions the toolchain enforces:

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
- **`ruff` is pinned exactly** (formatter output stability), in
  `[dependency-groups] dev` and in `.pre-commit-config.yaml`, which have to
  move together; bump only with a reason. Line length 80; indent width 2.
- **The rule set is 33 families wide and the repository is clean of it.** A
  family goes in unless it contradicts a convention documented here, and an
  exception is either a `# noqa` at the site with its reason or a
  `per-file-ignores` entry where one reason covers the whole file. `ALL` is
  not used: `E111` alone objects 2969 times to the indent width set above,
  `D4xx` re-litigates what `numpydoc` validates, and `COM812` fights the
  formatter. Two mechanics to know before running it: `--select` with `--fix`
  makes every `noqa` for a *non-selected* rule read as unused and **deletes
  it**, so narrow a run with `--statistics` or `--output-format`, never with
  `--select --fix`; and `B009` / `B010` rewrite a `getattr` / `setattr` whose
  whole purpose is to keep `ty` from rejecting the call a test means to make
  at runtime.
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

1. `AGENTS.md` (this file) — invariants and boundaries. Start with
   [the layers](#the-layers-and-which-way-they-depend) and
   [the four levels](#the-four-levels-of-one-construction): which layer you are
   in and which level you are on determine most of what follows, and
   [the rule index](#where-each-cross-cutting-rule-is-written-down) says
   where each cross-cutting decision is stated.
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
- **A submodule bump additionally runs the numerics suites**, which differ by
  submodule. For `lib/vinecopulib` and `lib/wdm`: `tests/test_torch_tll_bicop.py`,
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
- **A signature says `ArrayT`; `Any` is what is left when even `Array` will not
  do.** `ArrayT` is **bounded** by `Array` (`core/protocols.py`), which names
  `shape`, `ndim`, `dtype`, `device`, `__getitem__` / `__setitem__`, the four
  arithmetic operators with their reflections and the three comparisons -- so a
  body that indexes or does arithmetic writes `ArrayT` and needs no hatch, and
  `Array` itself is for a value only ever computed on rather than handed back.
  The bound is written out member by member rather than aliased to the array
  API standard's own, whose bound is `__array_namespace__` and which
  `torch.Tensor` therefore fails.

  A **subclass declaring `Tensor` where the base says `ArrayT` is not a
  narrowing** and needs no hatch: every torch class is parameterized
  (`TorchVinecop(VinecopBase[torch.Tensor])`, and likewise for the other
  four), so `Tensor` *is* the substituted signature. Three `_*_batched` hooks
  carried an `Any` on that mistaken reasoning and no longer do.

  What still takes `Any` is one of four things, each stated at the site: a
  `*args` / `**kwargs` it only forwards; a name resolved at attribute access
  (a module `__getattr__`, the deprecation shim); a foreign object the package
  cannot name (SciPy under a 3.11 floor, the ecosystem adapters in
  `core/_margins.py`, and `declared_eval`, whose return type is open because
  a margin may answer in another array type than it was handed); and a hook
  whose real type needs a structure that does not exist yet --
  `_build_batched` / `_ensure_batched` want a protocol naming what a batched
  vine offers, and `place` / `_prep` / `reference_array` want
  `PlacementMixin` to be `Generic[ArrayT]`, which it is not. Those last are
  the standing list; the rest are permanent. An `Any` in a *signature*
  erases the type for every caller and is published contract text, while one
  in a body is confined to an expression -- which is why a body hands its
  result back through `cast("ArrayT", ...)`, as `core/bicop_base.py` and
  `core/bicop_independence.py` do.

  Where a whole file's `Any` is one reason -- OpenTURNS having no types at all
  -- it goes in `per-file-ignores` with that reason stated once; everywhere
  else it is a `# noqa: ANN401` at the site, and `RUF100` fails the build when
  one goes stale.
- **`__init__.py` files use explicit `__all__`** to define the public
  surface; ruff's per-file ignore (`F403`/`F405`) covers the
  re-export pattern. No wildcard re-exports elsewhere.
- **Tests import from public namespaces** (`from
  pyvinecopulib.sklearn import VineDensity`), not deep internals.
  Underscore-prefixed modules are off limits to tests, with two
  carve-outs. The innermost numeric kernels in `torch/_bicop_fit_tll.py`
  (`_win_smoother`, `_ace`) are reached directly, because what they guarantee
  is not observable through the public surface at the precision that matters —
  a leaking per-lane freeze moves a vine's pdf by less than the arithmetic
  noise a batched fit has to tolerate, and is unmistakable one call in. Import
  inside the test function, as those do, so the module stays out of collection
  for a torch-free run. And the three plot modules
  (`core/_bicop_plot.py`, `core/_vinecop_plot.py`, `core/_margin_plot.py`)
  are called directly because the public surface reaches them only through
  `.plot()`, which draws and returns nothing: their grids, marks and limits
  are checkable at the function and nowhere above it.
- **Generated files stay generated.** `docstr.hpp` and every
  `__init__.pyi` are produced by `scripts/generate_docstring.py` and
  `scripts/generate_stubs.py` respectively. Do not hand-edit; do not
  commit. If a docstring or stub is wrong, fix the C++ source or the
  binding code, then rebuild.
- **Underscore-prefixed modules are internal.** Move helpers into a
  leading-underscore file inside the subpackage that uses them rather than
  exposing them, and name it for the level it serves (`core/_bicop_plot.py`,
  `torch/_bicop_fit_tll.py`).
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
- **One letter, one meaning, in every shape annotation.** `n` is the number of
  observations, `d` the dimension, `p` the number of exogenous covariates, and
  `k` a count of a *subset* of the variables -- the discrete ones in the
  `(n, d + k)` copula layout, the conditioners in `u_cond`'s `(n, k)`. So a
  covariate matrix is always `(n, p)`, which is what `validate_covariates`'
  own error message says; it was documented `(n, k)` at 78 sites, colliding
  with both of `k`'s other uses in the same files.
- **One spelling for a variable type, and one name for it: `var_type` /
  `var_types`, valued `"c"` / `"d"` / `"zi"`.** That is what `MarginLike`
  answers, what a pair copula, a vine and `VinedistBase.from_data` declare,
  what `Kde1d` and `TorchKde1d` take and report, and what the sklearn
  estimators record in `schema_["var_types"]`. Nothing in this package emits
  anything else.

  The long names -- `"continuous"`, `"discrete"`, `"zero-inflated"` -- are
  upstream's: `kde1d`'s `VarType` enum and the `as_str` that renders it. They
  reach this package in two places and leave it in neither. Every constructor
  **accepts** them, because upstream's own parser does and refusing what it
  takes would be gratuitous; and a stored `Kde1d` payload carries one in its
  `"type"` field, which is a format rather than an API and is normalized on
  the way in. There used to be a `_VAR_TYPE_OF` translation table copied into
  two tier-2 modules that cannot import each other, and a `Kde1d.type`
  property beside `Kde1d.var_type` answering the same question differently.

- **Two covariate-forwarding rules, and they are not interchangeable**
  (`core/_covariates.py`). `pair_eval` forwards `x` to a pair copula
  **whenever there is one**, and *raises* first when the pair declares
  `supports_covariates = False`. `declared_eval` forwards to a margin or a
  whole copula **only when it declares the flag**, and *skips* silently when
  it does not, because one distribution may hold conditional and
  unconditional parts side by side -- a per-column choice the caller made,
  not an accident. What enforces the second instead is that the object
  refuses covariates **nothing** reads (`VinedistBase._check_covariates`).
  At fit time there is no skipping: `reject_covariates` refuses outright.

  Both now read the same flag, and the difference that matters is
  raise-versus-skip, not read-versus-ignore. Collapsing the first into the
  second -- making the *forwarding itself* conditional on the flag -- would
  make a forgotten one
  a **silent unconditional fit**, which is the one outcome neither rule may
  produce. Raising on an absent declaration has no such failure mode, and the
  forwarding itself stays unconditional for anything that passes. The
  signature cannot serve as the declaration here even though `ty` makes every
  `BicopBase` subclass write `x`: `inspect.signature(Bicop.pdf)` is
  `(*args, **kwargs)`, so at runtime nothing can read it, and a vine
  evaluating with covariates over pairs that take none used to fail as
  nanobind's overload dump rather than by name.
- **A mixin at the placement position defines only ordinary private methods.**
  `core/_placement.PlacementMixin`, `QrngUniformMixin` and
  `torch/_placement.TensorPlacementMixin` all land **ahead of
  `torch.nn.Module`** in the MRO of every torch class, because the canonical
  bases precede it. So anything such a mixin defines that `nn.Module` also
  defines would shadow it silently -- `to`, `cpu`, `cuda`, `state_dict`,
  `load_state_dict`, `_apply`, `parameters`, `buffers`, `forward`, and every
  dunder. `VinecopBase.__getstate__` reaches `nn.Module.__getstate__` through
  that same chain, so a mixin defining `__getstate__` or `__reduce__` breaks
  every torch pickling test. Note `__repr__` is *already* shadowed this way
  through the bases: `repr(TorchTllBicop.from_data(u))` is `'TorchTllBicop()'`
  rather than the submodule tree a torch user expects. That is the failure a
  `__repr__` mixin would institutionalize, and the reason not to write one.

  Ordering is part of the contract, not a detail: a torch class lists
  `TensorPlacementMixin` **before** its canonical base, or the base's
  `PlacementMixin` linearizes first and `_prep` resolves to the array-API
  inference instead of the tensor one.

- **One spelling for a parameter's type, and it says which.** Three forms,
  each with a rule: `<T>, or None, optional` for a parameter that accepts
  `None`; `<T>, default=X` where the default is a real value; and
  `<T>, or None` for a *return* that may be `None`, never with `optional`,
  which means nothing in a `Returns` block. `PR04` is enabled, so a parameter
  with no type at all fails the docs check rather than rendering without one --
  which is how `Vinecop.from_structure` shipped four untyped parameters.
  Documented classes are linked from a type field automatically: `docs/conf.py`
  derives `numpydoc_xref_aliases` from the same `_CLASS_MODULE` table
  `process_cross_references` uses, so a name cannot be a link in prose and
  dead text in the table below it. `numpydoc_xref_ignore = "all"` is what
  keeps that safe -- numpydoc leaves every unrecognized token alone, so the
  shape grammar stays plain text instead of becoming unresolvable references
  under `-W`.

- **A capability flag exists where a consumer reads it, and nowhere else.**
  `supports_covariates` is declared on `BicopBase`, `MarginBase` and
  `VinecopBase` -- `pair_eval` reads the first and `declared_eval` the other
  two -- and it means the same thing at each: *this part accepts `x`*.
  `supports_weights` is declared on the same three, and means *a fit of this
  honors `controls.weights`*. Both are absent from `VinedistBase`, where the
  aggregate questions are derived from the parts instead
  (`margin_class` and `vinecop_class`), so there is nothing to declare.
  Adding a flag with no reader is the thing `supported_var_types` was deleted
  for, and `supports_controls` after it -- `controls_class is None` says
  "reads no controls" and carries the type besides.
- **A protocol declares what hosting a part needs, and supplies its own
  defaults.** `VinecopBase.bicop_class` is declared
  `type[BicopLike[Any]] | None`, so `BicopLike` is by declaration "what a vine
  may host" -- which is more than the six evaluating members. `cdf` and
  `with_var_types` are what a **discrete** edge needs, `flip` what *selection*
  needs, and `supports_covariates` / `supports_weights` what a consumer must
  read before handing over a matrix or weights. All are members.

  A member an implementation cannot serve gets a **default in the protocol
  body**, and that default is *the fallback the consumer used to apply* --
  `logpdf` is `safe_log(pdf)`, `var_type` is `"c"`, a capability flag is
  `False`, and `flip` / `cdf` / `sample` raise naming the class. So the
  `getattr` probes disappear from the consumers, and `isinstance` means
  "hostable" rather than "evaluable". `bicop_base.flip_of` and
  `continuous_of` existed only because `flip` and `with_var_types` were not
  members, so `ty` would not let a `BicopLike`-typed value call them; both are
  gone and their call sites read `pair.flip()` / `pair.with_var_types()`.

  The cost, accepted rather than overlooked: **a Protocol's default reaches only
  a nominal subclass.** `runtime_checkable` `isinstance` compares member
  *names*, so a duck-typed class must define every member itself. The
  documented way to implement a contract directly is therefore to *inherit*
  it -- `class Mine(BicopLike[np.ndarray])` -- and write the abstract members;
  everything else arrives with a default. A six-member duck type no longer
  satisfies `BicopLike`, which is the point: it could not be selected on,
  could not sit on a discrete edge, and said `True` anyway.
- **American English** in code, comments, documentation, commit messages,
  and changelog entries: *behavior*, *normalize*, *serialize*, *finalize*,
  *center*, *modeling*, *honored*, *color*. There is no legacy exemption.
  `codespell` enforces it in `make lint` and in pre-commit, using its
  `en-GB_to_en-US` dictionary, so the rule covers every British spelling
  rather than a list this repository happened to drift on.
  A banned word or phrase is occasionally the right one: wrap those lines in
  `codespell:ignore-begin` / `-end`, which both `codespell` and
  `tests/test_prose.py` skip. Use that form rather than the single-line
  `codespell:ignore <word>`, which does not match a hyphenated entry and
  reports the marker itself. The marker goes inside whatever comment the file
  uses -- `#`, `//`, or an HTML `<!-- ... -->` in Markdown; neither half
  anchors it to the end of a line, so the trailing ` -->` is fine. It catches
  ordinary typos in the same pass. Configuration -- the skip list for
  generated and vendored files, and the domain vocabulary it would otherwise
  flag -- lives in `[tool.codespell]` in `pyproject.toml`; add a word there
  only when it really is a term of art, never to silence a real
  misspelling.

- **Write plainly, and there is a list.** `.codespell-prose.txt` names the
  words this repository does not use, each with its replacement and the
  reason. It is a codespell dictionary, wired in through
  `[tool.codespell] dictionary`, so `make lint` and pre-commit both enforce
  it; every entry carries a reason after a comma, which is what makes
  codespell **report a suggestion without auto-fixing** -- the right
  replacement depends on the sentence. Two mechanics worth knowing before
  editing that setting: the file must come *first*
  (`".codespell-prose.txt,-"`), because a value beginning with `-` is read as
  another flag, and the trailing `-` is what preserves the `en-GB_to_en-US`
  builtins the rule above depends on.

  The list is the source of truth and is not repeated here -- a rule that
  quotes its own banned words fails itself, which is how this entry was first
  written. The substitutions mostly converge on vocabulary the repository
  already prefers rather than coining any; the largest of them replaced a term
  used 42 times with one already used 91 times for the same idea, which
  `core/_placement.py` had been using both of, two lines apart.

  codespell tokenizes, and cannot tell a docstring from an identifier, so
  `tests/test_prose.py` carries what it structurally cannot: the multi-word
  phrases, and the exemptions -- a `Tensor` parameter named for a mask in
  `torch/_bicop_fit_tll.py`, the ordinary noun in "test harness", and one phrase
  quoting upstream PyTorch's own support tier.

- **Never say "compiled" or "C++" in a docstring.** A user of
  `pyvinecopulib.core.Bicop` does not know — and does not need to know — that
  it is bound from C++, so the words are implementation detail leaking into
  rendered API text. Name the class instead: not "the compiled pair-copula
  controls" but "controls for a ``Bicop`` fit"; where a *contrast* with the
  PyTorch layer is the point, contrast the classes (``Vinecop`` versus
  ``TorchVinecop``) or say "core" rather than "compiled". Inline `#` comments
  are exempt — they address the next person to edit the code, for whom the
  distinction is required — and so are the internal `torch/_*.py` fidelity
  modules, where "reproduces the compiled `kde1d`" *is* the documented
  contract. Note `torch.compile` is a different sense of the word and stays:
  `compile_cascades` is actually about compilation. Same reasoning as the
  entry below, applied to vocabulary.

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
  the sites where libclang actually cannot disambiguate (they are
  enumerated in `src/include/vinecop/class.hpp`).
- **Doxygen upstream, numpydoc downstream.** The generator translates one
  into the other. Never "fix" an upstream `//!` comment into numpydoc.
- **One argument order on every estimator: the observations, then `controls`,
  then keyword-only whatever the object cannot infer.** `fit`, `select` and
  `from_data` take `(data, controls)` positionally on all four bases, on the
  torch lane and on the compiled `Bicop` / `Vinecop` -- and `structure`,
  `matrix`, `var_types`, `margins`, `margin_controls`, `names`, `x`, `weights`
  and the callback hooks are keyword-only. `Vinecop.from_data` used to take
  `controls` *fifth*, behind `structure`, so the call a user carries over from
  `fit` bound a controls object as a structure: a `TypeError` from the binding
  and an `AttributeError` naming `dim` from the Python lane. The rule is worth
  more than the two characters it costs at a call site, and the changelog
  claimed it before the code did. There is **no** exception: `Kde1d` used to be
  one -- its second positional was `weights`, on the grounds that it took no
  controls at all -- and `#339` gave it `FitControlsKde1d` and the same
  `(y, controls)` every other margin reads. What the exception actually cost
  was a validator (`reject_array_controls`) whose whole job was catching the
  carried-over `kde.fit(y, w)` at six call sites.
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
  Where two forms actually are one operation, prefer a single method that
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
- **Do not say "backend".** It named the sklearn adapter layer that `#339`
  deleted, and it was always the wrong word for the NumPy-vs-PyTorch array
  distinction. Say "array namespace", "array library", "array-agnostic", or
  just "NumPy or PyTorch"; for the two sides of the port, "lane".
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
  only allowed suppression is the short `nitpick_ignore_regex` list
  (upstream-C++ getter-name mismatches, `BicopFamily` value aliases,
  scikit-learn-generated methods); prefer fixing a reference over extending it.

### Maintaining this file

If a coding agent or reviewer keeps repeating the same correction —
about a convention this repo enforces — update `AGENTS.md` rather
than relying on tribal knowledge. Do not add ephemeral, user-specific,
or machine-local preferences here. The `CHANGELOG.md` is the place for
release-by-release context; this file is for invariants.

## Module boundaries

### The layers, and which way they depend

Four tiers. An import may point **down** a tier, never up, and within a tier
only where this list says so. `tests/test_import_surface.py` pins the whole
edge set, so crossing a layer is an edit to a declared table with
the reason written beside it — not something a stray import can do quietly.

```text
  __init__.py                                     the root re-export surface
      |
      v
  margins      torch      sklearn                 tier 2: may need an extra
      |          |           |
      +----------+-----------+---> core           tier 1: NumPy only
                                    |             (with families and utils)
                                    v
                            pyvinecopulib_ext     tier 0: the binding
                                    |
                                    v
                     lib/{vinecopulib,wdm,kde1d}   upstream C++
```

- **Tier 1 depends on no optional extra and on nothing above it.** That is
  what makes `import pyvinecopulib` work with nothing but NumPy installed,
  and most of the rules below follow from it. **No** import reaches up into
  `margins`, at module scope or inside a function. There used to be one, and
  what it resolved is why it is gone: `margins="parametric"` was a string
  `core` accepted and had to turn into `SciPyMargin`, a class it can *name*
  but not *contain* because that needs the SciPy extra. `#339` deleted the
  argument rather than the deferral — a caller names the class as a
  `margin_class`, which is their own import. A new such import needs an
  argument of the same kind, a `core` API whose contract names a class behind
  an extra, and not merely the same shape: everything an extension point can
  carry is registered by the module that owns the class instead.
- **Within tier 2 there is exactly one edge.** `sklearn` imports `margins`
  at module scope (both need no extra of `sklearn`'s own). It does **not**
  import `torch`, at module scope or inside a function: the PyTorch lane is
  opted into by handing an estimator `distribution=TorchVinedist`, so the
  caller's own import is the signal and no name in `sklearn` mentions torch.
  That is why the parameter takes a class rather than a string. `margins` and
  `torch` import neither of the other two.
- **The callables the binding looks up by name live under `core`.**
  `Bicop.plot`, `Vinecop.plot` and `Kde1d.plot` are bound to
  `core/_bicop_plot.py`, `core/_vinecop_plot.py` and `core/_margin_plot.py`,
  resolved by module path at call time — which is why `tests/test_plots.py` is
  what proves a repoint, and why a wrong path fails only when the method runs.
  Housing them beside the extension instead put them *above* `core` while
  `core` imported back into them, a cycle two declared edges recorded rather
  than forbade. Within tier 1 the one edge is `utils` -> `core`, for the
  SciPy-free normal and exponential scales in `core/_normal.py` that both
  `core/_bicop_plot.py` and `utils/_pair_plots.py` draw on.
- **The extras stay out of `__all__`.** They are reachable through the
  top-level `__getattr__` only, because `from pyvinecopulib import *`
  resolves every name in `__all__` and would otherwise make PyTorch a hard
  dependency of the one import form beginners reach for first.

### The four levels of one construction

A margin, a pair copula, a vine and a vine distribution are the same
construction four times: a `runtime_checkable` Protocol naming what a
*consumer* needs, a canonical base supplying everything derivable from it,
and a short list of members a subclass owes. Know which level you are writing
and
most of the rest is determined.

| Base | Protocol requires | Abstract — no evaluating without it | Reports its own absence — only fitting needs it | Names its parts as |
|---|---|---|---|---|
| `MarginBase` | `pdf`, `cdf`, `icdf` | `pdf`, `cdf` | `fit` | — |
| `BicopBase` | `pdf`, `hfunc1/2`, `hinv1/2`, `sample` | `_pdf_raw`, `_hfunc1_raw`, `_hfunc2_raw` | `fit`; `_flip_raw` and `_cdf_raw` to host the pair in *selection* or to declare it *discrete* | — |
| `VinecopBase` | `pdf`, `cdf`, `rosenblatt`, `inverse_rosenblatt`, `sample`, `structure` | `get_pair_copula` | `set_pair_copulas` | `bicop_class` |
| `VinedistBase` | the ten above plus `logpdf`, `loglik`, `margins`, `vinecop`, `copula_layout` | *(none)* | `_coerce_fit_data` | `vinecop_class`, `margin_class` |

Three things in that table have a reason. Do not undo them:

- **The two middle columns are different mechanisms on purpose.** A member is
  `@abstractmethod` when the object cannot be *evaluated* without it, and a
  stub that raises when it is needed only to *fit* — so a vine that merely
  hosts pairs, or an immutable one, is still a valid subclass and says so at
  the one call it cannot serve. The rule is stated once, under
  *"`get_pair_copula` reads, `set_pair_copulas` writes"*.
- **The protocol is still narrower than the base, and now by less.** It
  carries what a consumer in this library asks of the part, each member with
  the default the consumer used to supply; the base adds the machinery --
  the dispatcher/leaf split, the input pipeline, the inherited inverses and
  samplers, the estimator surface. What a protocol may *not* gain is a member
  no consumer reads.
- **`VinedistBase` has no abstract member at all.** A vine distribution is
  determined by its two halves, so nothing has to be declared to evaluate
  one; naming the part classes is what makes it *fittable*, and `_fit_copula`
  reports a `vinecop_class` of `None` rather than the base pretending it
  could fit one.

### Where each cross-cutting rule is written down

The decisions that span layers live next to the layer that motivated them
rather than in one chapter. This index is the map; the rule is stated once,
where the link points, and nowhere else.

| The decision | Stated in |
|---|---|
| The three input steps — placement / layout / domain — and why covariates are placed but never trimmed | `### pyvinecopulib.core`, *"One input pipeline, three separable steps"* |
| Which of the two covariate-forwarding rules applies to a callee | `### Coding conventions`, *"Two covariate-forwarding rules"* |
| Whether a capability flag may exist at all | `### Coding conventions`, *"A capability flag exists where a consumer reads it"* |
| The argument order every estimator method takes | `### Coding conventions`, *"One argument order"* |
| What `fit` / `select` / `from_data` each mean, on all four bases | `## Extension points`, *"Fitting has one shape across all four bases"* |
| What naming a part class buys, and what `None` means | `## Extension points`, *"Declare the parts, inherit the fitting"* |
| Which hook a subclass must implement to evaluate vs. to fit | `## Extension points`, *"`get_pair_copula` reads, `set_pair_copulas` writes"* |
| Which modules may carry a leading underscore | `### pyvinecopulib.margins`, *"A margin class is named for the ecosystem"* |
| Why a fitted slot's conditioning order is state, not a reading of the matrix | `### pyvinecopulib.core`, *"A selected slot's conditioning order"* |
| What may break, and what needs a deprecation alias | `### Stability tiers` |
| Which suites actually run where, and what silently skips | `### Which CI leg covers what` |

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

**`src/include/` may define a C++ type of its own, and one does.**
`src/include/kde1d/controls.hpp` declares `FitControlsKde1d`, which exists
nowhere upstream: `lib/kde1d`'s own API takes the kernel knobs as constructor
arguments and needs no struct. What the struct is for is **Python API shape** --
one controls object per class, so `Kde1d` reads `(y, controls)` like every
other estimator here -- which is this repo's concern and not upstream's.
That is the test for a second one: binding-side shape, not behavior. Anything
a C++ caller would want belongs upstream, where it also gets upstream's
tests.

### `pyvinecopulib_ext` (the nanobind extension)

- Single binding module compiled from `src/pyvinecopulib_ext.cpp` +
  `src/include/**`.
- `nanobind` (>= 2.7) is the binding system, **not** pybind11. The
  pybind11 → nanobind switch changed the API on purpose;
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
  `controls=`) and `"random_unweighted"`.
- `FitControlsVinecop.conditioning_set` (property + pickled, not a
  positional ctor arg) drives conditioning-aware selection — the fitted
  order ends with the given 1-based variables so
  `Vinecop.sample_conditional` / `reorient` can condition on them.
- `RVineStructure.from_trees(d, trees)` is the **faithful** inverse of
  `RVineStructure.get_trees()` (identity diagonal policy — each edge's
  `conditioned[0]` on the diagonal — so `from_trees(s.dim, s.get_trees()) == s`).
  Upstream `Vinecop.select` finalizes with the *same* (flip-free)
  convention, so `VinecopBase.select` assembles its selected trees through this
  same `from_trees` and matches the compiled selector's matrix exactly —
  one diagonal convention throughout.
- **Array-agnostic abstraction layer** (pure Python; `core` imports
  without PyTorch). The extension point for custom (e.g. neural,
  conditional) pair copulas and vines:
  - `BicopLike[ArrayT]` / `VinecopLike[ArrayT]` (`protocols.py`) —
    generic, `runtime_checkable` protocols mirroring the `Bicop` /
    `Vinecop` evaluation surface on any array namespace (NumPy or
    PyTorch). `Bicop`, `Vinecop` and `Kde1d` satisfy them **structurally**,
    which `ty` checks and `tests/test_core_protocols.py` pins, so downstream
    code can type against the contract instead of a concrete class. What makes
    that possible is that the protocols describe the *unconditional* surface:
    a protocol may ask for **fewer** parameters than an implementation
    provides, never more, so declaring the conditioning matrix `x` there would
    put the three classes that model no covariates outside their own
    contracts. `x` lives on the four bases instead, keyword-only, and a
    subclass that reads covariates declares it on the primitives that do.
  - **One input pipeline, three separable steps, one owner per level.** Every
    layer does the same three things to an incoming array, and they are kept
    apart because they do not always apply together:
    **placement** (`_placement.py`'s `place`, reached through the `_prep(a)`
    hook the four bases inherit from its `PlacementMixin`) puts the values on
    the namespace, dtype and device the object evaluates on; **layout** (the
    `_layout` hook, one per level -- a two-column check on `BicopBase`, a
    one-dimensional check on `MarginBase`, and `collapse_data`'s
    `var_types`-dependent widths on `VinecopBase`) says which shapes are
    admissible; **domain** (`_trim.py`'s `trim`) clamps copula arguments into
    the open unit square at the working precision. `trim` is not
    entry-only: two of its nine call sites are pipeline entry and seven clamp
    an h-function or distribution-function value the cascades produced, which
    is the same domain question asked on the way out. The
    composites that apply all three to a copula argument are `_prep_args` —
    `BicopBase._prep_args(u)`, `MarginBase._prep_args(y, name)` (placement
    plus the single-column layout -- a margin's argument is on the data scale,
    so it is never clamped) and `VinecopBase._prep_args(u, name, *,
    values_only)`. What forces the split is that **exogenous covariates are
    placed but never trimmed**: they are arbitrary reals, not copula
    arguments, and `prepare_covariates(onto, x, n)` is the composite applying
    exactly those two steps -- called at every entry point that takes an `x`,
    including the static fit engines, where an *array* is its own placement
    reference. Placing `x` is not cosmetic: a non-simplified vine concatenates
    it with the conditioning columns it gathered, so a NumPy `x` handed to a
    PyTorch vine has to be brought across before they can meet. Placement is *inferred* from the arrays an object already
    holds, so hosting a subclass on PyTorch requires writing none of it —
    override `_prep` only where those arrays live somewhere the inference
    misses, or where finding them again per evaluation costs more than naming
    them (`TorchTllBicop` names its grid; `torch/_placement.py` reads a
    module's own tensors for the classes that hold submodules, 12.1 us against
    `reference_array`'s 64.9). The one array a base manufactures from nothing is
    `BicopBase.plot`'s evaluation grid, which is why that is the one place the
    hook is required rather than a convenience.
    Since the inference is the whole contract, what it reads has to be right,
    and **every** search for a reference array ranks candidates the same way:
    a **floating-point** array wins. An object may hold an index table or a
    count buffer, and adopting `int64` from it placed every copula argument at
    zero — a wrong answer, not a failure. `reference_array` keeps an integer
    array as its fallback, whose dtype `place` then does *not* adopt, since one
    still names a namespace and a device; the torch lane's `reference_tensor`
    has no fallback at all, because each of its callers has a floating default
    of its own and an integer dtype is the one answer none of them can use.
    Two searches over one object that rank differently is what put a margin's
    `_prep` on `float64` while its own sampler drew in `int64`, so a new
    placement site adopts the rule rather than restating the first-hit loop.
    The return trip is `to_numpy`, in the same module and shared for the same
    reason: `np.asarray` alone raises on a tensor that requires grad and again
    on one held on an accelerator, so everything outside the array namespace
    that reads a value -- the criterion binding, the three plots, the sklearn
    estimator boundary -- goes through the one walk rather than a fourth copy
    of it.

    **Placement is exported; the rest of the pipeline is not, and the
    inference is best-effort.** `place`, `reference_array`, `to_numpy`,
    `prepare_covariates` and `covariate_row` are named in
    `pyvinecopulib.core.extend` — never in `core` itself, which the two
    namespaces being disjoint is tested for. What earns them the place is that
    `_prep` is a real override with a documented recipe, and that an extension
    which wrote its own covariate check instead ran a second, 1-d-accepting
    contract on the same object. `trim` is *not* among them: it is the domain
    step of a composite (`_prep_args`) that no subclass overrides, so nothing
    documented reaches it. And placement has a **third** answer besides a
    namespace and a
    failure: an object holding no array has nothing to infer from, so `place`
    returns the values untouched. That is right for a part that computes in
    whatever namespace it is handed and silently wrong for one that does not — a
    torch class that is no `nn.Module` at all and keeps its device as a handle
    rather than as a tensor is the case that hits it. Hence
    `reference_array(obj) is None` as the check and an overridden `_prep` as the
    fix, both stated on the hook; and hence *not* a loud `place`, which would
    refuse the functional part the `None` was written for. `TensorPlacementMixin`
    is no answer to it either: it reads a module's registered tensors, so it
    needs the `nn.Module` such a class does not have.
  - `BicopBase` (`bicop_base.py`) / `VinecopBase` (`vinecop_base.py`) —
    canonical partial implementations to subclass. A `BicopBase`
    subclass defines `pdf` / `hfunc1` / `hfunc2` and inherits `hinv1` /
    `hinv2` (bisection), `sample`, `loglik`, `plot` — which takes an optional
    single-row `x`, since a conditional pair's density is a different surface
    at every covariate value and a 2-d plot shows one slice (`flip` — needed
    only to host the pair in structure *selection* — defaults to
    raising); a `VinecopBase` subclass defines the one hook
    `get_pair_copula` and inherits the whole tree-by-tree cascade plus
    the public `fit` and `select`. `select` is an
    exact port of `Vinecop`'s Dissmann / Wilson structure selection
    (same matrix encoding, selection-time pairs reused via `flip`, no
    re-fit; parity is a hard guarantee). `threshold` acts twice there, and
    both halves are ported: it deprioritizes an edge in the spanning tree,
    *and* a surviving edge below it holds `IndependenceBicop` instead of a
    fit (`tools_select.ipp` `fit_or_reuse_pair_copula`). Porting only the
    weight is a silent divergence, and the default `threshold=0.0` hides
    it: nothing is below zero there, so every test that does not set it
    sees the two agree. `TorchTllBicop` / `TorchVinecop` are the torch
    subclasses.
  - **A selected slot's conditioning order is fitted state, not a reading of
    the matrix.** `select` fits each edge in the orientation the search built
    it in and reorients it at finalization with `flip`, which swaps the pair's
    two arguments and leaves its conditioning columns alone. So on a swapped
    slot the C1 order the finalized matrix names is the *other* endpoint's
    chain — the same conditioning set in a different order, measured at 19 of
    272 slots (7%) — and gathering `u_D` in it evaluates a conditional pair on
    a permutation of what it was estimated on. `_select_parts` therefore
    returns the order each pair was fitted on as a third value, `_set_cond_order`
    installs it, and `_cond_positions` answers with it where there is one. Three
    consequences: `fit` **drops** it (it fits along the structure's own order,
    so an order left over from a `select` is a claim about pairs that are gone);
    a simplified vine never gathers `u_D`, so all of this is inert there; and
    the order is carried as **variable labels**, not natural-order columns, so
    it survives the relabeling `conditioning_set=` performs. Do not "simplify"
    this back to `struct_array` — the two agree on 93% of slots, which is
    exactly enough for a spot check to pass.
  - **Variable types are state on the pair copula, not a wrapper around it.**
    `BicopBase` carries `var_types` (defaulting to `("c", "c")` as a *class*
    attribute, so a subclass writing its own constructor is still continuous)
    and a copy-returning `with_var_types`, exactly as `AbstractBicop` carries
    `var_types_` and `Bicop::with_var_types` returns a copy. The types are
    state because they enter the **fit**: a `tll` grid fitted on a discrete
    edge goes through `find_latent_sample` and is a different estimate, not
    the same copula viewed differently.
    `pdf` / `cdf` / `hfunc1` / `hfunc2` / `hinv1` / `hinv2` are therefore
    concrete **dispatchers** on the base, and a subclass writes the continuous
    leaves `_pdf_raw` / `_hfunc1_raw` / `_hfunc2_raw` (plus `_cdf_raw`,
    `_hinv1_raw`, `_hinv2_raw`, `_flip_raw` where it has them) — `_raw` being
    `AbstractBicop`'s own name for the primitive that ignores `var_types`. A
    leaf always sees two continuous columns, already placed, checked and
    clamped, so it never calls `_prep_args` itself.
    The vine still **owns** the layout: `_bind_vine(..., var_types=)` declares
    which variables have atoms, `pair_var_types(tree, edge)` derives each
    slot's types from the structure alone, the cascades hand a four-column
    `[u1, u2, u1^-, u2^-]` argument to a slot whose types include `"d"`, and
    the fit engines declare those types on each pair as it is fitted
    (`_declared`), so a `fit_edge` callback never has to. `continuous_of` is
    the one reading in the other direction, for the inverse Rosenblatt cascade
    and the density plot, both of which evaluate every pair continuously.
    **A discrete edge needs a pair that overrides `with_var_types`** — the
    protocol's own default answers only for a pair declaring no atoms, and
    refuses by name otherwise, which is what a `BicopBase` subclass and a
    `Bicop` both supply.
    The parity test that binds is the **normalization identity**
    `Σ_atoms c(u₁,u₂)·(u₁ − u₁⁻) = 1`: the quotients telescope, so it holds
    exactly and needs no reference implementation and no tolerance argument.
    It is what established that the Python quotients were right and the
    compiled `tll` pair was wrong (fixed upstream in vinecopulib#739 and
    pinned since). Parametrize pair-level parity over **every** family, not a
    representative couple — covering only `gaussian` and `clayton` is why that
    class of defect stayed invisible on both sides for as long as it did. Note
    the identity **cannot** catch a cache regression: it telescopes to the four
    corners, so it reads `1 − 2e-10` for a correct density and for a 38%-wrong
    one alike.
    **An atom's probability is read, not differenced, wherever the pair can
    read it** — the choice `AbstractBicop` makes per family, mirrored here.
    `BicopBase.rect_prob` / `cond_interval_prob` are ordinary methods with the
    four-corner and two-term defaults, built on the *leaves* rather than the
    public members, so the dispatchers cannot recur through them;
    `TorchTllBicop` overrides both onto its grid, as `KernelBicop` overrides
    them onto `InterpolationGrid`. Differencing amplifies an absolute error by
    `4/(w₁w₂)` in the atom widths and reading the mass by `1/w₂` alone, which
    at the widths a vine's inner trees reach is the difference between `8.5e-8`
    and `7.0e-14` on the torch↔`Vinecop` cascade comparison. `rect_mass` /
    `cond_interval_mass` stay `InterpolationGrid2D`'s own names, as upstream
    keeps them `InterpolationGrid`'s; nothing outside the torch grid reaches
    for them.
  - `sample_conditional` / `reorient` (`_vinecop_reorient.py`) — conditional
    sampling and
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
    `NonSimplifiedContext` (`vinecop_context.py`) — the per-edge policy that
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
  `BicopLike`. **`x` means exogenous covariates everywhere in the Python API,
  with no exception.** `Kde1d` was one: `lib/kde1d` names the observations `x`
  and `icdf`'s argument — a probability — `x` too, and the binding carried both
  through. `#339` renamed the three *fitting* verbs, whose signatures it writes
  itself, to `(y, controls, *, var_type, support, weights)`. The evaluation
  verbs (`pdf`, `cdf`, `logpdf`, `cdf_left`, `loglik`, `icdf`) still read `x` /
  `x`, and renaming them is an **upstream** change: their docstrings are lifted
  verbatim from `lib/kde1d`, so a local rename would document one name and bind
  another. Fix it there and bump the pin; do not patch `src/include/**`.
  Nothing in the library depends on either spelling, by design rather than by
  luck: the protocol makes the observations **positional-only** for exactly
  this reason, and `declared_eval` calls every margin method positionally, so
  the argument name is never used as a keyword. `pdf` means *the density with respect to the margin's own
  reference measure* — a Lebesgue density for a continuous margin, a
  probability mass at an atom — which is what makes
  `log f(x) = log c(u) + Σ_j log pdf_j(x_j)` hold verbatim for
  continuous, discrete and mixed margins with no branch in the
  likelihood path. `MarginBase` (`margin_base.py`) needs only `pdf` /
  `cdf` and supplies `icdf` (bisection), `logpdf`, `cdf_left`, `loglik`,
  `sample`, `plot`, `var_type`, `support`, `is_fitted`, the `nobs` /
  `n_parameters` a criterion penalizes against, and a raising `fit`.
  `plot` draws the density or the distribution function of any of the three
  variable types, on the `BicopBase.plot` pattern -- the grid is manufactured
  from nothing, so it is placed through `_prep`, and one covariate row is a
  slice of a conditional margin rather than the whole of it. It refuses an `x`
  a margin does not declare, because forwarding to a margin is by flag and
  `declared_eval` *skips* rather than raises: the alternative is the
  unconditional curve under a conditional-looking call.
  `{pdf, cdf, icdf}` are the only **abstract** members; `var_type` ∈
  `{"c","d","zi"}`, `support`, `cdf_left`, `logpdf`, `sample` and the two
  capability flags are members carrying the default a consumer used to apply
  — `"c"`, unbounded, `cdf`, `safe_log(pdf)`, a raise, and `False`. So a
  margin with nothing else to say still composes, and one that has something
  to say overrides it.
- **`Vinedist`** (`vinedist.py`) composes a `VinecopLike` with one
  margin per variable. It owns the copula-scale layout: `copula_data`
  builds the compact `(n, d + k)` matrix from `cdf` and `cdf_left`,
  clamps once, and checks `cdf_left <= cdf`, so callers never `hstack`
  a left-limit block by hand. `logpdf` sums logs rather than
  accumulating a product, because **both** terms carry a scale a product
  loses: the marginal densities underflow at a `d = 50` product, and so does
  the copula factor, itself a product of up to `d(d-1)/2` pair densities. Each
  half is therefore read in log space where it offers one — `logpdf` is an
  optional capability on `VinecopLike` and on `MarginLike` alike, read with
  `getattr` at this one consumer, with `safe_log` of the density as the
  fallback. `sample` is
  `marginal_icdf(copula.sample(n, ...))`, so it inherits the copula's
  quasi-random and seeding options and never calls a margin's own
  sampler. `sample_conditional` is the same sandwich around the
  copula's, with one difference the caller feels: a discrete
  conditioner needs no left-limit column, since it is derived from that
  variable's own margin. Both scales resolve an omitted
  `conditioning_set` through one `infer_conditioning_set`, so the
  column-to-variable rule cannot drift between them. Every method also takes optional exogenous covariates `x`,
  forwarded to each margin that declares `supports_covariates` and to a
  copula that declares it too. Two hooks keep the array namespace
  coherent: `_prep` (identity here, `torch.as_tensor` on
  `TorchVinedist`) coerces one input array onto the parts' namespace, so a
  caller may hand the type they have; and `copula_data` / `marginal_cdf` /
  `marginal_icdf` take `xp` from
  the *columns the margins returned*, never from the input — a margin may
  legitimately answer in another array type, and stacking that through the
  input's namespace either raises or silently detaches.

### `pyvinecopulib.margins`

The two ecosystem adapters, kept out of `core` because they are the only part
that needs an extra. The **contract internals live in `core`**, which own the
half a `Vinedist` fit runs on: `MarginLike` / `MarginBase`, `FitControlsMargin`
(`core/margin_controls.py`) and everything in `core/_margins.py` -- the two
registries and the `margin_controls=` resolution. None of those needs
SciPy -- they
import stdlib, NumPy and `core` -- and putting them here had `core` reaching
*up* a layer at ten sites, three of them into a private module of a package
above it, all deferred to hide the cycle. `pyvinecopulib.margins` re-exports
them, so its documented surface is unchanged and it stays where a user looks
for margins. No `core` -> `margins` import remains, deferred or otherwise.

Three groups:

- **Built-in margins** — `Kde1d` *is* the default margin, needing no
  wrapper; it takes `xmin` / `xmax` so a bounded variable is not fitted
  past its support — the sklearn estimators fill those in from a
  categorical's declared levels, since otherwise the density grid is
  padded past the data. What a bound *means* differs by variable type, and
  `docs/concepts.rst` (`concepts-kde-margins`) is where that is written down:
  for a discrete variable it is the **integer support**, so the bound and the
  data must both be integers and the fitted grid runs half a unit wider at each
  end. A categorical whose levels are not integers is therefore refused, by
  name, at fit time rather than by a bare `invalid_argument` from C++. Then
  `SciPyMargin`: named a family, `fit` estimates it; unnamed, `select` fits
  every admissible candidate and *becomes* the best one. Selection is a
  **method on the margin**, not a wrapper class — the shape `Bicop.select`
  has always had, and the reason there is no `MarginSelector`: if the wrapper
  were the right shape the library would want a `BicopSelector` and a
  `VinecopSelector` too. **Naming a family is the choice**, so `select` on a
  named margin reduces to `fit` rather than replacing what the caller asked
  for; an unnamed `SciPyMargin()` is the signal to search, and `family_set` is
  the explicit request to re-search a named one. A narrowed search anchors
  each family exactly as the curated search would (`_anchoring_group`), or the
  two estimate different parameter counts for the same family and their
  criteria stop being comparable.

  **A margin class is named for the ecosystem whose families it wraps** --
  `SciPyMargin` in `margins/scipy.py`, `OpenTURNSMargin` in
  `margins/openturns.py`. Neither module is underscore-prefixed, because each
  *is* an import path a user may reasonably reach for -- both are named for an
  ecosystem and behind its extra; the same-named modules do not shadow
  the real packages, since Python 3 resolves `import scipy` absolutely.

  **The underscore describes the module, not the names it exports.** It says
  "not an import path": `core/protocols.py`, `bicop_base.py`, `vinedist.py`,
  `margin_controls.py` and `independence.py` carry no underscore because each
  is one public thing, while `core/_vinecop_discrete.py`, `core/_margins.py`,
  `core/_placement.py` and `core/_covariates.py` keep theirs
  even though `as_margin`, `resolve_margin_controls`, the
  `margin_*_json` helpers and the five `core.extend` names they hold are
  public -- the
  internal layout helpers, the two registry tables, the per-ecosystem
  predicates, the specification shapes, the two mixins and the two forwarding
  rules are the bulk of those files, and the public names are reached
  through `core` or `margins`.
  Do not resolve a mismatch here by renaming a mixed module; resolve it by
  asking whether the module is something to import from.
- **Coercion** — `as_margin(obj)` is idempotent and routes **every**
  margin `Vinedist` receives, so a discrete SciPy object cannot slip
  past on a bare `pdf` (in SciPy's new API `pdf` is `+∞` at an atom;
  the mass is `pmf`). **`core` holds the two registries and names no
  ecosystem.** An adapter or a JSON reader is registered by the module that
  owns the class it produces — `margins/scipy.py`, `margins/openturns.py`,
  `torch/kde1d.py`, `torch/distribution_margin.py` — through the same
  `register_margin_adapter` / `register_margin_json` hooks a third party uses,
  so the first-party margins exercise the documented extension point rather
  than a private table beside it. Putting those tables in `core` instead is
  what forced `core` to name a class from every extra, and a `core` -> `torch`
  edge that `tests/test_import_surface.py` refuses outright. The one exception
  is `Kde1d`, which `core` owns and can therefore name.
- **Which class a variable gets is the distribution's `margin_class`**, not a
  per-call argument. `#339` deleted `margins=` from `VinedistBase` and from the
  estimators: one class per distribution is what `vinecop_class` already meant
  on the copula half, and a caller who wants another margin subclasses or
  composes rather than re-specifying it at every call. `_default_margins` is
  the hook a lane overrides where its margins need something no controls object
  carries — on the torch lane, the device and dtype the *copula* controls name.
  Margins still follow the library's **construct-then-`fit`** pattern (`fit`
  returns `self`), so one class is both the specification and the fitted
  object, and a margin with no `fit` of its own is *fixed*: the loop leaves it
  as it is instead of substituting one it can estimate.
- **Configuration** — `FitControlsMargin` is the marginal half of a `Vinedist`
  fit, and `resolve_margin_controls` expands `margin_controls=` four ways: one
  object broadcast per variable, a length-`d` sequence, or a mapping keyed by
  position or by name. It carries `family_set` / `selection_criterion` /
  `on_failure` and nothing else: **the declaration is not fit configuration**.
  What a variable *is* — its `var_type` and its `support` — travels keyword-only
  beside the controls, exactly as `var_types` does on `Bicop.from_data`, one
  entry per variable on `VinedistBase.from_data(y, controls, *, var_types=,
  supports=)`. That split is what retired `MarginBase.declare` and the three
  helpers that existed to move a declaration out of a controls object and into
  a constructor. A declaration is still a **default**, not an instruction — a
  margin the caller constructed keeps what it was built with — and it is what
  makes a bounded `Kde1d` reachable without naming a class. A margin that
  cannot honor a `family_set` **refuses** it rather than fitting one family and
  looking like it chose; whether controls are forwarded at all is the declared
  `controls_class`, because nanobind reports every bound signature as
  `(*args, **kwargs)` and introspection cannot answer it.
- **The kernel knobs are `FitControlsKde1d`**, a separate type rather than
  fields on `FitControlsMargin`: `multiplier` / `bandwidth` / `degree` /
  `grid_size` / `boundary_repair` configure a kernel density, and a family
  search has no use for any of them. `ControlsLike` already unifies the two --
  it asks only for `to_dict()`, which is how one signature accepts both
  `FitControlsVinecop` and `FitControlsTorchVinecop` -- so there is no
  hierarchy to build and no field either type has to ignore.

Conventions that bind: the fit is **two-step (IFM)** — margins first,
then the copula on the resulting pseudo-observations — never fit all of
SciPy (an unfiltered sweep ranks `vonmises` above the true `gamma` because its
reported support lies), and never silently skip a failed candidate: every
rejection is reported with its reason, and a column where everything fails
**raises**, naming each family and its cause. `on_failure="fallback"`
substitutes `Kde1d` with one warning instead -- available, but not the
default, because answering a parametric request nonparametrically is the
same class of silent downgrade the weights contract already refuses. The
substitution happens in `VinedistBase._fit_margin`, not in the margin: a
`SciPyMargin` would have to stop being parametric to make it, so the decision
belongs to whatever chooses which margin a column gets.

There is **no structured selection report**. The copula layer's
answer to the same question is `show_trace` printing to stdout, and margins
inventing a second, structured mechanism is what made the two layers
asymmetric in the first place. When diagnostics are designed, design both
layers at once.

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

- Re-exports `to_pseudo_obs`, `wdm`, `find_latent_sample`,
  `sobol`, `ghalton`, `sample_uniform` (all C++) plus the
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
  `needed_hfunc1` / `needed_hfunc2`) raise on a tree above `trunc_lvl`, where
  indexing the triangular array without a bounds check **segfaulted**. Fixed
  upstream in vinecopulib#756 and pinned since; the binding only translates the
  exception.
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

#### Which distribution an estimator fits

Estimators take **`distribution=`**, a `VinedistBase` subclass, and fit it:
`Vinedist` by default, `pyvinecopulib.torch.TorchVinedist` for the PyTorch
lane. That one class names both halves of the model, so naming it is the
whole lane choice -- `vinecop_class` says which vine is fitted, `margin_class`
and `_default_margins` say which margin an unaddressed column gets, and the
class itself is what `distribution_` ends up being. The copula's fit
configuration is `controls=` and a pre-specified structure is `structure=`;
both are stored verbatim and copied before the estimator writes to them.

There is **no backend layer**, and the word "backend" is not used anywhere in
this repository. `#218`'s `VinecopBackend` / `TorchVinecopBackend` /
`resolve_backend` / `_BackendVinecop` and the copy-on-write `with_*`
derivations are all gone: once `VinecopLike` declared `structure` and both
lanes accepted `num_threads`, every one of them restated something the
distribution classes already declare. `sklearn.clone` deep-copies a
non-estimator parameter, so it supplies the copy-on-write the `with_*`
builders were written for, and `clone(est).set_params(structure=...)` is the
ensembling idiom.

Two consequences worth knowing. **`sklearn` imports `torch` nowhere at all** --
the opt-in is the caller importing `TorchVinedist` to name it, which is why
the parameter takes a class rather than a string, and why
`tests/test_import_surface.py` declares no `sklearn -> torch` edge. And
**evaluation is the vine's own**: estimators call `self._vine.pdf(u,
num_threads=...)` directly and convert at their public boundary with
`core.extend.to_numpy(..., dtype=float)`, because `np.asarray` alone raises on
a tensor that requires grad or lives on an accelerator.

`TorchVinedist` lifts a `Kde1d` margin to `TorchKde1d` rather than refusing it,
while still refusing anything else that is no `nn.Module`: the lift is exact
and core's own class, whereas a SciPy margin would detach the graph silently.
It is what lets the `on_failure="fallback"` substitution land on this lane.

#### Estimator conventions (scikit-learn developer guide)

`#218` aligned every estimator with the
[scikit-learn third-party-estimator developer guide](https://scikit-learn.org/stable/developers/develop.html).
The rules that bind:

- `__init__` performs **no** validation — it stores arguments
  verbatim. All checks live in `_parameter_constraints` +
  `_validate_params()`, called at the top of `fit()`.
- The lane parameter is `distribution=` and takes a **class**, never a
  string shortcut: naming `TorchVinedist` is what imports torch, which a
  string could not be.
- Fitted attributes follow sklearn naming with a trailing underscore:
  `feature_names_in_`, `schema_`, `random_state_` (resolved RNG),
  `distribution_class_` and `controls_` (both resolved at fit time).
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

- `TorchTllBicop` / `TorchVinecop` — evaluators.
  - `TorchTllBicop` is a density on a grid; constructors:
    `TorchTllBicop(grid_points, values, cache_integrals=True, ...)`,
    `TorchTllBicop.from_bicop(cop, ...)` (lift a C++ `Bicop`),
    `TorchTllBicop.from_data(u, controls=None, ...)` (fit the grid).
  - `TorchVinecop` mirrors `pv.Vinecop`'s `pdf` / `cdf` /
    `rosenblatt` / `inverse_rosenblatt` / `sample` signatures.
- `TorchDistributionMargin` / `TorchVinedist` — the marginal and joint halves.
  - `TorchDistributionMargin` is a `MarginBase[Tensor]` that is *also* an
    `nn.Module`: `torch.distributions.Distribution` has no
    `.to(device)` and contributes nothing to `state_dict` as a plain
    attribute, so the parameters are registered and the distribution is
    **rebuilt per call** — the same shape `TorchTllBicop` uses for its
    grid. `TorchDistributionMargin.from_distribution(factory, parameters=...)` is
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
    `grid_points` / `values` / `type` / `prob0` and the declared bounds are the
    whole of what the compiled `pdf` / `cdf` / `icdf` read -- the bounds joined
    that list when kde1d#37 made them the discrete support. The grid is a **buffer**, not a
    parameter — the density is fitted, not learned — so optimizing it is the
    opt-in `values.requires_grad_(True)`, the `TorchTllBicop` precedent.
    `icdf` reproduces the C++ inversion exactly — a bracketed Newton within
    the cell holding the requested mass, bisecting where the density is flat,
    with the C++ early exit reproduced as a frozen-once-converged mask — and
    then reattaches an exact gradient by one Newton step, so the correction
    does not move the value while `dq/dp` and `dq/d values` are right. The
    residual is written in units of mass, so the total mass carries its share
    of `dq/d values`. This is the one parity claim in the port that is a
    tolerance rather than an equality, and has to be: the compiled
    quantile is not portable to a few ULPs -- rebuilding kde1d with
    `-march=native` alone moves it 19 -- so no port can equal every build of
    it. The correction is conditional on
    grad being enabled, not on the grid being learned — a fitted fixed grid
    still has to differentiate the quantile in `p`. One of `Kde1d`'s names could not
    be reused: `loglik` is the contract's method, so the fitted value is read
    as `loglik()`. `type` is `nn.Module`'s dtype cast, which is no longer a
    collision -- the declaration is `var_type` on both classes.
  - `torch/_margin_kde1d_interp.py` ports `kde1d`'s `InterpolationGrid`, and its
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
    `TorchTllBicop.cache_integrals` and for a stated reason — here the cached
    quantity would be an `O(m)` vector shared by the batch, not an `O(m^2)`
    integral per query.
- Discrete variables are declared with `var_types` on `TorchVinecop`'s three
  constructors. The stored pair copulas stay continuous interpolation grids and
  `get_pair_copula` hands out `with_var_types(...)` of one, so `state_dict` /
  `.to()` / pickling see only real `nn.Module` parameters and a continuous
  slot copies nothing. `TorchTllBicop.from_data`
  takes the four-column layout and reuses the compiled `find_latent_sample`,
  which is what `TllBicop::fit` now consumes for a discrete edge; the jittered
  ranks only seed the bandwidth. A discrete torch vine refuses the **batched
  fast path**, whose stacked levels gather `(N, n, 2)` where a discrete edge
  needs `(N, n, 4)` and the parallel left-limit cascade behind it, and whose
  `DELTA_MIN` split is per row. It does *not* refuse the **integral cache**: the prefix tables
  reconstruct the integral exactly, so a discrete edge can difference them and
  `cache_integrals` resolves the same way it does for a continuous vine.
- `FitControlsTorchBicop` / `FitControlsTorchVinecop` — fit-time
  dataclasses. Notable knobs:
  - `compile_fit` — off by default; fuses the bandwidth search's per-pass
    body with `torch.compile`. The pass is 39 launches over tensors whose
    arithmetic stays invisible even at `n = 30000`, so fusing it is worth
    2.0x end to end at `n = 500`, tapering to 1.04x at 30000 and never
    below 1. Off because the first call compiles for ~5 s and a second lane
    count for ~8 more — about 500 pair fits to break even, so a loss for one
    vine and a win for a few dozen. Compiled with `dynamic=True`
    on purpose: a `d = 20` vine presents 19 lane counts, static shapes
    would exceed `cache_size_limit` (8) and drop the widest levels back to
    eager in silence, and raising that limit from a library would reach
    every other compiled function in the process.
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
  - `rect_prob` / `cond_interval_prob` — the two `BicopBase` hooks, overridden
    onto `InterpolationGrid2D.rect_mass` / `cond_interval_mass`, and available
    in **both** cache modes: they read the density grid, not the prefix tables.
    A four-corner `cdf` difference turns an absolute error `ε` into
    `≈4ε/(w₁w₂)` in the atom widths; the rectangle amplifies by `1/w₂` alone,
    since only its `λ(b₂) − λ(a₂)` term cancels and that multiplies a term of
    order `w₁`. Measured on a `1.2e-4`-wide rectangle: `2.9e-12` against
    `8.7e-9`. That bound is earned only by expanding the `λ` difference over
    the common denominator, so it cancels against `b₂ − a₂` rather than against
    one — formed as `λ(b₂) − λ(a₂)` it swamps the first term on a low-mass
    rectangle, and the total is formed additively as `m_below + m_strip` for
    the same reason. `values >= 0` is a constructor precondition precisely
    because the nonnegative-weight bound depends on it. Note it is the
    **probability**, not the density's mass: `cdf` renormalizes each grid line
    by its own total, so the two differ, and a discrete edge is defined against
    the distribution function. The conditional one is **not** clamped into the
    open unit interval, unlike an h-function value, which is what makes the
    masses of a partition sum to one.
  - `compile` — runs the batched cascades through `torch.compile`, on CUDA
    replayed as a CUDA graph. Off by default: the first call at each input
    shape pays tens of seconds of Inductor, so it is worth it for a cascade
    called repeatedly and not for a single evaluation. Torch caps how many
    variants of one code object it will hold (`cache_size_limit`, 8), and a
    vine is a variant — so a process compiling more vines than that falls
    back to eager. The `batched` flag is **not** a control: it is resolved
    per device on each call, and overridable per call.
  - `batched_fit` — fits a whole tree level in one call instead of edge at a
    time, through the optional `fit_level` hook on `VinecopBase.fit` /
    `.select` (`TorchTllBicop.from_data_batched` is the pair-level entry point,
    taking `(P, n, 2)`). Resolved per device like the cascade's `batched`.
    The hook does not know what the pairs are for — `P`
    independent pairs on shared rows — so several vines' levels concatenate
    into the same axis as readily as one vine's. A level carrying a discrete
    edge or a conditioning context cannot stack and stays per-edge.
  - `device`, `dtype` — propagate to every tensor on construction;
    fitted modules respect `.to(device)` afterwards.
  - `trunc_lvl`, `tree_criterion`, `threshold`, `tree_algorithm`, `seeds`
    — `FitControlsTorchVinecop` only; the structure-selection knobs used
    when `TorchVinecop.from_data` is called with `structure=None`.
    Selection runs through `VinecopBase.select`, so it stays on the array
    namespace rather than round-tripping a compiled `pv.Vinecop`.
- `InterpolationGrid2D` (`torch/_bicop_interp.py`) — the 2-d bilinear grid
  backing `TorchTllBicop`; **internal** (not re-exported). Margin
  normalization uses Sinkhorn iterations to drive marginals to uniform.

### Top-level `pyvinecopulib`

`src/pyvinecopulib/__init__.py` is *not* empty: it re-exports the core
class surface so the long-standing
`from pyvinecopulib import Bicop, Vinecop, to_pseudo_obs` pattern keeps
working. Specifically:

- `__all__` is **16** names: the five copula classes (`Bicop`, `Vinecop` and
  the three `*VineStructure`s), `BicopFamily`, the three `FitControls*`,
  `Vinedist`, `to_pseudo_obs`, `__version__`, and the four subpackages that
  import with no extra (`core`, `families`, `margins`, `utils`). `Kde1d` is
  *not* among them -- it is reached as `pyvinecopulib.core.Kde1d`, or through
  the deprecation shim, never by star-import. The two that need one -- `sklearn`, `torch` -- are
  **out** of it: `from pyvinecopulib import *` resolves every name in
  `__all__`, so listing them made both extras hard requirements of the one
  import form a beginner reaches for first. Both stay reachable by attribute
  access and by `import pyvinecopulib.<name>`, and `__dir__` still names
  them, because discovery and star-binding are different questions.
- `__getattr__` provides two things: lazy import of `sklearn` (the
  extra is only triggered on `import pyvinecopulib.sklearn` or
  attribute access) and a deprecation shim for the 35 pre-#207
  top-level names (every family constant, every utility function from
  `utils`). Each access emits a `DeprecationWarning` pointing at the
  new canonical path.
- For new code, **always import from the canonical subpackage**
  (`from pyvinecopulib.families import gaussian`,
  `from pyvinecopulib.core import Kde1d`), not the top-level alias.

`pyvinecopulib.torch` and `pyvinecopulib.sklearn` resolve as attributes and
as submodule imports, but neither is in `__all__` (see above).

### Internal: `_deprecations`

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

- **`pyvinecopulib.core`** — `Bicop`, `Vinecop`, `Kde1d`, `RVineStructure`,
  `CVineStructure`, `DVineStructure`, `BicopFamily`, `FitControlsBicop`,
  `FitControlsVinecop`, `FitControlsMargin`, `FitControlsKde1d`; plus the
  array-agnostic
  abstraction layer
  `BicopLike`, `VinecopLike`, `BicopBase`, `VinecopBase`, `ControlsLike`,
  `IndependenceBicop`, `ConditioningContext`,
  `SimplifiedContext`, `NonSimplifiedContext`; plus the marginal layer
  `MarginLike`, `MarginBase` and the joint object with its contract and base,
  `Vinedist`, `VinedistLike`, `VinedistBase`; plus the margin serialization
  helpers `margin_from_json`, `margin_to_json`, `register_margin_json`; plus
  `ArrayT`, the type variable those signatures are written in.
- **`pyvinecopulib.core.extend`** — **six** names, and the count is the
  point: what an extension cannot be written without, rather than everything
  it might find useful. `place` and `reference_array`, the fix and the check
  of the `_prep` override; `prepare_covariates`, the composite every entry
  point taking an `x` routes through, with `covariate_row` for the single-row
  reading of one; `to_numpy`, the return trip; and `NotBatchable`, which a
  `_build_batched` override raises to decline the grid fast path. Reached as
  `pyvinecopulib.core.extend`, and kept out of `core`'s own namespace
  because using pyvinecopulib needs
  none of it. The four canonical bases and their protocols stay in `core`:
  `README.md` tells users to subclass them, so they are part of the surface
  rather than machinery behind it.
    - **An export is a promise, so the burden of proof is on keeping a name.**
      It earns one by being needed to implement a documented hook *correctly* —
      not because internal code calls it, and not because a private hook's
      docstring names it: none of `_prep`, `_prep_args`, `_layout`,
      `_build_batched` or `_coerce_fit_data` is rendered by Sphinx, so prose
      inside one is source-only. `trim` and `collapse_data` failed that test
      outright — `_prep_args` and `_layout` are overridden **nowhere** in
      `src/` and appear in no recipe below — and so did the three validators,
      `continuous_view`, the two fit-callback aliases and the model codec. All
      are private, and a subclass writes its own refusal message or a plain
      `def` instead. Note the aliases are named without backticks on purpose:
      the test below reads a backticked capitalized name here as a claim that
      the module exports it.
- **`pyvinecopulib.families`** — `BicopFamily` enum; per-family
  constants (`indep`, `gaussian`, `student`, `clayton`, `gumbel`,
  `frank`, `joe`, `bb1`, `bb6`, `bb7`, `bb8`, `tawn`, `tll`); group
  lists (`all`, `parametric`, `nonparametric`, `one_par`, `two_par`,
  `three_par`, `elliptical`, `archimedean`, `extreme_value`, `bb`,
  `rotationless`, `lt`, `ut`, `itau`, `analytic_derivs`).
- **`pyvinecopulib.utils`** — `to_pseudo_obs`, `wdm`,
  `find_latent_sample`, `sobol`, `ghalton`, `sample_uniform`,
  `pairs_copula_data`.
- **`pyvinecopulib.margins`** — `SciPyMargin`, `OpenTURNSMargin`,
  `FitControlsMargin`, `as_margin`, `register_margin_adapter`,
  `resolve_margin_controls`.
- **`pyvinecopulib.sklearn`** — `VineDensity`, `VineRegressor`.
- **`pyvinecopulib.torch`** — `TorchTllBicop`, `TorchVinecop`, `TorchKde1d`,
  `TorchDistributionMargin`, `TorchVinedist`, `FitControlsTorchBicop`,
  `FitControlsTorchVinecop`.

Top-level `pyvinecopulib` re-exports the ten classes named above and
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
- `batched_fit=True` ↔ `batched_fit=False`: numerically equivalent, and
  **not** bit-identical on any device — unlike `inverse_rosenblatt` above,
  which is. Batching changes how many elements the bandwidth search's `pow`
  is handed, and torch selects elementwise kernels by element count
  (vectorized past `2 * Vectorized<double>::size()`: 8 on AVX2, 4 on NEON),
  so a machine that vectorizes sooner diverges where another does not. Do
  not pin a fit comparison at `atol=rtol=0` on the strength of one machine
  agreeing; the exact claims available are that a lane's answer is
  independent of *which* lanes it travelled with (at a fixed shape) and that
  the selected structure matches, the tree criterion reading ranks rather
  than last bits.
- `sklearn.base.clone()` round-trip: every estimator clones cleanly
  with all `__init__` parameters preserved verbatim.

## Extension points

- **New copula families.** Add the C++ implementation upstream in
  `lib/vinecopulib`, bump the submodule, then extend
  `src/include/bicop/family.hpp` + `src/pyvinecopulib/families/__init__.py`
  to bind the new tag and group it appropriately.
- **Custom pair copulas / vines (`pyvinecopulib.core`).** Subclass
  `BicopBase` (define `_pdf_raw` / `_hfunc1_raw` / `_hfunc2_raw`) for a custom
  pair copula, and host it by subclassing `VinecopBase` (define the one hook
  `get_pair_copula`); both run on NumPy or PyTorch and inherit the full
  evaluation surface. Implement `BicopLike` / `VinecopLike` directly for
  an immutable / functional part. To put that pair on a **discrete**
  edge, add a `_cdf_raw` and return `pair.with_var_types(self.pair_var_types(t,
  e))` from `get_pair_copula`; the vine supplies the left-limit columns and
  the fit engines declare the types on each pair as they fit it, so a
  `fit_edge` callback does not have to. A pair implementing `BicopLike`
  directly carries no types and cannot sit on a discrete edge. For a
  **non-simplified / conditional** vine, pass a `NonSimplifiedContext` and drive
  `VinecopBase.fit` with a `fit_edge` callback. To condition on a subset of
  variables, implement `flip` as well and use `sample_conditional` /
  `select(conditioning_set=)`; see
  `examples/10_extending_pyvinecopulib.ipynb`. `TorchTllBicop` /
  `TorchVinecop` are the reference torch subclasses.

  Four things a subclass author hits that the signatures do not say, each
  reported by a third-party subclass built against this contract:

  - **Write leaves, never override a dispatcher.** The public members carry
    no per-call keywords, so a subclass with its own knobs is tempted to
    override `pdf` or `cdf` to keep one. That silently orphans the leaf:
    every internal caller goes through the public member (`loglik` ->
    `self.pdf`, the plots through `pair_eval`, the vine cascade likewise), so
    the mandatory `_pdf_raw` exists only to satisfy the ABC and never runs.
    Put evaluation knobs on the `ControlsLike` instead.
  - **Any pair that may sit in a conditional vine declares `x` on every
    leaf**, whether or not it reads one. `pair_eval` forwards unconditionally
    -- `method(*args) if x is None else method(*args, x=x)`, no introspection
    -- so a leaf without `x` raises `TypeError` at the first covariate call.
    `IndependenceBicop`'s `del x` is what that looks like.
  - **`_prep_args` is about clamp width, not rejection.** Overriding it on a
    *pair* installs no refusal on the cascade: `VinecopBase._prep_args` already
    trims, and every public vine entry point calls it before the cascade
    builds anything, so a pair only ever sees a clamped interior value.
    Override it where the subclass's domain step differs from the library's --
    `trim_bounds` is ~1e-10 at float64, which is about ten units of feature
    space for a model fitted on `logit(u)`. Test the hook, not a density: a
    density only notices where it varies.
  - **`_hinv1_raw` / `_hinv2_raw` are optional, and the default is expensive.**
    The inherited bisection of `_hfunc1_raw` costs `O(iterations)` in *model*
    calls, not in arithmetic, so a pair backed by a learned model or one that
    inverts in closed form should override both.
- **Fitting has one shape across all four bases.** `MarginBase`, `BicopBase`,
  `VinecopBase` and `VinedistBase` each expose `fit(...) -> Self` (mutates in
  place, returns the object) and `from_data(...) -> cls` (constructs one), and
  `BicopBase` / `VinecopBase` add `select(...) -> Self` where there is
  something to select — a family, a structure. `Bicop` / `Vinecop` return the
  object from `fit` / `select` too, so the idiom is the same whichever class
  you hold. Configuration travels as a `ControlsLike` — anything with
  `to_dict()`, which is how one signature accepts both `FitControlsVinecop`
  and `FitControlsTorchVinecop` — while data and callbacks stay explicit
  arguments.
- **`fit_edge` is keyword-only, and `VinecopBase.from_data` is the plain
  factory.** Keyword-only because a subclass that ships its own pair fitter
  (`TorchVinecop`) must be able to override the factory *compatibly*, and a
  required second positional would make that impossible; `ty` enforces it.
  For the same reason `from_data` carries only `fit_edge` / `controls` /
  `structure` / `var_types`: a **non-simplified** or covariate-driven fit is
  built the other way round — construct the vine with its
  `ConditioningContext`, then call `fit(u, fit_edge=..., x=...)`. `VinecopBase`
  keeps the array-agnostic engines behind `_fit_parts` / `_select_parts`, which
  return the loose parts a factory assembles, because `from_data` needs them
  before an object exists. Those two names are `staticmethod` bindings of
  `core/_vinecop_fit_engines.py`'s `fit_parts` / `select_parts`: both are module
  functions
  — neither reads `self` or `cls` — and they sat in the class body only for
  namespacing, 718 lines of a 2677-line class. The class keeps the names
  because that is what external drivers reach them through.
- **Declare the parts, inherit the fitting.** A base knows how to fit once it
  knows what its parts *are*: `VinecopBase.bicop_class` names the pair copula
  it fits, and `VinedistBase.vinecop_class` / `margin_class` name its two
  halves. Naming them is what makes `from_data` work with no callback — a part
  class is itself a fitter, since `fit` / `from_data` exist on every base — and
  what lets `VinecopBase.select` refuse a pair copula without `flip` *before*
  it reads the data. `None` means the object only ever hosts parts it is
  handed, and fitting then needs an explicit `fit_edge`.

  **`controls_class` is the fourth declaration of the same kind**, on all four
  bases and on the bound `Bicop` / `Vinecop` / `Kde1d` (set from the binding,
  the way `Kde1d` sets `supports_weights`): the fit configuration this class
  reads, or `None` where it reads none. It answers *what does `controls=None`
  mean here* once per class, so a consumer writes
  `cls.vinecop_class.controls_class()` instead of naming a `FitControls*` of
  its own — which seven sites used to do, three of them in one file. On a
  margin it also replaces a boolean: `controls_class is None` **is** "reads no
  controls", and carries the type besides, so `supports_controls` is gone.

  It is a declaration, not a contract. It is read with `getattr` like every
  other optional capability, and it is **not** on the protocols: those are the
  evaluation surface, and `fit` / `select` / `from_data` are not on them
  either. `ControlsLike` stays what it is — anything with `to_dict()` — so a
  consumer still reads the settings it owns and refuses what it cannot honor.
- **`get_pair_copula` reads, `set_pair_copulas` writes, and only the first is
  abstract.** Reading is required to *evaluate*, which every subclass does;
  writing is required only to *fit*, so a vine that merely hosts pairs — or an
  immutable one — stays valid without it, and `set_pair_copulas` reports its
  own absence instead. Same rule as `MarginBase.fit`. Neither has an
  underscore-prefixed twin: a public wrapper over a private hook bought
  nothing and is gone. What `set_pair_copulas` *does* owe is invalidation —
  it is the one place the pairs change without the structure changing — and
  `_invalidate_batched` is the name for it, so an implementation writes a call
  rather than an assignment to `_batched`. A lane holding further derived state
  overrides that hook and calls `super()`, which is how `TorchVinecop` clears
  the compiled cascades alongside the stacked grids.
- **A vine's controls *are* pair controls.** `FitControlsVinecop` derives from
  `FitControlsBicop` — the binding declares the C++ inheritance, as it does for
  `CVineStructure` — and `FitControlsTorchVinecop` from
  `FitControlsTorchBicop`. So one controls object configures both halves of a
  vine fit: a vine reads the settings it owns and the rest reach its pair
  copulas unchanged, with no nested object and no accessor. That is also how
  observation weights get to both halves.

- **Custom margins (`pyvinecopulib.core`).** Subclass `MarginBase` and
  define `pdf` / `cdf`; `icdf`, `logpdf`, `cdf_left`, `loglik`,
  `sample`, `plot`, `support`, `nobs` and `n_parameters` come with it --
  a fit records the first two into the `_nobs` / `_n_free` slots the base
  owns, so the criteria work without a subclass restating them. Add
  `fit(y, controls=None, *, var_type=None, support=None, x=None,
  weights=None) -> Self` to make it an estimator, or leave it out for a fixed
  margin — whether the class overrides `MarginBase.fit` is what a vine
  distribution's margin loop dispatches on, `is_fitted` being the wrong test
  (a fitted margin that *is* refittable must be re-estimated by `fit` and
  `select`). Override
  `cdf_left` whenever the family has an exact left limit (Poisson's
  `gammaincc(k, μ)`, a categorical's `cumsum(probs)[k-1]`): the derived
  `cdf(x) - pdf(x)` cancels in the right tail and `cdf(x - 1)` is
  meaningless off an integer lattice.
- **Custom vine distributions (`pyvinecopulib.core`).** Subclass
  `VinedistBase` and you have the whole data-scale surface immediately —
  evaluation needs no hook, because a vine distribution is determined by its
  two halves. Making it *fittable* is mostly a declaration: name
  `vinecop_class` and `margin_class` and `from_data` runs the two-step (IFM)
  estimator itself, in the base, leaving `_coerce_fit_data` the only real hook
  because the torch lane resolves a device and dtype before any part exists.
  **Deserialization is the same declaration.** `from_json` decodes, checks the
  payload's version and checks that its `kind` names this class, then hands a
  plain mapping to `_from_payload`, whose default rebuilds `vinecop_class`
  from its own `from_json` and every margin through the registry -- so
  `Vinedist` overrides nothing. Only `_from_payload` is a hook, and only for a
  copula that reads no JSON of its own, which is what `TorchVinedist` reports
  by naming `TorchVinecop`. Keeping the decode in the base is what lets the
  codec stay private: the payload is not plain JSON -- a non-finite float
  travels as a string, JSON having no spelling for one -- so an override
  reaching for the standard library's `json.loads` would read `-inf` back as
  that string, silently. `_json.read_payload` is the one place those three
  checks live, shared with `margin_from_json`.
  The one hook a lane normally overrides is `_copula_controls`, the single
  lane-specific step in the copula estimate: `Vinedist` writes `weights` into a
  copy of the controls there and `TorchVinedist` pins the device and dtype the
  margins resolved. It is read from **both** copula paths, which is what keeps
  them configured identically. `supports_weighted_copula` is `False` on the
  base for that reason: the inherited `_copula_controls` cannot weight the
  copula, so a lane declares the capability together with the override that
  honors it. Declaring one a lane cannot honor is what produces a half-applied
  fit, so the request is refused up front instead. The margin half has one
  hook of the same kind, `_adopt_margin`: the identity on the base, and on the
  torch lane the lift a `Kde1d` gets before it can join the others. It is
  reached only where the base *substitutes* a margin mid-fit
  (`on_failure="fallback"`), which is before `_bind_dist` runs and therefore
  before the lift that happens there -- and the copula-scale layout is
  assembled from the columns the margins return, so a substitute on another
  namespace would not stack with the rest.
- **The margin half fits the same way the copula half does.**
  `_default_margins` builds one margin per variable -- the hook a lane
  overrides to place them -- and `_fit_margin` calls that margin's own `fit` or
  `select` as `estimator(y, controls, ...)`, beside `_reestimate_copula`'s
  `estimator(u, resolved, ...)`. What it adds around that call is the three
  decisions a margin cannot make about itself: whether it is *fixed* (it
  overrides neither verb, so there is nothing to re-estimate), whether a
  request it cannot honor was made (`weights` it cannot use, a `family_set` it
  cannot search), and which margin the column gets when its own estimator
  raises. Two tempting shortcuts for the first are both wrong: `is_fitted`,
  because a fitted `SciPyMargin` must be re-estimated, and `fit` alone,
  because a margin that chooses between *kinds* of model overrides `select`
  and leaves `fit` raising. Whether the verb asked for is one the margin has
  is the verb's own to report.
- **`from_data` constructs the copula; `fit` and `select` re-estimate the one
  already held.** `_fit_copula` builds `vinecop_class` and is the *construction*
  path only; `fit` / `select` go through `_reestimate_copula`, which calls the
  held copula's own `fit` / `select`. That is what makes a hosted
  `VinecopLike` — a `VinecopBase` subclass a caller composed the distribution
  from — keep its class, its identity and its pair-copula types across a refit,
  exactly as a margin keeps its family. Building a fresh `vinecop_class` in
  `fit` silently replaced the caller's vine with the default one, which is the
  defect that established the split; a copula with no `fit` now reports that
  instead of being swapped. One thing neither verb may change is **which
  variables have atoms**: the held copula carries its own `var_types` and has
  no way to be told new ones, so a margin that gained or lost atoms in the
  refit is refused by name, pointing at `from_data`, rather than reaching the
  copula as a layout of the wrong width.
- **`supports_fit_covariates` is a lane-level "anything at all", not "both
  halves".** A conditional `Vinedist` is one whose *margins* read `x`: the
  compiled `Vinecop` models no covariates and takes no `x` argument, so the
  copula half is never conditional there and `x` reaches only the margins that
  declare it — the same per-part rule `declared_eval` applies at evaluation,
  and `_fit_copula` forwards to a copula class only when *it* declares
  `supports_covariates`. What enforces that is the object-level refusal:
  the flag says whether anything on the lane is fitted on covariates, and when
  nothing is, the request is refused rather than answered unconditionally.
  `Vinedist` (NumPy + compiled `Vinecop`) and `TorchVinedist` are the two
  reference subclasses; implement `VinedistLike` directly for an immutable /
  functional distribution.

- **Another ecosystem's distributions (`pyvinecopulib.margins`).** Call
  `register_margin_adapter(predicate, adapter)` — from a package or a
  notebook cell — rather than editing `core/_margins.py`. That is what keeps
  OpenTURNS, TFP and NumPyro out of `core` while remaining usable, and
  what lets `as_margin` stay the single funnel every margin passes
  through.
- **A new lane for the sklearn estimators.** There is nothing to subclass:
  write a `VinedistBase` subclass and hand it to an estimator as
  `distribution=`. Naming `vinecop_class` is what makes its vine fittable,
  `margin_class` (or a `_default_margins` override, where the margins need a
  placement) is what an unaddressed column gets, and `_copula_controls` is
  where a lane pins anything the margins resolved. The estimators add
  nothing of their own: they call `vinecop_class.from_data(...)` and then
  construct the distribution.
- **A different torch pair-copula fitter.** There is no method registry to
  add to any more -- `FitControlsTorchBicop.method` was a one-value enum
  nothing dispatched on, and it is gone. `TorchTllBicop.from_data` fits a TLL
  grid, full stop. A different estimator is a different **class**: subclass
  `BicopBase` **and** `torch.nn.Module` -- `BicopBase` is array-agnostic, so
  the second base is what supplies `state_dict` / `.to(device)` and lets a
  vine own the pair as a child, and its `__init__` has to be called
  explicitly, as `TorchTllBicop` does -- give it `fit` / `from_data`, and host
  it by naming it as a vine's `bicop_class` or by passing `fit_edge`. `torch/_bicop_fit_tll.py` is the reference for the kernel
  itself, and `examples/10_extending_pyvinecopulib.ipynb` for the hosting.
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
