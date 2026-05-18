# Changelog

## Unreleased

### `pyvinecopulib.sklearn`: backend system + scikit-learn-developer-guide cleanup

The sklearn module gains a thin public backend layer plus a focused
sweep of fixes against the
[scikit-learn third-party-estimator developer guide](https://scikit-learn.org/stable/developers/develop.html).
Both pieces land in the same PR so the rename touches each
constructor exactly once.

**Backend system** — new module `pyvinecopulib.sklearn.backends`:

- `VinecopBackend` (default) wraps `pyvinecopulib.Vinecop`; holds an
  optional `FitControlsVinecop` and an optional structure.
- `TorchVinecopBackend` wraps `pyvinecopulib.torch.TorchVinecop`; holds
  an optional `FitControlsTorchVinecop` (new — see torch entry below)
  and an optional structure. Constructing this class is the explicit
  opt-in signal that PyTorch is required.
- `VinecopLike` is a `runtime_checkable` Protocol describing the
  shared post-fit surface — `pdf(u, num_threads=, ...)`,
  `cdf(u, N, num_threads, seeds, ...)`,
  `simulate(n, qrng, num_threads, seeds, ...)`, plus a `structure`
  attribute. After PR1's `TorchVinecop` API alignment both
  `pv.Vinecop` and `pv.torch.TorchVinecop` satisfy the Protocol
  structurally.
- Estimators take a single `backend=` keyword; no string shortcuts
  (`backend="cpp"` / `"torch"` are no longer accepted — pass instances
  instead).

**Breaking sklearn-side API changes** (pre-release surface):

- `VineDensity(controls=..., structure=...)` →
  `VineDensity(backend=VinecopBackend(controls=..., structure=...))`.
- `VineRegressor` gains a real `normalize_weights=True` `__init__`
  parameter; the previous post-init `_normalize_weights` attribute is
  gone. Forests pass `normalize_weights=False` through `base_params`,
  and `sklearn.base.clone()` now preserves it.
- `seed`-style kwargs renamed to canonical `random_state` across the
  forest classes and on `VineDensity.sample` / `VineDensity.cdf` /
  `VineForestDensity.cdf`. The legacy `seeds=` and `seed=` kwargs
  are removed without a deprecation alias.

**scikit-learn-developer-guide fixes**:

- `__init__` now performs no validation on any estimator; per-class
  `_parameter_constraints` dicts + `_validate_params()` calls at the
  top of `fit()` replace the hand-rolled `isinstance`/`raise` blocks.
- DataFrame inputs now populate the canonical `feature_names_in_`
  (previously `_used_columns`).
- `schema_` (with trailing underscore) replaces `_schema` so
  `sklearn.base.clone()` round-trips correctly.
- `random_state_` (resolved RNG) and `backend_` (resolved backend
  pinned at fit time) are set as standard fitted attributes.
- Every post-fit method now calls
  `sklearn.utils.validation.check_is_fitted`.
- `base_params` is defensively copied at fit time so callers can't
  mutate the forest's view of it.

**Testing**:

- New `tests/test_sklearn_backends.py` exercises Protocol conformance,
  capability flags, `with_*` immutability, lazy torch isolation, clone
  round-trip, `feature_names_in_`, and cross-backend parity.
- New `tests/test_sklearn_check_estimator.py` runs
  `parametrize_with_checks` on every estimator. Genuine opt-outs
  (sparse inputs, 1-D / 1-feature degeneracies, complex inputs,
  density-log score in `Pipeline.score`, quantile-stacked regressor
  output) and known WIP items (memmap inputs, pickling round-trip,
  permutation invariance, `__sklearn_tags__` upgrades, …) are
  enumerated in the file as `expected_failed_checks` with one-line
  rationales each. The TODO-marked entries are intended to shrink in
  follow-up PRs.

**Documentation**:

- New Sphinx page `docs/concepts.rst` introduces Sklar's theorem,
  pair-copula construction / R-vines, and the default Transformed
  Local Likelihood (TLL) family in a ~5-minute read. Linked from the
  index toctree.
- The `pyvinecopulib.sklearn` module docstring now front-loads the
  "what does the default backend do" answer (no PyTorch required)
  and links to the concepts page.
- `pyvinecopulib.sklearn.backends` gains a "which backend should I
  pick?" subsection summarising the C++/torch trade-off (threading
  vs GPU vs autograd vs family-set coverage).
- `VineDensity` / `VineRegressor` class docstrings get short
  "use-the-torch-backend" examples and a "See also" block linking at
  the backend module and `pv.Vinecop`.
- `pyvinecopulib.torch` module docstring spells out TLL
  (*Transformed Local Likelihood*, Geenens 2014; Nagler 2018),
  motivates GPU / autograd / pipeline use-cases, and adds bibliographic
  references for VDC (Safaai 2026), the lazy-cascade design (Cheng
  et al. 2025), and TLL.
- `TorchBicop` / `TorchVinecop` class docstrings drop the
  `KernelBicop` aside in favour of plain-language descriptions and
  pick up reciprocal "See also" cross-refs at `pv.Bicop` / `pv.Vinecop`.
- `torch/_fit_vdc.py`: the IPFP / replicate-pad justification moves
  from the module docstring into an inline comment next to the
  function it actually documents; the module docstring is now a
  short user-facing summary anchored on the Safaai (2026) citation.
- README's "What are vine copulas?" expanded into a one-paragraph
  primer; a new "Optional backends" subsection introduces the
  sklearn and torch surfaces with one example each.

**Hygiene**:

- Removed `if TYPE_CHECKING:` guards in `torch/bicop.py` and
  `torch/_fit_vdc.py` — the gated import (`FitControlsTorchBicop`)
  was never a real cycle, so it now lives at module top.

### `pyvinecopulib.torch`: `FitControlsTorchVinecop`

`TorchVinecop.from_data` now accepts a
`controls=FitControlsTorchVinecop(...)` argument mirroring
`pv.Vinecop.from_data(controls=...)`. The dataclass bundles:

- a nested `FitControlsTorchBicop` (per-pair fit controls);
- vine-level placement / precision knobs (`cache_integrals`, `device`,
  `dtype`);
- runtime cascade knobs (`impl`, `batched`).

The previous inline kwargs on `from_data` (`grid_size`, `mult`,
`cache_integrals`, `grid_type`, `device`, `dtype`) keep working for
one cycle on `feature/torch-bicop` with a `DeprecationWarning`, then
will be removed before the next stable release.

### `pyvinecopulib.torch`: pluggable bicop fitters via `FitControlsTorchBicop`

`TorchBicop.from_data` now dispatches on a `FitControlsTorchBicop`
dataclass (mirroring `pv.FitControlsBicop`), opening the door to
alternative bicop fitters behind a single API. Two methods ship today:

- `method="tll"` (default) — the existing pure-torch TLL constant fit,
  unchanged in behaviour and still matching the C++ TLL fit to machine
  precision.
- `method="vdc"` — the pretrained amortized vine-copula estimator
  introduced by Safaai, H. (2026), *Amortized Vine Copulas for
  High-Dimensional Density and Information Estimation*,
  [arXiv:2604.20568](https://arxiv.org/abs/2604.20568). Reference
  implementation:
  [KempnerInstitute/vine-denoising-copula](https://github.com/KempnerInstitute/vine-denoising-copula).
  Not on PyPI yet; the `[vdc]` extra resolves it from GitHub via
  `[tool.uv.sources]` (`uv sync --extra vdc`). Plain pip users
  install with
  `pip install "vine-denoising-copula @ git+https://github.com/KempnerInstitute/vine-denoising-copula"`.
  The resulting `TorchBicop` reuses the standard interpolation-grid
  evaluation chain, so `pdf` / `cdf` / `hfunc` / `hinv` / `simulate`
  are identical to the TLL path.
  > **Upstream status**: vdc 0.1.0 on `main` ships an incomplete wheel —
  > `vdc.load_pretrained_model` references `vdc.inference` /
  > `vdc.vine.copula_diffusion` (not packaged). We ship a `sys.modules`
  > shim in `pyvinecopulib.torch._fit_vdc._install_upstream_shims` that
  > injects the two missing submodules with trivial stubs (the real
  > `scatter_to_hist` from `vdc.data.hist` is loaded directly via
  > `importlib.util`; `sample_density_grid` is stubbed since it's only
  > used in the diffusion-checkpoint path; `DiffusionCopulaModel` is
  > stubbed since it's import-referenced but never instantiated during
  > inference). The shim is installed at the first
  > `_load_bundle(...)` call and becomes a no-op once upstream restores
  > the missing subpackages.
  >
  > **Released-checkpoint accuracy caveat**: on the standard
  > Gaussian/Clayton precision bench (m=64, n ∈ {500, 2000, 10000}), the
  > `vdc-denoiser-m64-v1` checkpoint produces 10–20× worse pdf IAE than
  > the pure-torch TLL fit, with massive density spikes at the
  > anti-diagonal corners even for iid uniform samples (mean |pdf - 1| ≈
  > 0.9). The integration is correct and ready to use, but for parametric
  > targets TLL remains the better choice today.

The signature change is **breaking** for callers who passed `grid_size`,
`mult`, or `grid_type` as keyword arguments to `from_data`: those now
live on `FitControlsTorchBicop(...)`. `cache_integrals`, `device`, and
`dtype` remain direct keyword arguments on `from_data`.

`InterpolationGrid2D.normalize_margins` additionally accepts an optional
`tol` for early-stop (default `None` preserves byte-for-byte parity with
the C++ TLL pipeline).

### `pyvinecopulib.torch`: align `TorchVinecop` / `TorchBicop` with their `pv.Vinecop` / `pv.Bicop` counterparts

The torch evaluators now mirror the post-fit surface of the C++ classes
so downstream code (and the upcoming sklearn backend layer) can treat
either backend uniformly:

- `TorchVinecop.from_structure(structure | matrix, pair_copulas, var_types)`
  — new class method matching `pv.Vinecop.from_structure`. When
  `pair_copulas` is empty, every edge is populated with an independence
  `TorchBicop`, yielding the independence copula on `d` variables.
- `TorchVinecop.simulate(n, qrng=False, num_threads=1, seeds=[])` —
  new convenience method matching `pv.Vinecop.simulate`. Internally
  draws pseudo-random uniforms (`torch.rand` seeded from `seeds[0]`) or
  quasi-random uniforms (via `pyvinecopulib.utils.simulate_uniform`)
  and pushes them through `inverse_rosenblatt`.
- `TorchVinecop.cdf(u, N=10000, qrng=True, num_threads=1, seeds=[])` —
  new method that estimates the joint CDF via quasi-Monte-Carlo,
  matching `pv.Vinecop.cdf` to within MC error. A `block_size` kwarg
  caps the peak `(block, N, d)` scratch tensor for large query
  matrices.
- `TorchVinecop.pdf` / `rosenblatt` / `inverse_rosenblatt` now accept a
  `num_threads` keyword (default `1`) for signature parity with
  `pv.Vinecop`; on the torch backend it is a documented no-op. For CPU
  intraop parallelism call `torch.set_num_threads(N)` globally — note
  that mutates global state and is unsafe with concurrent workers.
- `TorchBicop.simulate(n, qrng=False, seeds=[])` — new method matching
  `pv.Bicop.simulate`. The previous `TorchBicop.sample(num_sample,
  seed, is_sobol)` is kept as a deprecated alias that forwards to
  `simulate(...)` and will be removed before the next stable release.

### CI

The notebook test and regenerate-notebooks jobs now install
`--extra torch` so that `examples/10_torch_backend.ipynb` can execute
under `nbmake`.

### `pyvinecopulib.sklearn`: VineRegressor and VineDensity

Two scikit-learn-compatible vine-copula estimators ship in the new
`pyvinecopulib.sklearn` submodule:

- `VineRegressor` — sklearn-style regressor that fits a vine copula to the
  joint distribution of `(X, y)` and predicts conditional means / quantiles.
- `VineDensity` — sklearn-style density estimator with `score_samples`,
  `score`, `sample`, `pdf`, and `cdf` methods.
- `VineForestDensity` / `VineForestRegressor` — ensembles of the above
  built via random search over vine structures and pruned by a model
  confidence set (MCS) on a held-out validation split. Predictions
  average across surviving structures. Structures are sampled either
  uniformly via Joe's algorithm (Joe, Cooke & Kurowicka 2011) or
  locally via Wilson's loop-erased random walk weighted by Kendall's
  τ (Wilson 1996); survivor selection uses a dual-split DA test
  adapted from the discrete-argmin-inference framework of Kim &
  Ramdas (2025), with Hansen, Lunde & Nason (2011) providing the
  foundational MCS definition.

Both follow the `BaseEstimator` / `RegressorMixin` / `DensityMixin` protocols
and handle mixed continuous/discrete inputs (DataFrame or ndarray). Class
docstrings include the full methodology — Sklar / pair-copula
factorization, Kde1d marginals, and the estimating-equation framework for
mean / quantile prediction — with references to Bedford & Cooke (2002),
Aas et al. (2009), Kraus & Czado (2017), and Nagler & Vatter (2024).

API surface tightened before v1:

- `VineRegressor.pdf` removed. Sklearn regressors don't expose density
  methods; users who need the joint or conditional density can call
  `pyvinecopulib.core.Vinecop.pdf` on the underlying fitted vine, or wait
  for the forest classes that surface ensemble-level log-likelihoods.
- `VineDensity.pdf(copula_only=...)` is now a real keyword argument (was
  documented but ignored).
- `VineDensity.cdf(X, N=10000, seeds=None)` added: returns the joint CDF
  via Monte-Carlo integration of the fitted vine copula.
- `schema` (both classes) and `normalize_weights` (regressor) are no
  longer `__init__` parameters. They remain settable via the
  `_schema` / `_normalize_weights` attributes for advanced / ensemble
  use; `clone()` won't preserve non-default values for these knobs.

The `[sklearn]` extra now pins `joblib>=1.3` and `scipy>=1.10` explicitly
(both are transitive sklearn dependencies today; pinning them
self-documents the forest requirement). Install via:

```bash
pip install pyvinecopulib[sklearn]
```

### Dependency changes

- `numpy>=2.0` is now a project-wide requirement (was `>=1.14`).
  `VineRegressor` needs `np.quantile(weights=...)` from NumPy 2.0; rather
  than pinning it only under `[sklearn]`, we bump everywhere — the wheel
  is now built and tested against the 2.x ABI consistently.
- `[sklearn]` extra: drops the redundant `numpy>=2.0` and adds
  `pandas>=2.0` (used by `VineBase.expand_factors` for DataFrame inputs).

## 1.0.0

### Breaking API changes in `pyvinecopulib`

The public API is now organized into four subpackages. Existing top-level
imports continue to work for backward compatibility, but the family constants
and most utilities emit a `DeprecationWarning` on access pointing at the new
canonical location. The deprecated aliases are scheduled for removal in 2.0.

| New canonical location | Names |
| --- | --- |
| `pyvinecopulib.core` | `Bicop`, `Vinecop`, `FitControlsBicop`, `FitControlsVinecop`, `RVineStructure`, `CVineStructure`, `DVineStructure`, `BicopFamily` |
| `pyvinecopulib.families` | `BicopFamily`, plus the constants (`indep`, `gaussian`, `student`, `clayton`, `gumbel`, `frank`, `joe`, `bb1`, `bb6`, `bb7`, `bb8`, `tawn`, `tll`) and the family groups (`all`, `parametric`, `nonparametric`, `one_par`, `two_par`, `three_par`, `elliptical`, `archimedean`, `extreme_value`, `bb`, `rotationless`, `lt`, `ut`, `itau`) |
| `pyvinecopulib.utils` | `Kde1d`, `wdm`, `to_pseudo_obs`, `pairs_copula_data`, `sobol`, `ghalton`, `simulate_uniform`, `benchmark` |
| `pyvinecopulib.sklearn` | (placeholder; estimator classes ship in a follow-up release) |

#### Kept at the top level indefinitely (no warning)

`Bicop`, `Vinecop`, `RVineStructure`, `CVineStructure`, `DVineStructure`,
`FitControlsBicop`, `FitControlsVinecop`, `BicopFamily`,
`to_pseudo_obs`, `__version__`.

#### Deprecated on access (warn, still resolves)

All family constants and groups (`indep`, `gaussian`, ..., `parametric`,
`nonparametric`, ..., `itau`); `Kde1d`, `wdm`, `sobol`, `ghalton`,
`simulate_uniform`, `benchmark`, `pairs_copula_data`. Use the canonical
subpackage path (`pyvinecopulib.families.<name>` or
`pyvinecopulib.utils.<name>`) to silence the warning.

#### `repr` and pickle now use canonical module paths

The C++ bindings now set `__module__` (and the hardcoded repr strings) per
subpackage: `Bicop.__module__ == "pyvinecopulib.core"`,
`BicopFamily.__module__ == "pyvinecopulib.families"`,
`Kde1d.__module__ == "pyvinecopulib.utils"`, etc. New pickles serialize via
the canonical path; pre-1.0 pickles still load via the deprecated top-level
alias (with a DeprecationWarning for symbols that fell off the kept list).

#### File renames

- `pyvinecopulib.pair_copuladata` (module) → `pyvinecopulib.utils.pairs_copula_data` (module). The function `pairs_copula_data` itself is unchanged.

### Build / packaging

- The `[sklearn]` optional extra was added (currently a placeholder; will
  install `scikit-learn>=1.4` and `numpy>=2.0` for the upcoming
  `pyvinecopulib.sklearn` estimators).
- Type stubs are now emitted per subpackage. The top-level stub also includes
  `def __getattr__(name: str) -> Any: ...` so static checkers don't flag access
  to deprecated names.
- `pytest` is configured to treat `pyvinecopulib`-originated
  `DeprecationWarning`s as errors so internal code stays on the canonical
  import paths.

## 0.7.6

This release was only used to pull the latest changes from `vinecopulib`'s `dev` branch, which includes some bug fixes.

## 0.7.5

### New features in `pyvinecopulib`

- Add support for 1d data with the `Kde1d` Python bindings (#189, #198)
- Add support for weighted dependence measures `wdm` (#194)
- Add an argument to control the nonparametric grid size in both `Kde1d` and `Bicop` (#191, #192)
- Release GIL for C++-only operations (#193)
- Add support for Python 3.14 (#200)

### Bug fixes in `pyvinecopulib`

- Improve repo health (#188)
- Improve pickling by making the state a dict instead of tuple for all classes (#190)

## 0.7.4

This release was yanked due to an issue with the wheel upload. Please use 0.7.5 instead.

## 0.7.3

### Breaking API changes in `pyvinecopulib`

- Rename the `mst_algorithm` argument in `FitControlsVinecop` to `tree_algorithm` to match the underlying C++ API (#178).
    - The accepted values are now `"mst_prim"`, `"mst_kruskal"` instead of `"prim"` and `"kruskal"`.
- Drop support for Python 3.8. Add support for Python 3.13 (#181).

### New features in `pyvinecopulib`

- Add support for random spanning trees in structure selection via `tree_algorithm="random_weighted"` and `"random_unweighted"` using Wilson's algorithm, enabling uniform or weight-proportional tree sampling (#178).
- Upgrade to `nanobind` v0.2.7, which removes the need for casting hacks and improves compatibility (#180).
- Introduce a `format()` method and improve `__str__` output for the `Vinecop` class (#144, #185).
- Add PEP 561-compatible stub files and full static type annotations across the package. Enable `mypy` checks in CI (#144).
- Add `.[dev]` and `.[examples]` extras in `pyproject.toml` for installing development and example dependencies (#181, #182, #183).
- Replace `flake8` with `ruff` for linting and add code formatting checks (#182).
- Execute and test Jupyter notebooks in CI, including Graphviz rendering (#181).

### Bug fixes in `pyvinecopulib`

- Fix non-deterministic structure selection in multithreaded environments by decoupling criterion computation from edge insertion (#178, [vinecopulib#640](https://github.com/vinecopulib/vinecopulib/pull/640)).
- Refactor unit test environment setup and cleanup; remove stale directories in tests (#183).
- Fix documentation rendering and improve docstrings across the package (#144, #185).

### Changes in `vinecopulib` version 0.7.3

These changes originate from the [`release 0.7.3 of vinecopulib`](https://github.com/vinecopulib/vinecopulib/releases/tag/v0.7.3), the C++ library which powers `pyvinecopulib`.

### BREAKING API CHANGES

- The `mst_algorithm` option to `FitControlsVinecop` has been renamed to `tree_algorithm` to
  allow for alternative spanning tree algorithms ([#637](https://github.com/vinecopulib/vinecopulib/pull/637)).
- `tree_algorithm`'s default value is now `"mst_prim"` instead of `"prim"`, and `"mst_kruskal"`
  replaces `"kruskal"` ([#637](https://github.com/vinecopulib/vinecopulib/pull/637)).
- The CMake option `VINECOPULIB_BUILD_SHARED_LIBS` has been changed to `VINECOPULIB_PRECOMPILED`
  to better reflect its purpose ([#641](https://github.com/vinecopulib/vinecopulib/pull/641)).

### NEW FEATURES

- Allow for random spanning trees as alternatives to the MST-based structure selection using
  `tree_algorithm` in `FitControlsVinecop` with `"random_weighted"` or `"random_unweighted"`
  ([#637](https://github.com/vinecopulib/vinecopulib/pull/637)).

### BUG FIXES

- Decouple edge insertion from criterion computation in `VinecopSelector` to fix randomness
  issues in structure selection when using multiple threads ([#640](https://github.com/vinecopulib/vinecopulib/pull/640))

### Changes in `vinecopulib` version 0.7.2

These changes originate from the [`release 0.7.2 of vinecopulib`](https://github.com/vinecopulib/vinecopulib/releases/tag/v0.7.2), the C++ library which powers `pyvinecopulib`.

### BUG FIXES

- More build system updates by @jschueller ([#633](https://github.com/vinecopulib/vinecopulib/pull/633))

- Fix deprecation warning in json header ([#634](https://github.com/vinecopulib/vinecopulib/pull/634))

- Fix TLL speed issues related to FFT ([#635](https://github.com/vinecopulib/vinecopulib/pull/635))

## 0.7.1

### New features in `pyvinecopulib`

- Add pickle support for all classes (#168)
- Add `allow_rotation` option to `FitControlsBicop` and `FitControlsVinecop` (#168)

### Bug fixes in `pyvinecopulib`

- Upgrade nanobind to allow for single row matrices (fix #169 and #170)

### Changes in `vinecopulib` version 0.7.1

These changes originate from the latest release of [`vinecopulib`](https://github.com/vinecopulib/vinecopulib/releases/tag/v0.7.1), the C++ library which powers `pyvinecopulib`.

#### NEW FEATURES

- Add `allow_rotation` option to `FitControlsBicop` and `FitControlsVinecop`
  to allow for the rotation of the pair copulas ([#628](https://github.com/vinecopulib/vinecopulib/pull/628)).
- Add a `FitControlsConfig` struct to create flexible and yet safe constructors
  for `FitControlsBicop` and `FitControlsVinecop` ([#629](https://github.com/vinecopulib/vinecopulib/pull/629)).

#### BUG FIXES

- Restrict parameter range for fitting Tawn copulas; fix handling of their
  shape/argument order ([#620](https://github.com/vinecopulib/vinecopulib/pull/620)).
- Compute and save loglik/nobs in `Vinecop::fit()` ([#623](https://github.com/vinecopulib/vinecopulib/pull/623))
- Disable unwanted compiler output related to BOOST_CONCEPT checks ([#624](https://github.com/vinecopulib/vinecopulib/pull/624))

## 0.7.0

This version introduces a switch to nanobind as a backend (#160): i.e., the C++ bindings, now use [nanobind](https://nanobind.readthedocs.io/) instead of [pybind11](https://pybind11.readthedocs.io/).
It allows for considerable performance improvements (~8x speedup in our latest benchmarks) and smaller binaries.

### Breaking API changes in `pyvinecopulib`

- Removal of the overloaded constructors:
    - For all classes, only one constructor is now available.
    The reason is that the overloaded constructors were un-Pythonic, error-prone, and could not be properly documented with Sphinx.
    They have been replaced by a single constructor for each class, along with factory `from_xzy` methods.
    - For the ``Bicop`` class:
        - ``Bicop.from_family()``: Instantiate from a family, rotation, parameters, and variable types.
        - ``Bicop.from_data()``: Instantiate from data, as well as optional controls and variable types.
        - ``Bicop.from_file()``: Instantiate from a file.
        - ``Bicop.from_json()``: Instantiate from a JSON-like string.
    - For the ``Vinecop`` class:
        - ``Vinecop.from_dimension()``: Instantiate an empty vine copula of a given dimension.
        - ``Vinecop.from_data()``: Instantiate from data, as well as an optional ``FitControlsVinecop``, an ``RVineStructure`` or matrix, and variable types.
        - ``Vinecop.from_structure()``: Instantiate from an ``RVineStructure`` or matrix, as well as optional pair-copulas and variable types.
        - ``Vinecop.from_file()``: Instantiate from a file.
        - ``Vinecop.from_json()``: Instantiate from a JSON-like string.
    - For the ``RVineStructure`` class:
        - ``RVineStructure.from_dimension()``: Instantiate a default structure of a given dimension and truncation level.
        - ``RVineStructure.from_order()``: Instantiate from an order vector.
        - ``RVineStructure.from_matrix()``: Instantiate from a matrix.
        - ``RVineStructure.from_file()``: Instantiate from a file.
        - ``RVineStructure.from_json()``: Instantiate from a JSON-like string.

### New features in `pyvinecopulib`

- Expose more structure methods to python (#157)
- Switch to nanobind as a backend (#160)
- New IO methods for `Bicop` and `Vinecop` classes to use JSON-like strings (#160)
- Extensive documentation revamp (#160)
- Adding a benchmark example (#160)
- Convertion of all examples to Jupyter notebooks (#160)

### Bug fixes in `pyvinecopulib`

- Install and test source distribution (#164)

### Changes in `vinecopulib`

These changes originate from the underlying C++ library, [`vinecopulib`](https://github.com/vinecopulib/vinecopulib), which powers `pyvinecopulib`.

#### New features

- Use analytical derivatives in discrete pdf/hfuncs ([#572](https://github.com/vinecopulib/vinecopulib/pull/572))
- Allow for alternative for `"prim"` vs `"kruskal"` in MST-based model selection ([#577](https://github.com/vinecopulib/vinecopulib/pull/577))
- Improve the dependencies install script to use it in other projects ([#576](https://github.com/vinecopulib/vinecopulib/pull/576))
- Add tawn copula ([#579](https://github.com/vinecopulib/vinecopulib/pull/579))
- Improve doc ([#580](https://github.com/vinecopulib/vinecopulib/pull/580), [#585](https://github.com/vinecopulib/vinecopulib/pull/585), [#607](https://github.com/vinecopulib/vinecopulib/pull/607))
- Allow for the discrete Rosenblatt transform ([#581](https://github.com/vinecopulib/vinecopulib/pull/581))
- Add `Vinecop::fit()` ([#584](https://github.com/vinecopulib/vinecopulib/pull/584))
- Improve `Bicop::str()` ([#588](https://github.com/vinecopulib/vinecopulib/pull/588)) and `Vinecop::str()` ([#589](https://github.com/vinecopulib/vinecopulib/pull/589))
- Properly handle discrete variables for the TLL family ([#597](https://github.com/vinecopulib/vinecopulib/pull/597))
- Weighted pseudo-observations ([#602](https://github.com/vinecopulib/vinecopulib/pull/602))
- Cross-platform random numbers and add seeds options to `to_pseudo_obs` ([#603](https://github.com/vinecopulib/vinecopulib/pull/603))
- Improve performance by
    - aligning with the `R` defaults (e.g., `BOOST_NO_AUTO_PTR`, `BOOST_ALLOW_DEPRECATED_HEADERS`, `BOOST_MATH_PROMOTE_DOUBLE_POLICY=false`, `std::string nonparametric_method = "constant"` for the TLL instead of `"quadratic"`, `-O3 -march=native` compiler flags) and add benchmarking example ([#592](https://github.com/vinecopulib/vinecopulib/pull/592), [#611](https://github.com/vinecopulib/vinecopulib/pull/611), [#613](https://github.com/vinecopulib/vinecopulib/pull/613)),
    - using `Eigen` element-wise operations instead of `boost` whenever possible ([#598](https://github.com/vinecopulib/vinecopulib/pull/598), [#612](https://github.com/vinecopulib/vinecopulib/pull/612)),
    - using binary search in the TLL for `get_indices` ([#613](https://github.com/vinecopulib/vinecopulib/pull/613)).

#### Bug fixes

- Improve stability in BB7 PDF ([#573](https://github.com/vinecopulib/vinecopulib/pull/573))
- Revamped CI/CD pipeline, tests discoverable by CTest, boost version on windows (([66cf8b0](https://github.com/vinecopulib/vinecopulib/commit/66cf8b0)))
- Fix ASAN issues ([#583](https://github.com/vinecopulib/vinecopulib/pull/583))
- Fix interface includes and other CMake issue ([#586](https://github.com/vinecopulib/vinecopulib/pull/586), [#599](https://github.com/vinecopulib/vinecopulib/pull/599), [#601](https://github.com/vinecopulib/vinecopulib/pull/601), [#608](https://github.com/vinecopulib/vinecopulib/pull/608)), by @jschueller
