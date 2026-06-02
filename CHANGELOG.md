# Changelog

## 1.0.0

### Breaking API changes in `pyvinecopulib`

- Reorganize the public API into the `core` / `families` / `utils` / `sklearn` subpackages (#207).
    - Top-level classes (`Bicop`, `Vinecop`, `RVineStructure`, `CVineStructure`, `DVineStructure`, `FitControlsBicop`, `FitControlsVinecop`, `BicopFamily`) and `to_pseudo_obs` are kept at the top level indefinitely.
    - Family constants / groups (`indep`, `gaussian`, …, `parametric`, …, `itau`), `Kde1d`, `wdm`, `sobol`, `ghalton`, `simulate_uniform`, `benchmark`, `pairs_copula_data` still resolve at the top level but emit a `DeprecationWarning` on access pointing at the canonical subpackage path. Aliases are scheduled for removal in 2.0.
    - `repr` and pickle now use canonical module paths (`Bicop.__module__ == "pyvinecopulib.core"`, `Kde1d.__module__ == "pyvinecopulib.utils"`, etc.); pre-1.0 pickles still load via the deprecated aliases.
- `pyvinecopulib.sklearn` estimators take a single `backend=` keyword (#218).
    - `VineDensity(controls=..., structure=..., seed=...)` → `VineDensity(backend=VinecopBackend(controls=..., structure=...), random_state=...)`. Same shape change for `VineRegressor` / `VineForestDensity` / `VineForestRegressor`. No string shortcuts; pass a `VinecopBackend` or `TorchVinecopBackend` instance.
    - `seed`-style kwargs renamed to `random_state` across the forests and on `VineDensity.sample` / `VineDensity.cdf` / `VineForestDensity.cdf`. Legacy `seed` / `seeds` kwargs removed without a deprecation alias.
    - `VineRegressor` gains a real `normalize_weights=True` `__init__` parameter; the previous post-init `_normalize_weights` attribute is gone.
- `TorchBicop.from_data` now dispatches on `controls=FitControlsTorchBicop(...)` (#217). Callers who passed `grid_size`, `mult`, or `grid_type` as keyword arguments must move them onto the dataclass. `cache_integrals`, `device`, and `dtype` remain direct kwargs on `from_data`.
- `TorchBicop.sample(num_sample, seed, is_sobol)` renamed to `TorchBicop.simulate(n, qrng=False, seeds=[])` for parity with `pv.Bicop.simulate` (#216).

### New features in `pyvinecopulib`

- Add `pyvinecopulib.sklearn` with scikit-learn-compatible `VineDensity` and `VineRegressor` estimators following the `BaseEstimator` / `DensityMixin` / `RegressorMixin` protocols, with mixed continuous/discrete input handling (DataFrame or ndarray) (#211).
- Add `VineForestDensity` and `VineForestRegressor`: ensembles of vine estimators sampled either uniformly (Joe's algorithm) or via Wilson's Kendall's-τ-weighted random walk, pruned by a model-confidence-set survivor test on a held-out split (#213).
- Add the `pyvinecopulib.torch` subpackage: pure-PyTorch `TorchBicop` / `TorchVinecop` evaluators (`nn.Module` subclasses) for GPU placement, autograd, and `nn.Module` pipeline composition. The torch cascade matches the C++ TLL fit to machine precision (#216).
- Add a sklearn-side backend layer (`pyvinecopulib.sklearn.backends`): `VinecopBackend` (default, C++) and `TorchVinecopBackend` (opt-in PyTorch) implement a shared `VinecopLike` protocol so the same estimator class routes through either backend (#218).
- Add the optional amortized `method="vdc"` pair-copula fit (Safaai 2026, [arXiv:2604.20568](https://arxiv.org/abs/2604.20568)) behind the `[vdc]` extra; reference implementation [KempnerInstitute/vine-denoising-copula](https://github.com/KempnerInstitute/vine-denoising-copula) (#217).
- Add `TorchVinecop.from_structure`, `TorchVinecop.simulate`, and `TorchVinecop.cdf` (quasi-MC) so the torch evaluators mirror their `pv.Vinecop` counterparts (#216).
- Flip torch defaults based on a bicop + vine benchmark sweep: `cache_integrals=True` everywhere (80–300× faster `cdf` / `hfunc` / `hinv` on cpu, 2–80× on cuda), `batched` resolves device-aware via `_default_batched()` (`True` on cuda, `False` on cpu) (#219).
- Add a tutorial-style `docs/concepts.rst` introducing Sklar's theorem, pair-copula construction, R-vines, and the TLL family in a ~5-minute read (#218).
- Use Sphinx autosummary on the four subpackage landing pages so module docstrings, classes, and free functions get their own indexed pages (#214).

### Build / packaging

- Migrate to `uv` + `scikit-build-core` for the editable / wheel build pipeline, with `[build-system].requires` mirroring the dev `[dependency-groups]` so `--no-build-isolation` works out of the box (#209).
- Replace `mypy` with `ty` (Astral's type checker, alpha) and enable strict checks against a Python 3.10 baseline; only `pyvinecopulib.pyvinecopulib_ext` is allowed as an unresolved import (#210).
- Refactor the build / docs / examples pipeline into a thin Makefile over `uv run` and rework `scripts/regenerate_notebooks.py` (#205).
- Add a `pyvinecopulib[sklearn]` extra (`scikit-learn>=1.4`, `pandas>=2.0`, `joblib>=1.3`, `scipy>=1.10`), a `pyvinecopulib[torch]` extra (`torch>=2.0`), and a `pyvinecopulib[vdc]` extra resolving `vine-denoising-copula` from GitHub via `[tool.uv.sources]` (#211, #216, #217).
- Treat `pyvinecopulib`-originated `DeprecationWarning`s as errors under pytest so internal call sites stay on the canonical import paths (#207).
- Install `--extra torch` in the notebook-test and regenerate-notebooks CI jobs so `examples/10_torch_backend.ipynb` executes under `nbmake` (#216).
- Fix osx / musllinux wheel builds: feed libclang the compiler's implicit system include dirs (plus `-ferror-limit=0` and the macOS SDK sysroot) so the C++ stdlib / intrinsic headers resolve, and abort `docstr.hpp` generation only on *fatal* libclang diagnostics (intrinsic-header `error`s are benign and no longer silently drop symbols).
- Migrate macOS CI to `macos-15` and re-add `macos-15-intel` (`MACOSX_DEPLOYMENT_TARGET=10.13` for nanobind's aligned `new`/`delete`); wheel matrix is now 5 platforms × 3 ABI (15 wheels).

### Bug fixes in `pyvinecopulib`

- Port the `integrate_2d` marginal-renormalisation fix to the torch backend (`InterpolationGrid2D.integrate_2d` and `integrate_2d_batched`) so `TorchBicop.cdf` enforces ``C(1, u_2) = u_2`` exactly, matching the post-vinecopulib#667 C++ CDF to machine precision on the on-the-fly path ([vinecopulib#667](https://github.com/vinecopulib/vinecopulib/pull/667)).

### Dependency changes

- `numpy>=2.0` is now a project-wide requirement (was `>=1.14`); `VineRegressor` needs `np.quantile(weights=...)` from NumPy 2.0 (#211).
- `[sklearn]` extra adds `pandas>=2.0` (used by `VineBase.expand_factors` for DataFrame inputs) (#211).

### Changes in `vinecopulib`

#### NEW FEATURES

- Early exit in vine selection when the structure is already a tree, avoiding redundant work in `select` ([vinecopulib#661](https://github.com/vinecopulib/vinecopulib/pull/661)).
- Per-family parameter / rotation / tail-dependence documentation on the `BicopFamily` enum members, surfaced through the Python `families` subpackage (#214, [vinecopulib#668](https://github.com/vinecopulib/vinecopulib/pull/668)).
- Numpydoc-compliant `//!` comments on every property getter / setter in the Python-binding surface, surfaced through the pyvinecopulib autosummary pages (#214, [vinecopulib#670](https://github.com/vinecopulib/vinecopulib/pull/670)).

#### BUG FIXES

- Fix `integrade_2d` numerical handling on the bivariate-copula CDF path ([vinecopulib#667](https://github.com/vinecopulib/vinecopulib/pull/667)).
- Fix a typo in `pdf_d_d` and tighten the threshold under which the analytical derivative is preferred over the finite-difference fallback ([vinecopulib#664](https://github.com/vinecopulib/vinecopulib/pull/664)).
- Drop an unused fit-controls option ([vinecopulib#662](https://github.com/vinecopulib/vinecopulib/pull/662)).
- Fix Doxygen typos surfaced by the `libclang` extractor ([vinecopulib#665](https://github.com/vinecopulib/vinecopulib/pull/665)).
- Bump CI off the deprecated `macos-13` runner ([vinecopulib#663](https://github.com/vinecopulib/vinecopulib/pull/663)).

### Changes in `kde1d`

- Numpydoc-compliant `//!` getter-return docstrings on `Kde1d`, surfaced through `pyvinecopulib.utils.Kde1d` (#214, [kde1d#26](https://github.com/vinecopulib/kde1d-cpp/pull/26)).
- Capitalize method first words, add a fit-summary, and fix the quantile description ([kde1d#27](https://github.com/vinecopulib/kde1d-cpp/pull/27)).

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
