# Changelog

## 0.8.0 (unreleased)

### Breaking API changes in `pyvinecopulib`

- Reorganize the public API into the `core` / `families` / `utils` / `sklearn` subpackages (#207).
    - Top-level classes (`Bicop`, `Vinecop`, `RVineStructure`, `CVineStructure`, `DVineStructure`, `FitControlsBicop`, `FitControlsVinecop`, `BicopFamily`) and `to_pseudo_obs` are kept at the top level indefinitely.
    - Family constants / groups (`indep`, `gaussian`, …, `parametric`, …, `itau`), `Kde1d`, `wdm`, `sobol`, `ghalton`, `simulate_uniform` (now `utils.sample_uniform`), `benchmark`, `pairs_copula_data` still resolve at the top level but emit a `DeprecationWarning` on access pointing at the canonical subpackage path. Aliases are scheduled for removal in the next major release.
    - `repr` and pickle now use canonical module paths (`Bicop.__module__ == "pyvinecopulib.core"`, `Kde1d.__module__ == "pyvinecopulib.utils"`, etc.); pickles from releases up to 0.7.x still load via the deprecated aliases.
- `pyvinecopulib.sklearn` estimators take a single `backend=` keyword (#218).
    - `VineDensity(controls=..., structure=..., seed=...)` → `VineDensity(backend=VinecopBackend(controls=..., structure=...), random_state=...)`. Same shape change for `VineRegressor`. No string shortcuts; pass a `VinecopBackend` or `TorchVinecopBackend` instance.
    - `seed`-style kwargs renamed to `random_state` on the estimator constructors and on `VineDensity.sample` / `VineDensity.cdf`. Legacy `seed` / `seeds` kwargs removed without a deprecation alias.
    - `VineRegressor` gains a real `normalize_weights=True` `__init__` parameter; the previous post-init `_normalize_weights` attribute is gone.
- `TorchBicop.from_data` now dispatches on `controls=FitControlsTorchBicop(...)` (#217). Callers who passed `grid_size`, `mult`, or `grid_type` as keyword arguments must move them onto the dataclass. `cache_integrals`, `device`, and `dtype` remain direct kwargs on `from_data`.
- `TorchBicop.sample` takes `(n, qrng=False, seeds=[])` instead of `(num_sample, seed, is_sobol)`, for parity with `pv.Bicop.sample` (#216).
- Remove `TorchBicop`'s `log_pdf` method. The vine `pdf` cascade now accumulates a product of per-edge `pdf` (rather than a log-sum-exp), so the pair-level log-density convenience method is no longer needed; it was never part of the `BicopLike` contract. Use `TorchBicop.pdf(u).log()` if you need a log density.
- `VineDensity` and `VineRegressor` put their scikit-learn mixin before `VineBase`, so `is_regressor` / `is_density_estimator` are now `True` and meta-estimators such as `StackingRegressor` accept them; previously they were rejected as "not a regressor" (#263).
- `VineRegressor.predict` keeps the sample axis for a one-row `X`: it returns shape `(1,)` (or `(1, k)` with quantiles) instead of collapsing to a 0-d scalar (#263).
- `BicopLike` / `BicopBase` take the conditioning matrix `x` as a keyword-only argument, matching `VinecopLike`. It shared a positional slot with the compiled `Bicop`'s per-row `parameters`, so hosting a `pv.Bicop` in a non-simplified `VinecopBase` silently evaluated the conditioning values as parameters; it now raises `TypeError` (#265).
- `FitControlsBicop` and `FitControlsVinecop` default `selection_criterion` to `"aic"` instead of `"bic"`, matching the C++ and R defaults; the same knob now selects the same model from all three languages. Pass `selection_criterion="bic"` explicitly to keep the previous selection (#251).
- Rename the sampling entry points to `sample`: `simulate` on `Bicop` / `Vinecop` / `RVineStructure` / `Kde1d` and `utils.simulate_uniform` stay as `DeprecationWarning` aliases, while the unreleased `simulate_conditional` → `sample_conditional` and `_simulate_uniform` → `_sample_uniform` hook get none (#297).
- Add `n_jobs` to `VineDensity` and `VineRegressor`, following the scikit-learn convention (`None` is one thread, `-1` every processor). It governs fitting and every evaluation, where the fit-time `num_threads` used to pin both: `VineRegressor.predict` on 2000 rows goes from 1.13 s to 0.26 s, and a 20-dimensional fit from 9.4 s to 1.8 s. The default stays serial, since a caller parallelizing over vines owns the parallelism. Results are bit-identical at any thread count (#297).
- Collapse the three hand-rolled Graphviz retry loops into one composite action on `nick-fields/retry`, with a per-manager retry shape: apt and brew fail slowly so they get few attempts and a long budget, Chocolatey fails fast so it gets more attempts and a short one, and each apt attempt is bounded from the inside by `sudo timeout` (#297).

- Bump `lib/vinecopulib` and `lib/wdm`, and port the two `tll` behavior changes the torch evaluator mirrors: the margin renormalization is the geometric mean of both sweep orders iterated to `1e-10` (`norm_times = 3` -> `norm_maxiter = 25`, renamed because the meaning changed from "exactly" to "at most"), and `hfunc` / `hinv` no longer floor the interpolated density at `1e-4`. Every `tll` fit moves, continuous and discrete alike, and `tll` is the default family (#305).
- `cache_integrals=True` changes meaning: the cached `cdf` / `hfunc*` are now exact rather than a bilinear surrogate with ~1e-3 mean error, and they carry an exact gradient in `values` as well as in `u`. `hinv*` no longer read a cache in either mode -- locating the bracketing cell needs the conditional cumulative along the whole free axis -- so they were off by up to 2e-2 and are now bit-identical to the on-the-fly path (#305).
### New features in `pyvinecopulib`

- Add backend-neutral pair-copula / vine contracts to `pyvinecopulib.core`: the generic, `runtime_checkable` protocols `BicopLike` / `VinecopLike` and `BicopBase`, their canonical array-backend-agnostic (NumPy / PyTorch, via `array_api_compat`) partial implementation. `BicopLike` mirrors the C++ `Bicop` evaluation surface exactly (`pdf` / `cdf` / `hfunc1` / `hfunc2` / `hinv1` / `hinv2`, with an optional trailing `x` conditioning matrix; no `dtype` / `device`, since an `nn.Module`-backed pair has no intrinsic precision/placement), so the C++ `Bicop` / `Vinecop` satisfy `BicopLike` / `VinecopLike` structurally (`isinstance` is `True`); the `x` conditioning matrix appears symmetrically on both (per-pair on `BicopLike`, per-call on `VinecopLike`). a custom pair need only implement `pdf` / `hfunc1` / `hfunc2` and inherits `hinv1` / `hinv2` (bisection), `sample` (inverse Rosenblatt of the pair, via a `_sample_uniform` RNG hook), `loglik` (`sum(log pdf)`, differentiable on autograd backends), `plot` (density contour / surface), and a structural `__repr__` from `BicopBase`. `sample` is part of the `BicopLike` contract (symmetric with `VinecopLike`; the C++ `Bicop` provides it). `TorchBicop` now subclasses `BicopBase[torch.Tensor]`. These are the extension point for custom (e.g. neural, conditional) pair copulas, and `pyvinecopulib.core` imports without PyTorch.
- Extract the vine evaluation cascades onto an array-agnostic `VinecopBase` in `pyvinecopulib.core` (NumPy / PyTorch): `pdf` / `cdf` / `rosenblatt` / `inverse_rosenblatt` / `sample`, plus `loglik` / `plot` / `__repr__` / `dim` / `get_pair_copula` and the public `fit` engine. Conditioning is threaded per edge through a pluggable `ConditioningContext` (`SimplifiedContext` default — external covariates only; `NonSimplifiedContext` also gathers the edge's conditioning-set values `u_D`): `fit(structure, u, fit_edge, *, context=, x=)` builds a **non-simplified / conditional** vine, assembling each pair's `x_e` in the fixed C1 column order (`u_D` in ascending conditioning-tree index, then external covariates `x`). A backend hosting arbitrary (e.g. scikit-style, conditional) pair copulas subclasses `VinecopBase` directly. `TorchVinecop` subclasses `VinecopBase[torch.Tensor]`, stays an `nn.Module` vine of `TorchBicop` pairs (GPU / autograd / `state_dict` honor the whole object), keeps the batched TLL fast path (falling back when a pair lacks the grid internals), and inherits `loglik` / `plot`. Behavior-preserving: the `pdf` cascade is a product of per-edge densities and the torch↔`Vinecop` (1e-10) / batched↔non-batched (1e-12) parity are unchanged.
- Evaluate the batched TLL cascade from the exact prefix tables (torch): with `cache_integrals=True` the per-tree-level batched `hfunc1` / `hfunc2` reconstruct the conditional distribution function in closed form from one `(N, m, m)` cumulative table per argument, replacing the stacked bilinear value caches. On CUDA the batched fast path runs ~6-18x faster than the non-batched cascade (`pdf` / `rosenblatt`, d = 10-20, n = 2000-10000; e.g. `pdf` at d = 20, n = 10000: ~388 -> 22 ms) and up to ~16-18x faster than single-threaded C++ (benchmark: `scripts/bench_torch_vinecop.py`) (#305).
- Add `pyvinecopulib.sklearn` with scikit-learn-compatible `VineDensity` and `VineRegressor` estimators following the `BaseEstimator` / `DensityMixin` / `RegressorMixin` protocols, with mixed continuous/discrete input handling (DataFrame or ndarray) (#211).
- Add the `pyvinecopulib.torch` subpackage: pure-PyTorch `TorchBicop` / `TorchVinecop` evaluators (`nn.Module` subclasses) for GPU placement, autograd, and `nn.Module` pipeline composition. The torch cascade matches the C++ TLL fit to machine precision (#216).
- Add a sklearn-side backend layer (`pyvinecopulib.sklearn.backends`): `VinecopBackend` (default, C++) and `TorchVinecopBackend` (opt-in PyTorch) route the same estimator class through either backend; both produce fitted vines that satisfy the canonical `pyvinecopulib.core.VinecopLike` protocol (#218).
- Add `TorchVinecop.from_structure`, `TorchVinecop.sample`, and `TorchVinecop.cdf` (quasi-MC) so the torch evaluators mirror their `pv.Vinecop` counterparts (#216).
- `TorchVinecop.from_data` now selects an R-vine structure automatically when `structure=None`, mirroring `pv.Vinecop.from_data`, while a supplied `structure` still fixes the skeleton and fits only the pair copulas. With the API gap closed, the two `pyvinecopulib.sklearn` backends collapse onto a single shared base: `VinecopBackend` and `TorchVinecopBackend` now differ only in the divergent members, the `name` / `supports_discrete` / `supports_cdf` / `supports_simulate` capability flags are gone (the continuous-only guard lives on `TorchVinecop.from_data`), and `resolve_backend` only defaults `None`. `pyvinecopulib.torch` no longer re-exports `InterpolationGrid2D` (it stays an internal detail of `TorchBicop`) (#241).
- Port vine structure selection into the array-agnostic `VinecopBase`, so `TorchVinecop.from_data(structure=None)` selects the R-vine **natively in torch** instead of round-tripping through `Vinecop`. `VinecopBase.select(u, fit_edge, *, trunc_lvl=, tree_criterion=, threshold=, tree_algorithm=, seeds=, to_numpy=)` is the array-agnostic (NumPy / PyTorch) Dissmann / Wilson selector, an **exact port** of `Vinecop`'s: it builds each tree's candidate graph under the proximity condition (same enumeration order), weights edges by `1 - |tau|` (Kendall's tau via `wdm`), keeps a spanning tree, fits the pair copulas to propagate h-functions, finalizes the chosen edges into an `RVineStructure` with the same peel and diagonal policy, and returns `(structure, pair_copulas)` — the selection-time pairs reoriented onto their finalized slots, so (like `Vinecop`) no re-fit happens. The fixed-structure fit engine is exposed as `VinecopBase.fit` (renamed from the unreleased `sequential_fit`), so `VinecopBase` now mirrors `Vinecop`'s `select` / `fit` by name — both differ only in being functional (they return the structure / fitted pairs rather than mutating the vine in place, since `VinecopBase` leaves pair storage to the subclass). Given identical pair fits, the selected matrix encoding, the reused pairs, and the resulting density are **identical** to `Vinecop.from_data`'s (bit-for-bit on the NumPy path; to TLL-fit precision, ~1e-13, on the torch path). Two scoped C++ primitives support it: `RVineStructure.from_trees(d, trees)` (reuses upstream's `RVineTrees` peel — the same finalization `Vinecop` uses) and an internal boost-based spanning-tree helper (prim / kruskal / Wilson). The `BicopLike` contract gains a `flip()` method (return the argument-swapped copula) used to reorient reused pairs: `Bicop.flip()` is now bound (returns a flipped copy), `TorchBicop.flip()` transposes its interpolation grid, and `BicopBase` ships a raising default (custom pairs only need it to be hosted in *selection*). The `structure_controls` field on `FitControlsTorchVinecop` (a nested `FitControlsVinecop`, added in the unreleased #241) is **replaced** by native fields `trunc_lvl` / `tree_criterion` / `threshold` / `tree_algorithm` / `seeds`; both sklearn backends now steer selection through these same-named fields, dropping the last backend-specific indirection. Structure selection is continuous-only and TLL-only on the torch path; automatic truncation / thresholding (`aic` / `bic` / `mbicv`) is not ported (`trunc_lvl` is a fixed cap) (#244).
- Flip torch defaults based on a bicop + vine benchmark sweep: `cache_integrals=True` everywhere (80–300× faster `cdf` / `hfunc` / `hinv` on cpu, 2–80× on cuda), `batched` resolves device-aware via `_default_batched()` (`True` on cuda, `False` on cpu) (#219).
- Add a tutorial-style `docs/concepts.rst` introducing Sklar's theorem, pair-copula construction, R-vines, and the TLL family in a ~5-minute read (#218).
- Use Sphinx autosummary on the four subpackage landing pages so module docstrings, classes, and free functions get their own indexed pages (#214).
- Add a custom `tree_criterion` for vine structure selection: set `FitControlsVinecop(tree_criterion="custom")` and supply `tree_criterion_function`, a callable `f(data, weights) -> float` mapping a two-column array of pair pseudo-observations (and observation weights) to a scalar edge weight (its absolute value is used). The callable round-trips through pickle when it is picklable (e.g. a module-level function); calls acquire the GIL, so they serialize under `num_threads > 1` ([vinecopulib#674](https://github.com/vinecopulib/vinecopulib/pull/674)).
- Add per-row-parameter overloads of `Bicop.pdf` / `cdf` / `hfunc1` / `hfunc2` / `hinv1` / `hinv2` / `loglik`: pass an `(n, p)` array of `parameters` (one row per observation, `p` family parameters each) plus an optional `num_threads` to evaluate the copula with a different parameter set per row in a single (optionally threaded) call, instead of reusing the object's stored parameters. Parametric families only ([vinecopulib#675](https://github.com/vinecopulib/vinecopulib/pull/675)).
- Expose the vine gradient/diagnostics surface on `Vinecop`: `pdf_full` (density plus, with `keep_all=True`, the per-edge densities and h-functions, returned as a dict of `[tree][edge]` nested arrays; the left-limit `hfunc*_sub` entries are only populated for models with discrete variables), `scores` (observation-wise score matrix), `hessian` (average Hessian), `hessian_full` (per-observation Hessians), and `scores_cov` (score covariance). Scores and Hessians are computed analytically (models with discrete variables fall back to finite differences); models with nonparametric pair copulas raise `RuntimeError`. Mirrors the R additions in [rvinecopulib#320](https://github.com/vinecopulib/rvinecopulib/pull/320) ([vinecopulib#679](https://github.com/vinecopulib/vinecopulib/pull/679), [vinecopulib#683](https://github.com/vinecopulib/vinecopulib/pull/683)).
- Port the closed-form TLL conditional-quantile inversion (vinecopulib#691) to the torch backend: `InterpolationGrid2D.inverse_integrate_1d` inverts the piecewise-quadratic conditional cdf exactly, replacing the vectorized ITP root-finder in `TorchBicop.hinv1` / `hinv2` (`cache_integrals=False`) and in the inverse-cache construction. The no-cache h-inverse now matches the C++ `Bicop.hinv*` to machine precision (~1e-15, was root-finder-limited) and is ~60–120× faster on cpu (1.4 ms vs 171 ms for 10k points; see `scripts/bench_torch_bicop_fit.py --mode hinv`); the internal `_util.solve_itp` helper is removed (a reference copy lives in the benchmark script).
- Add binary CBOR model persistence: `Bicop.to_file` / `Vinecop.to_file` / `RVineStructure.to_file` and the matching `from_file` factories select CBOR automatically when the filename ends in `.cbor`; all other filenames keep reading / writing JSON text. `to_json` / `from_json` (and pickling, which serializes via JSON strings) are unchanged ([vinecopulib#684](https://github.com/vinecopulib/vinecopulib/pull/684)).
- Add `Vinecop.gradient` (the observation-average of `scores`, mirroring how `hessian` averages `hessian_full`) and `Vinecop.scores_full` (the scores together with the per-edge derivative caches behind them, returned as a dict of `[tree][edge]` nested arrays with a `keep_all` option mirroring `pdf_full`; caches are only populated on the analytic path) ([vinecopulib#683](https://github.com/vinecopulib/vinecopulib/pull/683)).
- Add analytic derivatives of the bivariate copula density and h-functions to `Bicop`: `pdf_deriv`, `pdf_deriv2`, `hfunc1_deriv`, `hfunc1_deriv2`, `hfunc2_deriv`, `hfunc2_deriv2`, `logpdf_deriv`, and `logpdf_deriv2`, each taking a `deriv` selector string (`"par1"`, `"par2"`, …, `"u1"`, `"u2"`; `"par"` is short for `"par1"`; second-order selectors concatenate two components in any order, a single component meaning twice, e.g. `"par1u1"` or `"u1"` ≡ `"u1u1"`). Closed forms are used for the families in `families.analytic_derivs` (currently every parametric family); rotations are handled internally, so derivatives are always w.r.t. the natural parameters and arguments of the rotated copula. Like `pdf` and friends, each method accepts optional per-row `parameters` (an `(n, p)` array) and `num_threads`. For `tll`, the argument gradients `pdf_deriv` / `logpdf_deriv` w.r.t. `"u1"` / `"u2"` return the exact slope of the interpolation grid; parameter selectors, second-order, and h-function derivatives raise ([vinecopulib#683](https://github.com/vinecopulib/vinecopulib/pull/683), [vinecopulib#687](https://github.com/vinecopulib/vinecopulib/pull/687), [vinecopulib#694](https://github.com/vinecopulib/vinecopulib/pull/694)).
- Add the `families.analytic_derivs` group: families with closed-form density / h-function derivatives ([vinecopulib#683](https://github.com/vinecopulib/vinecopulib/pull/683), [vinecopulib#687](https://github.com/vinecopulib/vinecopulib/pull/687)).
- Add tail dependence and Blomqvist's beta to `Bicop`: the read-only `taildep` property (a 2×2 matrix collecting the tail dependence coefficients in the four corners of the unit square; NaN for `tll`) and `beta` property (Blomqvist's β = 4·C(0.5, 0.5) − 1, all families), plus `parameters_to_taildep` / `parameters_to_beta` for evaluation at arbitrary parameters, analogous to `parameters_to_tau` ([vinecopulib#682](https://github.com/vinecopulib/vinecopulib/pull/682)).
- Add `RVineStructure.from_struct_array(order, struct_array, natural_order=False, check=True)` (build a structure from an order vector and a `[tree][edge]` nested-list structure array) and `RVineStructure.get_struct_array(natural_order=False)` / `Vinecop.get_struct_array(natural_order=False)` (the full structure array as a nested list, complementing the per-entry `struct_array(tree, edge)`), on top of the `TriangularArray` conversions added in [vinecopulib#680](https://github.com/vinecopulib/vinecopulib/pull/680).
- Add the log-likelihood score family to `Bicop` (parametric, continuous families): `scores` (the `(n, p)` per-observation score matrix), `gradient` (its observation-average), `hessian` (the `(p, p)` average Hessian), `hessian_full` (a length-`n` list of `(p, p)` per-observation Hessians), `scores_cov` (the `(p, p)` score covariance), and `scores_full` (the scores bundled in a dict, for parity with `Vinecop.scores_full`). Each also accepts an optional per-row `parameters` `(n, p)` array and `num_threads`, mirroring the per-row evaluation overloads ([vinecopulib#699](https://github.com/vinecopulib/vinecopulib/pull/699)).
- Add per-observation-parameter overloads to `Vinecop.pdf` / `pdf_full` / `loglik` / `scores` / `scores_full` / `gradient` / `hessian` / `hessian_full` / `scores_cov`: pass an `(n, npars)` `parameters` array (one full-vine parameter vector per observation, columns in the `(tree, edge, parameter)` order of `scores`) to evaluate at per-observation parameters instead of the fitted ones. The argument is appended to each signature, so existing positional calls are unchanged; continuous, all-parametric models only ([vinecopulib#699](https://github.com/vinecopulib/vinecopulib/pull/699)).
- Add `Vinecop.sample_conditional(u_cond, qrng=False, num_threads=1, seeds=[])`: sample from the conditional distribution of a subset of variables given fixed values of the rest. The conditioning variables are the last `k = u_cond.shape[1]` of the vine order (their columns in the output reproduce `u_cond`); discrete conditioning variables take an extra left-limit column. Built on the Rosenblatt / inverse-Rosenblatt cascade ([vinecopulib#696](https://github.com/vinecopulib/vinecopulib/pull/696)).
- Add conditioning-aware vine structure selection: `FitControlsVinecop.conditioning_set` (a list of 1-based variable labels; a settable property, also pickled) makes `Vinecop.select` place the conditioning set's own optimal sub-vine at the tail of the order, so `sample_conditional` can condition on it. `Vinecop.reorient(conditioning_set)` relabels an already-fitted vine to an equivalent one whose order tail equals a given set without refitting — value-preserving, so `pdf` / `loglik` are invariant ([vinecopulib#697](https://github.com/vinecopulib/vinecopulib/pull/697)).
- Add the list-of-trees round-trip for structures (upstream's shared `RVineTrees` primitive): `Vinecop.get_trees()` returns the fitted vine as nested `[tree][edge]` lists of `{"conditioned": (a, b), "conditioning": [...], "pair_copula": Bicop}` dicts (carrying the fitted pair copulas), `RVineStructure.get_trees()` returns the bare `(a, b, conditioning)` decomposition, and `RVineStructure.from_trees(d, trees)` is its **faithful** inverse (identity `conditioned[0]` diagonal policy, so `RVineStructure.from_trees(s.dim, s.get_trees()) == s`); `RVineStructure` also gains `__eq__`. Upstream `Vinecop.select` finalizes with the same (flip-free) diagonal convention, so `VinecopBase.select` assembles its selected trees through this same `from_trees` and matches the compiled selector's matrix byte-for-byte — one diagonal convention throughout ([vinecopulib#698](https://github.com/vinecopulib/vinecopulib/pull/698), [vinecopulib#702](https://github.com/vinecopulib/vinecopulib/pull/702)).
- `RVineStructure.from_trees` also accepts the `{"conditioned": …, "conditioning": …}` edge mappings that `Vinecop.get_trees()` returns (ignoring their `"pair_copula"`), so a fitted vine's decomposition round-trips back to a structure without being rewritten as triples first (#251).
- Add `FitControlsVinecop.from_bicop_controls(controls, ...)`, which takes a `FitControlsBicop` in place of the eleven inherited pair-copula arguments, and the `bicop_controls` property that reads and replaces them as a group (#251).
- Add `Bicop.family_name` (the family's display name) and `Bicop.as_continuous()` (a copy with `var_types` reset to continuous) (#251).
- Add `RVineStructure.get_min_array()`, `get_needed_hfunc1()` and `get_needed_hfunc2()`, returning the whole `[tree][edge]` triangular array that the existing per-entry accessors index into (#251).
- `Vinecop.pair_copulas` is now settable, so pair copulas can be replaced as a group on a fitted vine; the assignment validates the nested shape against the structure (#251).

### Build / packaging

- Keep the persistent `build-dir` out of build-isolated installs. `uv sync` built the project in a throwaway environment and baked its temporary `ninja` path into `CMakeCache.txt`; once uv deleted that environment, `editable.rebuild` invoked a binary that no longer existed and every later import died inside cmake. `make sync` now passes `--no-install-project`, the Makefile exports `UV_NO_SYNC` so no target's `uv run` silently re-does it (which also stops `make check` reinstalling the project on every invocation), and the `build_sdist` job matches (#305).
- Migrate to `uv` + `scikit-build-core` for the editable / wheel build pipeline, with `[build-system].requires` mirroring the dev `[dependency-groups]` so `--no-build-isolation` works out of the box (#209).
- Replace `mypy` with `ty` (Astral's type checker, alpha) and enable strict checks against a Python 3.10 baseline; only `pyvinecopulib.pyvinecopulib_ext` is allowed as an unresolved import (#210).
- Add `bandit` security linting to `make check` and pre-commit (scanning `src/pyvinecopulib` + `scripts`). The previously-unused `[tool.bandit]` config silently excluded everything (`"lib"` matched `pyvinecopulib`); the config is fixed and the surfaced findings resolved.
- Refactor the build / docs / examples pipeline into a thin Makefile over `uv run` and rework `scripts/regenerate_notebooks.py` (#205).
- Add a `pyvinecopulib[sklearn]` extra (`scikit-learn>=1.4`, `pandas>=2.0`, `joblib>=1.3`, `scipy>=1.10`) and a `pyvinecopulib[torch]` extra (`torch>=2.0`) (#211, #216).
- Treat `pyvinecopulib`-originated `DeprecationWarning`s as errors under pytest so internal call sites stay on the canonical import paths (#207).
- Install `--extra torch` in the notebook-test and regenerate-notebooks CI jobs so `examples/09_torch_backend.ipynb` executes under `nbmake` (#216).
- Fix osx / musllinux wheel builds: feed libclang the compiler's implicit system include dirs (plus `-ferror-limit=0` and the macOS SDK sysroot) so the C++ stdlib / intrinsic headers resolve, and abort `docstr.hpp` generation only on *fatal* libclang diagnostics (intrinsic-header `error`s are benign and no longer silently drop symbols).
- Fix Windows wheel builds against the VS 2026 / MSVC 14.51 STL: rewrite `__builtin_verbose_trap` to `__builtin_trap()` during `docstr.hpp` generation. The newer STL emits this Clang-18 builtin from core allocator / string / `call_once` headers, which the pinned PyPI libclang (≤18) can't resolve; because those headers reach nearly every translation unit, the broken parse crashed symbol extraction. Also narrow a numpy-scalar union in the `Kde1d` plot helper so `make check` (`ty`) passes (#224).
- Migrate macOS CI to `macos-15` and re-add `macos-15-intel` (`MACOSX_DEPLOYMENT_TARGET=10.13` for nanobind's aligned `new`/`delete`); wheel matrix is now 5 platforms × 3 ABI (15 wheels).
- Fix the source build against Eigen 5.x (e.g. conda-forge `eigen>=5`): its CMake package exposes the include path only through the `Eigen3::Eigen` target and no longer sets the legacy `EIGEN3_INCLUDE_DIR`, leaving `docstr.hpp` generation (and the compile) unable to find `<Eigen/Dense>`. Derive `EIGEN3_INCLUDE_DIR` from the target's `INTERFACE_INCLUDE_DIRECTORIES` when unset.
- Build the Read the Docs site from source instead of the published PyPI wheel: enable the `lib/` submodules, install header-only Boost / Eigen via apt (Ubuntu's `libboost-dev` / `libeigen3-dev` clear vinecopulib's `>=1.56` / `3.4` floors), and `pip install .[doc,examples,sklearn,torch]`. The wheel-install shortcut only rendered the last release's API, so RTD's `dev` build failed once `dev` moved ahead of PyPI (missing `numpydoc`, then the entire `core` abstraction layer plus the `sklearn` / `torch` subpackages `conf.py` documents).
- Make the Sphinx docs cross-reference-clean and enforce it going forward: enable `nitpicky` so `make docs` (`-W`, run by the `verify_docs_build` CI job) fails on any unresolved reference. The bulk of previously-dead links are fixed by giving the autosummary class template's *Attributes* block a `:toctree:` (so each property gets a page for numpydoc's member-summary links to resolve); backticked class names are rewritten to fully-qualified targets in `process_cross_references` so they link from any page; `intersphinx` mappings (numpy / torch / python) resolve external types; nanobind's `numpy.ndarray[dtype=…]` signature annotations are collapsed to `numpy.ndarray`; and private-method references are reworded to plain literals. A short `nitpick_ignore_regex` covers the few irreducible refs (upstream-C++ getter-name mismatches, `BicopFamily` value aliases, scikit-learn-generated methods).

- Replace the torch integral cache with three `(m, m)` cumulative-trapezoid prefix tables, from which `cdf` / `hfunc*` read their value in closed form. The reconstruction is exact -- a bilinear interpolant is piecewise linear along a grid line, so its integral is piecewise linear across cells -- so the cache now costs nothing in accuracy: on a `d = 10` vine of `tll` pairs it agrees with the compiled `Vinecop` to 1.9e-13 on `pdf`, where the bilinear cache reached 2.6 relative on a single `hfunc1`. Three cumulative sums build it, against the previous O(m^4) pass (#305).
- Add `TorchBicop.rect_mass(a1, b1, a2, b2)`, the exact probability of a rectangle, available in both cache modes. It is the value a four-corner `cdf` difference defines, arranged so that almost none of it cancels: that difference amplifies an absolute error by `~4/(w1 w2)` in the widths, where `rect_mass` amplifies by `1/w2` alone. Against exact rational truth on a `1.2e-4`-wide rectangle, `2.9e-12` against `8.7e-9`. `values >= 0` is now a constructor precondition, since the nonnegative-weight bound depends on it (#305).
- Expose `utils.find_latent_sample(u, b, niter=3)`: recovers a continuous sample from interval-censored copula data, the transform a nonparametric fit on discrete margins runs on. The draw is deterministic and invariant to argument order (#305).
- `utils.wdm` accepts Chatterjee's xi (`"chatterjee"` / `"cxi"` / `"xi"`), the one asymmetric measure it offers -- it measures how far `y` is a function of `x`. `FitControlsVinecop.tree_criterion` accepts it too, spelled `"cxi"`, and symmetrizes it as `max(xi12, xi21)` because a spanning-tree edge weight has to be symmetric ([wdm#15](https://github.com/tnagler/wdm/pull/15), [vinecopulib#754](https://github.com/vinecopulib/vinecopulib/pull/754)).
- `Vinecop.reorient` and the `conditioning_set` argument of `rosenblatt` / `inverse_rosenblatt` accept a truncated model instead of throwing: the trees above the truncation are independence, so the peel has nothing to move there ([vinecopulib#743](https://github.com/vinecopulib/vinecopulib/pull/743), [vinecopulib#752](https://github.com/vinecopulib/vinecopulib/pull/752)).
- `utils.wdm` takes `seeds`, the tie-breaking seeds Chatterjee's xi applies to `x`. Ordering tied predictors by the response would manufacture dependence, so they are broken at random -- but from a fixed default, so xi is a function of its arguments: it previously drew from `std::random_device` and returned a different number on every call for tied data. Untied predictors never construct the generator, so continuous data is unaffected ([wdm#26](https://github.com/tnagler/wdm/pull/26)).
- `utils.wdm` raises on weights whose sum is not finite and positive, where it returned `NaN`. There is no weighted measure to compute from weights that carry no mass, and the `NaN` propagated silently into whatever consumed it ([wdm#21](https://github.com/tnagler/wdm/pull/21)).
### Bug fixes in `pyvinecopulib`

- Stop the bindings from holding a live instance per default argument. `Bicop.from_data` / `fit` / `select` and `Vinecop.from_data` / `fit` / `select` took a constructed `FitControls*` as their default, which nanobind keeps alive in the function record for the module's lifetime and reports at interpreter shutdown as "leaked instances ... likely caused by a reference counting issue in the binding code". They now take `None` and construct inside, so `controls` accepts `None` as well as a `FitControls*` and the rendered default is unchanged (#305).
- Make the unit-interval trim survive the working precision. It was the literal `1 - 1e-10`, which rounds to exactly `1.0` in `float32`, so on the advertised `device="cuda", dtype=torch.float32` configuration the clamp admitted the value it exists to exclude: `hfunc1` / `hfunc2` / `cdf` returned exactly `1.0` at the boundary, and `ndtri` of that is an infinity for whatever refits on it. `core._trim.trim_bounds` derives the bound from `finfo(dtype).eps` and returns the historical pair unchanged at `float64`, so every parity gate is unmoved (#305).
- Scale `_ace`'s outer tolerance with the dtype. `2e-15` is upstream's `tools_stats::ace` value and is seven orders below `float32` eps, so a `float32` fit iterated to `outer_iter_max` against its own rounding noise instead of converging (observed at rho = -0.7, n = 200). `float64` keeps the literal (#305).
- `TorchBicop.sample(seeds=...)` no longer calls `torch.manual_seed`, which reseeded the global CPU generator and every CUDA one as a side effect of seeding a single pair copula. It uses a device-local `torch.Generator`, matching `TorchVinecop._sample_uniform` (#305).
- `VinecopBase.select` no longer defaults its host transfer to `np.asarray`, which raises on a GPU tensor -- so the documented `VinecopBase` extension point could not select a structure from device-resident data without passing `to_numpy` by hand (#305).
- `BicopBase.plot` moves the density off the device before reshaping it, so plotting a CUDA-resident `TorchBicop` works instead of raising (#305).
- Give `InterpolationGrid2D._interval_weights`' zero branch the grid's dtype and device; it was a bare `torch.zeros(())`, i.e. a cpu `float32` scalar mixed into a `torch.where` over `float64` device tensors (#305).
- Keep `TorchVinecop`'s batched cache out of `state_dict()`, so a checkpoint taken after a `batched=True` call no longer carries derived buffers that `load_state_dict` then rejects on a fresh model (#264).
- Reject non-finite `X` / `y` in the `pyvinecopulib.sklearn` estimators instead of passing them to `Kde1d`, which reads NaN as a segmentation fault (#263).
- Port the `integrate_2d` marginal-renormalization fix to the torch backend (`InterpolationGrid2D.integrate_2d` and `integrate_2d_batched`) so `TorchBicop.cdf` enforces ``C(1, u_2) = u_2`` exactly, matching the post-vinecopulib#667 C++ CDF to machine precision on the on-the-fly path ([vinecopulib#667](https://github.com/vinecopulib/vinecopulib/pull/667)).
- Declare `BicopFamily` enum members as class attributes in the generated type stubs so the documented `pv.BicopFamily.clayton` access pattern passes static type checking (`ty` / pyright / mypy), not only the module-level constants (#223).
- Declare base classes in the generated type stubs, so `DVineStructure` and `CVineStructure` are recognized as `RVineStructure` subclasses. Passing either where an `RVineStructure` is expected -- the documented way to build a C- or D-vine -- no longer fails static type checking, although it always worked at runtime (#268).
- Drop `@torch.no_grad()` from the torch backend's cached-grid lookup and from `TorchBicop.hinv1` / `hinv2`. With `cache_integrals=True` the cached `cdf` / `hfunc*` / `hinv*` were gradient-dead, so `TorchVinecop.pdf`'s gradient in `u` was off by 111% of the largest entry on the non-batched path -- while the batched path, which never went through the cached lookup, was right. Both now agree with central differences to 1e-10, and the batched state is re-baked when grad tracking changes under it (#305).
- Re-bake `TorchVinecop`'s batched state when grad tracking changes under it. The bake copies each pair's grid, and `requires_grad_` mutates a flag in place, so calling `pdf(batched=True)` before `values.requires_grad_(True)` left a detached copy behind and `backward()` raised torch's generic "does not require grad". Device moves already invalidated; this is the same for the other thing a caller changes after fitting. Costs two flag reads per pair, 0.37% of a batched `pdf` at `d = 12`, `n = 5000` (#302).

### Dependency changes

- `numpy>=2.0` is now a project-wide requirement (was `>=1.14`); `VineRegressor` needs `np.quantile(weights=...)` from NumPy 2.0 (#211).
- `[sklearn]` extra adds `pandas>=2.0` (used by `VineBase.expand_factors` for DataFrame inputs) (#211).

### Changes in `vinecopulib`

#### BREAKING API CHANGES

- Require C++17, CMake 3.14 and Boost 1.75, and put `-march=native` behind `VINECOPULIB_NATIVE_ARCH` so the default release build is redistributable ([vinecopulib#711](https://github.com/vinecopulib/vinecopulib/pull/711)). pyvinecopulib wheels now set `-march=x86-64-v3` explicitly; editable builds take the plain baseline (#250).
- Remove `Vinecop::select_all`, `Vinecop::select_families` and the `*_truncation_level` accessors, deprecated since 0.3.1 ([vinecopulib#718](https://github.com/vinecopulib/vinecopulib/pull/718)). None were reachable from Python.
- The umbrella header no longer disables Boost's concept assertions ([vinecopulib#714](https://github.com/vinecopulib/vinecopulib/pull/714)).

#### BEHAVIOR CHANGES

- Kendall's τ of `bb6`, `bb7`, `bb8` and `tawn` changes: four numerical defects in `parameters_to_tau` are fixed, the worst of which returned about `1e-11` where the true value was `0.33` ([vinecopulib#713](https://github.com/vinecopulib/vinecopulib/pull/713)). `Bicop.tau`, `parameters_to_tau`, `str()` and family selection move for those four families.
- Setting `FitControlsVinecop.tree_criterion_function` while `tree_criterion != "custom"` now raises instead of being silently ignored ([vinecopulib#722](https://github.com/vinecopulib/vinecopulib/pull/722)).
- `Bicop.loglik` and `Bicop.fit` validate the column count and raise on a wrong one, where they previously read past the data ([vinecopulib#729](https://github.com/vinecopulib/vinecopulib/pull/729)).
- A `Vinecop` built from a full structure with no pair copulas treats the omitted ones as independence instead of indexing an empty store ([vinecopulib#729](https://github.com/vinecopulib/vinecopulib/pull/729)).

#### NEW FEATURES

- `Bicop.from_data` accepts discrete data. It took a statically two-column matrix, so the `n x (2 + k)` and `n x 4` layouts a discrete pair needs could not be passed at all ([vinecopulib#729](https://github.com/vinecopulib/vinecopulib/pull/729)).
- `Vinecop.sample_conditional` accepts an explicit `conditioning_set` and both discrete `u_cond` layouts — expanded `(n, 2k)` and compact `(n, k + k_d)` ([vinecopulib#729](https://github.com/vinecopulib/vinecopulib/pull/729)). Note that the two forms map `u_cond`'s columns differently: implicitly by the vine's order tail, explicitly by the given set.
- `Vinecop.rosenblatt` and `Vinecop.inverse_rosenblatt` accept an explicit `conditioning_set`, holding the named variables fixed instead of the ones at the tail of the vine order, and evaluating through a reoriented non-owning view rather than copying pair copulas ([vinecopulib#715](https://github.com/vinecopulib/vinecopulib/pull/715)). The argument is keyword-only, so `rosenblatt(u, 4)` still means `num_threads=4`.
- `Bicop.sample` draws one observation per parameter set: pass `parameters` (an `(n, p)` array) instead of `n`, whose row count then fixes the sample size ([vinecopulib#719](https://github.com/vinecopulib/vinecopulib/pull/719)). `parameters` is keyword-only, so the positional `sample(n, qrng, seeds)` keeps its meaning, and rejects a 1-d array, which would otherwise be read as a single row. `num_threads` applies only to this form.
- Early exit in vine selection when the structure is already a tree, avoiding redundant work in `select` ([vinecopulib#661](https://github.com/vinecopulib/vinecopulib/pull/661)).
- Per-family parameter / rotation / tail-dependence documentation on the `BicopFamily` enum members, surfaced through the Python `families` subpackage (#214, [vinecopulib#668](https://github.com/vinecopulib/vinecopulib/pull/668)).
- Numpydoc-compliant `//!` comments on every property getter / setter in the Python-binding surface, surfaced through the pyvinecopulib autosummary pages (#214, [vinecopulib#670](https://github.com/vinecopulib/vinecopulib/pull/670)).
- Allow a custom `tree_criterion` function for vine structure selection ([vinecopulib#674](https://github.com/vinecopulib/vinecopulib/pull/674)).
- Per-row parameter evaluation for parametric bivariate copulas: `Bicop::pdf` / `cdf` / `hfunc1` / `hfunc2` / `hinv1` / `hinv2` / `loglik` gain an overload taking an `n×p` parameter matrix (one set per observation) plus an optional thread count, backed by an eval-core refactor to a single parameter-aware leaf per family ([vinecopulib#675](https://github.com/vinecopulib/vinecopulib/pull/675)).
- Improve start parameters when fitting pair copulas on discrete data ([vinecopulib#677](https://github.com/vinecopulib/vinecopulib/pull/677)).
- `TriangularArray<T>` owns its conversions: `to_json()`, a JSON constructor, and `to_list()` (nested rows), used by the pyvinecopulib bindings for the structure-array and per-edge outputs ([vinecopulib#680](https://github.com/vinecopulib/vinecopulib/pull/680)).
- Analytic derivatives of the bivariate copula density and h-functions with respect to parameters and arguments for all parametric families, making `Vinecop` scores and Hessians analytic (exact and faster); models with nonparametric pair copulas are now rejected by `scores` / `hessian` ([vinecopulib#683](https://github.com/vinecopulib/vinecopulib/pull/683), [vinecopulib#687](https://github.com/vinecopulib/vinecopulib/pull/687)). The TLL family additionally gains the exact argument gradient of its density ([vinecopulib#694](https://github.com/vinecopulib/vinecopulib/pull/694)).
- Rename the `Vinecop` Hessian API: `hessian_avg` → `hessian` (average) and `hessian` → `hessian_full` (per-observation), mirrored by the Python bindings ([vinecopulib#679](https://github.com/vinecopulib/vinecopulib/pull/679)).

#### PERFORMANCE

- Speed up the bicop evaluation engine: vectorized closed-form pdf / h-function / CDF leaves and derivative-cascade allocation hygiene. Evaluations may shift by ≤ 1e-12 and fitted parameters by ≤ 1e-8 ([vinecopulib#681](https://github.com/vinecopulib/vinecopulib/pull/681)).
- Speed up `Vinecop` evaluation and structure selection: no per-edge copula copies in `inverse_rosenblatt`, in-place data collapsing, parallel allocation-free Monte-Carlo `cdf`, and selection fast paths. `pdf_full` no longer duplicates continuous h-functions into the `hfunc*_sub` buffers (they are now empty for models without discrete variables) ([vinecopulib#692](https://github.com/vinecopulib/vinecopulib/pull/692)).
- Speed up TLL fitting and evaluation: fused conditional-cdf interpolation and closed-form inversion of the conditional cdf. TLL `hinv1` / `hinv2` (hence `sample` / `inverse_rosenblatt` on TLL vines) may shift by ≤ 1e-9 ([vinecopulib#691](https://github.com/vinecopulib/vinecopulib/pull/691)).
- Speed up `tools_stats`: SIMD `qnorm`, leaner bivariate normal / t kernels, faster pseudo-observations and quasi-random fills. Gaussian / Student evaluations may shift at the ≤ 1e-12 level ([vinecopulib#690](https://github.com/vinecopulib/vinecopulib/pull/690)).
- Faster shared Eigen / thread / integration primitives; the relaxed integration tolerance (1e-12 → 1e-9) may shift Kendall's τ of the integrated families (BB6 / BB7 / BB8, Tawn) by up to ~1e-7 ([vinecopulib#689](https://github.com/vinecopulib/vinecopulib/pull/689)).
- Compute the Student t score's df-only terms once per call, speeding up `Vinecop` scores / Hessians and Student t maximum-likelihood fits ([vinecopulib#693](https://github.com/vinecopulib/vinecopulib/pull/693)).

#### BUG FIXES

- Fix an out-of-bounds read in the mixed discrete/continuous h-inverse (`hinv2` for `{d, c}` copulas with a numeric inverse: Frank, BB1/6/7/8, Tawn) ([vinecopulib#675](https://github.com/vinecopulib/vinecopulib/pull/675)).
- Add a missing thread-related include ([vinecopulib#676](https://github.com/vinecopulib/vinecopulib/pull/676)).
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
- Conversion of all examples to Jupyter notebooks (#160)

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
