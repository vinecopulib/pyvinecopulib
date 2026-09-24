# Changelog

## 1.0.1 (unreleased)

### Bug fixes in `pyvinecopulib`

- Evaluate a discrete or mixed `Bicop` or `Vinecop` read from JSON, a file or a pickle with its variable types, not as continuous (#352, [vinecopulib#789](https://github.com/vinecopulib/vinecopulib/pull/789)).
- Fix a crash in discrete fits and `find_latent_sample` when built against Eigen 5 (#352, [vinecopulib#792](https://github.com/vinecopulib/vinecopulib/pull/792)).
- Keep a `Kde1d` density positive between the smallest and largest observation; a wide gap, or ties at `degree=2`, gave exactly zero (#352, [kde1d#42](https://github.com/vinecopulib/kde1d-cpp/pull/42)).
- Select the same `Kde1d` bandwidth on every build; on tied data the plug-in estimate was `NaN` on some builds only (#352, [kde1d#42](https://github.com/vinecopulib/kde1d-cpp/pull/42)).
- Make a `degree=1` `Kde1d` fit scale equivariant, which changes every `degree=1` estimate (#352, [kde1d#42](https://github.com/vinecopulib/kde1d-cpp/pull/42)).
- Break ties the same way on every platform in `to_pseudo_obs(..., "random")`, Chatterjee's xi and every `tll` fit, and keep input order for `"first"`; tied-data results change once (#352, [wdm#28](https://github.com/tnagler/wdm/pull/28)).
- Compute Hoeffding's D correctly on tied data, in `wdm(..., "hoeffding")` and the `"hoeffd"` tree criterion; it could leave the unit interval (#352, [wdm#28](https://github.com/tnagler/wdm/pull/28)).
- Break ties in `TorchTllBicop.from_data` and `from_data_batched` as `Bicop.from_data` does, at random from a fixed seed rather than by row order (#351).

### Build / packaging

- Version an open cycle as `X.Y.Z.devN`, enforced by `scripts/check_version.py` (#352).
- Find conda-forge's `libclang-<major>.dll` on Windows (#352).
- Test against Eigen 5 in CI (#352).

### Dependency changes

- Bump `lib/vinecopulib` past `v1.0.0`, `lib/kde1d` to `v1.2.3` and `lib/wdm` past `v0.3.0` (#352).

## 1.0.0 (2026-09-18)

This release ships the whole vine-modeling stack rather than the copula half of
it. `Vinedist` and the new `pyvinecopulib.margins` layer put Sklar's theorem on
the data scale; `pyvinecopulib.torch` runs the evaluation cascade on the GPU
under autograd; `pyvinecopulib.sklearn` adds scikit-learn-compatible estimators;
and the array-agnostic layer in `pyvinecopulib.core` hosts custom pair copulas
and non-simplified vines. The public API is reorganized into the `core`,
`families`, `utils` and `margins` subpackages, with the pre-reorganization
top-level names kept as warning aliases until 2.0.

It also advances all three vendored C++ libraries, so nearly every `tll` and
`Kde1d` fit moves. Read the migration guide in `docs/migrating-to-1.0.md` before upgrading from
0.7.x.

### Breaking API changes in `pyvinecopulib`

- Use one argument order on every estimator: the observations first, then `controls`, then keyword-only whatever the object cannot infer. `Vinecop.from_data` took `controls` fifth, behind `structure`, so the call carried over from `fit` bound a controls object as a structure (#326).
    - `pv.Vinecop.from_data(u, structure, matrix, var_types, controls)` -> `pv.Vinecop.from_data(u, controls, structure=..., matrix=..., var_types=...)`
    - `pv.Bicop.from_data(u, controls, var_types)` -> `pv.Bicop.from_data(u, controls, var_types=...)`
    - `pv.core.Kde1d.fit(x, weights)` -> `pv.core.Kde1d.fit(y, FitControlsKde1d(weights=...))`, and likewise `select`; the kernel knobs and the weights are all `FitControlsKde1d`, and the observations are `y`, because `x` means exogenous covariates everywhere in the Python API (#339)
    - likewise on `BicopBase`, `VinecopBase`, `VinedistBase` and their torch subclasses; `MarginBase` already read this way

- Reorganize the public API into the `core`, `families` and `utils` subpackages; the family constants, `Kde1d` and the utility functions still resolve at the top level, but warn on access and are removed in 2.0 (#207, #292).
    - `pyvinecopulib.gaussian` -> `pyvinecopulib.families.gaussian`, and likewise every family constant and group
    - `pyvinecopulib.Kde1d` -> `pyvinecopulib.core.Kde1d`
    - `pyvinecopulib.wdm` -> `pyvinecopulib.utils.wdm`, and likewise `sobol`, `ghalton`, `benchmark` and `pairs_copula_data`
- Default `selection_criterion` to `"aic"` instead of `"bic"` on `FitControlsBicop` and `FitControlsVinecop`, matching the C++ and R defaults; pass `selection_criterion="bic"` to keep the previous selection (#251).
- Move every `tll` fit, and so every unqualified `Bicop.from_data` / `Vinecop.from_data` that selects one: the interpolation grid's margins are balanced across both sweep orders and iterated to convergence, and `hfunc` / `hinv` no longer floor the interpolated density at `1e-4` (#305, [vinecopulib#751](https://github.com/vinecopulib/vinecopulib/pull/751)).
- Move every `Kde1d` fit: the plug-in bandwidth shifts by `O(1/n)`, and a bounded fit shifts again under the new endpoint treatment (#312, [kde1d#29](https://github.com/vinecopulib/kde1d-cpp/pull/29)).
- Fit a finite `Kde1d` bound with a local-linear boundary estimator by default; pass `boundary_repair=False` for the previous transformed fit (#312, [kde1d#36](https://github.com/vinecopulib/kde1d-cpp/pull/36)).
- Read a discrete `Kde1d`'s `xmin` / `xmax` as the *integer support*, so both the bounds and the data must be integers and the fitted grid runs half a unit wider at each end (#312, [kde1d#37](https://github.com/vinecopulib/kde1d-cpp/pull/37)).
- Label a selected `RVineStructure` with the conditioned variable on the diagonal, so `get_matrix()` and `order` differ for the same model; densities and log-likelihoods do not (#251, [vinecopulib#702](https://github.com/vinecopulib/vinecopulib/pull/702)).
- Rename the sampling entry points to `sample`, keeping the released spellings as warning aliases (#297).
    - `Bicop.simulate(n)` -> `Bicop.sample(n)`, and likewise on `Vinecop`, `RVineStructure` and `Kde1d`
    - `pyvinecopulib.simulate_uniform` -> `pyvinecopulib.utils.sample_uniform`
- Return `self` from `Bicop.fit` / `.select` and `Vinecop.fit` / `.select` instead of `None`, so every estimator in the package composes the same way (#326).
- Make `parameters`, `num_threads`, `seeds` and `randomize_discrete` keyword-only on the evaluation methods of `Bicop` and `Vinecop`: each class was internally consistent and the two disagreed, `Bicop` reading `(u, parameters, num_threads)` where `Vinecop` read `(u, num_threads, parameters)`, so a call carried between them bound a parameter matrix as a thread count. `Bicop(family, rotation, parameters)`, `from_family` and `parameters_to_tau` keep theirs positional, where `parameters` specifies the model rather than the call. Keyword calls are unaffected (#345).
    - `cop.pdf(u, pars)` -> `cop.pdf(u, parameters=pars)`
    - `cop.sample(n, False, [1, 2])` -> `cop.sample(n, False, seeds=[1, 2])`
- Name the observations `u` on `Bicop` / `Vinecop`'s `fit` / `select` / `from_data`, which is what every evaluation method on the same class already called them. Renamed upstream rather than in the binding, whose docstrings are lifted verbatim ([vinecopulib#781](https://github.com/vinecopulib/vinecopulib/pull/781), #345).
    - `cop.fit(data=u)` -> `cop.fit(u=u)`, and likewise on `select` / `from_data`
- Drop `Vinecop.fit`'s `num_threads` argument, which duplicated `FitControlsBicop.num_threads` (#326).
    - `vine.fit(u, controls, num_threads=4)` -> `vine.fit(u, FitControlsBicop(num_threads=4))`
- Rename `Kde1d`'s `quantile` to `icdf`, the name modern SciPy and `torch.distributions` use for the inverse distribution function, with no alias (#292).
- Rename `Kde1d`'s `edf` to `npars`, the name `Bicop` and `Vinecop` already answer for the same quantity, with no alias. It was bound twice -- `edf` and `n_parameters` were the same getter -- so one class offered two live names for one number; `npars` is now the only spelling, on the compiled classes, on `MarginBase` and on the four contracts (#345).
    - `kde.edf` -> `kde.npars`
    - `kde.quantile(p)` -> `kde.icdf(p)`
- Serialize every margin the same way: `to_json` returns JSON **text** on `SciPyMargin`, `TorchKde1d` and `TorchDistributionMargin` as it already did on `Bicop` / `Vinecop` / `Kde1d` / `Vinedist`, and the reader is `from_json(json)` rather than `from_json_payload(payload)`. A margin's own `to_json` builds its mapping and hands it to `margin_json`, which stamps `kind` and `version` and runs the shared codec -- so a non-finite parameter still travels as a tagged string and a subclass still reads back. `register_margin_json` takes a `Callable[[str], Any]` (#345).
    - `margin.to_json()` returns `str`, not `dict`
    - `Cls.from_json_payload(payload)` -> `Cls.from_json(json)`
- Rename `Kde1d`'s `type` to `var_type`, on the property and on every constructor, and spell a variable's type `"c"` / `"d"` / `"zi"` throughout -- the one spelling `var_types` already used on `Bicop` and `Vinecop`. Both constructors still accept the long names, so only the attribute is a break; nothing in the package emits one (#339).
    - `pv.core.Kde1d(type="continuous")` -> `pv.core.Kde1d(var_type="c")`, and likewise on `from_params` / `from_grid`
    - `kde.type` -> `kde.var_type`, which answers `"c"` / `"d"` / `"zi"` where `type` answered `"continuous"` / `"discrete"` / `"zero_inflated"`
- Rename `Bicop.plot`'s first parameter from `type` to `plot_type`, which is what its own rendered documentation has always called it -- so the call the docs describe raised -- and what `BicopBase.plot` takes, so one call carries between the two (#330).
    - `cop.plot(type="contour")` -> `cop.plot(plot_type="contour")`
- Make `Kde1d.loglik` a method taking optional data rather than a property, matching `Bicop.loglik` and `Vinecop.loglik` (#292).
    - `kde.loglik` -> `kde.loglik()`
- Serve `repr` and pickle from the canonical module paths (`Bicop`'s `__module__` is now `pyvinecopulib.core`); pickles written by 0.7.x still load through the deprecated aliases (#207).
- Drop Python 3.9 and 3.10: `requires-python` is now `>=3.11` (#207, #292).

### New features in `pyvinecopulib`

#### Vine distributions and margins

- Add `Vinedist`, Sklar's theorem as an object: any vine copula plus one margin per variable, with `pdf` / `logpdf` / `cdf` / `loglik` / `sample` / `sample_conditional` / `rosenblatt` / `inverse_rosenblatt` on the data scale, for continuous, discrete and mixed margins alike (#292).
- Add `pyvinecopulib.margins`, the univariate half of a vine distribution, so fitting a vine copula to pseudo-observations becomes one configuration of a vine distribution rather than a separate workflow (#292).
- `Vinedist.from_data(y, controls, *, margin_controls=, var_types=, supports=)` runs the two-step (IFM) estimator, margins first and the copula on the pseudo-observations they produce. Which *class* each margin is comes from the distribution's `margin_class`, beside the `vinecop_class` naming its other half (#292, #339).
- `Kde1d` *is* the default margin, with no wrapper, and gains `var_type` / `support` / `cdf_left` / `logpdf` / `n_parameters` / `family_name` / `is_fitted` plus `fit` / `select` / `from_data` that read like every other estimator's and return `self` (#292, #339).
- Add `FitControlsKde1d`, the kernel knobs (`multiplier`, `bandwidth`, `degree`, `grid_size`, `boundary_repair`) plus `weights` as a controls object, so a `Kde1d` margin is configured per variable through `margin_controls=` like any other (#339).
- Observation weights ride in the controls at every level, as they already did on `FitControlsBicop` / `FitControlsVinecop`: `FitControlsKde1d` and `FitControlsMargin` carry a `weights` field too, and no `fit` / `select` / `from_data` in the package takes a `weights=` argument. Nothing copies weights from one controls object to another, which is what lets a vine distribution weight its margins and its copula differently -- or weight one and leave the other alone. `supports_weights` says whether a class honors `controls.weights`, and one that declares `False` refuses them rather than fitting them away (#339).
- Add parametric margins behind a new `pyvinecopulib[scipy]` extra: `SciPyMargin` wraps one `scipy.stats` family, or selects one from a curated candidate set by AIC / BIC / AICc (#292, #326).
- Add `as_margin`, which presents another ecosystem's distribution object as a margin, and `register_margin_adapter` for one it does not know (#292).
- Add `Vinedist.margin_summary()`, one row per variable naming the margin it ended up with, and `None` for any field a margin declines -- so a margin from another ecosystem contributes what it has (#292, #334).
- `Kde1d.plot` draws the distribution function too, with `kind="cdf"` (#330).

#### The array-agnostic extension layer

- Give the four canonical bases one fitting surface: `fit` returns `self`, `from_data` constructs, and `BicopBase` / `VinecopBase` add `select` where there is a family or a structure to choose (#326).
- Declare the parts instead of hooking them: `VinecopBase.bicop_class` names the pair copula a vine fits, and `VinedistBase.vinecop_class` / `margin_class` its two halves, so `from_data` needs no callback -- a part class is itself a fitter. `controls_class` names the fit configuration each reads, so what `controls=None` means is answered once per class rather than at every site that had to construct one (#326, #339).
    - naming the pair class also lets `select` refuse one without `flip` before it reads the data, instead of after fitting an edge
- Derive `FitControlsVinecop` from `FitControlsBicop` and `FitControlsTorchVinecop` from `FitControlsTorchBicop`, so one controls object configures both halves of a vine fit: a vine reads the settings it owns and the rest reach its pair copulas unchanged, `isinstance` included (#251, #326).
- Weight the array-agnostic tree criterion, so `VinecopBase` selection honors `controls.weights` and reproduces `Vinecop.select` exactly (#326).
- Support `tree_criterion="custom"` on the array-agnostic selector, which reaches the same criterion binding `Vinecop.select` does and receives the exogenous covariates `x` by keyword (#326).
- Refuse a PyTorch copula or margin in `Vinedist`, mirroring `TorchVinedist`'s refusal of a core `Vinecop`: either mix silently loses gradients and device placement (#326).
- Add `ControlsLike`, the fit-configuration contract — anything with `to_dict()` — and `to_dict()` on all four `FitControls*` classes. A consumer reads the settings it owns and refuses one it cannot honor rather than dropping it: the array-agnostic selector rejects `select_trunc_lvl`, `select_threshold`, `select_families` and `show_trace`, which it does not implement and cannot delegate to a pair-copula fit (#326).
- Add `VinedistLike` and `VinedistBase`, so a custom vine distribution has a contract and a base to subclass rather than the concrete class; `Vinedist` is now the NumPy subclass, and `TorchVinedist` derives from the base instead of from it (#326).
- `VinedistBase.from_data` owns the two-step (IFM) estimator for both array namespaces. A subclass that names its two part classes writes no fitting code at all: the one hook it may still need is `_coerce_fit_data`, where a lane puts the observations on its own array namespace (#326, #339).
- A vine distribution declares no capabilities of its own; it asks its parts. Weighted controls handed to a part that declares `supports_weights = False` are refused rather than fitted away, and covariates are refused by `from_data` when neither the margins it is about to build nor its `vinecop_class` reads any. So composing two parts that honor weights is all it takes to honor them, and there is no flag on the distribution to keep in step with the parts (#326, #339).
- Add `Vinedist.copula_layout`, the copula-scale layout a consumer of the fitted distribution needs (#326).
- Add `select` to all four bases, defaulting to `fit` where there is nothing to choose -- upstream's own `select_families = false` equivalence -- so `from_data` is `cls().select(...)` throughout, as `Bicop`'s data constructor has always been (#326).
- Add `VinedistBase.fit` and `.select`: `fit` re-estimates both halves along the structure and families it holds, `select` lets both change shape. Each margin is asked for its own `fit` or `select`, so it keeps its class and its family exactly as a hosted copula keeps its pairs; a margin with no estimator of its own is left as it was built. `fit` refuses a `family_set` it cannot honor rather than dropping it, and neither verb may change which variables have atoms -- the held copula carries its own `var_types`, so that is `from_data`'s to decide (#326, #339).
- Add `TorchTllBicop.fit`, which refits the density grid in place keeping the module's device, dtype and cache mode (#326).
- Select a margin's family with `select` on the margin class that owns the families, the shape `Bicop.select` has: after the call the margin *is* the winning family. Naming a family is itself the choice, so `select` on a named margin reduces to `fit` and `family_set` is how a caller asks for the search back (#292, #326).
- Add `FitControlsMargin` and `margin_controls=`, the marginal half of a vine-distribution fit -- `family_set`, `selection_criterion`, `on_failure`, `weights` -- resolved per variable four ways: one object broadcast, a length-`d` sequence, or a mapping keyed by position or name (#326, #339).
- Declare a variable rather than configure it: `var_types=` and `supports=` are keyword-only on `VinedistBase.from_data` and on every margin's `fit` / `select` / `from_data`, one entry per variable, exactly as `var_types` is on `Bicop.from_data`. So one call bounds the variables whose bounds are known and leaves the rest alone, without a controls object carrying a claim about the data (#326, #339).
- Add `MarginBase.aic` / `.bic` / `.aicc`, which `Bicop` and `Vinecop` have always had and margins did not (#326).
- Plot a custom object from its base: `BicopBase.plot`, `VinecopBase.plot` and `MarginBase.plot` draw what `Bicop.plot`, `Vinecop.plot` and `Kde1d.plot` draw, so a hand-written pair copula, vine or margin is inspected like a fitted one (#327, #330).
    - `BicopBase.plot` and `MarginBase.plot` take an optional single-row `x`, since a conditional object is a different surface at every covariate value and a plot shows one slice; both refuse an `x` their object reads no covariates from (#327, #330)
    - `MarginBase.plot` covers all three variable types, drawing a curve, marks on an integer support, or both for a zero-inflated variable (#330)
- Host any pair copula in a `TorchVinecop`, not only a grid one: a `BicopBase` that is also an `nn.Module`, with learnable parameters, composes and trains. `examples/10_extending_pyvinecopulib.ipynb` walks through it (#326).
- Re-export `FitControlsMargin` from `core` and the top level, so all three controls classes sit together where a caller looks for them (#326).
- Remove `benchmark`, a timing harness with no caller or test -- the one 0.7.6 top-level name this release does not keep, even as a warning alias (#326).
- Add the array-agnostic contracts `BicopLike`, `VinecopLike`, `MarginLike` and `VinedistLike` to `pyvinecopulib.core`, which `Bicop`, `Vinecop`, `Kde1d` and `Vinedist` satisfy structurally, so downstream code can type against the contract instead of a concrete class. Each declares everything the library may ask of that part -- the evaluation cascade, plus what hosting it needs: `cdf` and `with_var_types` for a **discrete** edge, `flip` for structure selection, and the `supports_covariates` / `supports_weights` a consumer reads before handing over a covariate matrix or weighted controls. A member an implementation cannot serve has a default in the protocol body, which is the fallback its consumer used to apply -- `logpdf` is `log(pdf)`, `var_type` is `"c"`, a capability flag is `False`, and `flip` / `cdf` / `sample` raise naming the class. So a class implementing a protocol directly is hostable, rather than passing `isinstance` and failing at the first thing the library needs (#236, #237, #265, #292, #326, #339).

- Add the canonical partial implementations `BicopBase`, `VinecopBase`, `MarginBase` and `VinedistBase`, which run on NumPy or PyTorch: a custom pair copula defines `pdf` / `hfunc1` / `hfunc2`, a custom vine the single `get_pair_copula` hook, a custom margin `pdf` / `cdf`, and a custom distribution nothing at all beyond its two halves (#236, #237, #292, #326).
    - all four take incoming arrays through one `_prep` hook whose default *infers* the namespace, dtype and device from the arrays they already hold, so a subclass on PyTorch writes no conversion code — including for `BicopBase.plot`'s evaluation grid, the one array a base manufactures itself (#327, #334)
    - `pyvinecopulib.core.extend` holds the six names writing that override needs: `place` and `reference_array`, `prepare_covariates` and `covariate_row`, `to_numpy`, and the `NotBatchable` a `_build_batched` override raises to decline the batched cascade (#334)
    - an object holding no array has nothing to infer from and gets its values back untouched; a torch class that is no `nn.Module` declares its `device` and `dtype` instead, through `TensorPlacementMixin` (#334)
- `VinecopBase` ships `fit`, `select` and `from_data`, an array-agnostic port of `Vinecop`'s Dissmann and Wilson structure selection whose selected matrix matches `Vinecop.select`'s exactly. All three place their observations, so a CUDA vine accepts NumPy input; the fitted pairs are stored through the `set_pair_copulas` hook, which drops the batched and compiled cascades through `_invalidate_batched`; and `fit_edge` / `fit_level` are keyword-only on all three (#237, #244, #317, #326, #328, #334).
- Build non-simplified / conditional vines with `ConditioningContext` and its `SimplifiedContext` (default) and `NonSimplifiedContext` policies (#237, #328).
    - `VinecopBase.select` and `.from_data` honor the context while they fit, so a selected conditional vine is estimated as the model it evaluates as; each slot carries the conditioning order its own pair was fitted on, which the finalizing `flip` would otherwise permute (#328)
    - every estimator on `BicopBase` and `VinecopBase` takes an optional `x`, and one that cannot condition on it refuses it rather than fitting the unconditional model (#328)
    - an `x` is placed onto the object's own namespace but never clamped, and a one-dimensional one is refused: the shape says nothing about which axis is which, so the message names both readings (#334)
- A pair copula carries its own `var_types`: `BicopBase.with_var_types` declares a slot discrete, and the base dispatches `pdf` / `cdf` / `hfunc*` / `hinv*` over the continuous `_pdf_raw` / `_hfunc1_raw` / `_hfunc2_raw` leaves a subclass writes, mirroring `AbstractBicop` (#306, #330, #339).
- Add `core.IndependenceBicop`, the pair a `VinecopBase.select` edge below `threshold` holds instead of a fit (#317, #330).
- Every `Vinedist` and margin method takes optional exogenous covariates as a keyword-only `x`, forwarded to each part that declares `supports_covariates` (#292, #328).

#### Discrete and conditional models

- Support discrete and mixed variables end to end on every layer: `Bicop.from_data` takes the discrete layouts, `VinecopBase` and `TorchVinecop` take `var_types`, and `Vinedist` builds the compact `(n, d + k)` layout itself (#257, #292, #306, [vinecopulib#729](https://github.com/vinecopulib/vinecopulib/pull/729)).
- Add `Vinecop.sample_conditional`, which holds a subset of the variables fixed and draws the rest, accepting an explicit `conditioning_set` and both discrete `u_cond` layouts (#246, #256, [vinecopulib#696](https://github.com/vinecopulib/vinecopulib/pull/696), [vinecopulib#729](https://github.com/vinecopulib/vinecopulib/pull/729)).
- Add conditioning-aware structure selection: `FitControlsVinecop.conditioning_set` puts the conditioning set at the tail of the fitted order, and `Vinecop.reorient` relabels an already-fitted vine onto a chosen tail without refitting (#246, [vinecopulib#697](https://github.com/vinecopulib/vinecopulib/pull/697)).
- `Vinecop.rosenblatt` and `Vinecop.inverse_rosenblatt` accept a keyword-only `conditioning_set` and a truncated model (#255, #306, [vinecopulib#715](https://github.com/vinecopulib/vinecopulib/pull/715), [vinecopulib#743](https://github.com/vinecopulib/vinecopulib/pull/743), [vinecopulib#752](https://github.com/vinecopulib/vinecopulib/pull/752)).
- The array-agnostic layer gets the same three, reproducing the compiled versions bit for bit: `VinecopBase.sample_conditional`, `reorient` and `select(conditioning_set=)`, plus the data-scale `Vinedist.sample_conditional` (#292, #306).
- Read an atom's probability rather than differencing four `cdf` values, where a pair copula can: `BicopBase.rect_prob` / `cond_interval_prob` are the two hooks, the mixed-discrete quotients call them, and `TorchTllBicop` overrides both onto its grid (#336, #339, [vinecopulib#771](https://github.com/vinecopulib/vinecopulib/pull/771), [vinecopulib#777](https://github.com/vinecopulib/vinecopulib/pull/777)).

#### Inference and diagnostics

- Add the analytic-inference surface to `Bicop` and `Vinecop`: `scores`, `gradient`,
  `hessian`, `scores_cov` and the per-observation `scores_full` / `hessian_full` /
  `pdf_full`, built on new derivative methods `pdf_deriv` / `pdf_deriv2` /
  `logpdf_deriv` / `logpdf_deriv2` / `hfunc{1,2}_deriv` / `hfunc{1,2}_deriv2` that are
  analytic for the families in `families.analytic_derivs` (#228, #231, #232, #246,
  [vinecopulib#679](https://github.com/vinecopulib/vinecopulib/pull/679), [vinecopulib#683](https://github.com/vinecopulib/vinecopulib/pull/683), [vinecopulib#687](https://github.com/vinecopulib/vinecopulib/pull/687),
  [vinecopulib#694](https://github.com/vinecopulib/vinecopulib/pull/694), [vinecopulib#699](https://github.com/vinecopulib/vinecopulib/pull/699)).
- Evaluate with per-observation parameters: `Bicop.pdf` / `cdf` / `hfunc*` / `hinv*` / `loglik` / `sample` take an `(n, p)` `parameters` array, and `Vinecop.pdf` / `pdf_full` / `loglik` and its score family an `(n, npars)` one (#227, #246, #254, [vinecopulib#675](https://github.com/vinecopulib/vinecopulib/pull/675), [vinecopulib#699](https://github.com/vinecopulib/vinecopulib/pull/699), [vinecopulib#719](https://github.com/vinecopulib/vinecopulib/pull/719)).
- Add tail dependence and Blomqvist's beta to `Bicop`: the `taildep` and `beta` properties plus `parameters_to_taildep` / `parameters_to_beta` (#230, [vinecopulib#682](https://github.com/vinecopulib/vinecopulib/pull/682)).
- Add `Vinecop.logpdf` and `VinecopBase.logpdf`, the per-observation log-density, plus a `"logpdf"` key on `Vinecop.pdf_full`. The accurate way to obtain it: a vine density is a product of up to `d(d-1)/2` pair densities, so `pdf` underflows to `0` well before the log-density stops being representable (#336, [vinecopulib#770](https://github.com/vinecopulib/vinecopulib/pull/770)).

#### PyTorch

- Add `pyvinecopulib.torch` behind a `pyvinecopulib[torch]` extra: `TorchTllBicop` and `TorchVinecop` are pure-PyTorch `nn.Module` evaluators mirroring `Bicop` / `Vinecop`, configured through the `FitControlsTorchBicop` / `FitControlsTorchVinecop` dataclasses (#216, #217, #225, #237, #244).
- Add the torch half of the vine-distribution stack: `TorchDistributionMargin` over any continuous `torch.distributions` family, and `TorchVinedist`, a vine distribution that is `nn.Module` throughout, so one `.to(device)` moves the whole joint distribution (#292).
- `TorchVinecop.from_data` selects an R-vine structure natively in torch when `structure=None`, instead of round-tripping through a compiled `Vinecop` (#241, #244).
- Add `TorchKde1d`, the torch marginal estimator, and the only torch margin that handles discrete and zero-inflated variables (#292, #312).
- The torch cascades gain a batched fast path that evaluates a whole tree level at once, resolved per device and 3-12x faster on CUDA than the per-edge cascade on
  `pdf` and `rosenblatt`; `FitControlsTorchVinecop.compile` additionally runs them through `torch.compile` (#219, #239, #264, #305, #307).

- `cache_integrals=True` is now the default and exact, with `cdf` and `hfunc*` read in closed form from cumulative-trapezoid prefix tables carrying an exact gradient, plus `TorchTllBicop.rect_prob` / `cond_interval_prob` for an atom's exact probability (#219, #305, #307, #336).
- Compute the `tll` window smoother in `O(n)` rather than `O(n**2)`: the window
  grows with the data, so at `n = 12000` it was 97% of a vine fit, and a
  fixed-structure CUDA fit at `d = 9` is 26x faster (#313, #315).

- Fit a whole tree level in one call: `FitControlsTorchVinecop.batched_fit` drives the new `fit_level` hook and `TorchTllBicop.from_data_batched`, making a vine fit 1.3-4.0x faster on CUDA, while the opt-in `compile_fit` fuses the `tll` bandwidth search (#316, #319).
- Invert the `tll` conditional cdf in closed form on the torch path, so `TorchTllBicop.hinv1` / `hinv2` match the compiled ones to machine precision and run 60-120x faster on CPU (#234, [vinecopulib#691](https://github.com/vinecopulib/vinecopulib/pull/691)).

#### scikit-learn

- Add `pyvinecopulib.sklearn` behind a `pyvinecopulib[sklearn]` extra: `VineDensity` and `VineRegressor`, scikit-learn-compatible estimators over mixed continuous and discrete input as a DataFrame or an ndarray (#211, #213, #263).

- Both sklearn estimators take `distribution=`, the `VinedistBase` subclass they fit: `Vinedist` by default, `TorchVinedist` to run the same pipeline under PyTorch. That one class names both halves of the model, so `controls=`, `structure=` and `margin_controls=` are all that go beside it (#218, #241, #339).
- Both sklearn estimators delegate the whole two-step fit to the distribution, handing it `margin_controls=` and the variable types and bounds `schema_` inferred, and publishing the result as `distribution_` alongside `schema_`, `structure_`, `margin_summary_`, `controls_`, `distribution_class_` and `random_state_` (#218, #292, #339).
- Add `n_jobs` to `VineDensity` and `VineRegressor`, governing fitting *and* every evaluation where the fit-time thread count used to pin both; results are bit-identical at any thread count (#297).
- `VineRegressor` accepts any continuous margin as its response, taking the `use_grid=True` quadrature on the probability scale so `n_nodes` fixes the number of probability levels (#292).
- `VineRegressor` publishes the two halves of its estimator: `conditional_weights(X)` gives the per-row copula weights and `y_nodes_` the response nodes they sit on, which is what a wrapper averaging several fitted vines works with and what `normalize_weights` changes -- a prediction is a ratio of them and is the same either way. `copula_marginal_density` is public for the same reason (#292, #345).

#### Structures, utilities and persistence

- Add the list-of-trees round-trip and the whole-array structure accessors: `Vinecop.get_trees()`, `RVineStructure.get_trees()` / `from_trees(d, trees)` /
  `from_struct_array(order, struct_array)` / `__eq__`, and `get_struct_array` / `get_min_array` / `get_needed_hfunc1` / `get_needed_hfunc2` (#246, #251, [vinecopulib#680](https://github.com/vinecopulib/vinecopulib/pull/680), [vinecopulib#698](https://github.com/vinecopulib/vinecopulib/pull/698), [vinecopulib#702](https://github.com/vinecopulib/vinecopulib/pull/702)).
- Supply your own edge weight for structure selection: set `FitControlsVinecop(tree_criterion="custom")` and pass `tree_criterion_function`, a callable `f(data, weights) -> float`; either without the other raises (#226, [vinecopulib#674](https://github.com/vinecopulib/vinecopulib/pull/674)).
- `utils.wdm` accepts Chatterjee's xi (`"chatterjee"` / `"cxi"` / `"xi"`) and a `seeds` argument for its predictor tie-breaking, and `FitControlsVinecop.tree_criterion` accepts it too, spelled `"cxi"` (#305, #312, [wdm#15](https://github.com/tnagler/wdm/pull/15), [wdm#26](https://github.com/tnagler/wdm/pull/26), [vinecopulib#754](https://github.com/vinecopulib/vinecopulib/pull/754)).
- Add binary CBOR model persistence: `to_file` / `from_file` on `Bicop`, `Vinecop` and `RVineStructure` select CBOR when the filename ends in `.cbor`, and JSON otherwise (#233, [vinecopulib#684](https://github.com/vinecopulib/vinecopulib/pull/684)).
- Persist the objects this release adds the way the copula classes already do:
  `Kde1d`, `Vinedist` and `SciPyMargin` gain `to_json` /
  `from_json`, and `Kde1d` / `Vinedist` also `to_file` / `from_file`, which write
  CBOR when the filename ends in `.cbor`. A margin type from outside the package
  joins in through `core.register_margin_json` (#320, #334).
    - a `VinedistBase` subclass inherits `from_json`: the base decodes, checks the payload's version and checks that its `kind` names the class being read, then rebuilds the halves from the declared `vinecop_class` — so reading back is a declaration, like fitting (#334)
- Expose `utils.find_latent_sample(u, b, niter=3)`, which recovers a continuous sample from interval-censored copula data (#305).
- Add `Kde1d.actual_grid_size`, the number of grid points a fit built (#312).
- Add `Bicop.family_name`, `Bicop.flip`, `Bicop.with_var_types()`, a settable `Vinecop.pair_copulas`, and `FitControlsVinecop.from_bicop_controls` with its `bicop_controls` property, which read and replace the inherited pair-copula settings as a group (#237, #251).

#### Documentation and examples

- Merge `03_vine_copulas_fit_sample.ipynb` into `02_vine_copulas.ipynb` and renumber, so the examples read bivariate copula, vine copula, vine distribution: `11_vine_distributions` becomes `03`, and `10_extending_pyvinecopulib` keeps its number. Fixes four cross-references that named the wrong notebook (#326).
- Document all of it: a `concepts.rst` primer on Sklar's theorem, pair-copula constructions, R-vines and the TLL family, a migration guide, and ten executed example notebooks (#211, #214, #216, #218, #240, #243, #246, #247, #252, #261, #262, #292, #326).

### Bug fixes in `pyvinecopulib`

- `Vinecop.loglik`, `aic`, `bic` and `mbicv` stay finite wherever the log-likelihood is representable, where the density was accumulated as a running product and underflowed to `0` below about `-745`; an observation carrying a `NaN` is left out of the total rather than making it `NaN`. `VinecopBase`, `TorchVinecop` and `BicopBase` follow, and `Vinedist.logpdf` reads the copula's log-density rather than logging its density (#335, #336, [vinecopulib#770](https://github.com/vinecopulib/vinecopulib/pull/770)).
- A multithreaded `Vinecop.fit` reproduces the serial one: concurrent edges read h-function columns their siblings were writing (#336, [vinecopulib#774](https://github.com/vinecopulib/vinecopulib/pull/774)).
- `Vinecop.pdf` and `logpdf` do not depend on how many rows they are handed, nor on `num_threads` (#336, [vinecopulib#779](https://github.com/vinecopulib/vinecopulib/pull/779)).
- Seeds a caller sets on `FitControlsVinecop` choose the vine structure, where a `random_weighted` fit through a sklearn estimator discarded them and reseeded from `random_state`. `FitControlsVinecop().seeds` is now empty rather than 20 materialized draws, so "the caller set seeds" is answerable; an unseeded fit still draws from entropy (#339, [vinecopulib#780](https://github.com/vinecopulib/vinecopulib/pull/780)).
- `VineRegressor`'s marginal copula quadrature runs in log space, where it summed densities that underflow at a grid node and then floored the total at `log(1e-10)`, scoring a row the model excludes as finite and wrong. Its other half moved to log space in #336 (#339).
- `Kde1d.fit`, `.select` and `.from_data` refuse inputs that leave no observation standing, where all-`NaN` or all-zero weights and all-`NaN` observations each terminated the process: a `NaN` observation, a `NaN` weight and a zero weight are drop markers, and the fit rescales by what survives them. An infinite or negative weight is refused too, the latter having silently fitted an unweighted density (#326, #336, [kde1d#40](https://github.com/vinecopulib/kde1d-cpp/pull/40)).
- Every discrete or mixed `tll` fit is corrected: the fit uses its latent sample rather than discarding it, and evaluation computes a real discrete density whose atom masses sum to one, where the midpoint density missed that sum by up
  to 10% and single cells by 40% (#306, [vinecopulib#739](https://github.com/vinecopulib/vinecopulib/pull/739)).
- `tree_algorithm="random_weighted"` no longer hangs: Wilson's walk could not leave a vertex whose incident weights were all zero, which happens for every edge at `n <= 10` and wherever a degenerate pair
  strands a node on a larger sample (#319, [vinecopulib#759](https://github.com/vinecopulib/vinecopulib/pull/759)).
- The per-entry structure accessors `RVineStructure.struct_array` / `min_array` / `needed_hfunc1` / `needed_hfunc2` raise on a tree above `trunc_lvl` instead of segfaulting the interpreter (#306, [vinecopulib#756](https://github.com/vinecopulib/vinecopulib/pull/756)).
- `Kde1d.fit` re-selects the bandwidth on a refit, where it fed the previous
  fit's bandwidth back into the selector (#292, [kde1d#28](https://github.com/vinecopulib/kde1d-cpp/pull/28)).
- Preserve `nobs`, the attained `loglik()` and the information criteria when pickling a fitted `Bicop` or `Vinecop`, which came back reporting the model as unfitted; pickles from earlier releases still load (#320).
- `Vinecop.aic()`, `bic()` and `mbicv()` called without data return the values
  attained at the fit, where `bic` / `mbicv` gave `-inf` and `aic` a value
  computed from no observations (#320).
- `Kde1d.from_grid` rejects a grid of fewer than four points -- the minimum its
  own constructors accept, so nothing constructible fails to unpickle, which used to raise `MemoryError` or hang the interpreter (#320).
- `Kde1d.icdf` returns an empty array for an empty input instead of terminating the process (#320).
- `Bicop.plot` evaluates a continuous copy of a discrete copula instead of retyping the caller's model in place (#320).
- `Kde1d.plot` draws a discrete margin on its own integer support, where rounding the fitted grid outwards added one level below the smallest observation and one above the largest; its continuous grid is monotone rather than overwritten at the first point (#330).
- `core`'s SciPy-free `norm_ppf` / `inv_erf` evaluate every polynomial coefficient, where `polevl` read `N` of them instead of cephes' `N + 1` and dropped each constant term: `norm_ppf(0.25)` was `-0.663271` against `-0.674490`, worst at the quartiles, which misplaced the axes of a normal-scale `Bicop.plot`. The test that would have caught it pinned the implementation's own output as a "known value" (#339).
- `utils.wdm` releases the GIL, so it no longer serializes against other Python threads (#320).
- `utils.wdm` raises on weights that are not finite and nonnegative with a positive sum, where it returned `NaN` (#305, [wdm#21](https://github.com/tnagler/wdm/pull/21)).
- `utils.to_pseudo_obs(..., ties_method="random")` breaks ties uniformly; the random offsets were one step out of phase, so a tie group of two always came out equal (#305, [wdm#17](https://github.com/tnagler/wdm/pull/17)).
- `Bicop.loglik` and `Bicop.fit` validate the column count instead of reading past the data (#306, [vinecopulib#729](https://github.com/vinecopulib/vinecopulib/pull/729)).
- `Bicop.tau` returns the limit for `joe` at `theta = 2`, a removable
  singularity where it was `NaN`, and keeps its digits nearby (#319,
  [vinecopulib#758](https://github.com/vinecopulib/vinecopulib/pull/758)).
- `Bicop.mbic` counts only the observations actually fitted, so on data containing `NaN` it no longer disagrees with `bic` (#251, [vinecopulib#710](https://github.com/vinecopulib/vinecopulib/pull/710)).
- `Vinecop.inverse_rosenblatt` is thread-safe (#251, [vinecopulib#712](https://github.com/vinecopulib/vinecopulib/pull/712)).
- A `Vinecop` built from a full structure with no pair copulas treats the omitted ones as independence instead of indexing an empty store (#306, [vinecopulib#729](https://github.com/vinecopulib/vinecopulib/pull/729), [vinecopulib#743](https://github.com/vinecopulib/vinecopulib/pull/743)).
- Declare the `BicopFamily` members and the `RVineStructure` base of `CVineStructure` / `DVineStructure` in the generated type stubs, so both documented patterns pass static type checking (#224, #268).

### Build / packaging

- Enforce the engineering conventions rather than documenting them: `ty` type-checks `src` and `tests` with ten off-by-default rules on, `ruff` runs 33 rule families so an unjustified `Any`, a swallowed exception or a stale `noqa` fails the build, `codespell` holds the American-English rule and a banned-word list, and `tests/test_prose.py` covers the multi-word phrases and identifiers a tokenizer cannot see (#210, #252, #262, #326, #339).
- Fail the documentation build on an unresolved cross-reference: Sphinx runs nitpicky with `-W`, `numpydoc` validation is a pre-commit check, and every code example is a notebook cell executed in CI rather than a paste (#214, #240, #243, #326).
- Resolve the Eigen include directory from the `Eigen3::Eigen` target, so a
  source build works against Eigen 5.x (#235).
- Build the x86-64 wheels for the x86-64-v3 baseline (AVX2 and FMA) instead of
  inheriting upstream's `-march=native` (#250).
- Refuse to import such a wheel on an x86-64 CPU without AVX2 and FMA, with an
  `ImportError` naming the missing feature rather than a `SIGILL`; a source
  build targets the plain baseline and is never refused (#320).
- Ship 10 wheels rather than 16: a cp311 wheel plus a cp312 ABI3 wheel for manylinux, musllinux, macOS x86-64, macOS arm64 and Windows, with macOS x86-64 returning on `macos-15-intel` (#220, #292).
- Require CMake 3.20 and a C++17 compiler for a source build, following upstream.
  Boost is looked for in config mode, `FindBoost` having been removed in CMake
  3.30, and then by a plain search for its headers -- which is all this package
  needs of it -- so an install with no `BoostConfig.cmake` works; point
  `Boost_INCLUDE_DIR` at the headers to skip both (#250, #336, #339, [vinecopulib#711](https://github.com/vinecopulib/vinecopulib/pull/711), [vinecopulib#774](https://github.com/vinecopulib/vinecopulib/pull/774)).
- Build the source distribution from an explicit allowlist, so a dirty worktree
  cannot leak build trees, caches or untracked notes into it (#320).
- Require, before the tag-triggered PyPI upload runs, that the tagged commit is an ancestor of `main`, that the documentation build passes, and that the version in `pyproject.toml`, `CHANGELOG.md`, `CITATION.cff` and `.zenodo.json` matches the tag (#245, #253, #267, #320).

- Move the build to `uv` and `scikit-build-core`, with `[build-system].requires` mirroring the development dependency groups so `--no-build-isolation` works out of the box; `[dev]` is now a `uv` dependency group rather than an installable extra (#205, #209, #210, #248, #249, #258, #321, #322, #323).

### Dependency changes

- Raise the NumPy floor to `numpy>=2.0`, up from `>=1.14` (#211).
- Add `array_api_compat>=1.7` as a runtime dependency, which the array-agnostic `core` layer resolves its array namespace through (#236).
- Add three extras: `[sklearn]` (`scikit-learn>=1.4`, `pandas>=2.0`), `[torch]` (`torch>=2.2`) and `[scipy]` (`scipy>=1.16`) (#211, #216, #292, #307).
- `[examples]` adds `xlrd>=2.0`, and `[doc]` takes version ranges instead of exact pins and adds `numpydoc` (#220, #259).
- Pin `lib/vinecopulib` to its 1.0.0 line (#229, #251, #305, #312, #319).
- Bump `lib/wdm` to `v0.3.0`, which is where Chatterjee's xi comes from (#305, #312).
- Bump `lib/kde1d` to `v1.2.2`, across thirteen pull requests (#220, #292, #312, #320, #336).

### Changes in `vinecopulib`

These changes originate from [`vinecopulib`](https://github.com/vinecopulib/vinecopulib), the C++ library which powers `pyvinecopulib`. The pin is a commit on its `main`, on the way to its own 1.0.0 release rather than a tagged version, so the entries below cite pull requests.

#### BREAKING API CHANGES

- Store the R-vine structure with the conditioned variable on the diagonal, so `get_matrix`, `order`, `get_struct_array` and edge orientation differ for the same model; densities and log-likelihoods do not ([#702](https://github.com/vinecopulib/vinecopulib/pull/702)).
- Require C++17, CMake 3.14 and Boost 1.75 -- narrowed to Graph, Math and Random -- and put `-march=native` behind `VINECOPULIB_NATIVE_ARCH`, so the default release build is redistributable ([#711](https://github.com/vinecopulib/vinecopulib/pull/711), [#714](https://github.com/vinecopulib/vinecopulib/pull/714)).
- Remove `Vinecop::select_all`, `Vinecop::select_families` and the `*_truncation_level` accessors, deprecated since 0.3.1 ([#718](https://github.com/vinecopulib/vinecopulib/pull/718)).
- Rename `Bicop::as_continuous()` to `with_var_types()`, which takes the variable types rather than assuming both continuous and so goes in either direction: the same fitted copula can be evaluated on a discrete or mixed edge without being refitted. The old spelling is the default argument ([#777](https://github.com/vinecopulib/vinecopulib/pull/777)).

#### BEHAVIOR CHANGES

- Every `tll` fit moves: the interpolation grid's margins are balanced across both sweep orders and iterated to a tolerance, `hfunc` and `hinv` no longer floor the interpolated density at `1e-4`, and the discrete latent sample no longer depends on which variable is passed first ([#751](https://github.com/vinecopulib/vinecopulib/pull/751)).
- Every discrete or mixed `tll` fit moves again: the fit uses its latent sample instead of the jittered ranks, and evaluation computes a real discrete density whose atom masses sum to one ([#739](https://github.com/vinecopulib/vinecopulib/pull/739)).
- `Vinecop::pdf()` is the exponential of a log-space sum rather than a running product of edge densities, so its values move at the `1e-15` level ([#770](https://github.com/vinecopulib/vinecopulib/pull/770)).
- A discrete `tll` pair reads each atom's probability off its interpolation grid instead of differencing four distribution values, moving the density by about `1e-9` and the log-likelihood by `9e-10`; the fourteen parametric families are bit-identical ([#771](https://github.com/vinecopulib/vinecopulib/pull/771)).
- Kendall's tau of `bb6`, `bb7`, `bb8` and `tawn` changes, the worst case having returned about `1e-11` where the true value was `0.33` ([#713](https://github.com/vinecopulib/vinecopulib/pull/713)).
- The density of a discrete/discrete pair moves by `O(atom width)` wherever an atom is narrower than `5e-5`, the collapsed argument now being the atom midpoint in both h-function evaluations ([#744](https://github.com/vinecopulib/vinecopulib/pull/744)).
- Maximum-likelihood estimates shift in the low digits: BOBYQA is replaced by Brent and BFGS, and `tawn` starts from a tau-based initial value ([#685](https://github.com/vinecopulib/vinecopulib/pull/685)).
- `Vinecop::fit` raises when the model has no pair copulas to fit instead of returning it unfitted, and a structure truncated at zero records the independence fit ([#752](https://github.com/vinecopulib/vinecopulib/pull/752)).

#### NEW FEATURES

- Add `Vinecop::logpdf()`, the per-observation log-density, and a `logpdf` field on `pdf_full()`; `loglik()` sums it ([#770](https://github.com/vinecopulib/vinecopulib/pull/770)).

- Add Chatterjee's xi as a `tree_criterion`, spelled `"cxi"` and symmetrized as `max(xi(x, y), xi(y, x))`; like `"hoeffd"` it picks up non-monotonic dependence, and unlike it, relationships that are not smooth ([#754](https://github.com/vinecopulib/vinecopulib/pull/754)).
- Add the analytic-inference surface: scores, gradient, Hessian and score covariance
  on `Bicop` and `Vinecop`, their per-observation `scores_full` / `hessian_full` /
  `pdf_full` companions, and the analytic density and h-function derivatives they are
  built on, for every parametric family plus `tll`'s ([#669](https://github.com/vinecopulib/vinecopulib/pull/669), [#683](https://github.com/vinecopulib/vinecopulib/pull/683), [#687](https://github.com/vinecopulib/vinecopulib/pull/687), [#694](https://github.com/vinecopulib/vinecopulib/pull/694), [#699](https://github.com/vinecopulib/vinecopulib/pull/699)).
- Add conditional simulation -- `simulate_conditional`, conditioning-aware structure selection, and `reorient` -- and extend all three to truncated models ([#696](https://github.com/vinecopulib/vinecopulib/pull/696), [#697](https://github.com/vinecopulib/vinecopulib/pull/697), [#752](https://github.com/vinecopulib/vinecopulib/pull/752)).
- Accept a conditioning set on `rosenblatt` and `inverse_rosenblatt`, evaluated through non-owning reoriented views rather than by copying pair copulas ([#715](https://github.com/vinecopulib/vinecopulib/pull/715)).
- Accept both the expanded `2d` and compact `d + k` discrete layouts consistently across the bivariate, vine, Rosenblatt and conditional-simulation surfaces ([#729](https://github.com/vinecopulib/vinecopulib/pull/729)).
- Evaluate a copula at one parameter set per observation, on `pdf`, `cdf`, the h-functions, `loglik`, `simulate` and the derivative methods ([#675](https://github.com/vinecopulib/vinecopulib/pull/675), [#699](https://github.com/vinecopulib/vinecopulib/pull/699), [#719](https://github.com/vinecopulib/vinecopulib/pull/719)).
- Add tail-dependence coefficients and Blomqvist's beta to `Bicop` ([#682](https://github.com/vinecopulib/vinecopulib/pull/682)).
- Allow a user-supplied `tree_criterion` function for structure selection, selected by `tree_criterion="custom"` ([#674](https://github.com/vinecopulib/vinecopulib/pull/674)).
- Add `RVineTrees`, an explicit edge-list view of a structure, with `Vinecop::get_trees` and an `RVineStructure` constructor from it ([#698](https://github.com/vinecopulib/vinecopulib/pull/698), [#680](https://github.com/vinecopulib/vinecopulib/pull/680)).
- Persist models as CBOR as well as JSON, selected by the file extension, and raise on an I/O failure instead of failing silently ([#684](https://github.com/vinecopulib/vinecopulib/pull/684)).
- Rewrite the `//!` comments the Python API documentation is generated from, including per-family parameter and tail-dependence text on `BicopFamily` and a numpydoc-compliant comment on every property ([#668](https://github.com/vinecopulib/vinecopulib/pull/668), [#670](https://github.com/vinecopulib/vinecopulib/pull/670), [#723](https://github.com/vinecopulib/vinecopulib/pull/723), [#733](https://github.com/vinecopulib/vinecopulib/pull/733)).

#### PERFORMANCE

- Speed up bivariate evaluation, `tll` fitting and evaluation, vine evaluation and structure selection, and the shared Eigen, thread and statistics primitives; evaluations may shift at the `1e-12` level and `tll` inverse h-functions by up to `1e-9` ([#681](https://github.com/vinecopulib/vinecopulib/pull/681), [#689](https://github.com/vinecopulib/vinecopulib/pull/689), [#690](https://github.com/vinecopulib/vinecopulib/pull/690), [#691](https://github.com/vinecopulib/vinecopulib/pull/691), [#692](https://github.com/vinecopulib/vinecopulib/pull/692), [#693](https://github.com/vinecopulib/vinecopulib/pull/693)).
- Normalize an interpolation grid's margins about four times faster, and peel `RVineTrees` faster with byte-identical results ([#751](https://github.com/vinecopulib/vinecopulib/pull/751), [#701](https://github.com/vinecopulib/vinecopulib/pull/701)).

#### BUG FIXES

- `Vinecop::loglik()` and the criteria built on it stay finite wherever the log-likelihood is representable, and leave out an observation carrying a `NaN` rather than returning `NaN` for the whole sample ([#770](https://github.com/vinecopulib/vinecopulib/pull/770)).
- A multithreaded `Vinecop::fit()` reproduces the serial one: concurrent edges read h-function columns their siblings were writing ([#774](https://github.com/vinecopulib/vinecopulib/pull/774)).
- `Vinecop::pdf()` and `logpdf()` do not depend on how many rows they are handed, nor on `num_threads`: past a thousand rows Eigen reaches its vectorized logarithm, which rounds differently from the scalar one ([#779](https://github.com/vinecopulib/vinecopulib/pull/779)).
- A `ThreadPool` stays usable after a job throws, where the stored exception was rethrown by every later `wait()` and the work queued in the meantime was cleared without running; the pool's own state is also read under its lock ([#764](https://github.com/vinecopulib/vinecopulib/pull/764), [#774](https://github.com/vinecopulib/vinecopulib/pull/774)).
- Widen `find_latent_sample`'s sweep counter, which wrapped to zero and never terminated for `niter > 65535`; the draw is unchanged ([#764](https://github.com/vinecopulib/vinecopulib/pull/764)).
- Validate user-supplied shapes before indexing them: a parameter matrix of the wrong
  shape on `parameters_to_tau` / `parameters_to_taildep` / `parameters_to_beta`, an
  empty `RVineStructure(mat)`, and a `Vinecop::set_var_types` vector shorter than the
  dimension each read past their own storage ([#761](https://github.com/vinecopulib/vinecopulib/pull/761)).
- Return the limit for Joe's Kendall's tau at `theta = 2`, a removable singularity where the value was `NaN` and nearby ones lost most of their significant digits ([#758](https://github.com/vinecopulib/vinecopulib/pull/758)).
- Floor the random-tree weights so Wilson's walk terminates: `tree_algorithm="random_weighted"` hung outright at `n <= 10`, and on larger samples wherever a degenerate pair left a node with only zero-weight edges ([#759](https://github.com/vinecopulib/vinecopulib/pull/759)).
- Bounds-check the per-entry structure accessors, which read past the trapezoid and segfaulted on any tree above a truncated structure's `trunc_lvl` ([#756](https://github.com/vinecopulib/vinecopulib/pull/756)).
- Read an empty pair-copula store as independence everywhere -- evaluation, `reorient`, `get_trees`, `truncate` and deserialization -- instead of indexing past its end ([#729](https://github.com/vinecopulib/vinecopulib/pull/729), [#743](https://github.com/vinecopulib/vinecopulib/pull/743)).
- Validate the column count in `Bicop::loglik` and `Bicop::fit`, which read past the data on a wrong one ([#729](https://github.com/vinecopulib/vinecopulib/pull/729)).
- Compute `Bicop::mbic` from the observations actually fitted, so on data containing `NaN` it no longer disagrees with `bic` ([#710](https://github.com/vinecopulib/vinecopulib/pull/710)).
- Make `inverse_rosenblatt` thread-safe ([#712](https://github.com/vinecopulib/vinecopulib/pull/712)).
- Evaluate a user-supplied `tree_criterion` function serially, since a Python callable is not required to be thread-safe ([#722](https://github.com/vinecopulib/vinecopulib/pull/722)).
- Fix the `tll` CDF and the two-dimensional integration behind it ([#667](https://github.com/vinecopulib/vinecopulib/pull/667)).
- Fix an out-of-bounds read in the mixed discrete/continuous h-inverse ([#675](https://github.com/vinecopulib/vinecopulib/pull/675)).
- Improve the starting parameters when fitting a pair copula to discrete data ([#677](https://github.com/vinecopulib/vinecopulib/pull/677)).
- Report `select_threshold` rather than `select_trunc_lvl` on the `"Select threshold"` line of `FitControlsVinecop`'s `str()` ([#735](https://github.com/vinecopulib/vinecopulib/pull/735)).

### Changes in `kde1d`

These changes originate from [`kde1d`](https://github.com/vinecopulib/kde1d-cpp), the 1-d kernel-density library behind `Kde1d`. The pin is its `v1.2.2` tag; the entries below cite the pull requests behind it.

#### BREAKING API CHANGES

- A discrete `Kde1d`'s `xmin` / `xmax` mean its integer support rather than grid endpoints: bounds and data must both be integers, the fitted grid runs from `xmin - 0.5` to `xmax + 0.5`, and the masses are normalized over the declared support ([#37](https://github.com/vinecopulib/kde1d-cpp/pull/37)).

#### BEHAVIOR CHANGES

- Numbers move for every fit: linear binning no longer drops the observation at the top of its range, so automatic bandwidths shift by `O(1/n)`, right-tail interpolation is centered at the grid endpoint, densities are truncated at finite bounds, and discrete CDF values are clamped to the unit interval ([#29](https://github.com/vinecopulib/kde1d-cpp/pull/29)).
- Fit a finite bound with a dedicated local-linear boundary estimator by default, given at least sixteen observations; pass `boundary_repair=False` to keep the transformed fit ([#36](https://github.com/vinecopulib/kde1d-cpp/pull/36)).
- Fit a one-sidedly bounded variable through a regularized fourth-root transform with a scale-adaptive offset instead of a log transform with a fixed `1e-5` one, so such a fit changes and is now equivariant under a change of units ([#34](https://github.com/vinecopulib/kde1d-cpp/pull/34), [#36](https://github.com/vinecopulib/kde1d-cpp/pull/36)).
- Apply the declared bounds to the continuous component of a zero-inflated variable, and accept a negative discrete support ([#37](https://github.com/vinecopulib/kde1d-cpp/pull/37)).

#### NEW FEATURES

- Restore a fitted estimator from an interpolation grid plus its state -- multiplier, requested and attained bandwidth, degree, grid size, boundary-repair flag, `edf` and `loglik` -- so a `Kde1d` survives a round trip complete rather than as a bare grid
  ([kde1d#39](https://github.com/vinecopulib/kde1d-cpp/pull/39)).

#### PERFORMANCE

- Speed up spline interpolation, integration and quantile inversion by caching coefficients and cumulative integrals and inverting a continuous quantile within its cell rather than by 35 global bisections, which also moves `icdf` in the last bits ([#32](https://github.com/vinecopulib/kde1d-cpp/pull/32), [#33](https://github.com/vinecopulib/kde1d-cpp/pull/33)).

#### BUG FIXES

- Refuse fit inputs that leave no observation standing, and weights that are infinite or negative: the rescaling divided by zero and the grid was built from an empty sample, which terminated the process ([#40](https://github.com/vinecopulib/kde1d-cpp/pull/40)).
- Avoid a zero-size `realloc` in `remove_nans`, which is deprecated in C17, undefined in C23 and an error under valgrind; reachable from a zero-inflated fit whose observations are all zero ([#41](https://github.com/vinecopulib/kde1d-cpp/pull/41)).
- Re-select the bandwidth when an estimator is refitted, where `fit` fed the previous fit's bandwidth back into the selector ([#28](https://github.com/vinecopulib/kde1d-cpp/pull/28)).
- Report defined `multiplier` / `bandwidth` / `degree` on a density built from a grid, where they were indeterminate ([#28](https://github.com/vinecopulib/kde1d-cpp/pull/28)).
- Rewrite the public `//!` comments this package lifts verbatim into its Python API documentation ([#26](https://github.com/vinecopulib/kde1d-cpp/pull/26), [#27](https://github.com/vinecopulib/kde1d-cpp/pull/27), [#38](https://github.com/vinecopulib/kde1d-cpp/pull/38)).

### Changes in `wdm`

These changes originate from [`wdm`](https://github.com/tnagler/wdm) `v0.3.0`, the weighted-dependence-measure library behind `utils.wdm`, up from `v0.2.6`.

#### NEW FEATURES

- Add Chatterjee's xi, the one asymmetric measure in the list, since it measures how far `y` is a function of `x` ([#15](https://github.com/tnagler/wdm/pull/15), [#18](https://github.com/tnagler/wdm/pull/18), [#19](https://github.com/tnagler/wdm/pull/19)).
- Break tied predictors in Chatterjee's xi from a fixed default seed rather than `std::random_device`, so xi is a function of its arguments, and take `seeds` to vary that ordering ([#26](https://github.com/tnagler/wdm/pull/26)).

#### BUG FIXES

- Break ties at random in `rank(x, weights, "random")`, whose offsets were one step out of phase, so a tie group of two always came out equal and larger groups lost a rank ([#17](https://github.com/tnagler/wdm/pull/17)).
- Raise on weights that are not finite and nonnegative, or whose sum is not finite and positive, where the measure returned `NaN` ([#21](https://github.com/tnagler/wdm/pull/21)).
- Apply Hoeffding's five-observation minimum to the `"hoeffd"` and `"d"` aliases too, which slipped through the `"hoeffding"`-only check ([#21](https://github.com/tnagler/wdm/pull/21)).
- Raise on a constant response and drop zero-weight observations in Chatterjee's xi ([#21](https://github.com/tnagler/wdm/pull/21)).

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
- Bind the new `allow_rotation` option on `FitControlsBicop` and
  `FitControlsVinecop` (#168,
  [vinecopulib#628](https://github.com/vinecopulib/vinecopulib/pull/628))

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
