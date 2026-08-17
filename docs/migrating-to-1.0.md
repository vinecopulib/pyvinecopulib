# Migrating to 1.0.0

This page lists everything in 1.0.0 that can make code written against
**0.7.6** -- the previous release -- behave differently, in the order you are
likely to hit it. There was no 0.8.0 release, so every change below lands in
one step for anyone upgrading from 0.7.x.

Read the first two sections before you upgrade: the wheels now require a
newer CPU baseline, and nearly every fitted number moves.

If you are coming from the C++ library, note that
[vinecopulib's own migration guide](https://github.com/vinecopulib/vinecopulib/blob/main/docs/migrating-to-1.0.md)
is mostly **not** about Python. Its central warning — a constructor argument
inserted in the middle, which compiles silently against old positional calls —
does not reach here, because every alternative constructor in this package is a
named factory (`Bicop.from_family`, `Vinecop.from_data`, …) rather than an
overload. The next two sections are the parts that do carry over.

## x86-64 wheels now require AVX2 and FMA

The distributed x86-64 wheels are compiled for the x86-64-v3 baseline, so they
do not run on a CPU or VM without AVX2 and FMA (roughly, Intel before Haswell
or AMD before Excavator). The wheel tags cannot express this, so `pip` will
install one anyway; `import pyvinecopulib` then raises an `ImportError` saying
so rather than dying on an illegal instruction.

If you are on such a machine, build from source, which targets the plain
x86-64 baseline:

```bash
pip install --no-binary pyvinecopulib pyvinecopulib
```

## Nearly every fitted number moves

The vendored C++ libraries advanced substantially, so a refit will not
reproduce a 0.7.6 model exactly. Nothing raises; the numbers differ.

- Every `tll` fit changes -- and `tll` is what family selection usually picks
  for a nonparametric pair. The interpolation grid's margins are now balanced
  across both sweep orders and iterated to convergence, and `hfunc` / `hinv`
  no longer floor the interpolated density.
- Every `Kde1d` fit changes: the plug-in bandwidth shifts, and a bounded fit
  shifts again under a new boundary estimator.
- A discrete `Kde1d`'s `xmin` / `xmax` are now read as the **integer
  support**, so both the bounds and the data must be integers.
- `FitControlsVinecop.threshold` now also replaces a surviving sub-threshold
  edge with an independence pair, rather than only deprioritizing it when the
  spanning tree is built.

Store a fitted model rather than a recipe if you need to reproduce results
across this boundary.

## Model selection now defaults to AIC

`FitControlsBicop` and `FitControlsVinecop` default `selection_criterion` to
`"aic"`. They previously defaulted to `"bic"` in Python while C++ and R already
used `"aic"`, so the same knob selected different models depending on which
language you called from.

This changes which families and structures get selected. Nothing raises; the
fitted model is simply different, usually slightly less parsimonious. To keep
the previous behavior, say so explicitly:

```python
controls = pv.FitControlsVinecop(selection_criterion="bic")
```

Saved models are unaffected — the criterion is a fit-time setting, not part of
a serialized model.

## Some arguments are keyword-only

`parameters` on `Bicop.sample`, and `conditioning_set` on
`Vinecop.rosenblatt` / `inverse_rosenblatt` / `sample_conditional`, are
keyword-only. This is what keeps the long-standing positional forms meaning
what they always meant: `rosenblatt(u, 4)` is still `num_threads=4`, and
`sample(1000, True)` is still `n=1000, qrng=True`.

You only need to change code that was already passing these by keyword, which
is to say: none.

## Sampling methods are now called `sample`

`sample` is what scikit-learn, `torch.distributions` and SciPy's modern
distribution classes all call this operation, so `pyvinecopulib` calls it that
too. The five names that shipped in 0.7.6 still work and emit a
`DeprecationWarning` naming their replacement:

| 0.7.6 | 1.0.0 |
| --- | --- |
| `Bicop.simulate` | `Bicop.sample` |
| `Vinecop.simulate` | `Vinecop.sample` |
| `RVineStructure.simulate` | `RVineStructure.sample` |
| `Kde1d.simulate` | `Kde1d.sample` |
| `pyvinecopulib.utils.simulate_uniform` | `pyvinecopulib.utils.sample_uniform` |

Arguments and return values are unchanged; only the spelling is. Like the
top-level aliases below, these are scheduled for removal in 2.0.

`Vinecop.simulate_conditional` is the exception: it is now
`Vinecop.sample_conditional` with **no** alias, because it never appeared in a
release, so only code written against a development build can be affected.

The RNG hook that `MarginBase`, `BicopBase` and `VinecopBase` draw their
uniforms through is renamed with them: `_simulate_uniform` becomes
`_sample_uniform`. A renamed hook is the one rename that cannot fail visibly on
its own — the base class simply stops calling the old name, so an override
under the old name would be ignored and the inherited default would raise as
though nothing had been overridden. A subclass that defines `_simulate_uniform`
and not `_sample_uniform` therefore raises `TypeError` at class-definition
time, naming both spellings.

Two `sample` conventions now coexist, deliberately. The `core` classes keep the
quasi-random arguments they always had — `sample(n, qrng=False, seeds=[])`,
where `seeds` is a list of `int`. The `pyvinecopulib.sklearn` estimators keep
`sample(n_samples, random_state)`, because that is the signature scikit-learn's
`check_estimator` requires. The method name is the same; the argument names
follow whichever ecosystem you are calling from.

## `Kde1d` moved, and two of its members were renamed

`Kde1d` now lives in `pyvinecopulib.core`, beside the margin contract it
satisfies, and two members changed name. The top-level `pyvinecopulib.Kde1d`
still resolves and warns; `utils.Kde1d` never shipped in a release and is gone.

| 0.7.6 | 1.0.0 |
| --- | --- |
| `kde.quantile(p)` | `kde.icdf(p)` — no alias |
| `kde.loglik` (property) | `kde.loglik()` — a method taking optional data |
| `pyvinecopulib.Kde1d` | `pyvinecopulib.core.Kde1d` |

`icdf` is the name modern SciPy and `torch.distributions` use for the inverse
distribution function, and `loglik()` now matches `Bicop.loglik` /
`Vinecop.loglik`.

## Python 3.10 is no longer supported

`requires-python` is `>=3.11`. 3.10 reaches end of life in October 2026.

## Top-level aliases still work, and still warn

`pyvinecopulib.gaussian`, `pyvinecopulib.Kde1d`, `pyvinecopulib.wdm` and the
other names that moved into subpackages continue to resolve at the top
level, emitting a `DeprecationWarning` that names the canonical path. They are
**not** removed in 1.0.0 — that release breaks enough on its own — but they are
scheduled for removal in 2.0. Move to the subpackage imports now:

```python
from pyvinecopulib.families import gaussian
from pyvinecopulib.core import Kde1d
from pyvinecopulib.utils import wdm
```

The eight core classes (`Bicop`, `Vinecop`, the three structures, the two
fit-controls classes, `BicopFamily`) and `to_pseudo_obs` stay at the top level
indefinitely.

## sklearn and Torch changes

The optional subpackages ship in the same distribution and have important 1.0
changes. sklearn estimators now take a single `backend=` object instead of
loose controls/structure/seed arguments; `seed` became `random_state`; and
`VineRegressor` keeps a sample axis for one-row predictions. Torch fitting now
uses `FitControlsTorchBicop` / `FitControlsTorchVinecop`, and `TorchBicop.sample`
uses the core-style `(n, qrng=False, seeds=[])` signature. See the complete
breaking-change inventory in `CHANGELOG.md` before upgrading either surface.

## What 1.0 does and does not promise

`core`, `families` and `utils` are stable, and follow semantic versioning from
here.

`margins`, `sklearn` and `torch` are **provisional in 1.x**: they ship for the
first time in this release, and their surfaces may change in a minor version.
The margin contract and the torch/C++ parity guarantees are the parts already
treated as load-bearing. Pin an exact version if you build on the rest.

## What did not change

- Serialized `Bicop` / `Vinecop` / `RVineStructure` models. JSON and CBOR
  files written by 0.7.x load unchanged, and so do pickles (through the
  deprecated aliases). Note that the classes 1.0 adds -- `Vinedist`, the
  `margins` types and the `torch` modules -- have no JSON/CBOR surface; use
  `pickle`, or `state_dict` for the torch modules.
- Every evaluation signature other than the keyword-only arguments above.
