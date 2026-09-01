# Migrating to 1.0.0

Most code needs no changes. This page lists everything in 1.0.0 that can
make working code behave differently, in the order you are likely to hit it.

If you are coming from the C++ library, note that
[vinecopulib's own migration guide](https://github.com/vinecopulib/vinecopulib/blob/main/docs/migrating-to-1.0.md)
is mostly **not** about Python. Its central warning — a constructor argument
inserted in the middle, which compiles silently against old positional calls —
does not reach here, because every alternative constructor in this package is a
named factory (`Bicop.from_family`, `Vinecop.from_data`, …) rather than an
overload. The next two sections are the parts that do carry over.

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
release. Only code written against a 0.8.0 development build can be affected.

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

## Top-level aliases still work, and still warn

`pyvinecopulib.gaussian`, `pyvinecopulib.Kde1d`, `pyvinecopulib.wdm` and the
other names moved into subpackages in 0.8.0 continue to resolve at the top
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

## What did not change

- Serialized models. JSON and CBOR files written by 0.8.x load unchanged.
- Every evaluation signature other than the keyword-only arguments above.
