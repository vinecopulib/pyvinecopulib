# Migrating to 1.0.0

Most code needs no changes. This page lists everything in 1.0.0 that can
make working code behave differently, in the order you are likely to hit it.

If you are coming from the C++ library, note that
[vinecopulib's own migration guide](https://github.com/vinecopulib/vinecopulib/blob/main/docs/migrating-to-1.0.md)
is mostly **not** about Python. Its central warning — a constructor argument
inserted in the middle, which compiles silently against old positional calls —
does not reach here, because every alternative constructor in this package is a
named factory (`Bicop.from_family`, `Vinecop.from_data`, …) rather than an
overload. The two items below are the parts that do carry over.

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

`parameters` on `Bicop.simulate`, and `conditioning_set` on
`Vinecop.rosenblatt` / `inverse_rosenblatt` / `simulate_conditional`, are
keyword-only. This is what keeps the long-standing positional forms meaning
what they always meant: `rosenblatt(u, 4)` is still `num_threads=4`, and
`simulate(1000, True)` is still `n=1000, qrng=True`.

You only need to change code that was already passing these by keyword, which
is to say: none.

## Top-level aliases still work, and still warn

`pyvinecopulib.gaussian`, `pyvinecopulib.Kde1d`, `pyvinecopulib.wdm` and the
other names moved into subpackages in 0.8.0 continue to resolve at the top
level, emitting a `DeprecationWarning` that names the canonical path. They are
**not** removed in 1.0.0 — that release breaks enough on its own — but they are
scheduled for removal in 2.0. Move to the subpackage imports now:

```python
from pyvinecopulib.families import gaussian
from pyvinecopulib.utils import Kde1d, wdm
```

The eight core classes (`Bicop`, `Vinecop`, the three structures, the two
fit-controls classes, `BicopFamily`) and `to_pseudo_obs` stay at the top level
indefinitely.

## What did not change

- Serialized models. JSON and CBOR files written by 0.8.x load unchanged.
- Every evaluation signature other than the keyword-only arguments above.
- `pyvinecopulib.sklearn` and `pyvinecopulib.torch`, which are on their own
  release cadence; see the changelog for those.
