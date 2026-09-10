"""The parts a subclass writes against, rather than the parts a caller uses.

Everything here is public and stable, and none of it is needed to *use*
pyvinecopulib. It is the short list an extension cannot be written without:
the placement step and the two names that diagnose and override it, the
covariate composite every entry point taking an ``x`` routes through, the
single-row reading of one, the return trip to NumPy, and the sentinel a
batched-evaluation override raises to decline.

The list is short on purpose. An export is a promise kept from 1.0.0 onward,
so a helper earns one by being needed to implement a documented hook
*correctly* -- not by being called internally, and not by appearing in the
docstring of a private one. The validators, the vine's layout step, the pair
unwrapper, the fit-callback aliases and the model codec were all reachable
that second way and are private again: a subclass that wants a refusal to read
like the library's own writes its own message, and one that fits a level
writes a plain ``def`` rather than importing an alias that cannot pin its
arity anyway.

Kept out of ``pyvinecopulib.core``'s own namespace because the two have
different audiences. `core` names ``Bicop``, ``Vinedist``,
``FitControlsVinecop`` and the four canonical bases -- what a user reaches for
-- and putting machinery beside them makes the surface harder to read for the
many to serve the few, puts it in every ``from pyvinecopulib.core import *``,
and lists it in the rendered API next to ``to_pseudo_obs``. Reached as
``pyvinecopulib.core.extend``, the way ``pyvinecopulib.sklearn.backends`` is.

The four canonical bases and their protocols stay in `core`: ``README.md``
tells users to subclass them "with the same confidence as on ``Vinecop``", so
they are a documented part of the surface rather than machinery behind it.
``ArrayT`` stays there too, being the type variable those signatures are
written in.
"""

from ._covariates import (
  covariate_row,
  prepare_covariates,
)
from ._placement import place, reference_array, to_numpy
from .vinecop_base import NotBatchable

__all__ = [
  "NotBatchable",
  "covariate_row",
  "place",
  "prepare_covariates",
  "reference_array",
  "to_numpy",
]
