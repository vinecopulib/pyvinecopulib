"""The parts a subclass writes against, rather than the parts a caller uses.

Everything here is public and stable, and none of it is needed to *use*
pyvinecopulib. It is what an extension needs: the three steps of the input
pipeline, the validators that make a refusal read the same as the library's
own, the sentinel a batched-evaluation override raises, the callback aliases
and the model codec.

Kept out of ``pyvinecopulib.core``'s own namespace because the two have
different audiences. `core` names ``Bicop``, ``Vinedist``,
``FitControlsVinecop`` and the four canonical bases -- what a user reaches for
-- and adding a dozen names like ``usable_observations`` beside them makes the
surface harder to read for the many to serve the few, puts them in every
``from pyvinecopulib.core import *``, and lists them in the rendered API next
to ``to_pseudo_obs``. Reached as ``pyvinecopulib.core.extend``, the way
``pyvinecopulib.sklearn.backends`` is.

The four canonical bases and their protocols stay in `core`: ``README.md``
tells users to subclass them "with the same confidence as on ``Vinecop``", so
they are a documented part of the surface rather than machinery behind it.
``ArrayT`` stays there too, being the type variable those signatures are
written in.
"""

from ._covariates import covariate_row, prepare_covariates
from ._json import MODEL_JSON_VERSION, loads as model_from_json
from ._placement import place, reference_array, to_numpy
from ._trim import trim
from ._validation import (
  reject_covariates,
  usable_observations,
  validate_weights,
)
from ._vinecop_discrete import collapse_data, continuous_view
from ._vinecop_fit_engines import FitEdge, FitLevel
from .vinecop_base import NotBatchable

__all__ = [
  "FitEdge",
  "FitLevel",
  "MODEL_JSON_VERSION",
  "NotBatchable",
  "collapse_data",
  "continuous_view",
  "covariate_row",
  "model_from_json",
  "place",
  "prepare_covariates",
  "reference_array",
  "reject_covariates",
  "to_numpy",
  "trim",
  "usable_observations",
  "validate_weights",
]
