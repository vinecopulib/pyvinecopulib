"""Core copula classes.

The core API exposes the bivariate copula model (:class:`Bicop`), the
vine copula model (:class:`Vinecop`), the three flavors of regular
vine structure used to describe how pair-copulas compose
(:class:`RVineStructure`, :class:`CVineStructure`,
:class:`DVineStructure`), and the fit-control objects that carry the
hyperparameters for both fits (:class:`FitControlsBicop` and
:class:`FitControlsVinecop`).

For higher-level scikit-learn-compatible wrappers around these
primitives (with ``fit`` / ``predict``-style interfaces, DataFrame
handling, and batched evaluation), see
:mod:`pyvinecopulib.sklearn`. The
``examples/02_vine_copulas.ipynb`` and
``examples/03_vine_copulas_fit_sample.ipynb`` notebooks walk through
end-to-end use of these classes.

Notes
-----
**When to use what.**

- *Bivariate dependence* — :class:`Bicop`. Pick a family explicitly
  via :meth:`Bicop.from_family` or fit one from data with
  :meth:`Bicop.from_data`; the family set considered during selection
  is controlled by :class:`FitControlsBicop` and the family tags from
  :mod:`pyvinecopulib.families`.
- *Multivariate dependence* — :class:`Vinecop`. The factory
  :meth:`Vinecop.from_data` runs the Dissmann algorithm by default
  (``tree_algorithm="mst_prim"`` on :class:`FitControlsVinecop`);
  passing ``"random_weighted"`` switches to Wilson-weighted random
  trees instead, and supplying a hand-crafted
  :class:`RVineStructure` (or one drawn from
  :meth:`RVineStructure.sample`) skips structure selection
  entirely.
- *C-vine / D-vine special cases* — :class:`CVineStructure` and
  :class:`DVineStructure` are the path-shaped and star-shaped
  specializations of :class:`RVineStructure` and can be passed
  anywhere an R-vine is accepted.

The :doc:`concepts page </concepts>` is a primer with formulas, the
:ref:`pair-copula construction <concepts-pcc>` factorization,
:ref:`available families <concepts-families>`, the
:ref:`two-step estimator <concepts-estimation>`, and
:ref:`structure selection <concepts-structure-selection>`.
"""

from ..pyvinecopulib_ext import (
  Bicop,
  BicopFamily,
  CVineStructure,
  DVineStructure,
  FitControlsBicop,
  FitControlsVinecop,
  Kde1d,
  RVineStructure,
  Vinecop,
)
from .._deprecations import _method_alias
from .bicop_base import BicopBase
from .context import (
  ConditioningContext,
  NonSimplifiedContext,
  SimplifiedContext,
)
from .margin_base import MarginBase
from .protocols import BicopLike, MarginLike, VinecopLike
from .vinedist import Vinedist
from .vinecop_base import VinecopBase

# Deprecated aliases for the pre-1.0 `simulate` spelling. The canonical name is
# `sample`, matching scikit-learn, `torch.distributions` and SciPy's modern
# distributions; these three shipped in 0.7.6, so they warn rather than vanish.
# `Vinecop.sample_conditional` has never been released and carries no alias.
Bicop.simulate = _method_alias(Bicop.sample, "simulate", "core.Bicop")
Kde1d.simulate = _method_alias(Kde1d.sample, "simulate", "core.Kde1d")
Vinecop.simulate = _method_alias(Vinecop.sample, "simulate", "core.Vinecop")
RVineStructure.simulate = staticmethod(
  _method_alias(RVineStructure.sample, "simulate", "core.RVineStructure")
)

__all__ = [
  "Bicop",
  "BicopBase",
  "BicopFamily",
  "BicopLike",
  "ConditioningContext",
  "CVineStructure",
  "DVineStructure",
  "FitControlsBicop",
  "FitControlsVinecop",
  "Kde1d",
  "MarginBase",
  "MarginLike",
  "NonSimplifiedContext",
  "RVineStructure",
  "SimplifiedContext",
  "Vinecop",
  "VinecopBase",
  "VinecopLike",
  "Vinedist",
]
