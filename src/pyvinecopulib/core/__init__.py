"""Core copula classes.

The core API exposes the bivariate copula model (:class:`Bicop`), the
vine copula model (:class:`Vinecop`), the three flavors of regular
vine structure used to describe how pair-copulas compose
(:class:`RVineStructure`, :class:`CVineStructure`,
:class:`DVineStructure`), and the fit-control objects that carry the
hyperparameters for both fits (:class:`FitControlsBicop` and
:class:`FitControlsVinecop`).

On top of that copula layer it exposes the pieces that turn a copula
into a distribution on the original scale: the univariate-margin
contract (:class:`MarginLike`, :class:`MarginBase`), the
boundary-corrected 1d kernel density estimator that serves as the
default margin (:class:`Kde1d`), and the joint object composing the two
halves of Sklar's theorem (:class:`Vinedist`).

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
- *A distribution, not just a copula* — :class:`Vinedist`. It pairs any
  :class:`VinecopLike` with one margin per variable and evaluates
  ``pdf`` / ``cdf`` / ``sample`` on the original scale;
  :meth:`Vinedist.from_data` fits the margins and an ``x``-independent compiled
  copula in one call. An externally conditional copula is fitted through the
  :meth:`VinecopBase.fit` extension seam and then composed with ``Vinedist``.
- *Margins* — :class:`Kde1d` is the default and needs no configuration
  beyond its variable type. Any object with ``pdf`` / ``cdf`` / ``icdf``
  works (:class:`MarginLike`); subclass :class:`MarginBase` to write
  one, or reach for the parametric families and family selection in
  :mod:`pyvinecopulib.margins`.
- *A custom pair copula on a discrete edge* — :class:`DiscretePair`.
  Wrap a continuous pair copula (anything with a ``cdf``) in the variable
  types :meth:`VinecopBase.pair_var_types` derives for its slot, and it
  supplies the mixed-discrete density and h-functions the cascades ask for.
- *An edge below a dependence threshold* — :class:`IndependencePair`.
  What :meth:`VinecopBase.select` places on a surviving edge whose criterion
  falls below ``threshold``, matching what the compiled selector leaves there:
  the edge is not fitted, and holds the independence copula.
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
from ._discrete import DiscretePair
from ._independence import IndependencePair
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
  "DiscretePair",
  "DVineStructure",
  "FitControlsBicop",
  "FitControlsVinecop",
  "IndependencePair",
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
