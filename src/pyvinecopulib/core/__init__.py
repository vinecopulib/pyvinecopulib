"""Core copula classes.

The core API exposes the bivariate copula model (:class:`Bicop`), the
vine copula model (:class:`Vinecop`), the three flavours of regular
vine structure used to describe how pair-copulas compose
(:class:`RVineStructure`, :class:`CVineStructure`,
:class:`DVineStructure`), and the fit-control objects that carry the
hyperparameters for both fits (:class:`FitControlsBicop` and
:class:`FitControlsVinecop`).

For higher-level scikit-learn-compatible wrappers around these
primitives (with ``fit`` / ``predict``-style interfaces, DataFrame
handling, and ensemble methods), see
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
  :meth:`RVineStructure.simulate`) skips structure selection
  entirely.
- *C-vine / D-vine special cases* — :class:`CVineStructure` and
  :class:`DVineStructure` are the path-shaped and star-shaped
  specialisations of :class:`RVineStructure` and can be passed
  anywhere an R-vine is accepted.

The :doc:`concepts page </concepts>` is a primer with formulas, the
:ref:`pair-copula construction <concepts-pcc>` factorisation,
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
  RVineStructure,
  Vinecop,
)
from .bicop_base import BicopBase
from ._context import (
  ConditioningContext,
  NonSimplifiedContext,
  SimplifiedContext,
)
from .protocols import BicopLike, VinecopLike
from ._vinecop_base import VinecopBase

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
  "NonSimplifiedContext",
  "RVineStructure",
  "SimplifiedContext",
  "Vinecop",
  "VinecopBase",
  "VinecopLike",
]
