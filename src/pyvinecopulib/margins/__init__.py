"""Univariate margins for vine distributions.

A margin is both a specification and a fitted object: hyperparameters go to the
constructor, and ``fit`` estimates the rest in place and returns ``self``, as
:class:`pyvinecopulib.core.Kde1d` and the scikit-learn estimators do. A margin
whose parameters are given at construction is already fitted.

- **Nonparametric** — ``Kde1d`` (the default; boundary-corrected,
  and handles continuous, discrete and zero-inflated data). It is also the
  only margin that takes observation weights.
- **Parametric** — :class:`ParametricMargin` wraps one SciPy family, and
  :class:`MarginSelector` chooses among a curated candidate set by AIC / BIC /
  AICc, keeps the winner on ``selected_`` and forwards evaluation to it. Both
  need the ``scipy`` extra. :class:`OpenTURNSMargin` and
  :class:`OpenTURNSSelector` are their OpenTURNS counterparts, over that
  library's families and its own ``FittingTest`` criteria; they need the
  ``openturns`` extra.
- **Interoperability** — :func:`as_margin` presents a distribution object from
  another ecosystem as a margin, and :func:`register_margin_adapter` teaches it
  about one it does not know. SciPy's modern and legacy distributions,
  ``torch.distributions`` and OpenTURNS' distributions are recognized out of
  the box.

- **Resolution** — :func:`resolve_margins` expands a ``margins=`` argument
  (an alias, one margin, a sequence, a mapping, or a callable) into one
  specification per variable.

Custom margins subclass :class:`pyvinecopulib.core.MarginBase`, which needs only
``pdf`` and ``cdf``.

Notes
-----
The :doc:`concepts page </concepts>` introduces the marginal side of Sklar's
theorem, and :class:`pyvinecopulib.core.MarginLike` states the contract these
all satisfy.
"""

from ._adapters import as_margin, register_margin_adapter

# Imported for its side effect as much as its names: it registers the OpenTURNS
# adapter with `as_margin`, and it imports OpenTURNS itself only when used.
from ._openturns import OpenTURNSMargin, OpenTURNSSelector
from ._resolve import resolve_margins
from .parametric import ParametricMargin
from .selection import MarginSelector

__all__ = [
  "MarginSelector",
  "OpenTURNSMargin",
  "OpenTURNSSelector",
  "ParametricMargin",
  "as_margin",
  "resolve_margins",
  "register_margin_adapter",
]
