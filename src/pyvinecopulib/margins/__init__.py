"""Univariate margins for vine distributions.

A margin is both a specification and a fitted object: hyperparameters go to the
constructor, and ``fit`` estimates the rest in place and returns ``self``, as
:class:`pyvinecopulib.utils.Kde1d` and the scikit-learn estimators do. A margin
whose parameters are given at construction is already fitted.

- **Nonparametric** — :class:`Kde1dMargin` (the default; boundary-corrected,
  and handles continuous, discrete and zero-inflated data) and
  :class:`EmpiricalMargin` (the rank transform behind
  :func:`pyvinecopulib.utils.to_pseudo_obs`, as a distribution).
- **Parametric** — :class:`ParametricMargin` wraps one SciPy family, and
  :class:`MarginSelector` chooses among a curated candidate set by AIC / BIC /
  AICc, keeps the winner on ``selected_`` and forwards evaluation to it. Both
  need the ``scipy`` extra.
- **Interoperability** — :func:`as_margin` presents a distribution object from
  another ecosystem as a margin, and :func:`register_margin_adapter` teaches it
  about one it does not know. SciPy's modern and legacy distributions and
  ``torch.distributions`` are recognized out of the box.

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
from ._resolve import resolve_margins
from .empirical import EmpiricalMargin
from .kde import Kde1dMargin
from .parametric import ParametricMargin
from .selection import MarginSelector

__all__ = [
  "EmpiricalMargin",
  "Kde1dMargin",
  "MarginSelector",
  "ParametricMargin",
  "as_margin",
  "resolve_margins",
  "register_margin_adapter",
]
