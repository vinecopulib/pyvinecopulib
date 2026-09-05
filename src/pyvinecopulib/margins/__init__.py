"""Univariate margins for vine distributions.

A margin is both a specification and a fitted object: hyperparameters go to the
constructor, and ``fit`` estimates the rest in place and returns ``self``, as
:class:`pyvinecopulib.core.Kde1d` and the scikit-learn estimators do. A margin
whose parameters are given at construction is already fitted.

- **Nonparametric** — ``Kde1d`` (the default; boundary-corrected,
  and handles continuous, discrete and zero-inflated data). It is also the
  only margin that takes observation weights.
- **Parametric** — :class:`SciPyMargin` wraps one SciPy family. Named,
  ``fit`` estimates its parameters; unnamed, ``select`` chooses the family too,
  from a curated candidate set scored by AIC / BIC / AICc. It needs the
  ``scipy`` extra. :class:`OpenTURNSMargin` is the OpenTURNS counterpart, over
  that library's families and its own ``FittingTest`` criteria; it needs the
  ``openturns`` extra.
- **Interoperability** — :func:`as_margin` presents a distribution object from
  another ecosystem as a margin, and :func:`register_margin_adapter` teaches it
  about one it does not know. SciPy's modern and legacy distributions,
  OpenTURNS distributions, and continuous ``torch.distributions`` families
  with an implemented ``cdf`` are recognized out of the box. Torch ``Normal``,
  ``Gamma`` and ``LogNormal`` are supported; discrete Torch distributions are
  rejected because they do not expose the left-limit cdf a vine needs.

- **Resolution** — :func:`resolve_margins` expands a ``margins=`` argument
  (an alias, one margin, a sequence, a mapping, or a callable) into one
  specification per variable, and :func:`resolve_margin_controls` does the same
  for ``margin_controls=``. The two are complementary: ``margins`` says which
  class each variable gets, and :class:`FitControlsMargin` says how to fit or
  select it — so one call can bound the two variables whose bounds are known
  and leave the rest alone.

Custom margins subclass :class:`pyvinecopulib.core.MarginBase`, which needs only
``pdf`` and ``cdf``; a custom *distribution* to put them in subclasses
:class:`pyvinecopulib.core.VinedistBase`.

Notes
-----
The :doc:`concepts page </concepts>` introduces the marginal side of Sklar's
theorem, and :class:`pyvinecopulib.core.MarginLike` states the contract these
all satisfy.
"""

from ._adapters import as_margin, register_margin_adapter

# Imported for its side effect as much as its names: it registers the OpenTURNS
# adapter with `as_margin`, and it imports OpenTURNS itself only when used.
from .openturns import OpenTURNSMargin
from .controls import FitControlsMargin
from ._resolve import resolve_margin_controls, resolve_margins
from .scipy import SciPyMargin

__all__ = [
  "FitControlsMargin",
  "OpenTURNSMargin",
  "SciPyMargin",
  "as_margin",
  "resolve_margin_controls",
  "resolve_margins",
  "register_margin_adapter",
]
