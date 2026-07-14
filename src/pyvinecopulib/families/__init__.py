"""Bivariate copula families.

Every bivariate copula model in pyvinecopulib (and every pair-copula
inside a vine) belongs to one of the families documented below. The
:class:`BicopFamily` enum holds the family tag; the module-level
constants ``indep``, ``gaussian``, ... are aliases for the enum members
so that callers can write ``family_set=[gaussian, student]`` without
the qualified prefix. The named family-group constants (``parametric``,
``elliptical``, ``itau``, ...) are pre-built lists you can pass to
:class:`pyvinecopulib.core.FitControlsBicop` /
:class:`pyvinecopulib.core.FitControlsVinecop` to constrain the
fitting search space.

See Also
--------
pyvinecopulib.core.Bicop : Bivariate copula model. Its ``family``
    attribute is a `BicopFamily` value.
pyvinecopulib.core.FitControlsBicop, pyvinecopulib.core.FitControlsVinecop
    Their ``family_set`` argument accepts any of the lists above
    (or a custom subset).

Notes
-----
**Available families**, grouped by mathematical class:

- *Independence* — ``indep``.
- *Elliptical* — ``gaussian``, ``student``.
- *Archimedean* (single-parameter) — ``clayton``, ``gumbel``,
  ``frank``, ``joe``.
- *Archimedean — BB family* (two-parameter extensions) — ``bb1``,
  ``bb6``, ``bb7``, ``bb8``.
- *Extreme-value* — ``tawn`` (three parameters); ``gumbel`` is also
  an extreme-value copula.
- *Nonparametric* — ``tll`` (transformation local-likelihood
  estimator).

See :class:`BicopFamily` for the per-member identifier.

**Family groups** — convenience lists of :class:`BicopFamily` values:

- ``all`` — every family above.
- ``parametric`` — every family except ``tll``.
- ``nonparametric`` — ``indep``, ``tll``.
- ``one_par`` — single-parameter parametric families
  (``gaussian``, ``clayton``, ``gumbel``, ``frank``, ``joe``).
- ``two_par`` — two-parameter parametric families
  (``student``, ``bb1``, ``bb6``, ``bb7``, ``bb8``).
- ``three_par`` — three-parameter parametric families (``tawn``).
- ``elliptical`` — ``gaussian``, ``student``.
- ``archimedean`` — ``clayton``, ``gumbel``, ``frank``, ``joe``,
  ``bb1``, ``bb6``, ``bb7``, ``bb8``.
- ``extreme_value`` — ``tawn``, ``gumbel``.
- ``bb`` — ``bb1``, ``bb6``, ``bb7``, ``bb8``.
- ``itau`` — families that support estimation by Kendall's-:math:`\\tau`
  inversion (``indep``, ``gaussian``, ``student``, ``clayton``,
  ``gumbel``, ``frank``, ``joe``).
- ``analytic_derivs`` — families with closed-form derivatives of the
  density and h-functions (currently every parametric family).
  ``Bicop.pdf_deriv`` and friends use these closed forms; parametric
  families outside the group fall back to central finite differences.
- ``lt`` — lower-tail dependent families (``clayton``, ``bb1``,
  ``bb7``, ``tawn``).
- ``ut`` — upper-tail dependent families (``gumbel``, ``joe``,
  ``bb1``, ``bb6``, ``bb7``, ``bb8``, ``tawn``).
- ``rotationless`` — families that already cover positive and negative
  dependence and therefore have no rotation
  (``indep``, ``gaussian``, ``student``, ``frank``, ``tll``).

The :doc:`concepts page </concepts>` provides a tabular overview of
every family above with parameter ranges, rotations, tail
dependence, and closed-form Kendall's :math:`\\tau`.
"""

from ..pyvinecopulib_ext import (
  BicopFamily,
  all,
  analytic_derivs,
  archimedean,
  bb,
  bb1,
  bb6,
  bb7,
  bb8,
  clayton,
  elliptical,
  extreme_value,
  frank,
  gaussian,
  gumbel,
  indep,
  itau,
  joe,
  lt,
  nonparametric,
  one_par,
  parametric,
  rotationless,
  student,
  tawn,
  three_par,
  tll,
  two_par,
  ut,
)

__all__ = [
  "BicopFamily",
  # Individual families
  "indep",
  "gaussian",
  "student",
  "clayton",
  "gumbel",
  "frank",
  "joe",
  "bb1",
  "bb6",
  "bb7",
  "bb8",
  "tawn",
  "tll",
  # Family groups
  "all",
  "parametric",
  "nonparametric",
  "one_par",
  "two_par",
  "three_par",
  "elliptical",
  "archimedean",
  "extreme_value",
  "bb",
  "rotationless",
  "lt",
  "ut",
  "itau",
  "analytic_derivs",
]
