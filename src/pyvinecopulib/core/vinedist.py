"""Vine distribution: a vine copula combined with univariate margins."""

from __future__ import annotations

import copy
from typing import Any, Generic, Optional, Sequence, cast

from array_api_compat import array_namespace

from .margin_base import _margin_eval, derive_cdf_left
from .protocols import ArrayT, MarginLike


def _named(spec: Any, name: Optional[str]) -> Any:
  """Label a selecting margin with the variable it is fitted to.

  A selector records the variable on each row of its report; without this the
  rows of a multi-variable report are indistinguishable.

  Parameters
  ----------
  spec : object
      A margin specification.
  name : str, or None
      The variable's name, or ``None`` when the data carry none.

  Returns
  -------
  object
      The specification, labeled. Only a nameless selector is copied and
      relabeled, so a specification the caller still holds is left untouched.
  """
  if name is None or not hasattr(spec, "report_"):
    return spec
  if getattr(spec, "name", "") is not None:
    return spec
  spec = copy.deepcopy(spec)
  spec.name = name
  return spec


__all__ = ["Vinedist"]

_TRIM_LO: float = 1e-10
_TRIM_HI: float = 1.0 - 1e-10


def _copula_eval(
  copula: Any, name: str, u: Any, x: Optional[Any], **kwargs: Any
) -> Any:
  """Call a copula method, forwarding ``x`` only when the copula reads it.

  The commonest half of a vine distribution is the compiled ``Vinecop``, whose
  second positional slot is per-row ``parameters`` rather than a conditioning
  matrix, so covariates must not reach it -- and a vine of unconditional pairs
  would refuse them one level further down. A copula declares that its pairs are
  conditional through ``supports_covariates``, exactly as a margin does.
  """
  method = getattr(copula, name)
  if x is None or not getattr(copula, "supports_covariates", False):
    return method(u, **kwargs)
  return method(u, x=x, **kwargs)


class Vinedist(Generic[ArrayT]):
  r"""A multivariate distribution built from a vine copula and its margins.

  By Sklar's theorem a joint distribution factorizes into a copula and its
  marginals, and this class is that factorization made into an object:

  .. math::

     F(\mathbf y) = C\bigl(F_1(y_1), \ldots, F_d(y_d)\bigr), \qquad
     \log f(\mathbf y) = \log c\bigl(F_1(y_1), \ldots, F_d(y_d)\bigr)
       + \sum_{j=1}^d \log f_j(y_j).

  The second identity holds for continuous, discrete and mixed margins alike,
  because each margin's ``pdf`` is a density with respect to its own reference
  measure and the copula factor is normalized by the marginal masses. Nothing
  in the evaluation path branches on the variable type.

    Every method takes the observations as ``y`` and optional exogenous
    covariates as a keyword-only ``x``, which are forwarded to both halves: to
    each margin that declares ``supports_covariates``, and to the copula when it
    declares ``supports_covariates`` too. Conditional margins under a copula that
    ignores covariates are the usual conditional vine, with the dependence
    structure held fixed across covariate values.

  Any :class:`VinecopLike` will do — the compiled
  :class:`~pyvinecopulib.core.Vinecop`, a :class:`VinecopBase` subclass, or the
  PyTorch evaluator — and any :class:`MarginLike`, including distributions from
  other libraries once passed through
  :func:`pyvinecopulib.margins.as_margin`.

  Parameters
  ----------
  copula : VinecopLike
      A fitted vine copula on ``[0, 1]^d``.
  margins : sequence of MarginLike, or MarginLike
      One fitted margin per variable. A single margin stands for every
      variable, which is the identical-margins case; the cascade evaluates it
      once per column, so a margin carrying one parameter *per variable* has to
      be passed as a sequence instead.

  Raises
  ------
  ValueError
      If the number of margins disagrees with the copula's dimension, or if the
      copula's ``var_types`` disagrees with what the margins declare.

  See Also
  --------
  pyvinecopulib.core.MarginLike : The margin contract.
  pyvinecopulib.core.VinecopLike : The copula contract.
  pyvinecopulib.margins.as_margin : Coerce a foreign distribution to a margin.

  Examples
  --------
  Combine a fitted copula with explicit parametric margins::

      import numpy as np, scipy.stats as st, pyvinecopulib as pv

      u = pv.utils.to_pseudo_obs(y)
      dist = pv.Vinedist(
        pv.Vinecop.from_data(u),
        margins=[st.Normal(mu=0.0, sigma=1.0), st.gamma(2.0), st.norm(0, 1)],
      )
      dist.logpdf(y)
      dist.sample(100, seeds=[1])
  """

  def __init__(
    self,
    copula: Any,
    margins: Sequence[Any] | Any,
  ) -> None:
    # `copula` is a `VinecopLike`, but typed `Any`: the compiled `Vinecop`
    # satisfies that contract nominally, not statically (its signatures spell
    # per-row `parameters` where the protocol spells the conditioning matrix
    # `x`), so narrowing here would reject the commonest call. Same reason the
    # sklearn backend layer returns `Any` from `fit_vine`.
    self._bind_dist(copula, margins)

  def _bind_dist(
    self,
    copula: Any,
    margins: Sequence[Any] | Any,
  ) -> None:
    """Install the copula and margins.

    The initialization seam a subclass calls once from its ``__init__``, after
    any framework base class has been initialized.

    Parameters
    ----------
    copula : VinecopLike
        The vine copula.
    margins : sequence of MarginLike, or MarginLike
        The margins.

    Returns
    -------
    None
    """
    from ..margins import as_margin

    d = int(
      getattr(copula, "dim", 0) or len(getattr(copula, "order", ()) or ())
    )
    if not d:
      structure = getattr(copula, "structure", None)
      d = int(getattr(structure, "dim", 0))
    if not d:
      raise ValueError("cannot determine the copula's dimension")

    if isinstance(margins, (list, tuple)):
      resolved = [as_margin(m) for m in margins]
    else:
      # One object standing for every variable. Shared rather than copied:
      # evaluation never mutates a fitted margin, and a copy per variable would
      # duplicate a kernel-density grid `d` times.
      resolved = [as_margin(margins)] * d
    if len(resolved) != d:
      raise ValueError(
        f"got {len(resolved)} margins for a {d}-dimensional copula"
      )

    self._copula = copula
    self._margins = tuple(resolved)
    self._var_types = self.copula_var_types(self._margins)

    declared = getattr(copula, "var_types", None)
    if declared is not None and list(declared) != self._var_types:
      raise ValueError(
        f"the copula declares var_types={list(declared)}, but the margins "
        f"imply {self._var_types}; refit the copula on the matching layout"
      )

  # --- structure ----------------------------------------------------------- #

  @property
  def copula(self) -> Any:
    """The vine copula.

    Returns
    -------
    VinecopLike
        The copula this distribution was built on, typed ``Any`` for the reason
        given on ``__init__``.
    """
    return self._copula

  @property
  def margins(self) -> tuple[MarginLike, ...]:
    """The margins, in variable order.

    Returns
    -------
    tuple of MarginLike
        One margin per variable.
    """
    return self._margins

  def selection_report(self) -> list[dict[str, Any]]:
    """Per-candidate family-selection rows, across every margin that selected.

    Margins that were given rather than selected contribute nothing, so an
    all-fixed or all-KDE distribution reports an empty list.

    Returns
    -------
    list of dict
        One row per candidate considered, in variable order. Each carries the
        variable, the family, its parameter count, the criteria, whether it was
        selected, and — for a candidate that was not fitted — why.
    """
    return [
      dict(row)
      for margin in self._margins
      for row in getattr(margin, "report_", ())
    ]

  @property
  def dim(self) -> int:
    """Number of variables.

    Returns
    -------
    int
        The dimension ``d``.
    """
    return len(self._margins)

  @property
  def var_types(self) -> list[str]:
    """Variable types as the copula sees them.

    Returns
    -------
    list of str
        ``"c"`` or ``"d"`` per variable; a zero-inflated margin is ``"d"``,
        since what the copula needs from it is the left limit.
    """
    return list(self._var_types)

  # --- marginal transforms ------------------------------------------------- #

  def _columns(self, y: ArrayT) -> tuple[Any, Any, int]:
    """Split ``y`` into columns and resolve its array namespace.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.

    Returns
    -------
    tuple
        The namespace, the array, and its row count.

    Raises
    ------
    ValueError
        If ``y`` does not have ``d`` columns.
    """
    ya: Any = y
    xp = array_namespace(ya)
    if ya.ndim != 2 or ya.shape[1] != self.dim:
      raise ValueError(
        f"y must have shape (n, {self.dim}); got {tuple(ya.shape)}"
      )
    return xp, ya, int(ya.shape[0])

  def marginal_cdf(self, y: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Apply each margin's ``cdf`` to its column.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, forwarded to every margin that reads them.

    Returns
    -------
    array, shape (n, d), dtype float
        Marginal distribution values, the copula-scale data.
    """
    self._check_covariates(x)
    _, ya, _ = self._columns(y)
    cols = [
      _margin_eval(m, "cdf", ya[:, j], x) for j, m in enumerate(self._margins)
    ]
    # The margins' namespace, not the input's: a torch copula hosting NumPy
    # margins is legal, and `torch.stack` cannot consume NumPy columns.
    xp = array_namespace(cols[0])
    return cast(ArrayT, xp.clip(xp.stack(cols, axis=-1), _TRIM_LO, _TRIM_HI))

  def marginal_icdf(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Apply each margin's ``icdf`` to its column.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Copula-scale values in ``[0, 1]^d``.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, forwarded to every margin that reads them.

    Returns
    -------
    array, shape (n, d), dtype float
        Observations on the original scale.
    """
    self._check_covariates(x)
    _, ua, _ = self._columns(u)
    cols = [
      _margin_eval(m, "icdf", ua[:, j], x) for j, m in enumerate(self._margins)
    ]
    xp = array_namespace(cols[0])
    return cast(ArrayT, xp.stack(cols, axis=-1))

  @staticmethod
  def copula_var_types(margins: Sequence[Any]) -> list[str]:
    """Variable types a copula must be fitted with to host these margins.

    Parameters
    ----------
    margins : sequence of MarginLike
        One margin per variable, in variable order. Foreign distribution
        objects are coerced with :func:`pyvinecopulib.margins.as_margin`, since
        a bare one carries no variable type of its own.

    Returns
    -------
    list of str
        ``"c"`` or ``"d"`` per variable; a zero-inflated margin is ``"d"``,
        since what the copula needs from it is the left limit.
    """
    from ..margins import as_margin

    return [
      "d" if getattr(as_margin(m), "var_type", "c") in ("d", "zi") else "c"
      for m in margins
    ]

  @staticmethod
  def copula_data(
    margins: Sequence[Any], y: Any, *, x: Optional[Any] = None
  ) -> Any:
    """Assemble the copula-scale data in the layout a copula expects.

    Continuous variables contribute one column each. A variable with atoms
    contributes a second, its left limit ``F(x^-)``, appended after the first
    block in variable order — the compact ``(n, d + k)`` layout. Users never
    build this by hand, which is the point, and neither does code that fits its
    own copula before handing it to ``Vinedist``.

    Parameters
    ----------
    margins : sequence of MarginLike
        One margin per variable, in variable order. Foreign distribution
        objects are coerced with :func:`pyvinecopulib.margins.as_margin`.
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, forwarded to every margin that reads them.

    Returns
    -------
    array, shape (n, d + k), dtype float
        The copula-scale data, clamped away from the unit square's boundary.

    Raises
    ------
    ValueError
        If ``y`` does not carry one column per margin, or if a margin reports a
        left limit above its own distribution function.
    """
    from ..margins import as_margin

    resolved = [as_margin(m) for m in margins]
    var_types = Vinedist.copula_var_types(resolved)
    ya: Any = y
    xp = array_namespace(ya)
    if ya.ndim != 2 or ya.shape[1] != len(resolved):
      raise ValueError(
        f"y must have shape (n, {len(resolved)}); got {tuple(ya.shape)}"
      )
    upper = [
      _margin_eval(m, "cdf", ya[:, j], x) for j, m in enumerate(resolved)
    ]
    lower = []
    for j, m in enumerate(resolved):
      if var_types[j] != "d":
        continue
      # A margin that declares atoms and carries no left limit of its own gets
      # the derived one: copying its `cdf` column would leave the pair copula a
      # zero-width rectangle, the marginal mass would stop canceling, and the
      # joint "density" would not integrate to one.
      left = getattr(m, "cdf_left", None)
      sub = (
        _margin_eval(m, "cdf_left", ya[:, j], x)
        if left is not None
        else derive_cdf_left(m, ya[:, j], x, var_types[j])
      )
      if bool(xp.any(sub > upper[j] + 1e-12)):
        raise ValueError(
          f"margin {j} reports cdf_left > cdf, which cannot happen for a "
          "distribution function; check its var_type and cdf_left"
        )
      lower.append(sub)
    block = xp.stack([*upper, *lower], axis=-1)
    return xp.clip(block, _TRIM_LO, _TRIM_HI)

  def _check_covariates(self, x: Optional[Any]) -> None:
    """Refuse covariates that neither half of this distribution reads.

    Evaluation ignores ``x`` per margin, which is what lets conditional and
    unconditional margins sit in one distribution. When *nothing* reads them the
    silence is indistinguishable from a conditional answer, so it is an error.

    Parameters
    ----------
    x : array, shape (n, k), or None
        The covariates the caller supplied.

    Raises
    ------
    ValueError
        If ``x`` is given and no margin and not the copula declares
        ``supports_covariates``.
    """
    if x is None:
      return
    readers = [
      getattr(m, "supports_covariates", False) for m in self._margins
    ] + [getattr(self._copula, "supports_covariates", False)]
    if not any(readers):
      raise ValueError(
        "covariates were given, but neither the margins nor the copula read "
        "them, so the result would be the unconditional one. Build the "
        "distribution from margins (or a copula) that declare "
        "supports_covariates, or drop x."
      )

  def _u_layout(self, y: ArrayT, x: Optional[ArrayT] = None) -> Any:
    """This distribution's copula-scale data for ``y``.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, forwarded to every margin that reads them.

    Returns
    -------
    array, shape (n, d + k), dtype float
        The copula-scale data.
    """
    self._check_covariates(x)
    return self.copula_data(self._margins, y, x=x)

  # --- evaluation ---------------------------------------------------------- #

  def logpdf(self, y: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Log-density of the joint distribution.

    Summed in log space rather than multiplied: for even a moderate ``d`` the
    product of marginal densities underflows long before the sum of their logs
    does, and the marginal term is the one carrying the scale.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, forwarded to every margin that reads them and to
        the copula.

    Returns
    -------
    array, shape (n,), dtype float
        Joint log-density values.
    """
    xp, ya, _ = self._columns(y)
    copula_term: Any = _copula_eval(
      self._copula, "pdf", cast(ArrayT, self._u_layout(y, x)), x
    )
    total = xp.log(xp.clip(copula_term, 1e-300, None))
    for j, m in enumerate(self._margins):
      if getattr(m, "logpdf", None) is not None:
        total = total + _margin_eval(m, "logpdf", ya[:, j], x)
      else:
        dens = _margin_eval(m, "pdf", ya[:, j], x)
        positive = dens > 0
        safe = xp.where(positive, dens, xp.ones_like(dens))
        total = total + xp.where(
          positive, xp.log(safe), xp.full_like(dens, float("-inf"))
        )
    return cast(ArrayT, total)

  def pdf(self, y: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Density of the joint distribution.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, forwarded to every margin that reads them and to
        the copula.

    Returns
    -------
    array, shape (n,), dtype float
        Joint density values, with respect to the product of the margins' own
        reference measures.
    """
    out: Any = self.logpdf(y, x=x)
    xp = array_namespace(out)
    return cast(ArrayT, xp.exp(out))

  def cdf(
    self, y: ArrayT, *, x: Optional[ArrayT] = None, **kwargs: Any
  ) -> ArrayT:
    """Distribution function of the joint distribution.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, forwarded to every margin that reads them and to
        the copula.
    **kwargs
        Forwarded to the copula's ``cdf`` (e.g. ``N``, ``seeds``), whose value
        is estimated by Monte-Carlo.

    Returns
    -------
    array, shape (n,), dtype float
        Joint distribution values in ``[0, 1]``.
    """
    # ``C(F_1(y_1), ..., F_d(y_d))`` needs no left limits, but a copula with
    # discrete variables accepts only the layouts that carry them, so hand it
    # the full layout and let it drop what it does not read.
    return _copula_eval(
      self._copula, "cdf", cast(ArrayT, self._u_layout(y, x)), x, **kwargs
    )

  def loglik(self, y: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Log-likelihood of the observations.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, forwarded to every margin that reads them and to
        the copula.

    Returns
    -------
    array, shape (), dtype float
        A 0-d array, so it stays differentiable on autograd backends.
    """
    terms: Any = self.logpdf(y, x=x)
    xp = array_namespace(terms)
    return cast(ArrayT, xp.sum(terms))

  def rosenblatt(
    self, y: ArrayT, *, x: Optional[ArrayT] = None, **kwargs: Any
  ) -> ArrayT:
    """Rosenblatt transform of observations to independent uniforms.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, forwarded to every margin that reads them and to
        the copula.
    **kwargs
        Forwarded to the copula's ``rosenblatt``.

    Returns
    -------
    array, shape (n, d), dtype float
        Independent uniforms.
    """
    return _copula_eval(
      self._copula,
      "rosenblatt",
      cast(ArrayT, self._u_layout(y, x)),
      x,
      **kwargs,
    )

  def inverse_rosenblatt(
    self, w: ArrayT, *, x: Optional[ArrayT] = None, **kwargs: Any
  ) -> ArrayT:
    """Inverse Rosenblatt transform, from independent uniforms to the data scale.

    Parameters
    ----------
    w : array, shape (n, d), dtype float
        Independent uniforms.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, forwarded to every margin that reads them and to
        the copula.
    **kwargs
        Forwarded to the copula's ``inverse_rosenblatt``.

    Returns
    -------
    array, shape (n, d), dtype float
        Observations on the original scale.
    """
    u = _copula_eval(self._copula, "inverse_rosenblatt", w, x, **kwargs)
    return self.marginal_icdf(u, x=x)

  def sample(
    self, n: int, *, x: Optional[ArrayT] = None, **kwargs: Any
  ) -> ArrayT:
    """Draw ``n`` samples from the joint distribution.

    Parameters
    ----------
    n : int
        Number of samples to draw. With covariates, ``x`` supplies one row per
        draw, so it must have ``n`` of them.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, forwarded to every margin that reads them and to
        the copula.
    **kwargs
        Forwarded to the copula's ``sample`` (e.g. ``qrng``, ``seeds``).

    Returns
    -------
    array, shape (n, d), dtype float
        Samples on the original scale.
    """
    u = _copula_eval(self._copula, "sample", n, x, **kwargs)
    return self.marginal_icdf(u, x=x)

  # --- construction from data ---------------------------------------------- #

  @classmethod
  def from_data(
    cls,
    y: Any,
    *,
    x: Optional[Any] = None,
    margins: Any = None,
    controls: Optional[Any] = None,
    structure: Optional[Any] = None,
    weights: Optional[Any] = None,
    names: Optional[Sequence[str]] = None,
  ) -> "Vinedist[Any]":
    """Fit margins and a vine copula to data, in that order.

    The two-step estimator: each margin is fitted from its own column, the data
    are transformed to the copula scale, and the copula is fitted on the result
    (Joe and Xu, 1996). Margins already fitted are left alone, so a fixed margin
    and one to estimate can be mixed freely.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates. Each margin that declares
        ``supports_covariates`` is fitted conditionally on them, and the copula
        is then fitted on the resulting conditional probability-integral
        transforms — the two-step estimator, with conditional margins.
    margins : object, optional
        What to use for each variable; see
        :func:`pyvinecopulib.margins.resolve_margins` for the accepted forms.
        The default is a kernel-density margin per variable.
    controls : FitControlsVinecop, or None, optional
        Copula fit controls.
    structure : RVineStructure, or None, optional
        A fixed vine structure; selected from the data when ``None``.
    weights : array, shape (n,), or None, optional
        Observation weights, applied to both the margins and the copula.
    names : sequence of str, or None, optional
        Variable names, so ``margins`` may be a mapping keyed by name.

    Returns
    -------
    Vinedist
        The fitted distribution.

    References
    ----------
    .. [1] Joe, H. and Xu, J. J. (1996). *The estimation method of inference
           functions for margins for multivariate models.* Technical Report
           166, Department of Statistics, University of British Columbia.
    """
    import numpy as np

    from ..margins import resolve_margins
    from ..margins._resolve import fit_margin
    from ..pyvinecopulib_ext import FitControlsVinecop, Vinecop

    data = np.asarray(y, dtype=float)
    if data.ndim != 2:
      raise ValueError(f"y must be two-dimensional; got {data.ndim} dimensions")
    d = data.shape[1]

    if names is None:
      # A DataFrame carries its own names, and `margins` is often keyed by them.
      # Duck-typed: pandas is an extra, not a dependency of `core`.
      columns = getattr(y, "columns", None)
      if columns is not None:
        names = [str(c) for c in columns]

    specs = resolve_margins(margins, d, names=names)
    if x is not None and not any(
      getattr(spec, "supports_covariates", False) for spec in specs
    ):
      raise ValueError(
        "covariates were given, but no margin reads them: every specification "
        "in `margins` is unconditional, so the fit would ignore `x` and return "
        "a model of f(y) labeled as one of f(y | x). Pass margins that declare "
        "supports_covariates, or drop x."
      )
    fitted = [
      fit_margin(
        _named(specs[j], names[j] if names is not None else None),
        data[:, j],
        x=x,
        weights=weights,
      )
      for j in range(d)
    ]

    # A copula needs its var_types up front, so both come from the fitted
    # margins before the copula exists; `_bind_dist` then re-derives and
    # cross-checks them.
    var_types = cls.copula_var_types(fitted)
    u = cls.copula_data(fitted, data, x=x)

    if controls is None:
      controls = FitControlsVinecop()
    if weights is not None and not len(controls.weights):
      # Copied first: the caller's controls object is theirs, and a `weights`
      # written into it would silently weight every later fit that reuses it.
      controls = copy.deepcopy(controls)
      controls.weights = np.asarray(weights, dtype=float).ravel()
    copula = Vinecop.from_data(
      data=u,
      structure=structure,
      var_types=var_types,
      controls=controls,
    )
    return cls(copula, list(fitted))

  def __repr__(self) -> str:
    """Return a structural representation.

    Returns
    -------
    str
        The dimension and the margin families.
    """
    families = ", ".join(
      str(getattr(m, "family_name", type(m).__name__)) for m in self._margins
    )
    return f"Vinedist(dim={self.dim}, margins=[{families}])"
