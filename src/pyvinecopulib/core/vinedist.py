"""Vine distribution: a vine copula combined with univariate margins."""

from __future__ import annotations

from typing import Any, Generic, Optional, Sequence, cast

from array_api_compat import array_namespace

from .protocols import ArrayT, MarginLike

__all__ = ["Vinedist"]

_TRIM_LO: float = 1e-10
_TRIM_HI: float = 1.0 - 1e-10


class Vinedist(Generic[ArrayT]):
  r"""A multivariate distribution built from a vine copula and its margins.

  By Sklar's theorem a joint distribution factorizes into a copula and its
  marginals, and this class is that factorization made into an object:

  .. math::

     F(\mathbf x) = C\bigl(F_1(x_1), \ldots, F_d(x_d)\bigr), \qquad
     \log f(\mathbf x) = \log c\bigl(F_1(x_1), \ldots, F_d(x_d)\bigr)
       + \sum_{j=1}^d \log f_j(x_j).

  The second identity holds for continuous, discrete and mixed margins alike,
  because each margin's ``pdf`` is a density with respect to its own reference
  measure and the copula factor is normalized by the marginal masses. Nothing
  in the evaluation path branches on the variable type.

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
      One fitted margin per variable. A single margin is accepted when it
      carries array-valued parameters and so already represents all ``d``.

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

      u = pv.utils.to_pseudo_obs(x)
      dist = pv.Vinedist(
        pv.Vinecop.from_data(u),
        margins=[st.Normal(mu=0.0, sigma=1.0), st.gamma(2.0), st.norm(0, 1)],
      )
      dist.logpdf(x)
      dist.simulate(100, seeds=[1])
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
      # One object standing for every variable (array-valued parameters).
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

  def _columns(self, x: ArrayT) -> tuple[Any, Any, int]:
    """Split ``x`` into columns and resolve its array namespace.

    Parameters
    ----------
    x : array, shape (n, d), dtype float
        Observations on the original scale.

    Returns
    -------
    tuple
        The namespace, the array, and its row count.

    Raises
    ------
    ValueError
        If ``x`` does not have ``d`` columns.
    """
    xa: Any = x
    xp = array_namespace(xa)
    if xa.ndim != 2 or xa.shape[1] != self.dim:
      raise ValueError(
        f"x must have shape (n, {self.dim}); got {tuple(xa.shape)}"
      )
    return xp, xa, int(xa.shape[0])

  def marginal_cdf(self, x: ArrayT) -> ArrayT:
    """Apply each margin's ``cdf`` to its column.

    Parameters
    ----------
    x : array, shape (n, d), dtype float
        Observations on the original scale.

    Returns
    -------
    array, shape (n, d), dtype float
        Marginal distribution values, the copula-scale data.
    """
    xp, xa, _ = self._columns(x)
    cols = [m.cdf(xa[:, j]) for j, m in enumerate(self._margins)]
    return cast(ArrayT, xp.clip(xp.stack(cols, axis=-1), _TRIM_LO, _TRIM_HI))

  def marginal_icdf(self, u: ArrayT) -> ArrayT:
    """Apply each margin's ``icdf`` to its column.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Copula-scale values in ``[0, 1]^d``.

    Returns
    -------
    array, shape (n, d), dtype float
        Observations on the original scale.
    """
    xp, ua, _ = self._columns(u)
    cols = [m.icdf(ua[:, j]) for j, m in enumerate(self._margins)]
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
  def copula_data(margins: Sequence[Any], x: Any) -> Any:
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
    x : array, shape (n, d), dtype float
        Observations on the original scale.

    Returns
    -------
    array, shape (n, d + k), dtype float
        The copula-scale data, clamped away from the unit square's boundary.

    Raises
    ------
    ValueError
        If ``x`` does not carry one column per margin, or if a margin reports a
        left limit above its own distribution function.
    """
    from ..margins import as_margin

    resolved = [as_margin(m) for m in margins]
    var_types = Vinedist.copula_var_types(resolved)
    xa: Any = x
    xp = array_namespace(xa)
    if xa.ndim != 2 or xa.shape[1] != len(resolved):
      raise ValueError(
        f"x must have shape (n, {len(resolved)}); got {tuple(xa.shape)}"
      )
    upper = [m.cdf(xa[:, j]) for j, m in enumerate(resolved)]
    lower = []
    for j, m in enumerate(resolved):
      if var_types[j] != "d":
        continue
      left = getattr(m, "cdf_left", None)
      sub = upper[j] if left is None else left(xa[:, j])
      if bool(xp.any(sub > upper[j] + 1e-12)):
        raise ValueError(
          f"margin {j} reports cdf_left > cdf, which cannot happen for a "
          "distribution function; check its var_type and cdf_left"
        )
      lower.append(sub)
    block = xp.stack([*upper, *lower], axis=-1)
    return xp.clip(block, _TRIM_LO, _TRIM_HI)

  def _u_layout(self, x: ArrayT) -> Any:
    """This distribution's copula-scale data for ``x``.

    Parameters
    ----------
    x : array, shape (n, d), dtype float
        Observations on the original scale.

    Returns
    -------
    array, shape (n, d + k), dtype float
        The copula-scale data.
    """
    return self.copula_data(self._margins, x)

  # --- evaluation ---------------------------------------------------------- #

  def logpdf(self, x: ArrayT) -> ArrayT:
    """Log-density of the joint distribution.

    Summed in log space rather than multiplied: for even a moderate ``d`` the
    product of marginal densities underflows long before the sum of their logs
    does, and the marginal term is the one carrying the scale.

    Parameters
    ----------
    x : array, shape (n, d), dtype float
        Observations on the original scale.

    Returns
    -------
    array, shape (n,), dtype float
        Joint log-density values.
    """
    xp, xa, _ = self._columns(x)
    copula_term: Any = self._copula.pdf(cast(ArrayT, self._u_layout(x)))
    total = xp.log(xp.clip(copula_term, 1e-300, None))
    for j, m in enumerate(self._margins):
      logpdf = getattr(m, "logpdf", None)
      if logpdf is not None:
        total = total + logpdf(xa[:, j])
      else:
        dens = m.pdf(xa[:, j])
        positive = dens > 0
        safe = xp.where(positive, dens, xp.ones_like(dens))
        total = total + xp.where(
          positive, xp.log(safe), xp.full_like(dens, float("-inf"))
        )
    return cast(ArrayT, total)

  def pdf(self, x: ArrayT) -> ArrayT:
    """Density of the joint distribution.

    Parameters
    ----------
    x : array, shape (n, d), dtype float
        Observations on the original scale.

    Returns
    -------
    array, shape (n,), dtype float
        Joint density values, with respect to the product of the margins' own
        reference measures.
    """
    out: Any = self.logpdf(x)
    xp = array_namespace(out)
    return cast(ArrayT, xp.exp(out))

  def cdf(self, x: ArrayT, **kwargs: Any) -> ArrayT:
    """Distribution function of the joint distribution.

    Parameters
    ----------
    x : array, shape (n, d), dtype float
        Observations on the original scale.
    **kwargs
        Forwarded to the copula's ``cdf`` (e.g. ``N``, ``seeds``), whose value
        is estimated by Monte-Carlo.

    Returns
    -------
    array, shape (n,), dtype float
        Joint distribution values in ``[0, 1]``.
    """
    # A distribution function needs only the marginal `F` values, and the
    # copula reads only those; but a copula with atoms validates the whole
    # layout, so give it the layout rather than the first block alone.
    return self._copula.cdf(cast(ArrayT, self._u_layout(x)), **kwargs)

  def loglik(self, x: ArrayT) -> ArrayT:
    """Log-likelihood of the observations.

    Parameters
    ----------
    x : array, shape (n, d), dtype float
        Observations on the original scale.

    Returns
    -------
    array, shape (), dtype float
        A 0-d array, so it stays differentiable on autograd backends.
    """
    terms: Any = self.logpdf(x)
    xp = array_namespace(terms)
    return cast(ArrayT, xp.sum(terms))

  def rosenblatt(self, x: ArrayT, **kwargs: Any) -> ArrayT:
    """Rosenblatt transform of observations to independent uniforms.

    Parameters
    ----------
    x : array, shape (n, d), dtype float
        Observations on the original scale.
    **kwargs
        Forwarded to the copula's ``rosenblatt``.

    Returns
    -------
    array, shape (n, d), dtype float
        Independent uniforms.
    """
    return self._copula.rosenblatt(cast(ArrayT, self._u_layout(x)), **kwargs)

  def inverse_rosenblatt(self, w: ArrayT, **kwargs: Any) -> ArrayT:
    """Inverse Rosenblatt transform, from independent uniforms to the data scale.

    Parameters
    ----------
    w : array, shape (n, d), dtype float
        Independent uniforms.
    **kwargs
        Forwarded to the copula's ``inverse_rosenblatt``.

    Returns
    -------
    array, shape (n, d), dtype float
        Observations on the original scale.
    """
    return self.marginal_icdf(self._copula.inverse_rosenblatt(w, **kwargs))

  def simulate(self, n: int, **kwargs: Any) -> ArrayT:
    """Draw ``n`` samples from the joint distribution.

    Parameters
    ----------
    n : int
        Number of samples to draw.
    **kwargs
        Forwarded to the copula's ``simulate`` (e.g. ``qrng``, ``seeds``).

    Returns
    -------
    array, shape (n, d), dtype float
        Samples on the original scale.
    """
    return self.marginal_icdf(self._copula.simulate(n, **kwargs))

  # --- construction from data ---------------------------------------------- #

  @classmethod
  def from_data(
    cls,
    x: Any,
    *,
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
    x : array, shape (n, d), dtype float
        Observations on the original scale.
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

    data = np.asarray(x, dtype=float)
    if data.ndim != 2:
      raise ValueError(f"x must be two-dimensional; got {data.ndim} dimensions")
    d = data.shape[1]

    specs = resolve_margins(margins, d, names=names)
    fitted = [
      fit_margin(specs[j], data[:, j], weights=weights) for j in range(d)
    ]

    # A copula needs its var_types up front, so both come from the fitted
    # margins before the copula exists; `_bind_dist` then re-derives and
    # cross-checks them.
    var_types = cls.copula_var_types(fitted)
    u = cls.copula_data(fitted, data)

    if controls is None:
      controls = FitControlsVinecop()
    if weights is not None and not len(controls.weights):
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
