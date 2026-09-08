"""Canonical vine distribution: the array-agnostic half of ``Vinedist``.

:class:`VinedistBase` is a :class:`~pyvinecopulib.core.VinecopLike` combined
with one :class:`~pyvinecopulib.core.MarginLike` per variable, evaluated on the
**data** scale rather than the copula scale: the joint density is Sklar's
factorization, and the Rosenblatt transforms and simulation are the copula's
own sandwiched between the marginal transforms. The cascades run on any array
namespace (NumPy or PyTorch), so a concrete subclass —
:class:`~pyvinecopulib.core.Vinedist` or
:class:`~pyvinecopulib.torch.TorchVinedist` — inherits the whole evaluation
surface as soon as its two halves are installed.

Two things live here rather than in a subclass. The copula-scale **layout** is
assembled in one place, :meth:`VinedistBase.copula_data`: one column per
variable, then a left-limit column for each variable with atoms, so nothing
downstream builds that block by hand. And the **two-step (IFM) estimator** is
here as well — :meth:`VinedistBase.from_data` for a fresh distribution,
:meth:`VinedistBase.fit` and :meth:`VinedistBase.select` for the parts an
object already holds — so the lanes cannot drift apart in what they validate,
or in which half of a fit a setting reaches. What each lane declares instead of
implementing is enumerated on :class:`VinedistBase` itself.

Arrays are typed ``Any`` inside the cascades per the
:mod:`pyvinecopulib.core` typing policy (the Array API namespace is untyped);
the public signatures carry ``ArrayT`` so a concrete subclass inherits precise
return types.
"""

from __future__ import annotations

from abc import ABC
from typing import Any, ClassVar, Optional, Self, Sequence, cast

from array_api_compat import array_namespace

from ._covariates import declared_eval, prepare
from .margin_base import derive_cdf_left, safe_log
from ._placement import PlacementMixin
from ._trim import trim
from ._validation import validate_covariates, validate_weights
from .protocols import (
  ArrayT,
  ControlsLike,
  MarginLike,
  VinecopLike,
  VinedistLike,
)
from .protocols import _VINEDIST_EXAMPLE


__all__ = ["VinedistBase"]


class VinedistBase(VinedistLike[ArrayT], PlacementMixin, ABC):
  r"""Canonical vine distribution: a vine copula and one margin per variable.

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

  Everything is read and returned on the **original** scale: :meth:`pdf` /
  :meth:`logpdf` / :meth:`cdf` / :meth:`loglik` / :meth:`rosenblatt` /
  :meth:`inverse_rosenblatt` / :meth:`sample` / :meth:`sample_conditional`
  take observations rather than copula arguments, and :meth:`copula_layout`
  hands back the copula-scale data they run on for a caller who wants it. Each
  of them also takes optional exogenous covariates as a keyword-only ``x``,
  forwarded to both halves: to every margin that declares
  ``supports_covariates``, and to the copula when it declares it too.
  Conditional margins under a copula that ignores covariates model ``Y | X``
  with the dependence held fixed across covariate values.

  Any :class:`VinecopLike` will do —
  :class:`~pyvinecopulib.core.Vinecop`, a :class:`VinecopBase` subclass, or the
  PyTorch evaluator — and any :class:`MarginLike`, including distributions from
  other libraries once passed through
  :func:`pyvinecopulib.margins.as_margin`.

  Parameters
  ----------
  vinecop : VinecopLike
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

  Notes
  -----
  Unlike :class:`~pyvinecopulib.core.BicopBase` and
  :class:`~pyvinecopulib.core.VinecopBase`, which each leave one evaluation
  hook for a subclass, this base needs **none**: a vine distribution is fully
  determined by its two halves, so every method above works the moment they are
  installed, and the class is directly usable.

  Fitting is the only namespace-specific half, and it too is implemented here,
  once: :meth:`from_data` fits a fresh distribution, while :meth:`fit` and
  :meth:`select` re-estimate the parts an object already holds. What a lane
  contributes is mostly *declared*:

  - ``vinecop_class`` and ``margin_class`` name the classes the two halves
    are. Naming them is what lets the estimators fit with no callback of any
    kind, a part class being its own fitter. Without a ``vinecop_class`` a
    lane cannot fit a copula at all unless it overrides ``_fit_copula``;
    without a ``margin_class`` an unaddressed variable falls back to a
    kernel-density margin.
  - ``_coerce_fit_data`` is the one real hook, and the only member fitting
    strictly requires: a lane may have to resolve a placement — a device and a
    dtype — before any part exists to read one from.
  - ``_fit_copula`` fits the copula half and defaults to ``vinecop_class``.
    Override it where a lane needs a step of its own:
    :class:`~pyvinecopulib.core.Vinedist` writes ``weights`` into a copy of the
    controls there.
  - ``_default_margins`` supplies the margin each variable gets when the caller
    named none, one ``margin_class`` per variable. Override it when those
    margins need an argument the class attribute cannot carry, as a placement
    is; override ``_default_margins`` alone when the margin class takes
    its variable type or its bounds at *construction*, as a kernel density
    does.
  - ``supports_weighted_copula`` and ``supports_fit_covariates`` declare what
    the lane can honor. Either one ``False`` makes the matching request an
    error up front, instead of a fit that applies it to one half only.

  See Also
  --------
  pyvinecopulib.core.Vinedist : The NumPy subclass, which fits a
      :class:`~pyvinecopulib.core.Vinecop`.
  pyvinecopulib.torch.TorchVinedist : The PyTorch subclass.
  pyvinecopulib.core.VinedistLike : The contract this implements.
  pyvinecopulib.core.MarginLike : The margin contract.
  pyvinecopulib.core.VinecopLike : The copula contract.
  pyvinecopulib.margins.as_margin : Coerce a foreign distribution to a margin.
  """

  def __init__(
    self,
    vinecop: Any,
    margins: Sequence[Any] | Any,
  ) -> None:
    # `copula` is a `VinecopLike`, but typed `Any`: the compiled `Vinecop`
    # satisfies that contract nominally, not statically (its signatures spell
    # per-row `parameters` where the protocol spells the conditioning matrix
    # `x`), so narrowing here would reject the commonest call. Same reason the
    # sklearn backend layer returns `Any` from `fit_vine`.
    self._bind_dist(vinecop, margins)

  def _bind_dist(
    self,
    vinecop: Any,
    margins: Sequence[Any] | Any,
  ) -> None:
    """Install the copula and margins.

    The initialization seam a subclass calls once from its ``__init__``, after
    any framework base class has been initialized.

    Parameters
    ----------
    vinecop : VinecopLike
        The vine copula.
    margins : sequence of MarginLike, or MarginLike
        The margins.

    Returns
    -------
    None
    """
    from ._adapters import as_margin

    d = int(
      getattr(vinecop, "dim", 0) or len(getattr(vinecop, "order", ()) or ())
    )
    if not d:
      structure = getattr(vinecop, "structure", None)
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

    self._vinecop = vinecop
    self._set_margins(resolved)
    self._var_types = self.copula_var_types(self._margins)

    declared = getattr(vinecop, "var_types", None)
    if declared is not None and list(declared) != self._var_types:
      raise ValueError(
        f"the copula declares var_types={list(declared)}, but the margins "
        f"imply {self._var_types}; refit the copula on the matching layout"
      )

  # --- structure ----------------------------------------------------------- #

  @property
  def vinecop(self) -> VinecopLike[ArrayT]:
    """The vine copula.

    Returns
    -------
    VinecopLike
        The copula this distribution was built on.
    """
    return cast(VinecopLike[ArrayT], self._vinecop)

  @property
  def margins(self) -> tuple[MarginLike[ArrayT], ...]:
    """The margins, in variable order.

    Returns
    -------
    tuple of MarginLike
        One margin per variable.
    """
    return self._margins

  def margin_summary(self) -> list[dict[str, Any]]:
    """One row per variable describing the margin that models it.

    Every field is optional, so a margin from another ecosystem contributes
    whatever it declares and ``None`` for the rest.

    Returns
    -------
    list of dict
        One row per variable in variable order, carrying its position, the
        margin's own ``name`` if it has one, the class, the family, the variable
        type as the copula sees it, the support, the number of parameters and
        the log-likelihood attained at the fit.
    """
    rows: list[dict[str, Any]] = []
    for j, margin in enumerate(self._margins):
      loglik = getattr(margin, "loglik", None)
      try:
        value = float(loglik()) if callable(loglik) else None
      except (RuntimeError, TypeError, ValueError, NotImplementedError):
        # `loglik()` with no data is only defined for a margin that was fitted
        # here; a fixed or foreign one has no fit to report.
        value = None
      rows.append(
        {
          "variable": j,
          "name": getattr(margin, "name", None),
          "margin": type(margin).__name__,
          "family": getattr(margin, "family_name", None),
          "var_type": self._var_types[j],
          "support": getattr(margin, "support", None),
          "n_parameters": getattr(margin, "n_parameters", None),
          "loglik": value,
        }
      )
    return rows

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

  # --- persistence --------------------------------------------------------- #

  def to_json(self) -> str:
    """Serialize to a JSON string.

    Both halves are stored: the copula through its own ``to_json``, and one
    payload per margin. A margin type this package does not ship must provide
    ``to_json`` and register a reader with
    :func:`~pyvinecopulib.core.register_margin_json`.

    Returns
    -------
    str
        A JSON string that :meth:`from_json` reads back.

    Raises
    ------
    TypeError
        If the copula or a margin cannot be serialized.

    See Also
    --------
    to_file : Write to a file, JSON or CBOR.
    """
    from ._serialization import MARGIN_JSON_VERSION, dumps, margin_to_json

    copula_to_json = getattr(self._vinecop, "to_json", None)
    if copula_to_json is None:
      raise TypeError(
        f"{type(self._vinecop).__name__} cannot be serialized: it has no "
        "`to_json`. Hold a copula that provides one, or persist this half "
        "the way its own library does."
      )
    return dumps(
      {
        "kind": type(self).__name__,
        "version": MARGIN_JSON_VERSION,
        "copula": copula_to_json(),
        "margins": [margin_to_json(m) for m in self._margins],
      }
    )

  def to_file(self, filename: str) -> None:
    """Write to a JSON file, or a CBOR file when the name ends in ``.cbor``.

    Parameters
    ----------
    filename : str
        Path to write.
    """
    from ._serialization import write_file

    write_file(filename, self.to_json())

  @classmethod
  def from_file(cls, filename: str) -> Self:
    """Instantiate from a JSON file, or a CBOR file by extension.

    Parameters
    ----------
    filename : str
        Path to read.

    Returns
    -------
    VinedistBase
        The distribution the file holds.

    See Also
    --------
    to_file : Write one out.
    """
    from ._serialization import read_file

    return cls.from_json(read_file(filename))

  # --- marginal transforms ------------------------------------------------- #

  def _set_margins(self, margins: Sequence[Any]) -> None:
    """Store the resolved margins on this distribution.

    A hook rather than a plain assignment because a lane may need its margins
    *registered* rather than merely held -- and because ``_bind_dist`` runs
    again on every refit, so the store has to be idempotent. The write
    counterpart of :attr:`margins`, and the same shape as
    :meth:`~pyvinecopulib.core.VinecopBase.set_pair_copulas` on the vine lane.
    """
    self._margins = tuple(margins)

  def _columns(self, y: ArrayT) -> tuple[Any, Any, int]:
    """Coerce ``y``, check its width, and resolve its array namespace.

    Raises
    ------
    ValueError
        If ``y`` does not have ``d`` columns.
    """
    ya: Any = self._prep(y)
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
    x : array, shape (n, p), or None, optional
        Exogenous covariates, forwarded to every margin that reads them.

    Returns
    -------
    array, shape (n, d), dtype float
        Marginal distribution values, on the array namespace the margins
        answer in -- which need not be the one ``y`` was passed on.
    """
    _, ya, n = self._columns(y)
    self._check_covariates(x, n)
    cols = [
      declared_eval(m, "cdf", ya[:, j], x) for j, m in enumerate(self._margins)
    ]
    # The margins' namespace, not the input's: a torch copula hosting NumPy
    # margins is legal, and `torch.stack` cannot consume NumPy columns.
    xp = array_namespace(cols[0])
    return cast(ArrayT, trim(xp, xp.stack(cols, axis=-1)))

  def marginal_icdf(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Apply each margin's ``icdf`` to its column.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Copula-scale values in ``[0, 1]^d``.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, forwarded to every margin that reads them.

    Returns
    -------
    array, shape (n, d), dtype float
        Observations on the original scale, on the array namespace the margins
        answer in.
    """
    _, ua, n = self._columns(u)
    self._check_covariates(x, n)
    cols = [
      declared_eval(m, "icdf", ua[:, j], x) for j, m in enumerate(self._margins)
    ]
    xp = array_namespace(cols[0])
    return cast(ArrayT, xp.stack(cols, axis=-1))

  @classmethod
  def copula_var_types(cls, margins: Sequence[Any]) -> list[str]:
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
    from ._adapters import as_margin

    return [
      "d" if getattr(as_margin(m), "var_type", "c") in ("d", "zi") else "c"
      for m in margins
    ]

  @classmethod
  def copula_data(
    cls, margins: Sequence[Any], y: Any, *, x: Optional[Any] = None
  ) -> Any:
    """Assemble the copula-scale data in the layout a copula expects.

    Continuous variables contribute one column each. A variable with atoms
    contributes a second, its left limit ``F(y^-)``, appended after the first
    block in variable order — the compact ``(n, d + k)`` layout. So no caller
    assembles a left-limit block by hand, which is the point — not even code
    that fits its own copula before handing it to ``Vinedist``.

    Parameters
    ----------
    margins : sequence of MarginLike
        One margin per variable, in variable order. Foreign distribution
        objects are coerced with :func:`pyvinecopulib.margins.as_margin`.
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, forwarded to every margin that reads them.

    Returns
    -------
    array, shape (n, d + k), dtype float
        The copula-scale data, clamped strictly inside ``(0, 1)``.

    Raises
    ------
    ValueError
        If ``y`` does not carry one column per margin, or if a margin reports a
        left limit above its own distribution function.
    """
    from ._adapters import as_margin

    resolved = [as_margin(m) for m in margins]
    var_types = cls.copula_var_types(resolved)
    ya: Any = y
    if ya.ndim != 2 or ya.shape[1] != len(resolved):
      raise ValueError(
        f"y must have shape (n, {len(resolved)}); got {tuple(ya.shape)}"
      )
    x = prepare(ya, x, int(ya.shape[0]))
    upper = [
      declared_eval(m, "cdf", ya[:, j], x) for j, m in enumerate(resolved)
    ]
    # The margins' namespace, not the input's, as `marginal_cdf` does: a margin
    # may legitimately return another array type than it was handed, and
    # stacking that through the input's namespace either raises or silently
    # detaches.
    xp = array_namespace(upper[0])
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
        declared_eval(m, "cdf_left", ya[:, j], x)
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
    return trim(xp, block)

  def _check_covariates(self, x: Optional[Any], n_rows: int) -> None:
    """Refuse covariates that neither half of this distribution reads.

    Evaluation ignores ``x`` per margin, which is what lets conditional and
    unconditional margins sit in one distribution. When *nothing* reads them the
    silence is indistinguishable from a conditional answer, so it is an error.

    Raises
    ------
    ValueError
        If ``x`` is not two-dimensional and row-aligned, or if no margin nor
        the copula declares ``supports_covariates``.
    """
    if x is None:
      return
    # Validated, not placed: every entry point forwards the caller's `x` on and
    # each part places it for itself, so placing here would be a tensor copy
    # per call whose result nothing reads.
    validate_covariates(x, n_rows)
    readers = [
      getattr(m, "supports_covariates", False) for m in self._margins
    ] + [getattr(self._vinecop, "supports_covariates", False)]
    if not any(readers):
      raise ValueError(
        "covariates were given, but neither the margins nor the copula read "
        "them, so the result would be the unconditional one. Build the "
        "distribution from margins (or a copula) that declare "
        "supports_covariates, or drop x."
      )

  def copula_layout(self, y: ArrayT, *, x: Optional[ArrayT] = None) -> Any:
    """This distribution's copula-scale data for ``y``.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, forwarded to every margin that reads them.

    Returns
    -------
    array, shape (n, d + k), dtype float
        The compact layout: one column per variable, then a left-limit column
        for each variable with atoms.

    See Also
    --------
    copula_data : The same assembly, over margins the caller supplies.
    """
    _, ya, n = self._columns(y)
    self._check_covariates(x, n)
    return self.copula_data(self._margins, ya, x=x)

  # --- evaluation ---------------------------------------------------------- #

  def logpdf(self, y: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Log-density of the joint distribution.

    The primitive, rather than the log of :meth:`pdf`: the marginal term is the
    one carrying the scale, and for even a moderate ``d`` the product of
    marginal densities underflows long before the sum of their logs does.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, forwarded to every margin that reads them and to
        the copula.

    Returns
    -------
    array, shape (n,), dtype float
        Joint log-density values.
    """
    _, ya, _ = self._columns(y)
    copula_term: Any = declared_eval(
      self._vinecop, "pdf", cast(ArrayT, self.copula_layout(y, x=x)), x
    )
    # The parts' namespace, not the input's, as `marginal_cdf` and
    # `copula_data` both do: a copula or a margin may legitimately answer in
    # another array type than it was handed -- a torch copula hosting NumPy
    # margins is legal -- and operating on that through the input's namespace
    # either raises or silently detaches. Each term is then coerced onto the
    # accumulator for the same reason: adding an ndarray to a tensor that
    # tracks grad sends NumPy looking for `__array__` and raises.
    xp = array_namespace(copula_term)
    total = xp.log(copula_term)
    for j, m in enumerate(self._margins):
      if getattr(m, "logpdf", None) is not None:
        term: Any = declared_eval(m, "logpdf", ya[:, j], x)
      else:
        term = safe_log(declared_eval(m, "pdf", ya[:, j], x))
      # Only when they actually differ: coercing a tensor that already
      # tracks grad through `asarray` warns about the flag it inherits.
      if array_namespace(term) is not xp:
        term = xp.asarray(term)
      total = total + term
    return cast(ArrayT, total)

  def pdf(self, y: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Density of the joint distribution.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, forwarded to every margin that reads them and to
        the copula.

    Returns
    -------
    array, shape (n,), dtype float
        Joint density values, with respect to the product of the margins' own
        reference measures -- a Lebesgue density where a margin is continuous,
        a probability at an atom.
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
    x : array, shape (n, p), or None, optional
        Exogenous covariates, forwarded to every margin that reads them and to
        the copula.
    **kwargs
        Forwarded to the copula's ``cdf`` (e.g. ``N``, ``seeds``), which
        estimates its value by Monte-Carlo.

    Returns
    -------
    array, shape (n,), dtype float
        Joint distribution values in ``[0, 1]``.

    Raises
    ------
    NotImplementedError
        If the copula reads ``x`` but does not implement a conditional CDF.
        In particular, ``VinecopBase``'s generic Monte-Carlo CDF supports only
        an unconditional copula.
    """
    # ``C(F_1(y_1), ..., F_d(y_d))`` needs no left limits, but a copula with
    # discrete variables accepts only the layouts that carry them, so hand it
    # the full layout and let it drop what it does not read.
    return declared_eval(
      self._vinecop,
      "cdf",
      cast(ArrayT, self.copula_layout(y, x=x)),
      x,
      **kwargs,
    )

  def loglik(self, y: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Log-likelihood of the observations.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, forwarded to every margin that reads them and to
        the copula.

    Returns
    -------
    array, shape (), dtype float
        A 0-d array, so it stays differentiable on an array namespace that
        carries gradients.
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
    x : array, shape (n, p), or None, optional
        Exogenous covariates, forwarded to every margin that reads them and to
        the copula.
    **kwargs
        Forwarded to the copula's ``rosenblatt``.

    Returns
    -------
    array, shape (n, d), dtype float
        Independent uniforms.
    """
    return declared_eval(
      self._vinecop,
      "rosenblatt",
      cast(ArrayT, self.copula_layout(y, x=x)),
      x,
      **kwargs,
    )

  def inverse_rosenblatt(
    self, w: ArrayT, *, x: Optional[ArrayT] = None, **kwargs: Any
  ) -> ArrayT:
    """Inverse Rosenblatt transform, from uniforms to the original scale.

    Parameters
    ----------
    w : array, shape (n, d), dtype float
        Independent uniforms.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, forwarded to every margin that reads them and to
        the copula.
    **kwargs
        Forwarded to the copula's ``inverse_rosenblatt``.

    Returns
    -------
    array, shape (n, d), dtype float
        Observations on the original scale.
    """
    _, wa, n = self._columns(w)
    self._check_covariates(x, n)
    u = declared_eval(
      self._vinecop, "inverse_rosenblatt", cast(ArrayT, wa), x, **kwargs
    )
    return self.marginal_icdf(u, x=x)

  def sample(
    self, n: int, *, x: Optional[ArrayT] = None, **kwargs: Any
  ) -> ArrayT:
    """Draw ``n`` samples from the joint distribution.

    The copula's own draw, mapped through each margin's ``icdf``: the copula's
    quasi-random and seeding options therefore apply, and no margin's own
    sampler is consulted.

    Parameters
    ----------
    n : int
        Number of samples to draw. With covariates, ``x`` supplies one row per
        draw, so it must have ``n`` of them.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, forwarded to every margin that reads them and to
        the copula.
    **kwargs
        Forwarded to the copula's ``sample`` (e.g. ``qrng``, ``seeds``).

    Returns
    -------
    array, shape (n, d), dtype float
        Samples on the original scale.
    """
    self._check_covariates(x, n)
    u = declared_eval(self._vinecop, "sample", n, x, **kwargs)
    return self.marginal_icdf(u, x=x)

  def sample_conditional(
    self,
    y_cond: ArrayT,
    *,
    conditioning_set: Optional[list[int]] = None,
    x: Optional[ArrayT] = None,
    **kwargs: Any,
  ) -> ArrayT:
    """Sample the remaining variables given fixed values of some of them.

    The data-scale counterpart of the copula's ``sample_conditional``: each
    row of ``y_cond`` is one conditioning point on the **original** scale, and
    the matching output row draws every variable from its distribution
    conditional on it. To draw many samples at one point, pass that point
    repeated over ``n`` rows.

    Parameters
    ----------
    y_cond : array, shape (n, k), dtype float
        Conditioning values on the original scale, one point per row. Column
        ``i`` is the value of conditioning variable ``i``. Unlike the copula
        scale, a discrete conditioner needs no left-limit column: it is derived
        from that variable's margin.
    conditioning_set : list of int or None, optional
        The 1-based variables to condition on, so column ``i`` of ``y_cond``
        is variable ``conditioning_set[i]``. ``None`` takes the last ``k``
        variables of the copula's sampling order, ``k`` being ``y_cond``'s
        width -- the same convention the copula scale uses.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, forwarded to every margin that reads them and to
        the copula.
    **kwargs
        Forwarded to the copula's ``sample_conditional`` (e.g. ``qrng``,
        ``seeds``, ``num_threads``).

    Returns
    -------
    array, shape (n, d), dtype float
        Samples on the original scale. The conditioning columns come back
        holding the values they were given, up to the margin's own round trip.
        With several discrete conditioners, later-drawn conditioning variables
        may land slightly outside their atom; the free-variable draws remain
        correct.

    Raises
    ------
    ValueError
        If ``y_cond`` is not two-dimensional, names more variables than there
        are, or its width matches no order tail when ``conditioning_set`` is
        ``None``; or if a named variable is out of range.
    """
    from .vinecop_base import infer_conditioning_set

    ya: Any = self._prep(y_cond)
    if ya.ndim != 2:
      raise ValueError(
        f"y_cond must be two-dimensional; got shape {tuple(ya.shape)}"
      )
    self._check_covariates(x, int(ya.shape[0]))
    k = ya.shape[1]
    if conditioning_set is None:
      # One rule for both scales. On the data scale a discrete conditioner
      # contributes no extra column, so the width *is* k.
      cond = infer_conditioning_set(
        [int(v) for v in self._vinecop.structure.order],
        ["c"] * self.dim,
        k,
      )
    else:
      cond = [int(v) for v in conditioning_set]
      if len(cond) != k:
        raise ValueError(
          f"conditioning_set names {len(cond)} variables but y_cond has {k} "
          "columns"
        )
      if any(v < 1 or v > self.dim for v in cond):
        raise ValueError(
          f"conditioning_set entries must be in 1, ..., {self.dim}; got {cond}"
        )
    # `copula_data` over the conditioning variables' margins only: one
    # column each, then the left limit of every discrete one.
    u_cond = self.copula_data([self._margins[v - 1] for v in cond], ya, x=x)
    u = declared_eval(
      self._vinecop,
      "sample_conditional",
      u_cond,
      x,
      conditioning_set=cond,
      **kwargs,
    )
    return self.marginal_icdf(u, x=x)

  # --- fitting ------------------------------------------------------------- #

  # The two part classes `from_data` fits; see the class docstring. Plain
  # comments, not `#:` ones: autosummary cannot generate a page for an
  # attribute whose value is a class, so a subclass that sets one would leave
  # the inherited entry dangling and fail the nitpicky docs build.
  vinecop_class: ClassVar[Optional[type]] = None
  margin_class: ClassVar[Optional[type]] = None

  #: Whether this lane's copula fitter accepts observation weights. Declared
  #: rather than discovered, so weights given to a lane that cannot apply them
  #: to *both* halves are refused up front instead of yielding a fit whose
  #: margins are weighted and whose copula is not. ``False`` here because the
  #: default ``_fit_copula`` cannot apply them -- it has only the part
  #: class and the caller's controls -- so a subclass declares ``True``
  #: together with the override that honors them, as
  #: :class:`~pyvinecopulib.core.Vinedist` does.
  supports_weighted_copula: bool = False

  #: Whether any margin on this array namespace can be fitted conditionally.
  #: When ``False`` the refusal is categorical, so the estimators say so rather
  #: than suggesting the caller pass margins that read covariates.
  supports_fit_covariates: bool = True

  @classmethod
  def _coerce_fit_data(
    cls,
    y: Any,
    weights: Optional[Any],
    controls: Optional[ControlsLike],
  ) -> tuple[Any, Any]:
    """Raise; override to put the fit inputs on this subclass's namespace.

    Called once before any part is fitted, by every estimator on the class.
    That is why it is a classmethod and not the instance-level ``_prep``: the
    placement a subclass fits on can come from ``controls`` (a device and a
    dtype, say), and on :meth:`from_data` there is no object yet to read it
    from.

    Parameters
    ----------
    y : object
        The caller's observations, in whatever form they passed — including a
        DataFrame.
    weights : object, or None
        The caller's observation weights, or ``None``.
    controls : ControlsLike, or None
        Fit configuration, which may carry the placement.

    Returns
    -------
    tuple
        The observations and the weights, both on this subclass's namespace so
        they cannot end up in different places.

    Raises
    ------
    NotImplementedError
        Always, unless a subclass provides the coercion.
    """
    raise NotImplementedError(
      f"{cls.__name__}._coerce_fit_data is not defined; implement it to "
      "fit this distribution from data, or compose an already-fitted copula "
      "and margins by construction."
    )

  @classmethod
  def _default_margins(
    cls,
    d: int,
    controls: Optional[ControlsLike],
    margin_controls: Optional[Sequence[Any]] = None,
  ) -> Optional[Sequence[Any]]:
    """The margin each variable gets when the caller named none.

    One ``margin_class`` per variable, which is why most subclasses declare
    that instead of overriding this. Override it when the margins need an
    argument the class attribute cannot carry -- a placement, or a family.

    The base ignores ``margin_controls``, which is right for a margin class
    that reads its controls at fit time. A subclass whose ``margin_class``
    takes its variable type or its bounds at *construction* overrides this --
    ``Vinedist`` does, since a kernel density fitted unbounded has already
    padded past the data by the time anything could tell it otherwise. That is what makes a
    bounded default reachable without naming a class.

    Parameters
    ----------
    d : int
        Number of variables.
    controls : ControlsLike, or None
        Copula fit configuration, which may carry a placement the margins
        share.
    margin_controls : sequence, or None, optional
        One marginal configuration per variable, already resolved, or ``None``
        when the caller named none.

    Returns
    -------
    sequence, or None
        One unfitted margin per variable, or ``None`` to accept
        :func:`~pyvinecopulib.margins.resolve_margins`' own default (a
        kernel-density margin throughout).
    """
    del controls
    if cls.margin_class is None:
      return None
    return [cls.margin_class() for _ in range(d)]

  @classmethod
  def _copula_controls(
    cls, controls: Optional[ControlsLike], u: Any, weights: Optional[Any]
  ) -> Any:
    """The controls the copula half is estimated with, for this lane.

    The one lane-specific step in the copula estimate, and the only hook a
    subclass normally overrides: :class:`~pyvinecopulib.core.Vinedist` writes
    ``weights`` into a copy of the controls here, and
    :class:`~pyvinecopulib.torch.TorchVinedist` pins the device and dtype the
    margins resolved. Reading it in one place is what lets ``from_data`` and
    ``fit`` share it -- one constructs the copula and the other re-estimates
    the one already held, and both must configure it identically.

    ``weights`` arrive explicitly rather than folded into ``controls``, so a
    lane whose fitter cannot use them is not silently handed them: declare
    ``supports_weighted_copula`` ``False`` instead and the estimators refuse
    the request before fitting anything.

    Parameters
    ----------
    controls : ControlsLike, or None
        What the caller passed.
    u : array, shape (n, d + k), dtype float
        The copula-scale layout, for a lane that reads its placement.
    weights : array, shape (n,), or None
        Observation weights.

    Returns
    -------
    ControlsLike, or None
        The controls to estimate with; ``controls`` unchanged by default.

    Raises
    ------
    NotImplementedError
        If weights are given and this class declares
        :attr:`supports_weighted_copula` without overriding this hook, which
        is the half-weighted fit the declaration exists to prevent.
    """
    del u
    if weights is not None:
      raise NotImplementedError(
        f"{cls.__name__} declares `supports_weighted_copula` but does not "
        "override `_copula_controls`, which is where weights reach the "
        "copula half. Override it to write them into the controls, or "
        "declare `supports_weighted_copula = False`."
      )
    return controls

  @classmethod
  def _fit_copula(
    cls,
    u: Any,
    *,
    var_types: list[str],
    controls: Optional[ControlsLike],
    structure: Optional[Any],
    weights: Optional[Any],
    x: Optional[Any] = None,
  ) -> Any:
    """Construct and fit the copula half on the pseudo-observations.

    The **construction** path, used by :meth:`from_data` when there is no
    copula yet. :meth:`fit` and :meth:`select` re-estimate the copula the
    object already holds instead -- see :meth:`_reestimate_copula` -- so a
    hosted :class:`~pyvinecopulib.core.VinecopLike` survives a refit.

    Defaults to ``vinecop_class``, so most subclasses declare that and
    override :meth:`_copula_controls` rather than this.

    Parameters
    ----------
    u : array, shape (n, d + k), dtype float
        The copula-scale layout the margins produced.
    var_types : list of str
        One ``"c"`` or ``"d"`` per variable.
    controls : ControlsLike, or None
        Fit configuration.
    structure : RVineStructure, or None
        A fixed structure, or ``None`` to select one from the data.
    weights : array, shape (n,), or None
        Observation weights.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, forwarded to ``vinecop_class.from_data`` only
        when that class declares ``supports_covariates`` -- the rule
        ``_covariates.declared_eval`` applies to a whole copula at evaluation,
        applied here at fitting. A ``Vinecop`` of ``Bicop`` pairs models none,
        so a conditional ``Vinedist`` is one whose *margins* read ``x``; what
        keeps that honest is the refusal one level up, where
        ``supports_fit_covariates`` says whether anything on the lane is fitted
        on covariates at all.

    Returns
    -------
    VinecopLike
        The fitted copula.

    Raises
    ------
    NotImplementedError
        If this class names no ``vinecop_class`` and does not override this
        hook, or if it declares :attr:`supports_weighted_copula` without
        overriding this hook to apply the weights.
    """
    if cls.vinecop_class is None:
      raise NotImplementedError(
        f"{cls.__name__} names no `vinecop_class` and does not override "
        "`_fit_copula`, so it cannot fit a copula. Set one, or compose an "
        "already-fitted copula and margins by construction."
      )
    resolved = cls._copula_controls(controls, u, weights)
    kwargs: dict[str, Any] = {"structure": structure, "controls": resolved}
    if var_types is not None:
      kwargs["var_types"] = var_types
    # Declared, not attempted: a copula class that models no covariates takes
    # no `x` argument at all, so forwarding one would raise rather than fit
    # something. Whether *anything* on the lane reads them is settled before
    # this, by `supports_fit_covariates`.
    if x is not None and getattr(
      cls.vinecop_class, "supports_covariates", False
    ):
      kwargs["x"] = x
    return cast("Any", cls.vinecop_class).from_data(u, **kwargs)

  def _reestimate_copula(
    self,
    u: Any,
    *,
    var_types: list[str],
    controls: Optional[ControlsLike],
    weights: Optional[Any],
    x: Optional[Any],
    verb: str,
  ) -> Any:
    """Re-estimate the copula this distribution already holds, in place.

    The counterpart of :meth:`_fit_copula` for :meth:`fit` and :meth:`select`,
    and the reason they mean what they say: the copula is asked to re-estimate
    *itself*, so a hosted :class:`~pyvinecopulib.core.VinecopLike` keeps its
    class and its pair-copula types across a refit -- exactly as a margin keeps
    its family. Constructing a fresh ``vinecop_class`` here would silently
    replace a caller's own vine with the default one.

    Parameters
    ----------
    u : array, shape (n, d + k), dtype float
        The copula-scale layout the margins produced.
    var_types : list of str
        One ``"c"`` or ``"d"`` per variable.
    controls : ControlsLike, or None
        Fit configuration.
    weights : array, shape (n,), or None
        Observation weights.
    x : array, shape (n, p), or None
        Exogenous covariates.
    verb : str
        ``"fit"`` to keep the copula's structure, ``"select"`` to re-select it.

    Raises
    ------
    NotImplementedError
        If the held copula has no estimator of that name -- an immutable or
        purely functional copula, which cannot be re-estimated in place.
    """
    copula = self._vinecop
    estimator = getattr(copula, verb, None)
    if estimator is None:
      raise NotImplementedError(
        f"{type(copula).__name__} has no `{verb}`, so this "
        f"{type(self).__name__} cannot re-estimate its copula half in place. "
        f"Use {type(self).__name__}.from_data to fit a fresh distribution."
      )
    resolved = type(self)._copula_controls(controls, u, weights)
    # Forwarded one at a time, so a copula whose estimator declares only what
    # it uses still works -- the idiom `BicopBase.select` follows. `var_types`
    # is left to the copula it already carries: a refit along itself cannot
    # change which variables have atoms.
    passed: dict[str, Any] = {}
    if x is not None and getattr(copula, "supports_covariates", False):
      passed["x"] = x
    estimator(u, resolved, **passed)
    return copula

  def fit(
    self,
    y: Any,
    /,
    controls: Optional[ControlsLike] = None,
    *,
    margin_controls: Optional[Any] = None,
    x: Optional[Any] = None,
    weights: Optional[Any] = None,
  ) -> Self:
    """Re-estimate both halves of this distribution, in place.

    The two-step estimator run over the parts this object already holds: each
    margin is re-estimated from its own column with the family it already has,
    and the copula's pairs along the structure it already has. That shape is
    what the call holds fixed -- the data-scale analog of
    :meth:`~pyvinecopulib.core.VinecopBase.fit`.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    controls : ControlsLike, or None, optional
        Copula fit configuration.
    margin_controls : object, or None, optional
        Marginal fit configuration, one object per variable: a single object
        broadcast to every variable, a length-``d`` sequence, or a mapping
        keyed by position -- or by name, for margins that name themselves --
        leaving the variables it does not address unconfigured. See
        :func:`pyvinecopulib.margins.resolve_margin_controls`.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, reaching each margin that declares
        ``supports_covariates``.
    weights : array, shape (n,), or None, optional
        Observation weights, applied to both halves.

    Returns
    -------
    VinedistBase
        ``self``, so the call chains.

    Notes
    -----
    Both halves are re-estimated **in place**, so a part held elsewhere shows
    the new fit and keeps its class -- including a
    :class:`~pyvinecopulib.core.VinecopLike` a caller composed the
    distribution from. The copula refits its own pair copulas along its own
    structure, so it holds the families it holds: a ``family_set`` in
    ``controls`` does not apply here, exactly as it does not to
    :meth:`~pyvinecopulib.core.Vinecop.fit`, and :meth:`select` is what
    chooses families. The margin half *refuses* the same request instead of
    ignoring it, because it can tell it was made --
    ``FitControlsMargin.family_set`` defaults to ``None`` where
    ``FitControlsVinecop.family_set`` defaults to every family, so only one of
    the two can distinguish "search these" from "the default".

    See Also
    --------
    select : Re-select the structure and the margin families as well.
    from_data : Fit a fresh distribution, choosing the margins.
    """
    return self._reestimate(
      y,
      controls,
      margin_controls,
      x=x,
      weights=weights,
      structure=self.vinecop.structure,
      verb="fit",
    )

  def select(
    self,
    y: Any,
    /,
    controls: Optional[ControlsLike] = None,
    *,
    margin_controls: Optional[Any] = None,
    x: Optional[Any] = None,
    weights: Optional[Any] = None,
  ) -> Self:
    """Re-select and re-estimate both halves, in place.

    As :meth:`fit`, except that both halves may change shape: each margin
    chooses a family from its own candidate set, and the copula re-selects its
    structure. The margin classes stay as they are -- to change those, fit a
    fresh distribution with :meth:`from_data`.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    controls : ControlsLike, or None, optional
        Copula fit configuration.
    margin_controls : object, or None, optional
        Marginal fit configuration, per variable as in :meth:`fit`; a declared
        family set is what bounds each margin's search.
    x : array, shape (n, p), or None, optional
        Exogenous covariates, reaching each margin that declares
        ``supports_covariates``.
    weights : array, shape (n,), or None, optional
        Observation weights, applied to both halves.

    Returns
    -------
    VinedistBase
        ``self``, so the call chains.

    See Also
    --------
    fit : Keep the structure and the margin families, re-estimate parameters.
    from_data : Fit a fresh distribution, choosing the margins.
    """
    return self._reestimate(
      y,
      controls,
      margin_controls,
      x=x,
      weights=weights,
      structure=None,
      verb="select",
    )

  def _reestimate(
    self,
    y: Any,
    controls: Optional[ControlsLike],
    margin_controls: Optional[Any],
    *,
    x: Optional[Any],
    weights: Optional[Any],
    structure: Optional[Any],
    verb: str,
  ) -> Self:
    """Re-estimate the held parts and rebind them, for ``fit`` and ``select``.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    controls : ControlsLike, or None
        Copula fit configuration.
    margin_controls : object, or None
        Marginal fit configuration.
    x : array, shape (n, p), or None
        Exogenous covariates.
    weights : array, shape (n,), or None
        Observation weights.
    structure : RVineStructure, or None
        The structure to fit along, or ``None`` to select one.
    verb : str
        ``"fit"`` or ``"select"``, the estimator asked of each margin.

    Raises
    ------
    ValueError
        If ``y`` is not two-dimensional, if its width does not match this
        distribution's dimension, or if weights are given that this lane
        cannot apply to the copula half.
    NotImplementedError
        If covariates are given and this lane cannot estimate both halves on
        them, or if the held copula has no estimator to re-run.
    """
    from ._resolve import resolve_margin_controls
    from ._resolve import fit_margin

    cls = type(self)
    data, weights = cls._coerce_fit_data(y, weights, controls)
    if data.ndim != 2:
      raise ValueError(f"y must be two-dimensional; got {data.ndim} dimensions")
    n, d = int(data.shape[0]), int(data.shape[1])
    if d != self.dim:
      raise ValueError(
        f"y has {d} columns but this {cls.__name__} has {self.dim} "
        "variables; "
        f"use {cls.__name__}.from_data to fit a distribution of another "
        "dimension"
      )
    x, weights = self._check_fit_inputs(x, weights, n, data, verb)

    # A mapping may be keyed by name, but the only labels this object has are
    # whatever its margins carry -- nothing here assigns one. So a name-keyed
    # mapping resolves only for margins that name themselves, and otherwise
    # `resolve_margin_controls` says to key by position.
    labels = [getattr(m, "name", None) for m in self._margins]
    named = None if any(label is None for label in labels) else labels
    per_variable = resolve_margin_controls(
      margin_controls, d, names=cast("Optional[Sequence[str]]", named)
    )
    margins = [
      fit_margin(
        margin,
        data[:, j],
        x=x,
        weights=weights,
        controls=per_variable[j],
        verb=verb,
        refit=True,
      )
      for j, margin in enumerate(self._margins)
    ]
    var_types = cls.copula_var_types(margins)
    u = cls.copula_data(margins, data, x=x)
    # The copula re-estimates *itself*, so a hosted `VinecopLike` keeps its
    # class -- `structure` is already the held one for `fit`, and `select`
    # re-selects it.
    del structure
    vinecop = self._reestimate_copula(
      u,
      var_types=var_types,
      controls=controls,
      weights=weights,
      x=x,
      verb=verb,
    )
    self._bind_dist(vinecop, margins)
    return self

  @classmethod
  def _check_fit_inputs(
    cls,
    x: Optional[Any],
    weights: Optional[Any],
    n: int,
    data: Any,
    verb: str,
  ) -> tuple[Optional[Any], Optional[Any]]:
    """Validate and place the covariates, and validate the weights.

    Returns both, because placing is half the job: the covariates a margin's
    estimator is handed have to live on the same namespace as the column it is
    handed, and returning only the weights left the placed value dying with
    this frame. The coerced observations are the reference for both: ``x`` is
    placed onto them and their first column types the weights. ``verb`` is the
    calling method's name, quoted in the error messages.

    Raises
    ------
    NotImplementedError
        If covariates are given and this lane declares it cannot estimate
        both halves on them.
    ValueError
        If weights are given that the copula half cannot apply.
    """
    if x is not None and not cls.supports_fit_covariates:
      raise NotImplementedError(
        f"{cls.__name__}.{verb} takes no covariates: no part on this array "
        f"namespace is fitted on them ({cls.__name__} declares "
        "`supports_fit_covariates = False`), so the fit would come back "
        "unconditional while the call suggested otherwise. Fit the parts "
        "yourself if only one half is conditional."
      )
    x = prepare(data, x, n)
    weights = validate_weights(weights, data[:, 0])
    if weights is not None and not cls.supports_weighted_copula:
      raise ValueError(
        f"{cls.__name__}.{verb} cannot weight the copula half, so the "
        "margins would be weighted and the copula would not -- which is not "
        "the weighted fit of anything. Drop `weights`, or fit the margins "
        "yourself and compose them with a copula you weighted."
      )
    return x, weights

  @classmethod
  def from_data(
    cls,
    y: Any,
    /,
    controls: Optional[ControlsLike] = None,
    *,
    margins: Any = None,
    margin_controls: Optional[Any] = None,
    structure: Optional[Any] = None,
    x: Optional[Any] = None,
    weights: Optional[Any] = None,
    names: Optional[Sequence[str]] = None,
  ) -> Self:
    """Fit margins and a vine copula to data, in that order.

    The two-step estimator: each margin is fitted from its own column, the
    data are transformed to the copula scale, and the copula is fitted on the
    result (Joe and Xu, 1996). Margins already fitted are left alone, so a
    fixed margin and one to estimate can be mixed freely.

    Unlike :meth:`fit` and :meth:`select`, which re-estimate the parts a
    distribution already holds, this is where the margins are *chosen*: what
    each variable gets comes from ``margins``, or from the lane's own default
    when the caller names none.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale. A DataFrame's column names are
        read as ``names`` when that is not given.
    controls : ControlsLike, or None, optional
        Copula fit configuration, in the form this subclass's fitter takes.
    margins : object, or None, optional
        What to use for each variable; see
        :func:`pyvinecopulib.margins.resolve_margins` for the accepted forms.
    margin_controls : object, or None, optional
        Marginal fit configuration, reaching every margin that is estimated
        here. The two halves are configured separately because they are fitted
        separately: ``margins`` says which class each variable gets, and this
        says how to fit or select it.
    structure : RVineStructure, or None, optional
        A fixed vine structure; selected from the data when ``None``.
    x : array, shape (n, p), or None, optional
        Exogenous covariates. Each margin that declares
        ``supports_covariates`` is fitted conditionally on them, and the copula
        is then fitted on the resulting conditional probability-integral
        transforms.
    weights : array, shape (n,), or None, optional
        Observation weights, applied to both halves. Refused when this
        subclass's copula fitter cannot apply them.
    names : sequence of str, or None, optional
        Variable names, so ``margins`` may be a mapping keyed by name.

    Returns
    -------
    VinedistBase
        The fitted distribution.

    Raises
    ------
    ValueError
        If ``y`` is not two-dimensional, or if weights are given that the
        copula half cannot apply.
    NotImplementedError
        If covariates are given and this lane cannot estimate both halves on
        them.

    See Also
    --------
    fit : Re-estimate the parts a distribution already holds.
    select : Re-estimate them, letting both halves change shape.

    References
    ----------
    .. [1] Joe, H. and Xu, J. J. (1996). *The estimation method of inference
           functions for margins for multivariate models.* Technical Report
           166, Department of Statistics, University of British Columbia.
    """
    from ._resolve import resolve_margin_controls, resolve_margins
    from ._resolve import fit_margin

    data, weights = cls._coerce_fit_data(y, weights, controls)
    if data.ndim != 2:
      raise ValueError(f"y must be two-dimensional; got {data.ndim} dimensions")
    n, d = int(data.shape[0]), int(data.shape[1])
    x, weights = cls._check_fit_inputs(x, weights, n, data, "from_data")

    if names is None:
      # A DataFrame carries its own names, and `margins` is often keyed by
      # them. Duck-typed: pandas is an extra, not a dependency of `core`.
      columns = getattr(y, "columns", None)
      if columns is not None:
        names = [str(c) for c in columns]

    per_variable = resolve_margin_controls(margin_controls, d, names=names)
    specs = resolve_margins(
      margins,
      d,
      names=names,
      default=lambda: cls._default_margins(d, controls, per_variable),
    )
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
        specs[j],
        data[:, j],
        x=x,
        weights=weights,
        controls=per_variable[j],
      )
      for j in range(d)
    ]

    # A copula needs its var_types up front, so both come from the fitted
    # margins before the copula exists; `_bind_dist` then re-derives and
    # cross-checks them.
    var_types = cls.copula_var_types(fitted)
    u = cls.copula_data(fitted, data, x=x)
    vinecop = cls._fit_copula(
      u,
      var_types=var_types,
      controls=controls,
      structure=structure,
      weights=weights,
      x=x,
    )
    return cls(vinecop, list(fitted))

  @classmethod
  def from_json(cls, json: str) -> Self:
    """Raise; override to read this subclass back from a JSON string.

    Deserialization has to name the concrete copula class it rebuilds, which
    only a subclass knows.

    Parameters
    ----------
    json : str
        A string produced by :meth:`to_json`.

    Returns
    -------
    VinedistBase
        The deserialized distribution — only when a subclass overrides this.

    Raises
    ------
    NotImplementedError
        Always, unless a subclass provides a reader.
    """
    del json
    raise NotImplementedError(
      f"{cls.__name__}.from_json is not defined: reading a distribution back "
      "requires naming the copula class to rebuild. Use "
      "pyvinecopulib.core.Vinedist for a core Vinecop, or persist this "
      "subclass the way its own library does."
    )

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
    return f"{type(self).__name__}(dim={self.dim}, margins=[{families}])"


VinedistBase.__doc__ = (VinedistBase.__doc__ or "") + _VINEDIST_EXAMPLE
