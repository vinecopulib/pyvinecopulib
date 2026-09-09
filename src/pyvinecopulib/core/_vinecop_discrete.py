"""Mixed-discrete evaluation for the vine cascades.

A discrete variable is described by two numbers, ``F(x)`` and its left limit
``F(x^-)``, and the copula quantities that are derivatives in a continuous
argument become difference quotients over the atom ``[F(x^-), F(x)]``. This
module owns that layer: the data-layout normalization vine input needs, the
per-edge variable types a vine's own types imply, and :class:`DiscreteBicop` --
a pair copula that reads a four-column argument and supplies the mixed-discrete
density and h-functions from the *continuous* ``pdf`` / ``cdf`` / ``hfunc1`` /
``hfunc2`` of the pair it wraps.

Keeping the difference quotients here rather than on the pair copulas is what
lets :class:`~pyvinecopulib.core.BicopLike` stay a purely continuous, two-column
contract: a custom pair copula gains discrete support by implementing ``cdf``,
and needs to know nothing else. A pair that handles the four-column layout
itself -- ``Bicop`` above all -- is hosted directly and never wrapped.

Internal: the vine cascades and the fit engines call these helpers, and
:class:`DiscreteBicop` is re-exported from :mod:`pyvinecopulib.core`.
"""

from __future__ import annotations

from types import ModuleType
from typing import Any, Optional, Protocol, cast

from array_api_compat import array_namespace

from ..pyvinecopulib_ext import RVineStructure
from ._rootfind import solve_increasing
from ._trim import trim

from ._covariates import pair_eval, prepare_covariates
from .bicop_base import BicopBase, flip_of
from .protocols import ArrayT, BicopLike

__all__ = ["DiscreteBicop"]


class _ContinuousBicop(Protocol[ArrayT]):
  """The unconditional evaluations :class:`DiscreteBicop` calls on its pair.

  ``cdf`` is what separates this from
  :class:`~pyvinecopulib.core.BicopLike`, where it is an optional capability:
  a difference quotient over an atom is built from it. The two inverses are
  here because a discrete argument's inverse falls back to them on the
  continuous side. ``sample`` is not, since ``DiscreteBicop`` draws through its
  own inverse Rosenblatt.

  This is also a surface ``Bicop`` satisfies structurally -- ``BicopLike`` it
  satisfies only nominally, its methods taking per-row ``parameters`` where the
  protocol takes a keyword-only ``x``. A conditioning matrix, when there is
  one, is forwarded dynamically (see ``pair_eval``), which is what makes a pair
  that cannot accept one fail loudly rather than silently.
  """

  def pdf(self, u: ArrayT) -> ArrayT: ...

  def cdf(self, u: ArrayT) -> ArrayT: ...

  def hfunc1(self, u: ArrayT) -> ArrayT: ...

  def hfunc2(self, u: ArrayT) -> ArrayT: ...

  def hinv1(self, u: ArrayT) -> ArrayT: ...

  def hinv2(self, u: ArrayT) -> ArrayT: ...


#: Atom width below which a difference quotient is numerically unstable and the
#: derivative is used instead (``AbstractBicop``'s threshold).
DELTA_MIN: float = 5e-5
#: The unit box the reference pair copula trims its arguments to.


def check_var_types(var_types: Optional[list[str]], d: int) -> tuple[str, ...]:
  """Normalize and validate a vine's per-variable types.

  Parameters
  ----------
  var_types : list of str, or None, optional
      Per-variable types, ``"c"`` or ``"d"``; ``None`` means all continuous.
  d : int
      Dimension the types must cover.

  Returns
  -------
  tuple of str
      The validated types, one per variable.

  Raises
  ------
  ValueError
      If the length is not ``d`` or an entry is outside ``{"c", "d"}``.
  """
  types = ("c",) * d if var_types is None else tuple(var_types)
  if len(types) != d:
    raise ValueError(f"var_types has {len(types)} entries, expected {d}")
  bad = [t for t in types if t not in ("c", "d")]
  if bad:
    raise ValueError(f"var_types entries must be 'c' or 'd'; got {bad[0]!r}")
  return types


def disc_cols(var_types: tuple[str, ...]) -> tuple[int, ...]:
  """Offsets of the left-limit columns within the compact layout's second block.

  Parameters
  ----------
  var_types : tuple of str
      Per-variable types.

  Returns
  -------
  tuple of int
      Variable ``i``'s left limit sits at column ``d + result[i]``; the entry is
      meaningless (``0``) at a continuous variable, which has no such column.
  """
  offsets, seen = [0] * len(var_types), 0
  for i, t in enumerate(var_types):
    if t == "d":
      offsets[i] = seen
      seen += 1
  return tuple(offsets)


def collapse_data(
  u: ArrayT,
  d: int,
  var_types: tuple[str, ...],
  name: str,
  *,
  values_only: bool = False,
) -> ArrayT:
  """Validate ``u``'s column layout and reduce it to the columns needed.

  Accepts the layouts ``Vinecop`` accepts: ``(n, d)`` for an all-continuous
  model, and the compact ``(n, d + k)`` or expanded ``(n, 2d)`` form when ``k``
  of the ``d`` variables are discrete. The plain ``(n, d)`` is rejected on a
  discrete model, because silently reusing each value as its own left limit
  would evaluate a continuous density.

  Parameters
  ----------
  u : array, shape (n, d), (n, d + k) or (n, 2d), dtype float
      The matrix to validate.
  d : int
      Dimension of the model.
  var_types : tuple of str
      Per-variable types.
  name : str
      Calling-method name, used only in the error message.
  values_only : bool, default=False
      Return just the ``d`` value columns, and accept a plain ``(n, d)`` input
      even on a discrete model. Set by the callers that never read a left limit.

  Returns
  -------
  array, shape (n, d + k) or (n, d), dtype float
      ``u`` in the compact layout, or its value block when ``values_only``.

  Raises
  ------
  ValueError
      If ``u`` is not 2-d or its column count matches no accepted layout.
  """
  ua: Any = u
  k = var_types.count("d")
  accepted = {d + k, 2 * d} | ({d} if values_only else set())
  if ua.ndim != 2 or int(ua.shape[1]) not in accepted:
    shapes = ", ".join(f"(n, {c})" for c in sorted(accepted))
    raise ValueError(
      f"{name}: u must have shape {shapes} for a vine with var_types="
      f"{list(var_types)}; got {tuple(ua.shape)}"
    )
  if values_only or k == 0:
    return cast("ArrayT", ua[:, :d])
  if int(ua.shape[1]) == d + k:
    return u
  # Expanded (n, 2d) -> compact (n, d + k): keep the left-limit columns of the
  # discrete variables only, in variable order.
  xp = array_namespace(ua)
  keep = [d + i for i, t in enumerate(var_types) if t == "d"]
  return cast("ArrayT", xp.concat([ua[:, :d], ua[:, keep]], axis=1))


def pair_var_types(
  structure: RVineStructure, var_types: tuple[str, ...]
) -> tuple[tuple[tuple[str, str], ...], ...]:
  """Per-edge variable types implied by a structure and its variable types.

  Port of ``Vinecop::set_var_types_internal``: tree 0 reads the natural-order
  variable types, and every later tree inherits from the pair that produced its
  input column -- ``hfunc2`` carries the first variable's type, ``hfunc1`` the
  second's. The types are therefore a function of the structure and never have
  to be stored on a pair copula.

  Parameters
  ----------
  structure : RVineStructure
      The vine structure.
  var_types : tuple of str
      Per-variable types, in variable order.

  Returns
  -------
  tuple of tuple of tuple of str
      ``result[tree][edge]`` is that edge's ``(type1, type2)``.
  """
  d = int(structure.dim)
  trunc_lvl = int(structure.trunc_lvl)
  order = [int(v) for v in structure.order]
  natural = [var_types[v - 1] for v in order]
  table: list[tuple[tuple[str, str], ...]] = []
  for tree in range(trunc_lvl):
    row: list[tuple[str, str]] = []
    for edge in range(d - tree - 1):
      m = int(structure.min_array(tree, edge))
      on_diagonal = m == int(structure.struct_array(tree, edge, True))
      if tree == 0:
        row.append((natural[edge], natural[m - 1]))
      else:
        prev = table[tree - 1]
        row.append((prev[edge][0], prev[m - 1][0 if on_diagonal else 1]))
    table.append(tuple(row))
  return tuple(table)


def seed_left_limits(
  u: ArrayT,
  d: int,
  order: tuple[int, ...],
  var_types: tuple[str, ...],
  offsets: tuple[int, ...],
  xp: ModuleType,
) -> Optional[ArrayT]:
  """Natural-order left limits read off the compact layout, or ``None``.

  ``None`` for an all-continuous model, which is what switches the whole
  left-limit cascade off. A continuous variable's column holds its own value: a
  pair only ever reads the left-limit column of a variable it declares discrete,
  and this keeps the four-column edge input well defined anyway.

  Parameters
  ----------
  u : array, shape (n, d + k), dtype float
      Prepared observations in the compact layout.
  d : int
      Dimension of the model.
  order : tuple of int
      The structure's variable order (1-based).
  var_types : tuple of str
      Per-variable types.
  offsets : tuple of int
      Left-limit column offsets, from :func:`disc_cols`.
  xp : module
      The array namespace of ``u``.

  Returns
  -------
  array, shape (n, d), dtype float, or None
      Left limits in natural order, or ``None`` when nothing is discrete.
  """
  if "d" not in var_types:
    return None
  ua: Any = u
  sub = xp.empty((ua.shape[0], d), dtype=ua.dtype, device=ua.device)
  for j in range(d):
    v = order[j] - 1
    sub[:, j] = ua[:, d + offsets[v] if var_types[v] == "d" else v]
  return cast("ArrayT", sub)


def edge_columns(
  structure: RVineStructure,
  pair_types: Optional[tuple[tuple[tuple[str, str], ...], ...]],
  tree: int,
  edge: int,
  hfunc1: ArrayT,
  hfunc2: ArrayT,
  hfunc1_sub: Optional[ArrayT],
  hfunc2_sub: Optional[ArrayT],
) -> tuple[ArrayT, ArrayT, Optional[tuple[ArrayT, ArrayT]], tuple[str, str]]:
  """Resolve one edge's pair-copula input columns and its variable types.

  ``m`` is the min-array entry: the natural-order index of the column finalized
  in a previous tree. The second pair input comes from ``hfunc2`` when ``m`` sits
  on the natural-order diagonal, else from ``hfunc1``
  (``class.ipp:1026-1034``). The left-limit pair is returned only when the edge
  has a discrete variable, and mirrors ``Bicop::format_data``: a continuous
  variable's left limit is its own value.

  Parameters
  ----------
  structure : RVineStructure
      The vine structure being walked.
  pair_types : tuple of tuple of tuple of str, or None, optional
      Per-edge types from :func:`pair_var_types`; ``None`` when all continuous.
  tree : int
      Tree index (``0``-based).
  edge : int
      Edge index within the tree (``0``-based).
  hfunc1, hfunc2 : array, shape (n, d), dtype float
      The h-function scratch matrices.
  hfunc1_sub, hfunc2_sub : array, shape (n, d), dtype float, or None, optional
      The left-limit scratch matrices; ``None`` when all continuous.

  Returns
  -------
  col0, col1 : array, shape (n,), dtype float
      The pair's two value inputs.
  subs : tuple of array, or None
      Their left limits, or ``None`` when the edge is fully continuous.
  types : tuple of str
      The edge's ``(type1, type2)``.
  """
  h1: Any = hfunc1
  h2: Any = hfunc2
  h2_sub: Any = hfunc2_sub
  m = int(structure.min_array(tree, edge))
  on_diagonal = m == int(structure.struct_array(tree, edge, True))
  col0 = cast("ArrayT", h2[:, edge])
  col1 = cast("ArrayT", h2[:, m - 1] if on_diagonal else h1[:, m - 1])
  types = ("c", "c") if pair_types is None else pair_types[tree][edge]
  if hfunc2_sub is None or "d" not in types:
    return col0, col1, None, types
  sub0 = cast("ArrayT", h2_sub[:, edge]) if types[0] == "d" else col0
  if types[1] != "d":
    sub1 = col1
  elif on_diagonal:
    sub1 = cast("ArrayT", h2_sub[:, m - 1])
  else:
    sub1 = cast("ArrayT", cast("Any", hfunc1_sub)[:, m - 1])
  return col0, col1, (sub0, sub1), types


def stack_edge(
  xp: ModuleType,
  col0: ArrayT,
  col1: ArrayT,
  subs: Optional[tuple[ArrayT, ArrayT]],
) -> ArrayT:
  """Assemble a pair-copula argument from its value columns and left limits.

  Parameters
  ----------
  xp : module
      The array namespace to build on.
  col0, col1 : array, shape (n,), dtype float
      The pair's two value inputs.
  subs : tuple of array, or None, optional
      Their left limits, or ``None`` for a fully continuous edge.

  Returns
  -------
  array, shape (n, 2) or (n, 4), dtype float
      ``[u1, u2]``, or ``[u1, u2, u1^-, u2^-]`` when left limits are given.
  """
  cols = [col0, col1] if subs is None else [col0, col1, *subs]
  return cast("ArrayT", xp.stack(cols, axis=-1))


def with_left_limit(u_e: ArrayT, arg: int) -> ArrayT:
  """A four-column edge input with one argument replaced by its left limit.

  The input the cascade needs for a left-limit h-function: conditioning on (or
  evaluating at) the lower end of a discrete argument's atom.

  Parameters
  ----------
  u_e : array, shape (n, 4), dtype float
      An edge input ``[u1, u2, u1^-, u2^-]``.
  arg : int
      Which argument to replace, ``0`` or ``1``.

  Returns
  -------
  array, shape (n, 4), dtype float
      ``u_e`` with column ``arg`` replaced by column ``2 + arg``.
  """
  edge: Any = u_e
  xp = array_namespace(edge)
  cols = [edge[:, 0], edge[:, 1], edge[:, 2], edge[:, 3]]
  cols[arg] = edge[:, 2 + arg]
  return cast("ArrayT", xp.stack(cols, axis=-1))


def continuous_view(pair: object) -> _ContinuousBicop[ArrayT]:
  """Return ``pair`` evaluated as a continuous copula, when it can be.

  A pair copula may carry its own variable types -- ``Bicop`` does -- in which
  case its two-column evaluation surface is unavailable. Such a pair is used
  through its optional ``as_continuous`` capability; one that does not advertise
  it is already continuous-only and is returned as is.

  Parameters
  ----------
  pair : BicopLike
      The pair copula to view.

  Returns
  -------
  BicopLike
      The continuous view, or ``pair`` itself.

  Raises
  ------
  ValueError
      If ``pair`` declares discrete variables but offers no continuous view.
  """
  view = getattr(pair, "as_continuous", None)
  if view is not None:
    return cast("_ContinuousBicop[ArrayT]", view())
  types = getattr(pair, "var_types", None)
  if types is not None and any(t != "c" for t in types):
    raise ValueError(
      f"{type(pair).__name__} declares var_types={list(types)} but has no "
      "as_continuous(); the inverse Rosenblatt cascade evaluates every pair "
      "copula as continuous, since it produces the values a left limit would "
      "be taken of. Host the continuous copula in a DiscreteBicop instead, or "
      "add as_continuous()."
    )
  return cast("_ContinuousBicop[ArrayT]", pair)


class DiscreteBicop(BicopBase[ArrayT]):
  """A continuous pair copula evaluated as a mixed-discrete one.

  Reads the four-column layout ``[F(u1), F(u2), F(u1^-), F(u2^-)]`` and returns
  the mixed-discrete density and h-functions, built from the wrapped pair's
  continuous ``pdf`` / ``cdf`` / ``hfunc1`` / ``hfunc2``: in a discrete argument
  a derivative becomes a difference quotient over the atom, and the derivative
  itself is used where the atom is narrower than ``5e-5``. A continuous
  argument's left-limit column is ignored -- it equals its value.

  This is how a pair copula that knows nothing about discreteness is hosted on a
  discrete edge of a :class:`~pyvinecopulib.core.VinecopBase`: wrap it in
  ``get_pair_copula`` with the types the vine derives for that slot,
  :meth:`~pyvinecopulib.core.VinecopBase.pair_var_types`. The h-functions of a
  discrete argument, and the density of a pair with two, are built from the
  copula's distribution function, so a pair hosted on a discrete edge must
  implement ``cdf`` -- the cascades evaluate h-functions at every edge.

  Parameters
  ----------
  pair : BicopLike
      The pair copula to wrap; only its two-column ``pdf`` / ``cdf`` / ``hfunc1``
      / ``hfunc2`` are used, through ``as_continuous()`` when it advertises one.
  var_types : tuple of str
      The edge's two variable types, ``"c"`` or ``"d"``.

  See Also
  --------
  pyvinecopulib.core.VinecopBase.pair_var_types : The types to wrap a slot with.
  pyvinecopulib.core.Bicop : The reference pair copula, discrete-aware itself.

  Notes
  -----
  ``examples/10_extending_pyvinecopulib.ipynb`` hosts a custom pair copula on a
  discrete edge and checks it against ``Vinecop``.

  The output matches ``Bicop``'s own discrete surface bit-for-bit, for every
  family including the nonparametric ``tll``, and satisfies
  ``sum_atoms c * (u1 - u1^-) == 1`` exactly. That identity is what pins the
  quotients: they telescope over the atoms, so no tolerance argument is needed
  to say which of two candidate surfaces is right.
  """

  #: The batched grid surface has no distribution-function lookup, which the
  #: discrete h-functions need, so a wrapped pair never takes that fast path.
  supports_batched: bool = False

  def __init__(
    self, pair: _ContinuousBicop[ArrayT], var_types: tuple[str, str]
  ) -> None:
    self._pair = continuous_view(pair)
    self.var_types = list(var_types)

  @property
  def var_types(self) -> list[str]:
    """The two variable types, ``"c"`` (continuous) or ``"d"`` (discrete).

    Returns
    -------
    list of str
        The wrapped edge's ``[type1, type2]``.
    """
    return list(self._var_types)

  @var_types.setter
  def var_types(self, value: list[str]) -> None:
    types = check_var_types(list(value), 2)
    self._var_types = types
    self._d1 = types[0] == "d"
    self._d2 = types[1] == "d"

  def as_continuous(self) -> BicopLike[ArrayT]:
    """The wrapped pair copula, evaluated without any left limits.

    Returns
    -------
    BicopLike
        The continuous pair copula this wraps.
    """
    return cast("BicopLike[ArrayT]", self._pair)

  def flip(self) -> "DiscreteBicop[ArrayT]":
    """The argument-swapped copula, with its variable types swapped too.

    Returns
    -------
    DiscreteBicop
        A wrapper around the flipped pair, reading ``[u2, u1, u2^-, u1^-]``.
    """
    flipped = flip_of(self._pair)
    return DiscreteBicop(flipped, (self._var_types[1], self._var_types[0]))

  def sample(
    self,
    n: int,
    *,
    x: Optional[ArrayT] = None,
    qrng: bool = False,
    seeds: Optional[list[int]] = None,
  ) -> ArrayT:
    """Draw from the wrapped continuous copula.

    Discreteness changes how a copula density is evaluated at atoms, not the
    continuous latent uniforms the copula samples. The wrapped pair therefore
    owns the RNG and inverse Rosenblatt transform.

    Parameters
    ----------
    n : int
        Number of samples.
    x : array, shape (n, p), or None, optional
        External covariates for a conditional draw.
    qrng : bool, default=False
        Draw quasi-random base uniforms instead of pseudo-random ones.
    seeds : list of int, or None, optional
        RNG seeds.

    Returns
    -------
    array, shape (n, 2), dtype float
        Samples in the unit square.
    """
    x = prepare_covariates(self, x, n)
    method = cast("BicopLike[ArrayT]", self._pair).sample
    draw_seeds = list(seeds) if seeds else []
    if x is None:
      return method(n, qrng=qrng, seeds=draw_seeds)
    return method(n, x=x, qrng=qrng, seeds=draw_seeds)

  def __repr__(self) -> str:
    return f"DiscreteBicop({self._pair!r}, var_types={list(self._var_types)})"

  # --- argument handling ------------------------------------------------ #
  # The helpers below compute on the columns rather than forwarding them --
  # differences, quotients, comparisons -- so they hold their arguments as
  # `Any`, which is what `protocols.py` says an unbounded `ArrayT` requires of
  # a body that does arithmetic. Every one is private to this class; the
  # module-level functions above, which only forward, are typed `ArrayT`.
  def _split(self, u: ArrayT) -> tuple[ModuleType, Any, Any, Any, Any]:
    """Namespace plus the two values and their left limits."""
    ua: Any = u
    expected = 4 if (self._d1 or self._d2) else 2
    if ua.ndim != 2 or int(ua.shape[1]) != expected:
      raise ValueError(
        f"u must have shape (n, {expected}) for var_types="
        f"{list(self._var_types)}; got {tuple(ua.shape)}"
      )
    xp = array_namespace(ua)
    # Trimmed before anything is subtracted, as ``Bicop::prep_for_abstract``
    # does: an h-function feeding the next tree may land exactly on 0 or 1, and
    # the atom width in the denominator has to be the width of the trimmed atom.
    ut: Any = trim(xp, u)
    u1, u2 = ut[:, 0], ut[:, 1]
    # A continuous argument's left limit is its own value
    # (``Bicop::format_data``), so the cascade's column for it is never read.
    return (
      xp,
      u1,
      u2,
      ut[:, 2] if self._d1 else u1,
      ut[:, 3] if self._d2 else u2,
    )

  def _pdf(self, xp: ModuleType, a: Any, b: Any, x: Optional[ArrayT]) -> Any:
    return pair_eval(self._pair.pdf, xp.stack([a, b], axis=-1), x)

  def _cdf(self, xp: ModuleType, a: Any, b: Any, x: Optional[ArrayT]) -> Any:
    return pair_eval(self._pair.cdf, xp.stack([a, b], axis=-1), x)

  def _h1(self, xp: ModuleType, a: Any, b: Any, x: Optional[ArrayT]) -> Any:
    return pair_eval(self._pair.hfunc1, xp.stack([a, b], axis=-1), x)

  def _h2(self, xp: ModuleType, a: Any, b: Any, x: Optional[ArrayT]) -> Any:
    return pair_eval(self._pair.hfunc2, xp.stack([a, b], axis=-1), x)

  @staticmethod
  def _take(value: Optional[Any], mask: Any) -> Optional[Any]:
    """Select rows from an optional conditioning matrix."""
    return None if value is None else value[mask]

  @staticmethod
  def _quotient(xp: ModuleType, num: Any, delta: Any, fallback: Any) -> Any:
    """``|num / delta|`` over a wide-enough atom, else ``|fallback|``."""
    wide = delta > DELTA_MIN
    safe = xp.where(wide, delta, xp.ones_like(delta))
    return xp.abs(xp.where(wide, num / safe, fallback))

  def _pdf_mixed(
    self,
    xp: ModuleType,
    u1: Any,
    u2: Any,
    u1m: Any,
    u2m: Any,
    x: Optional[ArrayT],
    *,
    discrete: int,
  ) -> Any:
    """Evaluate only the quotient or derivative each row requires."""
    delta = xp.abs((u1 - u1m) if discrete == 1 else (u2 - u2m))
    wide = delta > DELTA_MIN
    out = xp.empty_like(delta)
    if bool(xp.any(wide)):
      x_wide = self._take(x, wide)
      if discrete == 1:
        num = self._h2(xp, u1[wide], u2[wide], x_wide) - self._h2(
          xp, u1m[wide], u2m[wide], x_wide
        )
      else:
        num = self._h1(xp, u1[wide], u2[wide], x_wide) - self._h1(
          xp, u1m[wide], u2m[wide], x_wide
        )
      out[wide] = num / delta[wide]
    narrow = ~wide
    if bool(xp.any(narrow)):
      out[narrow] = self._pdf(
        xp,
        0.5 * (u1[narrow] + u1m[narrow]),
        0.5 * (u2[narrow] + u2m[narrow]),
        self._take(x, narrow),
      )
    return xp.abs(out)

  def _rect(
    self,
    xp: ModuleType,
    a1: Any,
    b1: Any,
    a2: Any,
    b2: Any,
    x: Optional[ArrayT],
  ) -> Any:
    """``P((a1, b1] x (a2, b2])`` as the four-corner difference.

    A pair that can compute the rectangle without the cancellation this carries
    -- ``TorchTllBicop.rect_mass`` does -- would be more accurate here, by 7.6x at
    a `1/8`-wide atom and far more at the widths the inner trees reach. It is
    not used: the density divides by the atom's area, and the
    discrete cascade then amplifies a 1e-15 pair-level difference to 8.5e-8 at
    the vine, which is a visible divergence from ``Vinecop``. Parity between
    ``TorchVinecop``'s cascade and ``Vinecop``'s is a documented guarantee, so
    this route stays the reference's, exactly. Taking the rectangle upstream is
    what would let both sides use it (vinecopulib#757).
    """
    # Summed in two pairs, as the compiled pair copula sums them: the grouping
    # is what makes the two agree to the last bit rather than to rounding.
    return (self._cdf(xp, b1, b2, x) + self._cdf(xp, a1, a2, x)) - (
      self._cdf(xp, a1, b2, x) + self._cdf(xp, b1, a2, x)
    )

  def _strip(
    self,
    xp: ModuleType,
    a1: Any,
    b1: Any,
    b2: Any,
    x: Optional[ArrayT],
    axis: int,
  ) -> Any:
    """``P((a1, b1] x (0, b2])`` for ``axis=1``, transposed for ``axis=2``.

    The rectangle anchored at the origin, which an h-function's numerator is.
    Its second pair of corners vanishes, so this is a two-term difference rather
    than evaluating ``cdf`` at zero. See :meth:`_rect` on why the exact route is
    not taken.
    """
    if axis == 1:
      return self._cdf(xp, b1, b2, x) - self._cdf(xp, a1, b2, x)
    return self._cdf(xp, b2, b1, x) - self._cdf(xp, b2, a1, x)

  # --- the mixed-discrete surface --------------------------------------- #
  def pdf(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Density with respect to each argument's own reference measure.

    A continuous argument contributes a derivative, a discrete one the
    probability of its atom, so a mixed pair gives a difference quotient and two
    discrete arguments the rectangle probability, each divided by the atom
    widths.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        ``[u1, u2]`` when both arguments are continuous, else
        ``[u1, u2, u1^-, u2^-]``.
    x : array, shape (n, p), or None, optional
        Conditioning variables, forwarded to the wrapped pair when given.

    Returns
    -------
    array, shape (n,), dtype float
        Density values.
    """
    xp, u1, u2, u1m, u2m = self._split(u)
    if self._d1 and self._d2:
      return cast("ArrayT", self._pdf_d_d(xp, u1, u2, u1m, u2m, x))
    if self._d1:
      return cast(
        "ArrayT", self._pdf_mixed(xp, u1, u2, u1m, u2m, x, discrete=1)
      )
    if self._d2:
      return cast(
        "ArrayT", self._pdf_mixed(xp, u1, u2, u1m, u2m, x, discrete=2)
      )
    return cast("ArrayT", self._pdf(xp, u1, u2, x))

  def _pdf_d_d(
    self,
    xp: ModuleType,
    u1: Any,
    u2: Any,
    u1m: Any,
    u2m: Any,
    x: Optional[ArrayT],
  ) -> Any:
    """Rectangle probability per unit area, with the degenerate fallbacks."""
    d1, d2 = xp.abs(u1 - u1m), xp.abs(u2 - u2m)
    m1, m2 = 0.5 * (u1 + u1m), 0.5 * (u2 + u2m)
    narrow1, narrow2 = d1 < DELTA_MIN, d2 < DELTA_MIN
    both = xp.where(d1 > d2, d1, d2) < DELTA_MIN
    only1 = narrow1 & ~both
    only2 = narrow2 & ~both
    wide = ~(both | only1 | only2)
    out = xp.empty_like(d1)
    if bool(xp.any(wide)):
      out[wide] = self._rect(
        xp,
        u1m[wide],
        u1[wide],
        u2m[wide],
        u2[wide],
        self._take(x, wide),
      ) / (d1[wide] * d2[wide])
    if bool(xp.any(only1)):
      x_only1 = self._take(x, only1)
      # A collapsed argument is held at the atom's midpoint in both terms.
      out[only1] = (
        self._h1(xp, m1[only1], u2[only1], x_only1)
        - self._h1(xp, m1[only1], u2m[only1], x_only1)
      ) / d2[only1]
    if bool(xp.any(only2)):
      x_only2 = self._take(x, only2)
      out[only2] = (
        self._h2(xp, u1[only2], m2[only2], x_only2)
        - self._h2(xp, u1m[only2], m2[only2], x_only2)
      ) / d1[only2]
    if bool(xp.any(both)):
      out[both] = self._pdf(xp, m1[both], m2[both], self._take(x, both))
    return xp.abs(out)

  def hfunc1(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """``P(U2 <= u2 | U1)``, conditioning on the atom when ``U1`` is discrete.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        See :meth:`pdf`.
    x : array, shape (n, p), or None, optional
        Conditioning variables, forwarded to the wrapped pair when given.

    Returns
    -------
    array, shape (n,), dtype float
        Conditional distribution values.
    """
    xp, u1, u2, u1m, _ = self._split(u)
    if not self._d1:
      return cast("ArrayT", self._h1(xp, u1, u2, x))
    # Conditioning on `u1^- < U1 <= u1` divides the rectangle probability by the
    # atom's width; the second argument enters at its value either way.
    return cast(
      "ArrayT",
      self._quotient(
        xp,
        self._strip(xp, u1m, u1, u2, x, axis=1),
        xp.abs(u1 - u1m),
        self._h1(xp, 0.5 * (u1 + u1m), u2, x),
      ),
    )

  def hfunc2(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """``P(U1 <= u1 | U2)``, conditioning on the atom when ``U2`` is discrete.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        See :meth:`pdf`.
    x : array, shape (n, p), or None, optional
        Conditioning variables, forwarded to the wrapped pair when given.

    Returns
    -------
    array, shape (n,), dtype float
        Conditional distribution values.
    """
    xp, u1, u2, _, u2m = self._split(u)
    if not self._d2:
      return cast("ArrayT", self._h2(xp, u1, u2, x))
    return cast(
      "ArrayT",
      self._quotient(
        xp,
        self._strip(xp, u2m, u2, u1, x, axis=2),
        xp.abs(u2 - u2m),
        self._h2(xp, u1, 0.5 * (u2 + u2m), x),
      ),
    )

  def cdf(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Distribution function, which the left limits do not enter.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        See :meth:`pdf`.
    x : array, shape (n, p), or None, optional
        Conditioning variables, forwarded to the wrapped pair when given.

    Returns
    -------
    array, shape (n,), dtype float
        Distribution values.
    """
    xp, u1, u2, _, _ = self._split(u)
    return cast("ArrayT", self._cdf(xp, u1, u2, x))

  def hinv1(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Inverse of :meth:`hfunc1` in its second argument.

    Analytic through the wrapped pair when the conditioning argument is
    continuous, and a monotone bisection of the mixed-discrete ``hfunc1``
    otherwise.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        Column ``0`` is the conditioning value and column ``1`` the level; see
        :meth:`pdf` for the layout.
    x : array, shape (n, p), or None, optional
        Conditioning variables, forwarded to the wrapped pair when given.

    Returns
    -------
    array, shape (n,), dtype float
        The inverse values.
    """
    xp, u1, p, u1m, _ = self._split(u)
    if not self._d1:
      return cast(
        "ArrayT", pair_eval(self._pair.hinv1, xp.stack([u1, p], axis=-1), x)
      )
    return cast(
      "ArrayT",
      solve_increasing(
        lambda v: self.hfunc1(
          cast("ArrayT", xp.stack([u1, v, u1m, v], axis=-1)), x=x
        ),
        p,
      ),
    )

  def hinv2(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Inverse of :meth:`hfunc2` in its first argument.

    Analytic through the wrapped pair when the conditioning argument is
    continuous, and a monotone bisection of the mixed-discrete ``hfunc2``
    otherwise.

    Parameters
    ----------
    u : array, shape (n, 2) or (n, 4), dtype float
        Column ``0`` is the level and column ``1`` the conditioning value; see
        :meth:`pdf` for the layout.
    x : array, shape (n, p), or None, optional
        Conditioning variables, forwarded to the wrapped pair when given.

    Returns
    -------
    array, shape (n,), dtype float
        The inverse values.
    """
    xp, p, u2, _, u2m = self._split(u)
    if not self._d2:
      return cast(
        "ArrayT", pair_eval(self._pair.hinv2, xp.stack([p, u2], axis=-1), x)
      )
    return cast(
      "ArrayT",
      solve_increasing(
        lambda v: self.hfunc2(
          cast("ArrayT", xp.stack([v, u2, v, u2m], axis=-1)), x=x
        ),
        p,
      ),
    )
