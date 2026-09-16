"""Array-agnostic contracts for the four modeling layers.

:class:`BicopLike`, :class:`VinecopLike`, :class:`MarginLike` and
:class:`VinedistLike` define what a pair copula, a vine, a univariate margin
and a data-scale vine distribution must provide, independent of the array
namespace (NumPy or PyTorch). They are the extension point for custom
implementations: implement the protocol — or, far more easily, subclass the
canonical base for that layer (:class:`~pyvinecopulib.core.BicopBase`,
:class:`~pyvinecopulib.core.VinecopBase`,
:class:`~pyvinecopulib.core.MarginBase`,
:class:`~pyvinecopulib.core.VinedistBase`), which fills in most or all of it —
and the object plugs into the rest of the library: it can be hosted in a vine,
composed into a distribution, or fitted by the sklearn estimators. The
reference implementations are :class:`~pyvinecopulib.core.Bicop`,
:class:`~pyvinecopulib.core.Vinecop`, :class:`~pyvinecopulib.core.Kde1d` and
:class:`~pyvinecopulib.core.Vinedist`, and their PyTorch counterparts in
:mod:`pyvinecopulib.torch`.

:class:`ControlsLike` is the odd one out: it describes *fit configuration*
rather than a model, and asks for a single ``to_dict``.

**Conditioning is not here.** These contracts describe the *unconditional*
surface — what every implementation serves, the reference ones included. A model
that depends on conditioning variables (conditioning-set values and/or external
covariates, row-aligned with the data) takes them as a keyword-only ``x``, which
is a **widening** of the contract rather than part of it, and which the four
canonical bases declare: subclass :class:`~pyvinecopulib.core.BicopBase` or
:class:`~pyvinecopulib.core.MarginBase` and add ``x`` to the members that read
one. In a vine, each pair's ``x`` is assembled per edge by a
:class:`~pyvinecopulib.core.ConditioningContext`.

Declaring ``x`` here instead would oblige every implementation to accept a
parameter most of them model nothing by, and would make the contract
unsatisfiable for :class:`~pyvinecopulib.core.Bicop`,
:class:`~pyvinecopulib.core.Vinecop` and :class:`~pyvinecopulib.core.Kde1d`,
which take none — the asymmetry that a protocol may ask for *less* than an
implementation provides, never more.

**Scope.** These describe the *evaluation* surface — enough to host a pair
copula in a vine, to consume a fitted vine, or to compose a distribution — not
the whole core API. :class:`~pyvinecopulib.core.Bicop` and
:class:`~pyvinecopulib.core.Vinecop` carry considerably more (score and
derivative families, per-row parameters, serialization, discrete data layouts)
that a custom implementation is not expected to provide. A performance knob is
likewise not part of a contract: ``num_threads`` lives on the classes that mean
something by it and not here, since a protocol parameter would oblige every
implementation to accept one. Accepting extra keyword arguments is a widening,
so a class that takes more than a protocol asks still satisfies it.

Every protocol here is ``runtime_checkable``, which compares member *names*
only, so ``isinstance(cop, BicopLike)`` reports that the names are present and
nothing about their signatures. The signatures nonetheless agree: ``Bicop``,
``Vinecop`` and ``Kde1d`` satisfy their contracts *structurally*, which ``ty``
checks and ``tests/test_core_protocols.py`` pins, so a consumer may type against
the contract instead of a concrete class.

**Typing.** :data:`ArrayT` is the array type an implementation commits to --
``numpy.ndarray`` or ``torch.Tensor`` -- carried on every signature here and on
everything that implements them, so a concrete class returns its own array
type: :class:`~pyvinecopulib.torch.TorchTllBicop` returns a ``torch.Tensor``,
not something vaguer.

It is bounded by ``Array``, which says what any array must provide: a
shape, indexing, arithmetic and comparison. Write ``ArrayT`` in a signature and
``Array`` where a value is only ever computed on.

``Namespace`` is the same idea for the module that operates on arrays --
conventionally ``xp`` -- and ``array_namespace`` resolves an array to it.
"""

from __future__ import annotations

from abc import abstractmethod
from collections.abc import Sequence
from typing import (
  Any,
  Protocol,
  Self,
  TypeVar,
  cast,
  runtime_checkable,
)

from array_api_compat import array_namespace as _array_namespace

from ..pyvinecopulib_ext import RVineStructure

# The right-hand side of an array operator, or an index: an array, a scalar, a
# slice, a mask, or a tuple of those. `Any` once here rather than at nineteen
# sites, and not a shrug: a protocol parameter is checked bivariantly, and the
# two reference array types accept unions that differ from each other, so any
# narrower spelling would exclude one of them.
_Operand = Any


@runtime_checkable
class Array(Protocol):
  """What an array provides: a shape, indexing, arithmetic and comparison.

  ``numpy.ndarray`` and ``torch.Tensor`` both satisfy it as they are -- there
  is nothing to inherit from and nothing to register. It is the bound on
  :data:`ArrayT`, so writing ``ArrayT`` in a signature already promises all of
  this, and ``Array`` itself is what to write where a value is only ever
  computed on rather than handed back.

  ``dtype`` and ``device`` are opaque: pass either back to the
  ``Namespace`` that produced the array rather than reading into it.
  """

  # `__array_namespace__`, `__eq__`, `T` / `mT`, `to` / `astype` and `__neg__`
  # are all left out: one of the two reference array types cannot satisfy each.
  # Use `xp.equal`, `xp.matrix_transpose` and `a * -1.0` instead.

  @property
  def shape(self) -> tuple[int, ...]: ...

  @property
  def ndim(self) -> int: ...

  @property
  def dtype(self) -> object: ...

  @property
  def device(self) -> object: ...

  def __getitem__(self, key: _Operand, /) -> Self: ...

  def __setitem__(self, key: _Operand, value: _Operand, /) -> None: ...

  def __add__(self, other: _Operand, /) -> Self: ...
  def __radd__(self, other: _Operand, /) -> Self: ...
  def __sub__(self, other: _Operand, /) -> Self: ...
  def __rsub__(self, other: _Operand, /) -> Self: ...
  def __mul__(self, other: _Operand, /) -> Self: ...
  def __rmul__(self, other: _Operand, /) -> Self: ...
  def __truediv__(self, other: _Operand, /) -> Self: ...
  def __rtruediv__(self, other: _Operand, /) -> Self: ...

  # A comparison answers in booleans, so it is a `BoolArray` rather than
  # `Self`: a boolean array is a different type from the float array compared,
  # and could not satisfy `Self` even in principle.
  def __lt__(self, other: _Operand, /) -> BoolArray: ...
  def __le__(self, other: _Operand, /) -> BoolArray: ...
  def __gt__(self, other: _Operand, /) -> BoolArray: ...

  def __bool__(self) -> bool: ...


@runtime_checkable
class BoolArray(Array, Protocol):
  """An array of booleans: what a comparison returns, and what indexes a mask.

  Everything ``Array`` provides, plus ``~``, ``&`` and ``|``.
  """

  # Those three live here rather than on `Array` because they raise on a float
  # array, which is what the bases are parameterized by.
  def __invert__(self) -> Self: ...
  def __and__(self, other: _Operand, /) -> Self: ...
  def __or__(self, other: _Operand, /) -> Self: ...


# PEP 695 syntax (``class BicopLike[ArrayT]``) needs 3.12 and a PEP 696
# ``default=Any`` needs 3.13 in the standard library; the floor here is 3.11,
# so neither is usable yet. Revisit both in one pass when it moves.
#: Array type an implementation commits to (``numpy.ndarray`` | ``torch.Tensor``).
ArrayT = TypeVar("ArrayT", bound=Array)

__all__ = [
  "Array",
  "ArrayT",
  "BicopLike",
  "BoolArray",
  "ControlsLike",
  "FInfo",
  "MarginLike",
  "Namespace",
  "VinecopLike",
  "VinedistLike",
  "array_namespace",
]


@runtime_checkable
class ControlsLike(Protocol):
  """Contract for a fit-configuration object.

  Controls carry *how* something is fitted — the family set, the selection
  criterion, the tree algorithm — and are handed to a ``fit`` / ``select`` /
  ``from_data`` call, which reads the settings it owns. The contract is a
  single method rather than a list of fields, because the settings themselves
  differ per layer and per array namespace: :class:`FitControlsBicop` and the
  PyTorch :class:`~pyvinecopulib.torch.FitControlsTorchBicop` have no field
  name in common, so a field-by-field contract would fit neither.

  Note this describes only the *shape*, not that any given consumer honors
  every setting. A consumer reads the keys it owns and passes the object on —
  a vine's structure selection reads ``tree_criterion`` while the per-edge
  ``family_set`` is the pair-copula fitter's business. A setting the consumer
  can neither honor nor delegate is refused rather than dropped.

  See Also
  --------
  pyvinecopulib.core.FitControlsBicop : Controls for a
      :class:`~pyvinecopulib.core.Bicop` fit.
  pyvinecopulib.core.FitControlsVinecop : Controls for a
      :class:`~pyvinecopulib.core.Vinecop` fit.
  """

  @abstractmethod
  def to_dict(self) -> dict[str, Any]:
    """Return the settings as a plain dictionary.

    Returns
    -------
    dict
        One entry per setting, keyed by the attribute name.
    """


@runtime_checkable
class FInfo(Protocol):
  """The floating-point limits of a dtype, as ``Namespace.finfo`` reports.

  Attributes
  ----------
  eps : float
      Smallest representable difference from one.
  tiny : float
      Smallest positive normal value.
  """

  eps: float
  tiny: float


class Namespace(Protocol[ArrayT]):
  """The module that operates on arrays -- ``xp`` by convention.

  Generic in the array type, so ``xp.stack(...)`` on a PyTorch model gives a
  ``torch.Tensor`` and on a NumPy one an ``ndarray``. Reach it through
  ``array_namespace`` rather than importing a particular array library,
  which is what lets one implementation run on either.

  ``dtype`` and ``device`` travel as opaque values: read one off an
  ``Array`` and pass it back here.
  """

  # -- construction ------------------------------------------------------ #
  def asarray(
    self,
    obj: object,
    /,
    *,
    dtype: object = None,
    device: object = None,
    copy: bool | None = None,
  ) -> ArrayT: ...

  def empty(
    self,
    shape: int | tuple[int, ...],
    /,
    *,
    dtype: object = None,
    device: object = None,
  ) -> ArrayT: ...

  def zeros(
    self,
    shape: int | tuple[int, ...],
    /,
    *,
    dtype: object = None,
    device: object = None,
  ) -> ArrayT: ...

  def full(
    self,
    shape: int | tuple[int, ...],
    fill_value: bool | float,
    /,
    *,
    dtype: object = None,
    device: object = None,
  ) -> ArrayT: ...

  def empty_like(
    self, x: ArrayT, /, *, dtype: object = None, device: object = None
  ) -> ArrayT: ...

  def zeros_like(
    self, x: ArrayT, /, *, dtype: object = None, device: object = None
  ) -> ArrayT: ...

  def ones_like(
    self, x: ArrayT, /, *, dtype: object = None, device: object = None
  ) -> ArrayT: ...

  def full_like(
    self,
    x: ArrayT,
    /,
    fill_value: bool | float,
    *,
    dtype: object = None,
    device: object = None,
  ) -> ArrayT: ...

  # -- shape ------------------------------------------------------------- #
  def stack(self, arrays: Sequence[ArrayT], /, *, axis: int = 0) -> ArrayT: ...

  def concat(
    self, arrays: Sequence[ArrayT], /, *, axis: int | None = 0
  ) -> ArrayT: ...

  def reshape(
    self, x: ArrayT, /, shape: tuple[int, ...], *, copy: bool | None = None
  ) -> ArrayT: ...

  def matrix_transpose(self, x: ArrayT, /) -> ArrayT: ...

  # -- elementwise ------------------------------------------------------- #
  def abs(self, x: ArrayT, /) -> ArrayT: ...

  def exp(self, x: ArrayT, /) -> ArrayT: ...

  def log(self, x: ArrayT, /) -> ArrayT: ...

  def round(self, x: ArrayT, /) -> ArrayT: ...

  def minimum(self, x1: ArrayT, x2: ArrayT, /) -> ArrayT: ...

  def maximum(self, x1: ArrayT, x2: ArrayT, /) -> ArrayT: ...

  def clip(
    self,
    x: ArrayT,
    /,
    # Mirrors the array API standard's own `clip` signature, where these are
    # the keyword names. Renaming them would describe a different function.
    min: float | ArrayT | None = None,  # noqa: A002
    max: float | ArrayT | None = None,  # noqa: A002
  ) -> ArrayT: ...

  # -- predicates, which answer in booleans ------------------------------ #
  def isnan(self, x: ArrayT, /) -> BoolArray: ...

  def isinf(self, x: ArrayT, /) -> BoolArray: ...

  def isfinite(self, x: ArrayT, /) -> BoolArray: ...

  def where(self, condition: Array, x1: ArrayT, x2: ArrayT, /) -> ArrayT: ...

  # `any` / `all` reduce to a zero-dimensional array, which every caller here
  # immediately passes to `bool()`.
  def any(
    self, x: Array, /, *, axis: int | tuple[int, ...] | None = None
  ) -> Array: ...

  def all(
    self, x: Array, /, *, axis: int | tuple[int, ...] | None = None
  ) -> Array: ...

  # -- reduction and dtype ----------------------------------------------- #
  def sum(
    self,
    x: ArrayT,
    /,
    *,
    axis: int | tuple[int, ...] | None = None,
    dtype: object = None,
  ) -> ArrayT: ...

  def mean(
    self, x: ArrayT, /, *, axis: int | tuple[int, ...] | None = None
  ) -> ArrayT: ...

  def astype(
    self, x: ArrayT, dtype: object, /, *, copy: bool = True
  ) -> ArrayT: ...

  def isdtype(self, dtype: object, kind: str | tuple[str, ...], /) -> bool: ...

  # The array API standard names this parameter `type`. Positional-only, so
  # no caller spells it, but the protocol has to describe the signature.
  def finfo(self, type: object, /) -> FInfo: ...  # noqa: A002


def array_namespace(*arrays: object) -> Namespace[Any]:
  """The namespace that operates on ``arrays``.

  ``array_api_compat``'s own resolver, typed. Loose in its arguments and its
  parameterization, because a caller may hand it a NumPy array and
  read a PyTorch one back out of the namespace it gets (``place`` does exactly
  that), so pinning the two together would be wrong.

  Parameters
  ----------
  *arrays : array
      One or more arrays, which must share a namespace.

  Returns
  -------
  Namespace
      The namespace, typed to the functions this package calls.
  """
  return cast("Namespace[Any]", _array_namespace(*arrays))


_VINEDIST_EXAMPLE = """

  Examples
  --------
  A distribution on the data scale from a fitted copula and explicit margins::

      import numpy as np, scipy.stats as st, pyvinecopulib as pv

      rng = np.random.default_rng(0)
      z = rng.normal(size=(300, 3))
      y = np.column_stack([z[:, 0], np.abs(z[:, 1]) + 0.5, z[:, 2]])

      u = pv.utils.to_pseudo_obs(y)
      dist = pv.Vinedist(
        pv.Vinecop.from_data(u),
        margins=[st.norm(0, 1), st.gamma(2.0), st.norm(0, 1)],
      )
      dist.logpdf(y)
      dist.sample(100, seeds=[1])

  Or fitted end to end, margins first and the copula on the pseudo-observations
  they produce::

      dist = pv.Vinedist.from_data(y)
      dist.margin_summary()
"""

# The worked examples are shared verbatim between each contract (``*Like``) and
# its canonical base (``*Base``, in bicop_base.py / vinecop_base.py) — defined
# once here and appended to both docstrings so the two never drift apart.
_BICOP_EXAMPLE = """

  Examples
  --------
  A minimal independence pair on NumPy — implement only the three primitives and
  inherit ``hinv1`` / ``hinv2`` / ``sample`` / ``loglik`` / ``plot`` /
  ``__repr__`` from :class:`~pyvinecopulib.core.BicopBase`::

      import numpy as np
      from pyvinecopulib.core import BicopBase

      class Independence(BicopBase[np.ndarray]):
        def _pdf_raw(self, u):
          return np.ones(u.shape[0])

        def _hfunc1_raw(self, u):
          return u[:, 1]

        def _hfunc2_raw(self, u):
          return u[:, 0]

        def _sample_uniform(self, n, qrng, seeds):
          return np.random.default_rng(seeds[0] if seeds else 0).uniform(
            size=(n, 2)
          )

      cop = Independence()
      cop.hinv1(np.array([[0.3, 0.7]]))   # -> array([0.7]) (numerical inverse)
"""

_VINECOP_EXAMPLE = """

  Examples
  --------
  Host copulas in a vine by subclassing
  :class:`~pyvinecopulib.core.VinecopBase`; the only required hook is
  ``get_pair_copula``. Under the default
  :class:`~pyvinecopulib.core.SimplifiedContext` it is a classic (unconditional,
  simplified) vine, and hosting :class:`~pyvinecopulib.core.Bicop` pairs
  reproduces :meth:`~pyvinecopulib.core.Vinecop.from_structure`::

      import numpy as np
      import pyvinecopulib as pv
      from pyvinecopulib.core import VinecopBase

      class ListVinecop(VinecopBase[np.ndarray]):
        def __init__(self, pairs, structure):
          self._pairs = pairs
          self._bind_vine(structure)          # SimplifiedContext by default

        def get_pair_copula(self, tree, edge):
          return self._pairs[tree][edge]

      struct = pv.RVineStructure.from_order([1, 2, 3])
      g = pv.families.gaussian
      pairs = [
        [pv.Bicop(family=g, parameters=np.array([[0.5]])),
         pv.Bicop(family=g, parameters=np.array([[0.4]]))],
        [pv.Bicop(family=g, parameters=np.array([[0.2]]))],
      ]
      vine = ListVinecop(pairs, struct)
      ref = pv.Vinecop.from_structure(structure=struct, pair_copulas=pairs)
      u = np.random.default_rng(0).uniform(size=(4, 3))
      np.allclose(vine.pdf(u), ref.pdf(u))    # -> True
"""


_MARGIN_EXAMPLE = """

  Examples
  --------
  A shifted exponential margin on NumPy — implement only the two primitives and
  inherit ``icdf`` / ``logpdf`` / ``cdf_left`` / ``loglik`` / ``__repr__`` from
  :class:`~pyvinecopulib.core.MarginBase`::

      import numpy as np
      from pyvinecopulib.core import MarginBase

      class ShiftedExp(MarginBase[np.ndarray]):
        def __init__(self, rate=1.0, shift=0.0):
          self.rate, self.shift = rate, shift

        @property
        def support(self):
          return (self.shift, float("inf"))

        def pdf(self, y):
          return self.rate * np.exp(-self.rate * (y - self.shift))

        def cdf(self, y):
          return 1.0 - np.exp(-self.rate * (y - self.shift))

      m = ShiftedExp(rate=2.0, shift=1.0)
      m.icdf(np.array([0.5]))   # -> array([1.34657359]) (numerical inverse)

  Neither primitive declares ``x``, which is what marks this margin
  unconditional: covariates are never forwarded to it. A conditional margin
  takes ``(self, y, *, x=None)`` instead and sets
  ``supports_covariates = True``, and then every inherited member forwards the
  covariates it was called with.
"""


@runtime_checkable
class BicopLike(Protocol[ArrayT]):
  """Contract for a bivariate (optionally conditional) pair copula.

  A pair copula maps pseudo-observations ``u`` of shape ``(n, 2)`` (in the unit
  square, clamped strictly inside it) to a density (``pdf``), the two
  conditional distributions
  ``hfunc1(u) = P(U2 <= u2 | U1 = u1)`` / ``hfunc2(u) = P(U1 <= u1 | U2 = u2)``
  and their inverses (``hinv1`` / ``hinv2``, inverting in the second / first
  argument), plus a sampler (``sample``). The distribution ``cdf`` is an
  optional capability rather than a member; see below. The optional ``x`` of
  shape ``(n, p)`` carries conditioning variables: a conditional copula reads
  them, an unconditional one ignores them.

  The easy way to satisfy this contract is to subclass
  :class:`~pyvinecopulib.core.BicopBase`, which supplies ``hinv1`` / ``hinv2``
  (numerical inversion) on top of ``pdf`` / ``hfunc1`` / ``hfunc2``; providing
  its ``_sample_uniform`` hook enables the inherited ``sample``.
  :class:`pyvinecopulib.core.Bicop` and
  :class:`pyvinecopulib.torch.TorchTllBicop` are the reference implementations.

  **What is required is the evaluation surface.** ``pdf`` / ``hfunc1`` /
  ``hfunc2`` / ``hinv1`` / ``hinv2`` / ``sample`` are all of it -- the
  surface :class:`pyvinecopulib.core.Bicop` presents, so that the contract is
  something a foreign pair copula can be typed against rather than a list of
  whichever methods the cascades happen to call today. (``sample`` is in it for
  that reason: no vine cascade asks a pair to sample, a vine drawing by inverse
  Rosenblatt, but a pair copula that cannot be drawn from is not one.)
  Everything else is an **optional capability**, read with ``getattr`` where it
  is needed, so a foreign object provides it only if it applies. ``cdf`` and
  ``flip`` are the two, and :class:`~pyvinecopulib.core.BicopBase` supplies
  both as raising stubs, which is where the message explaining each lives -- so
  a subclass gets a good error and a foreign object simply omits them. A third,
  ``supports_batched``, is a plain declaration: whether a vine may stack this
  pair into its stacked grid cascade, which reads an interpolation grid off
  each pair. Absent means it may not, and
  :class:`~pyvinecopulib.core.BicopBase` declares it ``False`` so the answer is
  findable rather than only discoverable by tripping the error.

  **Covariates are a widening, not a member.** There is no
  ``supports_covariates`` flag here, unlike on a margin or a whole copula,
  because for a pair copula the *signature* decides: a pair that reads
  conditioning variables declares a keyword-only ``x`` on the members that read
  one, and ``pair_eval`` forwards a matrix whenever there is one -- so a pair
  that models none, :class:`~pyvinecopulib.core.Bicop` above all, raises rather
  than quietly answering unconditionally. The contract above is the
  unconditional surface every pair serves; declaring ``x`` in it would put
  ``Bicop`` outside its own contract.

  See Also
  --------
  pyvinecopulib.core.BicopBase : Canonical partial implementation to subclass.
  pyvinecopulib.core.Bicop : The reference pair copula.
  VinecopLike : The vine-level evaluator contract.
  """

  @abstractmethod
  def pdf(self, u: ArrayT) -> ArrayT:
    """Pair-copula density ``c(u)`` at each observation.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pair pseudo-observations in the unit square.

    Returns
    -------
    array, shape (n,), dtype float
        Density values.
    """

  @abstractmethod
  def hfunc1(self, u: ArrayT) -> ArrayT:
    r"""First h-function ``P(U2 <= u2 | U1 = u1)``.

    .. math::

       h_1(u_1, u_2) = \\mathbb{P}(U_2 \\le u_2 \\mid U_1 = u_1).

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pair pseudo-observations; conditions on the first column.

    Returns
    -------
    array, shape (n,), dtype float
        Conditional distribution values in ``[0, 1]``.
    """

  @abstractmethod
  def hfunc2(self, u: ArrayT) -> ArrayT:
    r"""Second h-function ``P(U1 <= u1 | U2 = u2)``.

    .. math::

       h_2(u_1, u_2) = \\mathbb{P}(U_1 \\le u_1 \\mid U_2 = u_2).

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pair pseudo-observations; conditions on the second column.

    Returns
    -------
    array, shape (n,), dtype float
        Conditional distribution values in ``[0, 1]``.
    """

  @abstractmethod
  def hinv1(self, u: ArrayT) -> ArrayT:
    """Inverse of :meth:`hfunc1` in its second argument.

    Solves ``hfunc1([u1, .], x) = u2`` for the second argument.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Column 0 is the conditioning value ``u1``; column 1 is the level to
        invert.

    Returns
    -------
    array, shape (n,), dtype float
        The inverted values in ``[0, 1]``.
    """

  @abstractmethod
  def hinv2(self, u: ArrayT) -> ArrayT:
    """Inverse of :meth:`hfunc2` in its first argument.

    Solves ``hfunc2([., u2], x) = u1`` for the first argument.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Column 0 is the level to invert; column 1 is the conditioning value
        ``u2``.

    Returns
    -------
    array, shape (n,), dtype float
        The inverted values in ``[0, 1]``.
    """

  @abstractmethod
  def sample(
    self,
    n: int,
    *,
    qrng: bool = False,
    seeds: list[int] | None = None,
  ) -> ArrayT:
    """Draw ``n`` samples from the pair copula.

    Parameters
    ----------
    n : int
        Number of samples to draw.
    qrng : bool, default=False
        Draw quasi-random base uniforms instead of pseudo-random ones.
    seeds : list of int, or None, optional
        RNG seeds.

    Returns
    -------
    array, shape (n, 2), dtype float
        Samples in the unit square.
    """


BicopLike.__doc__ = (BicopLike.__doc__ or "") + _BICOP_EXAMPLE


@runtime_checkable
class VinecopLike(Protocol[ArrayT]):
  """Contract for a post-fit vine-copula evaluator.

  A vine evaluator exposes the joint ``pdf`` / ``cdf``, the ``rosenblatt`` and
  ``inverse_rosenblatt`` transforms, and a ``sample`` sampler, plus a
  ``structure`` attribute (the :class:`~pyvinecopulib.core.RVineStructure` it was
  built on). The optional keyword-only ``x`` covariate matrix (row-aligned with
  ``u``) is the conditional extension — left ``None`` for the usual
  (unconditional) case. Subclass :class:`~pyvinecopulib.core.VinecopBase` to get
  every method from a small set of hooks;
  :class:`pyvinecopulib.core.Vinecop` and
  :class:`pyvinecopulib.torch.TorchVinecop` are the reference implementations.

  Everything past that surface is an **optional capability**, read with
  ``getattr`` where it is needed. ``logpdf`` is the one: the joint log-density,
  which a vine distribution prefers over the logarithm of the density because
  the density is a product of up to ``d (d - 1) / 2`` pair densities and
  underflows on a deep or strongly dependent model. Both reference
  implementations supply it.

  See Also
  --------
  pyvinecopulib.core.VinecopBase : Canonical partial implementation to subclass.
  pyvinecopulib.core.Vinecop : The reference vine.
  BicopLike : The pair-copula contract.
  """

  @property
  @abstractmethod
  def structure(self) -> RVineStructure:
    """The R-vine structure the model was built on.

    Declared read-only, which is the weaker requirement: a vine that exposes a
    settable attribute satisfies it just as one with a getter alone does.

    Returns
    -------
    RVineStructure
        The structure the pair copulas are indexed by.
    """

  @abstractmethod
  def pdf(self, u: ArrayT) -> ArrayT:
    """Joint vine-copula density ``c(u_1, ..., u_d)`` at each observation.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Pseudo-observations in ``[0, 1]^d``.

    Returns
    -------
    array, shape (n,), dtype float
        Joint density values.
    """

  @abstractmethod
  def cdf(
    self,
    u: ArrayT,
    *,
    N: int = 10000,
    seeds: list[int] | None = None,
  ) -> ArrayT:
    """Joint vine-copula distribution ``C(u)`` via Monte-Carlo.

    Parameters
    ----------
    u : array, shape (m, d), dtype float
        Query points in ``[0, 1]^d``.
    N : int, default=10000
        Number of Monte-Carlo samples.
    seeds : list of int, or None, optional
        RNG seeds.

    Returns
    -------
    array, shape (m,), dtype float
        Distribution values in ``[0, 1]``.
    """

  @abstractmethod
  def rosenblatt(self, u: ArrayT) -> ArrayT:
    """Rosenblatt transform: dependent uniforms to independent uniforms.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Pseudo-observations in ``[0, 1]^d``.

    Returns
    -------
    array, shape (n, d), dtype float
        Independent uniforms.
    """

  @abstractmethod
  def inverse_rosenblatt(self, u: ArrayT) -> ArrayT:
    """Inverse Rosenblatt transform: independent uniforms to dependent uniforms.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Independent uniforms in ``[0, 1]^d``.

    Returns
    -------
    array, shape (n, d), dtype float
        Dependent uniforms.
    """

  @abstractmethod
  def sample(
    self,
    n: int,
    *,
    qrng: bool = False,
    seeds: list[int] | None = None,
  ) -> ArrayT:
    """Draw ``n`` samples from the fitted vine copula.

    Parameters
    ----------
    n : int
        Number of samples to draw.
    qrng : bool, default=False
        Draw quasi-random base uniforms instead of pseudo-random ones.
    seeds : list of int, or None, optional
        RNG seeds.

    Returns
    -------
    array, shape (n, d), dtype float
        Samples in ``[0, 1]^d``.
    """


VinecopLike.__doc__ = (VinecopLike.__doc__ or "") + _VINECOP_EXAMPLE


@runtime_checkable
class MarginLike(Protocol[ArrayT]):
  """Contract for a fitted univariate margin.

  A margin maps observations ``y`` on the original scale to a density
  (``pdf``), a distribution (``cdf``), and the inverse distribution
  (``icdf``) — enough to turn a vine copula into a full multivariate
  distribution and back.

  Every member also accepts optional exogenous covariates ``x``, so a margin
  may be *conditional*: ``f(y | x)`` rather than ``f(y)``. As on
  :class:`BicopLike`, ``x`` is keyword-only and the observations are
  positional-only — a margin from another library names its first parameter
  whatever it likes (:class:`pyvinecopulib.core.Kde1d` calls it ``x``, which is
  exactly the collision this shape avoids), and a covariate matrix must never
  bind to some unrelated positional slot. A margin that does not model
  covariates ignores ``x``, as one vine distribution may mix conditional and
  unconditional margins.

  ``pdf`` is the density with respect to **the margin's own reference
  measure**: a Lebesgue density for a continuous margin, a probability mass at
  an atom, and whichever applies pointwise for a mixed (e.g. zero-inflated)
  one. That is what makes the Sklar factorization hold unchanged in all three
  cases, with no branch on the variable type::

      log f(x) = log c(F_1(x_1), ..., F_d(x_d)) + sum_j log pdf_j(x_j)

  Notes
  -----
  Discreteness, the left-limit cdf, log-densities, sampling and support are
  **optional capabilities** rather than members of this contract, so that
  objects from other ecosystems can satisfy it. Consumers discover them with
  ``getattr(margin, name, None)``:

  - ``var_type`` — ``"c"``, ``"d"`` or ``"zi"``; absent means ``"c"``.
  - ``cdf_left`` — ``F(y^-)``, the left limit a margin with atoms needs;
    absent means it coincides with ``cdf``.
  - ``logpdf`` — absent means ``log(pdf)``.
  - ``sample`` — absent means ``icdf`` of uniforms.
  - ``support`` — a ``(lo, hi)`` pair; absent means unbounded. It describes the
    margin as a whole, not one conditional slice: a margin whose support moves
    with ``x`` overrides ``icdf`` instead.
  - ``supports_covariates`` — whether ``x`` is read rather than ignored; absent
    means it is ignored, and a consumer then omits it entirely.

  See Also
  --------
  pyvinecopulib.core.MarginBase : Canonical partial implementation to subclass.
  pyvinecopulib.core.Kde1d : The default nonparametric margin.
  BicopLike : The pair-copula contract.
  """

  @abstractmethod
  def pdf(self, y: ArrayT, /) -> ArrayT:
    """Density of the margin with respect to its own reference measure.

    Parameters
    ----------
    y : array, shape (n,), dtype float
        Observations on the original scale.

    Returns
    -------
    array, shape (n,), dtype float
        A density, a probability mass, or a mixture of the two, depending on
        the margin.
    """

  @abstractmethod
  def cdf(self, y: ArrayT, /) -> ArrayT:
    """Distribution function ``F(y)``, right-continuous.

    Parameters
    ----------
    y : array, shape (n,), dtype float
        Observations on the original scale.

    Returns
    -------
    array, shape (n,), dtype float
        Distribution values in ``[0, 1]``.
    """

  @abstractmethod
  def icdf(self, p: ArrayT, /) -> ArrayT:
    """Inverse distribution function ``inf{y : F(y) >= p}``.

    Parameters
    ----------
    p : array, shape (n,), dtype float
        Probabilities in ``[0, 1]``.

    Returns
    -------
    array, shape (n,), dtype float
        Quantiles on the original scale.
    """


MarginLike.__doc__ = (MarginLike.__doc__ or "") + _MARGIN_EXAMPLE


@runtime_checkable
class VinedistLike(Protocol[ArrayT]):
  """Contract for a vine distribution on the data scale.

  Where :class:`VinecopLike` evaluates on the copula scale ``[0, 1]^d``, this
  evaluates on the original scale: it is a copula combined with one margin per
  variable, i.e. Sklar's theorem as an object. The surface mirrors the copula
  one — ``pdf`` / ``cdf`` / ``rosenblatt`` / ``inverse_rosenblatt`` / ``sample``
  — plus ``logpdf`` (the primitive, summed in log space) and ``loglik``, and it
  reads ``y`` rather than ``u``. Subclass
  :class:`~pyvinecopulib.core.VinedistBase` to get all of it from the two
  halves; :class:`pyvinecopulib.core.Vinedist` and
  :class:`pyvinecopulib.torch.TorchVinedist` are the reference
  implementations.

  ``cdf`` / ``rosenblatt`` / ``inverse_rosenblatt`` / ``sample`` keep a
  ``**kwargs`` that :class:`VinecopLike` spells out as ``N`` / ``seeds``: a
  vine distribution forwards whatever options the copula it holds accepts, and
  those differ between
  :class:`~pyvinecopulib.core.Vinecop` and
  :class:`~pyvinecopulib.torch.TorchVinecop`.

  Notes
  -----
  Discreteness, conditioning and the fit-time reports are **optional
  capabilities** rather than members of this contract, discovered with
  ``getattr``: ``dim``, ``var_types``, ``sample_conditional`` and
  ``margin_summary`` -- the last of which
  :class:`pyvinecopulib.sklearn.VineDensity` reads to publish
  ``margin_summary_``, so a distribution class without one is
  told so by name.

  A distribution declares no ``supports_covariates`` of its own. It reads the
  flag on the parts it holds -- every margin, and the copula -- and refuses an
  ``x`` that none of them reads; there is nothing above it to read a flag of
  its own, and a capability nothing consumes is one that goes stale.
  Serialization is likewise out of scope, as it is for the other contracts.

  See Also
  --------
  pyvinecopulib.core.VinedistBase : Canonical partial implementation.
  pyvinecopulib.core.Vinedist : The reference vine distribution.
  VinecopLike : The copula half's contract.
  MarginLike : The marginal half's contract.
  """

  vinecop: VinecopLike[ArrayT]
  margins: Sequence[MarginLike[ArrayT]]

  @abstractmethod
  def logpdf(self, y: ArrayT) -> ArrayT:
    """Joint log-density at each observation.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.

    Returns
    -------
    array, shape (n,), dtype float
        Joint log-density values.
    """

  @abstractmethod
  def pdf(self, y: ArrayT) -> ArrayT:
    """Joint density at each observation.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.

    Returns
    -------
    array, shape (n,), dtype float
        Joint density values.
    """

  @abstractmethod
  def loglik(self, y: ArrayT) -> ArrayT:
    """Log-likelihood of the observations.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.

    Returns
    -------
    array, shape (), dtype float
        The summed log-density, kept zero-dimensional so it stays
        differentiable under autograd.
    """

  @abstractmethod
  def copula_layout(self, y: ArrayT) -> ArrayT:
    """This distribution's copula-scale data for ``y``.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.

    Returns
    -------
    array, shape (n, d + k), dtype float
        The compact copula-scale layout: one column per variable, followed by
        a left-limit column for each discrete variable.
    """

  @abstractmethod
  def cdf(
    self,
    y: ArrayT,
    **kwargs: Any,  # noqa: ANN401 - forwarded to the implementation
  ) -> ArrayT:
    """Joint distribution function at each observation.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    **kwargs : Any
        Forwarded to the copula's ``cdf``.

    Returns
    -------
    array, shape (n,), dtype float
        Distribution values in ``[0, 1]``.
    """

  @abstractmethod
  def rosenblatt(
    self,
    y: ArrayT,
    **kwargs: Any,  # noqa: ANN401 - forwarded to the implementation
  ) -> ArrayT:
    """Rosenblatt transform: observations to independent uniforms.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    **kwargs : Any
        Forwarded to the copula's ``rosenblatt``.

    Returns
    -------
    array, shape (n, d), dtype float
        Independent uniforms.
    """

  @abstractmethod
  def inverse_rosenblatt(
    self,
    w: ArrayT,
    **kwargs: Any,  # noqa: ANN401 - forwarded to the implementation
  ) -> ArrayT:
    """Inverse Rosenblatt transform: independent uniforms to observations.

    Parameters
    ----------
    w : array, shape (n, d), dtype float
        Independent uniforms in ``[0, 1]^d``.
    **kwargs : Any
        Forwarded to the copula's ``inverse_rosenblatt``.

    Returns
    -------
    array, shape (n, d), dtype float
        Observations on the original scale.
    """

  @abstractmethod
  def sample(
    self,
    n: int,
    **kwargs: Any,  # noqa: ANN401 - forwarded to the implementation
  ) -> ArrayT:
    """Draw observations on the original scale.

    Parameters
    ----------
    n : int
        Number of observations.
    **kwargs : Any
        Forwarded to the copula's ``sample``.

    Returns
    -------
    array, shape (n, d), dtype float
        The drawn observations.
    """


VinedistLike.__doc__ = (VinedistLike.__doc__ or "") + _VINEDIST_EXAMPLE
