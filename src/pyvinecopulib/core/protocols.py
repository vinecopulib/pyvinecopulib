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
composed into a distribution, or consumed by the sklearn backend layer. The
reference implementations are :class:`~pyvinecopulib.core.Bicop`,
:class:`~pyvinecopulib.core.Vinecop`, :class:`~pyvinecopulib.core.Kde1d` and
:class:`~pyvinecopulib.core.Vinedist`, and their PyTorch counterparts in
:mod:`pyvinecopulib.torch`.

:class:`ControlsLike` is the odd one out: it describes *fit configuration*
rather than a model, and asks for a single ``to_dict``.

**Conditioning.** Every method carries an optional **keyword-only** ``x`` — the
conditioning variables the model depends on (conditioning-set values and/or
external covariates), row-aligned with the data. Unconditional models leave it
``None``; a conditional pair copula reads it. In a vine, each pair's ``x`` is
assembled per edge by a :class:`~pyvinecopulib.core.ConditioningContext`.

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
only: ``isinstance(cop, BicopLike)`` reports that the names are present, not
that the signatures agree.

**Typing.** :data:`ArrayT` is an unbounded ``TypeVar`` carried only on these
public signatures, so a concrete implementation (e.g.
:class:`~pyvinecopulib.torch.TorchBicop`) inherits precise ``torch.Tensor``
return types. The numeric implementations in
:mod:`~pyvinecopulib.core.bicop_base` operate on arrays as ``Any`` (the Array API
namespace ``array_api_compat`` is itself untyped).
"""

from __future__ import annotations

from abc import abstractmethod
from typing import (
  TYPE_CHECKING,
  Any,
  Optional,
  Protocol,
  TypeVar,
  runtime_checkable,
)

if TYPE_CHECKING:
  from ..pyvinecopulib_ext import RVineStructure

#: Array type an implementation commits to (``numpy.ndarray`` | ``torch.Tensor``).
ArrayT = TypeVar("ArrayT")

__all__ = [
  "ArrayT",
  "BicopLike",
  "ControlsLike",
  "MarginLike",
  "VinecopLike",
  "VinedistLike",
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

  Passing it on works because a vine's controls *are* pair controls:
  :class:`FitControlsVinecop` derives from :class:`FitControlsBicop`, and
  :class:`~pyvinecopulib.torch.FitControlsTorchVinecop` from
  :class:`~pyvinecopulib.torch.FitControlsTorchBicop`. So one object configures
  both halves of a vine fit, which is also how observation weights reach the
  pair copulas.

  See Also
  --------
  pyvinecopulib.core.FitControlsBicop : Controls for a
      :class:`~pyvinecopulib.core.Bicop` fit.
  pyvinecopulib.core.FitControlsVinecop : Controls for a
      :class:`~pyvinecopulib.core.Vinecop` fit.
  """

  @abstractmethod
  def to_dict(self) -> dict:
    """Return the settings as a plain dictionary.

    Returns
    -------
    dict
        One entry per setting, keyed by the attribute name.
    """


_VINEDIST_EXAMPLE = """

  Examples
  --------
  A distribution on the data scale from a fitted copula and explicit margins::

      import numpy as np, scipy.stats as st, pyvinecopulib as pv

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
        def pdf(self, u, *, x=None):
          return np.ones(u.shape[0])

        def hfunc1(self, u, *, x=None):
          return u[:, 1]

        def hfunc2(self, u, *, x=None):
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
  square, clamped to ``[1e-10, 1 - 1e-10]``) to a density (``pdf``), a
  distribution (``cdf``), the two conditional distributions
  ``hfunc1(u) = P(U2 <= u2 | U1 = u1)`` / ``hfunc2(u) = P(U1 <= u1 | U2 = u2)``
  and their inverses (``hinv1`` / ``hinv2``, inverting in the second / first
  argument), plus a sampler (``sample``). The optional ``x`` of shape
  ``(n, k)`` carries conditioning variables: a conditional copula reads them, an
  unconditional one ignores them.

  The easy way to satisfy this contract is to subclass
  :class:`~pyvinecopulib.core.BicopBase`, which supplies ``hinv1`` / ``hinv2``
  (numerical inversion) on top of ``pdf`` / ``hfunc1`` / ``hfunc2``; providing
  its ``_sample_uniform`` hook enables the inherited ``sample``.
  :class:`pyvinecopulib.core.Bicop` and
  :class:`pyvinecopulib.torch.TorchBicop` are the reference implementations.

  See Also
  --------
  pyvinecopulib.core.BicopBase : Canonical partial implementation to subclass.
  pyvinecopulib.core.Bicop : The reference pair copula.
  VinecopLike : The vine-level evaluator contract.
  """

  @abstractmethod
  def pdf(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Pair-copula density ``c(u)`` at each observation.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pair pseudo-observations in the unit square.
    x : array, shape (n, k), or None, optional
        Conditioning variables; ignored by an unconditional copula.

    Returns
    -------
    array, shape (n,), dtype float
        Density values.
    """

  @abstractmethod
  def cdf(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Pair-copula distribution ``C(u)`` at each observation.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pair pseudo-observations in the unit square.
    x : array, shape (n, k), or None, optional
        Conditioning variables; ignored by an unconditional copula.

    Returns
    -------
    array, shape (n,), dtype float
        Distribution values in ``[0, 1]``.
    """

  @abstractmethod
  def hfunc1(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """First h-function ``P(U2 <= u2 | U1 = u1)``.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pair pseudo-observations; conditions on the first column.
    x : array, shape (n, k), or None, optional
        Conditioning variables; ignored by an unconditional copula.

    Returns
    -------
    array, shape (n,), dtype float
        Conditional distribution values in ``[0, 1]``.
    """

  @abstractmethod
  def hfunc2(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Second h-function ``P(U1 <= u1 | U2 = u2)``.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Pair pseudo-observations; conditions on the second column.
    x : array, shape (n, k), or None, optional
        Conditioning variables; ignored by an unconditional copula.

    Returns
    -------
    array, shape (n,), dtype float
        Conditional distribution values in ``[0, 1]``.
    """

  @abstractmethod
  def hinv1(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Inverse of :meth:`hfunc1` in its second argument.

    Solves ``hfunc1([u1, .], x) = u2`` for the second argument.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Column 0 is the conditioning value ``u1``; column 1 is the level to
        invert.
    x : array, shape (n, k), or None, optional
        Conditioning variables; ignored by an unconditional copula.

    Returns
    -------
    array, shape (n,), dtype float
        The inverted values in ``[0, 1]``.
    """

  @abstractmethod
  def hinv2(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Inverse of :meth:`hfunc2` in its first argument.

    Solves ``hfunc2([., u2], x) = u1`` for the first argument.

    Parameters
    ----------
    u : array, shape (n, 2), dtype float
        Column 0 is the level to invert; column 1 is the conditioning value
        ``u2``.
    x : array, shape (n, k), or None, optional
        Conditioning variables; ignored by an unconditional copula.

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
    x: Optional[ArrayT] = None,
    qrng: bool = False,
    seeds: Optional[list[int]] = None,
  ) -> ArrayT:
    """Draw ``n`` samples from the pair copula.

    Parameters
    ----------
    n : int
        Number of samples to draw.
    x : array, shape (n, k), or None, optional
        Conditioning variables (one row per sample) for a conditional draw.
    qrng : bool, default=False
        Draw quasi-random base uniforms instead of pseudo-random ones.
    seeds : list of int, or None, optional
        RNG seeds.

    Returns
    -------
    array, shape (n, 2), dtype float
        Samples in the unit square.
    """

  @abstractmethod
  def flip(self) -> "BicopLike[ArrayT]":
    """Return the pair copula with its two arguments swapped.

    The flipped copula satisfies ``c'(u1, u2) = c(u2, u1)``, with the two
    h-functions (and their inverses) exchanged accordingly. Structure
    selection (:meth:`~pyvinecopulib.core.VinecopBase.select`) uses it to
    reorient a fitted pair onto its slot in the finalized vine, mirroring
    :class:`~pyvinecopulib.core.Vinecop`.

    Returns
    -------
    BicopLike
        The argument-swapped pair copula; the object itself is left
        unchanged.
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

  See Also
  --------
  pyvinecopulib.core.VinecopBase : Canonical partial implementation to subclass.
  pyvinecopulib.core.Vinecop : The reference vine.
  BicopLike : The pair-copula contract.
  """

  structure: RVineStructure

  @abstractmethod
  def pdf(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Joint vine-copula density ``c(u_1, ..., u_d)`` at each observation.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Pseudo-observations in ``[0, 1]^d``.
    x : array, shape (n, p), or None, optional
        External covariates threaded to each pair copula. A simplified vine may
        still depend on them; ``None`` means there are no external covariates.

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
    x: Optional[ArrayT] = None,
    N: int = 10000,
    seeds: Optional[list[int]] = None,
  ) -> ArrayT:
    """Joint vine-copula distribution ``C(u)`` via Monte-Carlo.

    Parameters
    ----------
    u : array, shape (m, d), dtype float
        Query points in ``[0, 1]^d``.
    x : array, shape (m, p), or None, optional
        External covariates, or ``None``. The canonical
        :class:`~pyvinecopulib.core.VinecopBase` Monte-Carlo implementation does
        not support a non-``None`` value; a conditional implementation may
        override that limitation.
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
  def rosenblatt(self, u: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Rosenblatt transform: dependent uniforms to independent uniforms.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Pseudo-observations in ``[0, 1]^d``.
    x : array, shape (n, p), or None, optional
        External covariates threaded to each pair copula, or ``None``.

    Returns
    -------
    array, shape (n, d), dtype float
        Independent uniforms.
    """

  @abstractmethod
  def inverse_rosenblatt(
    self, u: ArrayT, *, x: Optional[ArrayT] = None
  ) -> ArrayT:
    """Inverse Rosenblatt transform: independent uniforms to dependent uniforms.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Independent uniforms in ``[0, 1]^d``.
    x : array, shape (n, p), or None, optional
        External covariates threaded to each pair copula, or ``None``.

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
    x: Optional[ArrayT] = None,
    qrng: bool = False,
    seeds: Optional[list[int]] = None,
  ) -> ArrayT:
    """Draw ``n`` samples from the fitted vine copula.

    Parameters
    ----------
    n : int
        Number of samples to draw.
    x : array, shape (n, p), or None, optional
        External covariates for a conditional draw (one row per sample), else
        ``None``.
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

  :class:`pyvinecopulib.core.Kde1d` satisfies this contract directly, in all
  three of its modes. Beware that a discrete distribution from another library
  may spell the mass ``pmf`` and use ``pdf`` for the (infinite) Lebesgue
  density; coerce foreign objects rather than relying on member names.

  The easy way to satisfy this contract is to subclass
  :class:`~pyvinecopulib.core.MarginBase`, which supplies ``icdf`` (numerical
  inversion of ``cdf`` over ``support``) on top of ``pdf`` / ``cdf``, plus the
  optional capabilities below.

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
  def pdf(self, y: ArrayT, /, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Density of the margin with respect to its own reference measure.

    Parameters
    ----------
    y : array, shape (n,), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, one row per observation. Ignored by a margin that
        does not model them.

    Returns
    -------
    array, shape (n,), dtype float
        A density, a probability mass, or a mixture of the two, depending on
        the margin.
    """

  @abstractmethod
  def cdf(self, y: ArrayT, /, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Distribution function ``F(y)``, right-continuous.

    Parameters
    ----------
    y : array, shape (n,), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, one row per observation.

    Returns
    -------
    array, shape (n,), dtype float
        Distribution values in ``[0, 1]``.
    """

  @abstractmethod
  def icdf(self, p: ArrayT, /, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Inverse distribution function ``inf{y : F(y) >= p}``.

    Parameters
    ----------
    p : array, shape (n,), dtype float
        Probabilities in ``[0, 1]``.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, one row per observation.

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
  ``getattr``: ``dim``, ``var_types``, ``sample_conditional``,
  ``supports_covariates`` and ``margin_summary``.
  Serialization is likewise out of scope, as it is for the other contracts.

  See Also
  --------
  pyvinecopulib.core.VinedistBase : Canonical partial implementation.
  pyvinecopulib.core.Vinedist : The reference vine distribution.
  VinecopLike : The copula half's contract.
  MarginLike : The marginal half's contract.
  """

  vinecop: object
  margins: object

  @abstractmethod
  def logpdf(self, y: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Joint log-density at each observation.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates, forwarded to each part that reads them.

    Returns
    -------
    array, shape (n,), dtype float
        Joint log-density values.
    """

  @abstractmethod
  def pdf(self, y: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Joint density at each observation.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates.

    Returns
    -------
    array, shape (n,), dtype float
        Joint density values.
    """

  @abstractmethod
  def loglik(self, y: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """Log-likelihood of the observations.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates.

    Returns
    -------
    array, shape (), dtype float
        The summed log-density, kept zero-dimensional so it stays
        differentiable on an autograd backend.
    """

  @abstractmethod
  def copula_layout(self, y: ArrayT, *, x: Optional[ArrayT] = None) -> ArrayT:
    """This distribution's copula-scale data for ``y``.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates.

    Returns
    -------
    array, shape (n, d + k), dtype float
        The compact copula-scale layout: one column per variable, followed by
        a left-limit column for each discrete variable.
    """

  @abstractmethod
  def cdf(
    self, y: ArrayT, *, x: Optional[ArrayT] = None, **kwargs: Any
  ) -> ArrayT:
    """Joint distribution function at each observation.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates.
    **kwargs : Any
        Forwarded to the copula's ``cdf``.

    Returns
    -------
    array, shape (n,), dtype float
        Distribution values in ``[0, 1]``.
    """

  @abstractmethod
  def rosenblatt(
    self, y: ArrayT, *, x: Optional[ArrayT] = None, **kwargs: Any
  ) -> ArrayT:
    """Rosenblatt transform: observations to independent uniforms.

    Parameters
    ----------
    y : array, shape (n, d), dtype float
        Observations on the original scale.
    x : array, shape (n, k), or None, optional
        Exogenous covariates.
    **kwargs : Any
        Forwarded to the copula's ``rosenblatt``.

    Returns
    -------
    array, shape (n, d), dtype float
        Independent uniforms.
    """

  @abstractmethod
  def inverse_rosenblatt(
    self, w: ArrayT, *, x: Optional[ArrayT] = None, **kwargs: Any
  ) -> ArrayT:
    """Inverse Rosenblatt transform: independent uniforms to observations.

    Parameters
    ----------
    w : array, shape (n, d), dtype float
        Independent uniforms in ``[0, 1]^d``.
    x : array, shape (n, k), or None, optional
        Exogenous covariates.
    **kwargs : Any
        Forwarded to the copula's ``inverse_rosenblatt``.

    Returns
    -------
    array, shape (n, d), dtype float
        Observations on the original scale.
    """

  @abstractmethod
  def sample(
    self, n: int, *, x: Optional[ArrayT] = None, **kwargs: Any
  ) -> ArrayT:
    """Draw observations on the original scale.

    Parameters
    ----------
    n : int
        Number of observations.
    x : array, shape (n, k), or None, optional
        Exogenous covariates.
    **kwargs : Any
        Forwarded to the copula's ``sample``.

    Returns
    -------
    array, shape (n, d), dtype float
        The drawn observations.
    """


VinedistLike.__doc__ = (VinedistLike.__doc__ or "") + _VINEDIST_EXAMPLE
