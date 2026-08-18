"""Array-agnostic pair-copula / vine / margin contracts.

:class:`BicopLike` and :class:`VinecopLike` define what a pair copula / vine
evaluator must provide, independent of the array namespace (NumPy or PyTorch).
They are the extension point for custom implementations: implement the protocol
— or, far more easily, subclass the canonical
:class:`~pyvinecopulib.core.BicopBase` / :class:`~pyvinecopulib.core.VinecopBase`,
which fill in most of it — and the object plugs into the rest of the library
(e.g. it can be hosted in a vine, or consumed by the sklearn backend layer). The
reference implementations are :class:`~pyvinecopulib.core.Bicop` /
:class:`~pyvinecopulib.core.Vinecop` (the default) and
:class:`~pyvinecopulib.torch.TorchBicop` /
:class:`~pyvinecopulib.torch.TorchVinecop`.

**Conditioning.** Every method carries an optional **keyword-only** ``x`` — the
conditioning variables the copula / vine depends on (conditioning-set values
and/or external covariates), row-aligned with ``u``. Unconditional models leave
it ``None``; a conditional pair copula reads it. In a vine, each pair's ``x`` is
assembled per edge by a :class:`~pyvinecopulib.core.ConditioningContext`.

**Scope.** These describe the *evaluation* surface — enough to host a pair
copula in a vine and to consume a fitted vine — not the whole compiled API.
:class:`~pyvinecopulib.core.Bicop` and :class:`~pyvinecopulib.core.Vinecop`
carry considerably more (score and derivative families, per-row parameters,
serialization, discrete data layouts) that a custom implementation is not
expected to provide. Both protocols are ``runtime_checkable``, which compares
method *names* only: ``isinstance(cop, BicopLike)`` reports that the names are
present, not that the signatures agree.

**Typing.** :data:`ArrayT` is an unbounded ``TypeVar`` carried only on these
public signatures, so a concrete implementation (e.g.
:class:`~pyvinecopulib.torch.TorchBicop`) inherits precise ``torch.Tensor``
return types. The numeric implementations in
:mod:`~pyvinecopulib.core.bicop_base` operate on arrays as ``Any`` (the Array API
namespace ``array_api_compat`` is itself untyped).
"""

from __future__ import annotations

from abc import abstractmethod
from typing import Optional, Protocol, TypeVar, runtime_checkable

#: Array type an implementation commits to (``numpy.ndarray`` | ``torch.Tensor``).
ArrayT = TypeVar("ArrayT")

__all__ = ["ArrayT", "BicopLike", "MarginLike", "VinecopLike"]


# The worked examples are shared verbatim between each contract (``*Like``) and
# its canonical base (``*Base``, in bicop_base.py / vinecop_base.py) — defined
# once here and appended to both docstrings so the two never drift apart.
_BICOP_EXAMPLE = """
  Examples
  --------
  A minimal independence pair on NumPy — implement only the three primitives and
  inherit ``hinv1`` / ``hinv2`` / ``simulate`` / ``loglik`` / ``plot`` /
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

      cop = Independence()
      cop.hinv1(np.array([[0.3, 0.7]]))   # -> array([0.7]) (numerical inverse)
"""

_VINECOP_EXAMPLE = """
  Examples
  --------
  Host copulas in a vine by subclassing
  :class:`~pyvinecopulib.core.VinecopBase`; the only required hook is
  ``_get_pair_copula``. Under the default
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

        def _get_pair_copula(self, tree, edge):
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

        def pdf(self, x):
          return self.rate * np.exp(-self.rate * (x - self.shift))

        def cdf(self, x):
          return 1.0 - np.exp(-self.rate * (x - self.shift))

      m = ShiftedExp(rate=2.0, shift=1.0)
      m.icdf(np.array([0.5]))   # -> array([1.34657359]) (numerical inverse)
"""


@runtime_checkable
class BicopLike(Protocol[ArrayT]):
  """Contract for a bivariate (optionally conditional) pair copula.

  A pair copula maps pseudo-observations ``u`` of shape ``(n, 2)`` (in the unit
  square, clamped to ``[1e-10, 1 - 1e-10]``) to a density (``pdf``), a
  distribution (``cdf``), the two conditional distributions
  ``hfunc1(u) = P(U2 <= u2 | U1 = u1)`` / ``hfunc2(u) = P(U1 <= u1 | U2 = u2)``
  and their inverses (``hinv1`` / ``hinv2``, inverting in the second / first
  argument), plus a sampler (``simulate``). The optional ``x`` of shape
  ``(n, k)`` carries conditioning variables: a conditional copula reads them, an
  unconditional one ignores them.

  The easy way to satisfy this contract is to subclass
  :class:`~pyvinecopulib.core.BicopBase`, which supplies ``hinv1`` / ``hinv2``
  (numerical inversion) and ``simulate`` on top of ``pdf`` / ``hfunc1`` /
  ``hfunc2``. :class:`pyvinecopulib.core.Bicop` and
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
  def simulate(
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
  ``inverse_rosenblatt`` transforms, and a ``simulate`` sampler, plus a
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

  structure: object

  @abstractmethod
  def pdf(
    self, u: ArrayT, *, x: Optional[ArrayT] = None, num_threads: int = 1
  ) -> ArrayT:
    """Joint vine-copula density ``c(u_1, ..., u_d)`` at each observation.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Pseudo-observations in ``[0, 1]^d``.
    x : array, shape (n, p), or None, optional
        External covariates threaded to each pair copula (non-simplified vines);
        ``None`` for the unconditional / simplified case.
    num_threads : int, default=1
        Accepted for parity with :meth:`pyvinecopulib.core.Vinecop.pdf`.

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
    num_threads: int = 1,
    seeds: Optional[list[int]] = None,
  ) -> ArrayT:
    """Joint vine-copula distribution ``C(u)`` via Monte-Carlo.

    Parameters
    ----------
    u : array, shape (m, d), dtype float
        Query points in ``[0, 1]^d``.
    x : array, shape (m, p), or None, optional
        External covariates (non-simplified vines), else ``None``.
    N : int, default=10000
        Number of Monte-Carlo samples.
    num_threads : int, default=1
        Accepted for parity with :meth:`pyvinecopulib.core.Vinecop.cdf`.
    seeds : list of int, or None, optional
        RNG seeds.

    Returns
    -------
    array, shape (m,), dtype float
        Distribution values in ``[0, 1]``.
    """

  @abstractmethod
  def rosenblatt(
    self, u: ArrayT, *, x: Optional[ArrayT] = None, num_threads: int = 1
  ) -> ArrayT:
    """Rosenblatt transform: dependent uniforms to independent uniforms.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Pseudo-observations in ``[0, 1]^d``.
    x : array, shape (n, p), or None, optional
        External covariates (non-simplified vines), else ``None``.
    num_threads : int, default=1
        Accepted for parity with :meth:`pyvinecopulib.core.Vinecop.rosenblatt`.

    Returns
    -------
    array, shape (n, d), dtype float
        Independent uniforms.
    """

  @abstractmethod
  def inverse_rosenblatt(
    self, u: ArrayT, *, x: Optional[ArrayT] = None, num_threads: int = 1
  ) -> ArrayT:
    """Inverse Rosenblatt transform: independent uniforms to dependent uniforms.

    Parameters
    ----------
    u : array, shape (n, d), dtype float
        Independent uniforms in ``[0, 1]^d``.
    x : array, shape (n, p), or None, optional
        External covariates (non-simplified vines), else ``None``.
    num_threads : int, default=1
        Accepted for parity with
        :meth:`pyvinecopulib.core.Vinecop.inverse_rosenblatt`.

    Returns
    -------
    array, shape (n, d), dtype float
        Dependent uniforms.
    """

  @abstractmethod
  def simulate(
    self,
    n: int,
    *,
    x: Optional[ArrayT] = None,
    qrng: bool = False,
    num_threads: int = 1,
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
    num_threads : int, default=1
        Accepted for parity with :meth:`pyvinecopulib.core.Vinecop.simulate`.
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

  A margin maps observations ``x`` on the original scale to a density
  (``pdf``), a distribution (``cdf``), and the inverse distribution
  (``icdf``) — enough to turn a vine copula into a full multivariate
  distribution and back.

  ``pdf`` is the density with respect to **the margin's own reference
  measure**: a Lebesgue density for a continuous margin, a probability mass at
  an atom, and whichever applies pointwise for a mixed (e.g. zero-inflated)
  one. That is what makes the Sklar factorization hold unchanged in all three
  cases, with no branch on the variable type::

      log f(x) = log c(F_1(x_1), ..., F_d(x_d)) + sum_j log pdf_j(x_j)

  :class:`pyvinecopulib.utils.Kde1d` satisfies this contract directly, in all
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
  - ``cdf_left`` — ``F(x^-)``, the left limit a margin with atoms needs;
    absent means it coincides with ``cdf``.
  - ``logpdf`` — absent means ``log(pdf)``.
  - ``simulate`` — absent means ``icdf`` of uniforms.
  - ``support`` — a ``(lo, hi)`` pair; absent means unbounded.

  See Also
  --------
  pyvinecopulib.core.MarginBase : Canonical partial implementation to subclass.
  pyvinecopulib.utils.Kde1d : The default nonparametric margin.
  BicopLike : The pair-copula contract.
  """

  @abstractmethod
  def pdf(self, x: ArrayT) -> ArrayT:
    """Density of the margin with respect to its own reference measure.

    Parameters
    ----------
    x : array, shape (n,), dtype float
        Observations on the original scale.

    Returns
    -------
    array, shape (n,), dtype float
        A density, a probability mass, or a mixture of the two, depending on
        the margin.
    """

  @abstractmethod
  def cdf(self, x: ArrayT) -> ArrayT:
    """Distribution function ``F(x)``, right-continuous.

    Parameters
    ----------
    x : array, shape (n,), dtype float
        Observations on the original scale.

    Returns
    -------
    array, shape (n,), dtype float
        Distribution values in ``[0, 1]``.
    """

  @abstractmethod
  def icdf(self, p: ArrayT) -> ArrayT:
    """Inverse distribution function ``inf{x : F(x) >= p}``.

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
