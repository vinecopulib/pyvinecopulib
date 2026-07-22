"""Backend-neutral pair-copula / vine contracts.

:class:`BicopLike` and :class:`VinecopLike` define what a pair copula / vine
evaluator must provide, independent of the array backend (NumPy or PyTorch).
They are the extension point for custom backends: implement the protocol — or,
far more easily, subclass the canonical :class:`~pyvinecopulib.core.BicopBase` /
:class:`~pyvinecopulib.core.VinecopBase`, which fill in most of it — and the
object plugs into the rest of the library (e.g. it can be hosted in a vine, or
consumed by the sklearn backend layer). The reference implementations are
:class:`pyvinecopulib.core.Bicop` / :class:`pyvinecopulib.core.Vinecop` (the
compiled default) and :class:`pyvinecopulib.torch.TorchBicop` /
:class:`pyvinecopulib.torch.TorchVinecop`.

**Conditioning.** Every method carries an optional trailing ``x`` — the
conditioning variables the copula / vine depends on (conditioning-set values
and/or external covariates), row-aligned with ``u``. Unconditional models leave
it ``None``; a conditional pair copula reads it. In a vine, each pair's ``x`` is
assembled per edge by a :class:`~pyvinecopulib.core.ConditioningContext`.

**Typing.** :data:`ArrayT` is an unbounded ``TypeVar`` carried only on these
public signatures, so a concrete backend (e.g. ``TorchBicop``) inherits precise
``torch.Tensor`` return types. The numeric implementations in
:mod:`~pyvinecopulib.core._base` operate on arrays as ``Any`` (the Array API
namespace ``array_api_compat`` is itself untyped).
"""

from __future__ import annotations

from abc import abstractmethod
from typing import Optional, Protocol, TypeVar, runtime_checkable

#: Array type a backend commits to (``numpy.ndarray`` | ``torch.Tensor`` | ...).
ArrayT = TypeVar("ArrayT")

__all__ = ["ArrayT", "BicopLike", "VinecopLike"]


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
  pyvinecopulib.core.Bicop : The compiled reference pair copula.
  VinecopLike : The vine-level evaluator contract.

  Examples
  --------
  A minimal independence pair on NumPy — implement only the three primitives and
  inherit ``hinv1`` / ``hinv2`` / ``simulate`` / ``loglik`` / ``plot`` from
  :class:`~pyvinecopulib.core.BicopBase`::

      import numpy as np
      from pyvinecopulib.core import BicopBase

      class Independence(BicopBase[np.ndarray]):
        def pdf(self, u, x=None):
          return np.ones(u.shape[0])

        def hfunc1(self, u, x=None):
          return u[:, 1]

        def hfunc2(self, u, x=None):
          return u[:, 0]

      cop = Independence()
      cop.hinv1(np.array([[0.3, 0.7]]))   # -> array([0.7]) (numerical inverse)
  """

  @abstractmethod
  def pdf(self, u: ArrayT, x: Optional[ArrayT] = None) -> ArrayT: ...

  @abstractmethod
  def cdf(self, u: ArrayT, x: Optional[ArrayT] = None) -> ArrayT: ...

  @abstractmethod
  def hfunc1(self, u: ArrayT, x: Optional[ArrayT] = None) -> ArrayT: ...

  @abstractmethod
  def hfunc2(self, u: ArrayT, x: Optional[ArrayT] = None) -> ArrayT: ...

  @abstractmethod
  def hinv1(self, u: ArrayT, x: Optional[ArrayT] = None) -> ArrayT: ...

  @abstractmethod
  def hinv2(self, u: ArrayT, x: Optional[ArrayT] = None) -> ArrayT: ...

  @abstractmethod
  def simulate(
    self,
    n: int,
    *,
    x: Optional[ArrayT] = None,
    qrng: bool = False,
    seeds: Optional[list[int]] = None,
  ) -> ArrayT: ...


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
  pyvinecopulib.core.Vinecop : The compiled reference vine.
  BicopLike : The pair-copula contract.
  """

  structure: object

  @abstractmethod
  def pdf(
    self, u: ArrayT, *, x: Optional[ArrayT] = None, num_threads: int = 1
  ) -> ArrayT: ...

  @abstractmethod
  def cdf(
    self,
    u: ArrayT,
    *,
    x: Optional[ArrayT] = None,
    N: int = 10000,
    seeds: Optional[list[int]] = None,
  ) -> ArrayT: ...

  @abstractmethod
  def rosenblatt(
    self, u: ArrayT, *, x: Optional[ArrayT] = None, num_threads: int = 1
  ) -> ArrayT: ...

  @abstractmethod
  def inverse_rosenblatt(
    self, u: ArrayT, *, x: Optional[ArrayT] = None, num_threads: int = 1
  ) -> ArrayT: ...

  @abstractmethod
  def simulate(
    self,
    n: int,
    *,
    x: Optional[ArrayT] = None,
    qrng: bool = False,
    seeds: Optional[list[int]] = None,
  ) -> ArrayT: ...
