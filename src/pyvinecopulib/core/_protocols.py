"""Backend-neutral pair-copula / vine Protocols (the extension contracts).

These generic, :func:`~typing.runtime_checkable` Protocols mirror the
nanobind-exposed :class:`~pyvinecopulib.Bicop` / :class:`~pyvinecopulib.Vinecop`
surfaces and are the extension point for custom pair-copula and vine backends
(torch is one; numpy another). Their members are ``@abstractmethod`` so the
Protocols double as subclassable bases for the canonical implementations in
:mod:`pyvinecopulib.core._base`.

The only extension over the C++ surface is an OPTIONAL trailing ``x`` giving the
conditioning variables (conditioning-set values and/or external covariates) a
copula depends on. It appears symmetrically on both Protocols: on each
:class:`BicopLike` method (per pair copula, assembled by the vine's context
policy) and on each :class:`VinecopLike` method (the per-call covariate matrix,
row-aligned with ``u``). Unconditional users ignore it; a simplified vine passes
``x=None`` down and each pair sees ``x=None``.

Typing policy for ``pyvinecopulib.core``: :data:`ArrayT` is an unbounded
``TypeVar`` used only on these public signatures (so a concrete backend such as
``TorchBicop`` inherits precise ``torch.Tensor`` returns). The neutral numeric
implementations in :mod:`~pyvinecopulib.core._base` /
:mod:`~pyvinecopulib.core._rootfind` operate on arrays as ``Any`` — the Array
API namespace (``array_api_compat``) is itself untyped, so array creation and
elementwise ops are ``Any``; parameters that must be indexed / have attributes
read are aliased to ``Any`` at point of use.
"""

from __future__ import annotations

from abc import abstractmethod
from typing import Optional, Protocol, TypeVar, runtime_checkable

#: Array type a backend commits to (``numpy.ndarray`` | ``torch.Tensor`` | ...).
ArrayT = TypeVar("ArrayT")

__all__ = ["ArrayT", "BicopLike", "VinecopLike"]


@runtime_checkable
class BicopLike(Protocol[ArrayT]):
  """Bivariate (optionally conditional) pair copula.

  Every method takes ``u`` of shape ``(n, 2)`` in the clamped domain
  ``[1e-10, 1 - 1e-10]`` and an OPTIONAL ``x`` of shape ``(n, k)`` giving the
  conditioning variables the copula depends on (conditioning-set values and/or
  external covariates). Unconditional copulas ignore ``x``. When the pair copula
  is hosted in a vine, ``x`` is assembled for it by the vine's context policy
  (a simplified vine passes ``x=None``, or forwards only external covariates).
  ``pdf`` is the density primitive (there is no ``log_pdf`` — mirroring the C++
  ``Bicop`` surface); ``dtype`` / ``device`` report the array precision and
  placement (replacing any interpolation-grid coupling).
  """

  @property
  @abstractmethod
  def dtype(self) -> object: ...

  @property
  @abstractmethod
  def device(self) -> object: ...

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


@runtime_checkable
class VinecopLike(Protocol[ArrayT]):
  """Post-fit vine evaluator; mirrors the nanobind ``Vinecop``.

  The optional keyword-only ``x`` covariate matrix (row-aligned with ``u``) is
  the conditional extension, kept symmetric with :class:`BicopLike`. It defaults
  to ``None``, so the unconditional ``Vinecop`` and any torch/numpy vine satisfy
  the Protocol whether or not they act on ``x``.
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
