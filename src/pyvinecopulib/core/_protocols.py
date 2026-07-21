"""Backend-neutral pair-copula / vine Protocols (the extension contracts).

These generic, :func:`~typing.runtime_checkable` Protocols mirror the
nanobind-exposed :class:`~pyvinecopulib.Bicop` / :class:`~pyvinecopulib.Vinecop`
surfaces and are the extension point for custom pair-copula and vine backends
(torch is one; numpy another). Their members are ``@abstractmethod`` so the
Protocols double as subclassable bases for the canonical implementations in
:mod:`pyvinecopulib.core._base`.

The only extension over the C++ surface is an OPTIONAL trailing ``cond`` (per
pair copula) — a conditioning matrix assembled by the vine's context policy —
and an optional per-call ``x`` covariate on the concrete vine (not on the
:class:`VinecopLike` Protocol, which mirrors the unconditional ``Vinecop``).

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
  ``[1e-10, 1 - 1e-10]`` and an OPTIONAL ``cond`` of shape ``(n, k)`` giving the
  conditioning variables the copula depends on (e.g. conditioning-set values
  and/or external covariates). Unconditional copulas ignore ``cond``. When the
  pair copula is hosted in a vine, ``cond`` is assembled for it by the vine's
  context policy (and a simplified vine simply passes ``cond=None``).
  ``dtype`` / ``device`` report the array precision and placement (replacing any
  interpolation-grid coupling).
  """

  @property
  @abstractmethod
  def dtype(self) -> object: ...

  @property
  @abstractmethod
  def device(self) -> object: ...

  @abstractmethod
  def log_pdf(self, u: ArrayT, cond: Optional[ArrayT] = None) -> ArrayT: ...

  @abstractmethod
  def pdf(self, u: ArrayT, cond: Optional[ArrayT] = None) -> ArrayT: ...

  @abstractmethod
  def cdf(self, u: ArrayT, cond: Optional[ArrayT] = None) -> ArrayT: ...

  @abstractmethod
  def hfunc1(self, u: ArrayT, cond: Optional[ArrayT] = None) -> ArrayT: ...

  @abstractmethod
  def hfunc2(self, u: ArrayT, cond: Optional[ArrayT] = None) -> ArrayT: ...

  @abstractmethod
  def hinv1(self, u: ArrayT, cond: Optional[ArrayT] = None) -> ArrayT: ...

  @abstractmethod
  def hinv2(self, u: ArrayT, cond: Optional[ArrayT] = None) -> ArrayT: ...


@runtime_checkable
class VinecopLike(Protocol[ArrayT]):
  """Post-fit vine evaluator; mirrors the nanobind ``Vinecop`` (unconditional).

  Concrete vines may add an optional per-call ``x`` covariate (the conditional
  extension); that is not part of this Protocol, which stays faithful to the
  ``Vinecop`` surface so both ``Vinecop`` and a torch/numpy vine satisfy it.
  """

  structure: object

  @abstractmethod
  def pdf(self, u: ArrayT, *, num_threads: int = 1) -> ArrayT: ...

  @abstractmethod
  def cdf(
    self, u: ArrayT, *, N: int = 10000, seeds: Optional[list[int]] = None
  ) -> ArrayT: ...

  @abstractmethod
  def rosenblatt(self, u: ArrayT, *, num_threads: int = 1) -> ArrayT: ...

  @abstractmethod
  def inverse_rosenblatt(
    self, u: ArrayT, *, num_threads: int = 1
  ) -> ArrayT: ...

  @abstractmethod
  def simulate(
    self, n: int, *, qrng: bool = False, seeds: Optional[list[int]] = None
  ) -> ArrayT: ...
