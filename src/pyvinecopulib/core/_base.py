"""Canonical partial implementations of the pair-copula / vine contracts.

:class:`BicopBase` is the canonical (array-backend-agnostic) implementation of
:class:`~pyvinecopulib.core._protocols.BicopLike`: a copula that supplies only
``log_pdf`` / ``hfunc1`` / ``hfunc2`` (plus the ``dtype`` / ``device`` seam)
inherits ``pdf`` (``exp ∘ log_pdf``) and ``hinv1`` / ``hinv2`` (numerical
inversion of the h-functions) for free; each is overridable when a native exact
form exists. (Standalone per-pair ``simulate`` is intentionally not provided —
it is not part of the :class:`BicopLike` contract and the vine samples via
``inverse_rosenblatt``, not per-pair sampling.)

Written against the Array API (:func:`array_api_compat.array_namespace`) so the
same code runs on numpy and torch. It is *not* portable to immutable backends
(the vine cascades use in-place scratch); such a backend implements
:class:`~pyvinecopulib.core._protocols.VinecopLike` directly. Per the
``pyvinecopulib.core`` typing policy, array values are handled as ``Any`` inside
these implementations (the generic ``ArrayT`` lives on the public signatures).
"""

from __future__ import annotations

from abc import ABC
from typing import Any, Optional, cast

from array_api_compat import array_namespace

from ._protocols import ArrayT, BicopLike
from ._rootfind import solve_increasing

__all__ = ["BicopBase"]


class BicopBase(BicopLike[ArrayT], ABC):
  """Canonical partial implementation of :class:`BicopLike` (numpy / torch).

  Subclasses must implement ``dtype`` / ``device`` and ``log_pdf`` / ``hfunc1``
  / ``hfunc2``; ``pdf`` / ``hinv1`` / ``hinv2`` have working defaults here
  (``cdf`` raises unless overridden). Set the class attribute
  ``supports_batched = True`` only if the concrete pair exposes the grid/cache
  internals the torch batched path needs.
  """

  #: Whether this pair can enter the torch grid-batched fast path.
  supports_batched: bool = False

  def pdf(self, u: ArrayT, cond: Optional[ArrayT] = None) -> ArrayT:
    logp = self.log_pdf(u, cond)
    xp = array_namespace(logp)
    return cast(ArrayT, xp.exp(logp))

  def hinv1(self, u: ArrayT, cond: Optional[ArrayT] = None) -> ArrayT:
    ua: Any = u
    xp = array_namespace(ua)
    u1, p = ua[:, 0], ua[:, 1]
    return cast(
      ArrayT,
      solve_increasing(
        lambda x: self.hfunc1(xp.stack([u1, x], axis=-1), cond), p
      ),
    )

  def hinv2(self, u: ArrayT, cond: Optional[ArrayT] = None) -> ArrayT:
    ua: Any = u
    xp = array_namespace(ua)
    p, u2 = ua[:, 0], ua[:, 1]
    return cast(
      ArrayT,
      solve_increasing(
        lambda x: self.hfunc2(xp.stack([x, u2], axis=-1), cond), p
      ),
    )

  def cdf(self, u: ArrayT, cond: Optional[ArrayT] = None) -> ArrayT:
    raise NotImplementedError(
      f"{type(self).__name__}.cdf is not defined; the vine cdf uses "
      "Monte-Carlo simulation and does not require a per-pair cdf."
    )
