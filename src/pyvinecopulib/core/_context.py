"""Conditioning-context policy for the vine cascades.

A :class:`ContextPolicy` decides, for a single pair-copula edge, the ``cond``
matrix the pair copula receives — assembled from the edge's conditioning-set
values ``u_D`` (gathered by the vine cascade) and an optional external covariate
matrix ``x`` (passed per call). It is the pluggable seam that turns the
otherwise simplified / unconditional cascades into non-simplified / conditional
ones.

**Column-order contract.** :meth:`ContextPolicy.edge_context` returns
``concat([u_D, x], axis=-1)`` with the ``u_D`` columns in the fixed order the
vine gathers them (ascending conditioning-tree index) and ``x`` appended last.
Conditional pair copulas consume ``cond`` positionally, so this order is a hard
contract that must be identical wherever ``cond`` is assembled (evaluation and
fitting alike).
"""

from __future__ import annotations

from typing import Optional, Protocol, runtime_checkable

from array_api_compat import array_namespace

from ._protocols import ArrayT

__all__ = ["ContextPolicy", "NonSimplifiedContext", "SimplifiedContext"]


@runtime_checkable
class ContextPolicy(Protocol[ArrayT]):
  """Per-edge conditioning-context assembler.

  Attributes
  ----------
  assembles_conditioning : bool
      Whether the cascade should gather the edge's conditioning-set values
      ``u_D`` before calling :meth:`edge_context`. When ``False`` the gather is
      skipped entirely (zero cost for simplified vines) and ``u_D`` is ``None``.
  """

  assembles_conditioning: bool

  def edge_context(
    self, *, u_D: Optional[ArrayT], x: Optional[ArrayT]
  ) -> Optional[ArrayT]:
    """Assemble the ``cond`` matrix for one edge.

    Parameters
    ----------
    u_D : array, shape (n, |D|), or None
        The edge's conditioning-set values (gathered by the cascade when
        :attr:`assembles_conditioning` is ``True``), else ``None``.
    x : array, shape (n, p), or None
        External covariates for this call, or ``None``.

    Returns
    -------
    array, shape (n, k), or None
        The concatenated conditioning matrix, or ``None`` for the
        unconditional (simplified) case.
    """
    ...


class SimplifiedContext:
  """Simplified vine: no conditioning, ``cond`` is always ``None``.

  The default policy; reproduces the classic simplified cascade at zero extra
  cost (the ``u_D`` gather is skipped).
  """

  assembles_conditioning: bool = False

  def edge_context(
    self, *, u_D: Optional[ArrayT] = None, x: Optional[ArrayT] = None
  ) -> None:
    return None


class NonSimplifiedContext:
  """Non-simplified vine: ``cond = concat([u_D, x])`` per edge.

  Parameters
  ----------
  use_conditioning : bool, default=True
      Whether to gather and include the edge's conditioning-set values ``u_D``.
      Set ``False`` for the covariate-only case (the pair copula depends on the
      external covariates ``x`` but not on ``u_D``), which lets the cascade skip
      the ``u_D`` gather.
  """

  def __init__(self, *, use_conditioning: bool = True) -> None:
    self.assembles_conditioning = use_conditioning

  def edge_context(
    self, *, u_D: Optional[ArrayT] = None, x: Optional[ArrayT] = None
  ) -> Optional[ArrayT]:
    parts = [t for t in (u_D, x) if t is not None]
    if not parts:
      return None
    if len(parts) == 1:
      return parts[0]
    return array_namespace(parts[0]).concat(parts, axis=-1)
