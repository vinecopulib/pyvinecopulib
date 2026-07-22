"""Conditioning-context policy for the vine cascades.

A :class:`ConditioningContext` decides, for a single pair-copula edge, the ``x`` matrix
(``x_e`` in the cascades) the pair copula receives — assembled from the edge's
conditioning-set values ``u_D`` (gathered by the vine cascade) and the optional
external covariate matrix ``x`` (passed per call). It is the pluggable seam that
turns the otherwise simplified cascades into non-simplified ones:

* :class:`SimplifiedContext` (default) forwards only the external covariates
  ``x`` — the pair-copula parameters may depend on ``x`` but **not** on the
  conditioning set ``u_D`` (the classic simplified vine when ``x`` is ``None``).
* :class:`NonSimplifiedContext` forwards ``concat([u_D, x])`` — the pair copula
  ``c_{a,b;D}`` depends on its conditioning-set values (and optionally ``x``).

**Column-order contract (C1).** :meth:`ConditioningContext.edge_context` returns
``concat([u_D, x], axis=-1)`` with the ``u_D`` columns in the fixed order the
vine gathers them (ascending conditioning-tree index) and ``x`` appended last.
Conditional pair copulas consume this matrix positionally, so the order is a
hard contract that must be identical wherever it is assembled (evaluation and
fitting alike).

Examples
--------
Pass a policy to a vine (or to :meth:`VinecopBase.sequential_fit`) to control how
each pair's conditioning matrix ``x_e`` is built::

    from pyvinecopulib.core import NonSimplifiedContext, SimplifiedContext

    SimplifiedContext().edge_context(u_D=None, x=None)   # -> None (unconditional)
    ctx = NonSimplifiedContext()
    ctx.assembles_conditioning                           # -> True (u_D gathered)
"""

from __future__ import annotations

from abc import abstractmethod
from typing import Optional, Protocol, runtime_checkable

from array_api_compat import array_namespace

from .protocols import ArrayT

__all__ = ["ConditioningContext", "NonSimplifiedContext", "SimplifiedContext"]


@runtime_checkable
class ConditioningContext(Protocol[ArrayT]):
  """Per-edge conditioning-context assembler.

  Attributes
  ----------
  assembles_conditioning : bool
      Whether the cascade should gather the edge's conditioning-set values
      ``u_D`` before calling :meth:`edge_context`. When ``False`` the gather is
      skipped entirely (zero cost for simplified vines) and ``u_D`` is ``None``.
  """

  assembles_conditioning: bool

  @abstractmethod
  def edge_context(
    self, *, u_D: Optional[ArrayT], x: Optional[ArrayT]
  ) -> Optional[ArrayT]:
    """Assemble the per-edge conditioning matrix ``x_e``.

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
        The assembled conditioning matrix, or ``None`` when the pair copula is
        unconditional for this edge.
    """
    ...


class SimplifiedContext(ConditioningContext[ArrayT]):
  """Simplified vine: forward only the external covariates ``x``.

  The default policy. The pair-copula parameters may depend on external
  covariates ``x`` but not on the conditioning set ``u_D``, so the ``u_D``
  gather is skipped (zero extra cost). With ``x=None`` this reproduces the
  classic simplified / unconditional cascade (each pair sees ``x_e = None``).
  """

  assembles_conditioning: bool = False

  def edge_context(
    self, *, u_D: Optional[ArrayT] = None, x: Optional[ArrayT] = None
  ) -> Optional[ArrayT]:
    """Return the external covariates ``x`` unchanged (``u_D`` is ignored).

    Parameters
    ----------
    u_D : array or None, optional
        Ignored (never gathered under a simplified context).
    x : array, shape (n, p), or None, optional
        External covariates for this call.

    Returns
    -------
    array, shape (n, p), or None
        ``x`` unchanged.
    """
    return x


class NonSimplifiedContext(ConditioningContext[ArrayT]):
  """Non-simplified vine: ``x_e = concat([u_D, x])`` per edge.

  Gathers the edge's conditioning-set values ``u_D`` (in the C1 column order)
  and appends the external covariates ``x`` when present, so the pair copula
  ``c_{a,b;D}`` conditions on both.
  """

  assembles_conditioning: bool = True

  def edge_context(
    self, *, u_D: Optional[ArrayT] = None, x: Optional[ArrayT] = None
  ) -> Optional[ArrayT]:
    """Concatenate ``[u_D, x]`` (C1 order); drop whichever is ``None``.

    Parameters
    ----------
    u_D : array, shape (n, |D|), or None, optional
        The edge's conditioning-set values, in ascending conditioning-tree
        order.
    x : array, shape (n, p), or None, optional
        External covariates for this call, appended last.

    Returns
    -------
    array, shape (n, |D| + p), or None
        ``concat([u_D, x])``, a single one of them when the other is ``None``,
        or ``None`` when both are ``None``.
    """
    parts = [t for t in (u_D, x) if t is not None]
    if not parts:
      return None
    if len(parts) == 1:
      return parts[0]
    return array_namespace(parts[0]).concat(parts, axis=-1)
