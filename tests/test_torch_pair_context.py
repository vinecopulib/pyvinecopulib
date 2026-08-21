"""Tests for the conditional / non-simplified vine path (on ``VinecopBase``).

The conditional / non-simplified capability lives on the array-agnostic
``VinecopBase`` (``ConditioningContext`` + ``x``-threaded cascades + the
``fit`` engine), not on ``TorchVinecop`` (which stays an
``nn.Module`` vine of ``TorchBicop`` pairs). These tests host a toy conditional
``GaussianBicop`` (correlation depends on ``x`` via a position-weighted link, so
it is genuinely non-simplified *and* sensitive to the C1 column order) in a
minimal ``VinecopBase`` subclass — the pattern a downstream package (e.g. npcc)
uses to build a conditional vine of scikit-style pairs. Covers:

* the master correctness property — ``inverse_rosenblatt(rosenblatt(u, x), x)
  ≈ u`` for a non-simplified vine;
* the C1 column-order contract is load-bearing (reversing ``u_D`` changes the
  density);
* the public ``VinecopBase.fit`` conditional-fit seam;
* array-namespace agnosticism (the same ``VinecopBase`` cascades match on numpy
  and torch);
* row alignment and the batched auto-fallback for a non-grid pair.
"""

from __future__ import annotations

from typing import Any, Optional

import numpy as np
import pytest
from array_api_compat import array_namespace

import pyvinecopulib as pv

torch = pytest.importorskip("torch")

from pyvinecopulib.core import (  # noqa: E402
  NonSimplifiedContext,
  SimplifiedContext,
  VinecopBase,
)

from tests.conftest import GaussianBicop, HostedVinecop  # noqa: E402

_SCALE = 0.6
_BASE_RHO = 0.3


def _pairs(d: int) -> list[list[GaussianBicop]]:
  """Nested ``[tree][edge]`` list of GaussianBicop pairs for a full d-vine."""
  return [
    [GaussianBicop(scale=_SCALE, base_rho=_BASE_RHO) for _ in range(d - 1 - t)]
    for t in range(d - 1)
  ]


def _vine(d: int, context) -> HostedVinecop:
  structure = pv.RVineStructure.from_order(list(range(1, d + 1)))
  return HostedVinecop(_pairs(d), structure, context=context)


def test_nonsimplified_bijection() -> None:
  """inverse_rosenblatt(rosenblatt(u, x), x) ≈ u for a conditional vine."""
  d, n, p = 4, 300, 2
  vine = _vine(d, NonSimplifiedContext())
  rng = np.random.default_rng(0)
  u = torch.as_tensor(rng.uniform(0.05, 0.95, (n, d)), dtype=torch.float64)
  x = torch.as_tensor(rng.standard_normal((n, p)), dtype=torch.float64)

  w = vine.rosenblatt(u, x=x)
  u_back = vine.inverse_rosenblatt(w, x=x)
  err = (u_back - u).abs().flatten()
  assert torch.quantile(err, 0.99).item() < 1e-8


def test_conditioning_changes_density() -> None:
  """A non-simplified vine's density differs from its simplified counterpart."""
  d, n, p = 4, 200, 2
  rng = np.random.default_rng(1)
  u = torch.as_tensor(rng.uniform(0.05, 0.95, (n, d)), dtype=torch.float64)
  x = torch.as_tensor(rng.standard_normal((n, p)), dtype=torch.float64)

  pdf_cond = _vine(d, NonSimplifiedContext()).pdf(u, x=x)
  pdf_simplified = _vine(d, SimplifiedContext()).pdf(u)  # x=None
  assert not torch.allclose(pdf_cond, pdf_simplified, atol=1e-6)


class _ReversedCondContext(NonSimplifiedContext):
  """NonSimplifiedContext that reverses the u_D column order (breaks C1)."""

  def edge_context(
    self, *, u_D: Optional[Any] = None, x: Optional[Any] = None
  ) -> Optional[Any]:
    if u_D is not None:
      xp = array_namespace(u_D)
      u_D = xp.flip(u_D, axis=-1)
    return NonSimplifiedContext.edge_context(self, u_D=u_D, x=x)


def test_c1_column_order_is_load_bearing() -> None:
  """Reversing the u_D column order changes the density (pins the C1 contract).

  At d=4 the single tree-2 edge conditions on ``|D| = 2`` variables; because
  ``GaussianBicop`` weights conditioning columns by position, reversing them
  changes ``rho`` and hence the density.
  """
  d, n = 4, 200
  rng = np.random.default_rng(2)
  # x=None so x_e is purely u_D (conditioning-set values).
  u = torch.as_tensor(rng.uniform(0.05, 0.95, (n, d)), dtype=torch.float64)

  pdf_c1 = _vine(d, NonSimplifiedContext()).pdf(u)
  pdf_rev = _vine(d, _ReversedCondContext()).pdf(u)
  assert not torch.allclose(pdf_c1, pdf_rev, atol=1e-6)


def test_fit_conditional_seam() -> None:
  """``VinecopBase.fit`` threads x_e (C1 widths) into a conditional fit.

  This is the public seam a downstream package drives to build a non-simplified
  vine: ``fit_edge(tree, edge, u_e, x_e) -> BicopLike`` receives ``x_e`` assembled
  per edge (conditioning-set values ``u_D`` in C1 order, then covariates ``x``).
  """
  d, n, p = 3, 400, 2
  rng = np.random.default_rng(3)
  u = torch.as_tensor(rng.uniform(0.05, 0.95, (n, d)), dtype=torch.float64)
  x = torch.as_tensor(rng.standard_normal((n, p)), dtype=torch.float64)
  structure = pv.RVineStructure.from_order([1, 2, 3])

  seen: list[tuple[int, int, Optional[int]]] = []

  def fit_edge(
    tree: int, edge: int, u_e: Any, x_e: Optional[Any]
  ) -> GaussianBicop:
    seen.append((tree, edge, None if x_e is None else x_e.shape[1]))
    return GaussianBicop(scale=_SCALE, base_rho=_BASE_RHO)

  pairs = VinecopBase.fit(
    structure, u, fit_edge, context=NonSimplifiedContext(), x=x
  )

  # Every edge saw an x_e in the C1 order: tree-0 edges have |D|=0 so width == p;
  # the tree-1 edge has |D|=1 so width == 1 + p.
  widths = {(t, e): w for (t, e, w) in seen}
  assert widths[(0, 0)] == p
  assert widths[(0, 1)] == p
  assert widths[(1, 0)] == 1 + p

  # The fitted pairs assemble into a working conditional vine.
  vine = HostedVinecop(pairs, structure, context=NonSimplifiedContext())
  u_back = vine.inverse_rosenblatt(vine.rosenblatt(u, x=x), x=x)
  err = (u_back - u).abs().flatten()
  assert torch.quantile(err, 0.99).item() < 1e-8


def test_numpy_and_torch_backends_match() -> None:
  """The same ``VinecopBase`` cascades agree on numpy and torch (agnostic)."""
  d, n, p = 4, 150, 2
  structure = pv.RVineStructure.from_order(list(range(1, d + 1)))
  rng = np.random.default_rng(7)
  u = rng.uniform(0.05, 0.95, (n, d))
  x = rng.standard_normal((n, p))

  np_vine = HostedVinecop(_pairs(d), structure, context=NonSimplifiedContext())
  torch_vine = HostedVinecop(
    _pairs(d), structure, context=NonSimplifiedContext()
  )
  u_t = torch.as_tensor(u, dtype=torch.float64)
  x_t = torch.as_tensor(x, dtype=torch.float64)

  np.testing.assert_allclose(
    np_vine.pdf(u, x=x),
    torch_vine.pdf(u_t, x=x_t).numpy(),
    atol=1e-10,
    rtol=1e-10,
  )
  np.testing.assert_allclose(
    np_vine.rosenblatt(u, x=x),
    torch_vine.rosenblatt(u_t, x=x_t).numpy(),
    atol=1e-10,
    rtol=1e-10,
  )


def test_row_permutation_alignment() -> None:
  """pdf(u[P], x[P]) == pdf(u, x)[P] — rows never leak across the cascade."""
  d, n, p = 4, 200, 2
  vine = _vine(d, NonSimplifiedContext())
  rng = np.random.default_rng(4)
  u = torch.as_tensor(rng.uniform(0.05, 0.95, (n, d)), dtype=torch.float64)
  x = torch.as_tensor(rng.standard_normal((n, p)), dtype=torch.float64)
  perm = torch.as_tensor(rng.permutation(n))

  base = vine.pdf(u, x=x)
  permuted = vine.pdf(u[perm], x=x[perm])
  torch.testing.assert_close(permuted, base[perm], atol=1e-12, rtol=0)


def test_batched_falls_back_for_non_grid_pair() -> None:
  """batched=True on a non-grid (GaussianBicop) vine falls back cleanly.

  ``HostedVinecop`` has no ``_build_batched`` override, so the base raises
  ``_NotBatchable`` and the dispatch layer transparently uses the non-batched
  cascade.
  """
  d, n = 4, 150
  vine = _vine(
    d, SimplifiedContext()
  )  # x=None so batched is otherwise eligible
  rng = np.random.default_rng(5)
  u = torch.as_tensor(rng.uniform(0.05, 0.95, (n, d)), dtype=torch.float64)

  out_batched = vine.pdf(u, batched=True)
  out_plain = vine.pdf(u, batched=False)
  torch.testing.assert_close(out_batched, out_plain, atol=1e-12, rtol=0)


def test_vinecopbase_loglik_matches_sum_log_pdf() -> None:
  """``VinecopBase.loglik`` equals ``sum(log(pdf))``."""
  d, n, p = 4, 120, 2
  vine = _vine(d, NonSimplifiedContext())
  rng = np.random.default_rng(6)
  u = torch.as_tensor(rng.uniform(0.05, 0.95, (n, d)), dtype=torch.float64)
  x = torch.as_tensor(rng.standard_normal((n, p)), dtype=torch.float64)
  ref = vine.pdf(u, x=x).clamp_min(1e-20).log().sum()
  torch.testing.assert_close(vine.loglik(u, x=x), ref, atol=1e-12, rtol=1e-12)


def test_vinecopbase_dim_and_repr() -> None:
  """The ``dim`` accessor and structural ``__repr__`` reflect the structure."""
  vine = _vine(4, SimplifiedContext())
  assert vine.dim == 4
  r = repr(vine)
  assert "dim=4" in r and "trunc_lvl=" in r and "order=" in r


def test_vinecopbase_plot_runs() -> None:
  """The inherited ``plot`` renders the tree structure without error (Agg)."""
  import matplotlib.pyplot as plt

  _vine(4, SimplifiedContext()).plot(layout="spring_layout")
  plt.close("all")
