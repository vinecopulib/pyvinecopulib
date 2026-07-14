"""Batched primitive equivalence against :mod:`pyvinecopulib.torch._interp`.

The batched primitives in :mod:`pyvinecopulib.torch._batched` are the
``(N, m, m)`` analogs of the operations in :mod:`._interp`. This file
asserts the two paths produce identical outputs when ``N`` independent
grids are fed through the batched version vs. an explicit Python loop
over per-pair :class:`InterpolationGrid2D` calls.
"""

from __future__ import annotations

import numpy as np
import pytest

import pyvinecopulib as pv

torch = pytest.importorskip("torch")

from pyvinecopulib.torch import TorchBicop  # noqa: E402
from pyvinecopulib.torch._batched import (  # noqa: E402
  _batched_cell_index,
  _build_cell_lookup,
  int_on_grid_batched,
  integrate_1d_batched,
  integrate_2d_batched,
  interpolate_batched,
)
from pyvinecopulib.torch._interp import InterpolationGrid2D  # noqa: E402


def _fit_tll_bicop(seed: int) -> pv.Bicop:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  u_fit = cop.simulate(2000, seeds=[seed, seed + 1, seed + 2])
  return pv.Bicop.from_data(
    u_fit,
    controls=pv.FitControlsBicop(family_set=[pv.families.tll], num_threads=1),
  )


def _eval_grid(n: int, seed: int) -> np.ndarray:
  rng = np.random.default_rng(seed)
  return rng.uniform(0.02, 0.98, size=(n, 2))


def test_interpolate_batched_matches_unbatched() -> None:
  """N=3 stacked grids vs three independent ``InterpolationGrid2D.interpolate``."""
  cops = [_fit_tll_bicop(seed) for seed in (10, 20, 30)]
  bcs = [TorchBicop.from_bicop(c) for c in cops]
  grid_points = bcs[0].interp_grid.grid_points
  values = torch.stack([bc.interp_grid.values for bc in bcs], dim=0)

  u_np = _eval_grid(500, seed=11)
  u_t = torch.from_numpy(u_np)
  u_batch = u_t.unsqueeze(0).expand(3, -1, -1).contiguous()

  out_batched = interpolate_batched(grid_points, values, u_batch).numpy()
  out_each = np.stack(
    [bc.interp_grid.interpolate(u_t).numpy() for bc in bcs], axis=0
  )
  np.testing.assert_allclose(out_batched, out_each, atol=1e-15, rtol=1e-15)


def test_int_on_grid_batched_matches_unbatched() -> None:
  """Batched trapezoidal integration vs the per-pair ``_int_on_grid``."""
  bc = TorchBicop.from_bicop(_fit_tll_bicop(42))
  grid_points = bc.interp_grid.grid_points
  m = grid_points.shape[0]
  rng = np.random.default_rng(7)
  N, n = 3, 200
  upr = torch.from_numpy(rng.uniform(0.05, 0.95, size=(N, n)))
  vals = torch.from_numpy(rng.uniform(0.1, 5.0, size=(N, n, m)))

  out = int_on_grid_batched(grid_points, upr, vals).numpy()
  ref = np.empty((N, n))
  for k in range(N):
    for ll in range(n):
      ref[k, ll] = bc.interp_grid._int_on_grid(
        upr[k, ll : ll + 1], vals[k, ll : ll + 1, :]
      ).item()
  np.testing.assert_allclose(out, ref, atol=1e-15, rtol=1e-15)


def test_integrate_1d_batched_matches_unbatched() -> None:
  """3 stacked grids vs three ``InterpolationGrid2D.integrate_1d`` calls."""
  cops = [_fit_tll_bicop(seed) for seed in (1, 2, 3)]
  bcs = [TorchBicop.from_bicop(c) for c in cops]
  grid_points = bcs[0].interp_grid.grid_points
  values = torch.stack([bc.interp_grid.values for bc in bcs], dim=0)

  u_np = _eval_grid(300, seed=15)
  u_t = torch.from_numpy(u_np)
  u_batch = u_t.unsqueeze(0).expand(3, -1, -1).contiguous()

  for cond_var in (1, 2):
    out = integrate_1d_batched(grid_points, values, u_batch, cond_var).numpy()
    ref = np.stack(
      [
        bc.interp_grid.integrate_1d(u_t, cond_var=cond_var).numpy()
        for bc in bcs
      ]
    )
    np.testing.assert_allclose(out, ref, atol=1e-13, rtol=1e-13)


def test_integrate_2d_batched_matches_unbatched() -> None:
  cops = [_fit_tll_bicop(seed) for seed in (4, 5, 6)]
  bcs = [TorchBicop.from_bicop(c) for c in cops]
  grid_points = bcs[0].interp_grid.grid_points
  values = torch.stack([bc.interp_grid.values for bc in bcs], dim=0)

  u_np = _eval_grid(300, seed=16)
  u_t = torch.from_numpy(u_np)
  u_batch = u_t.unsqueeze(0).expand(3, -1, -1).contiguous()

  out = integrate_2d_batched(grid_points, values, u_batch).numpy()
  ref = np.stack([bc.interp_grid.integrate_2d(u_t).numpy() for bc in bcs])
  np.testing.assert_allclose(out, ref, atol=1e-13, rtol=1e-13)


@pytest.mark.parametrize("m", [30, 64])
def test_bucket_cell_index_matches_searchsorted(m: int) -> None:
  """The vinecopulib#691 bucket path returns exactly the same cell index as
  the ``searchsorted`` fallback, on the Phi-spaced normal grid and on an
  arbitrary ascending grid (incl. exact-knot queries)."""
  normal = InterpolationGrid2D.make_grid_points("normal", m)
  # An arbitrary ascending grid on [0, 1] (endpoints pinned like the ctor).
  rng = np.random.default_rng(m)
  interior = np.sort(rng.uniform(0.01, 0.99, size=m - 2))
  arbitrary = torch.tensor([0.0, *interior.tolist(), 1.0], dtype=torch.float64)

  for grid_points in (normal, arbitrary):
    cell_lookup, max_advance = _build_cell_lookup(grid_points)
    # Random queries plus every exact knot (the boundary case where the
    # right=False convention and the guarded advance interact).
    u_rand = torch.from_numpy(rng.uniform(0.0, 1.0, size=(4, 500)))
    u = torch.cat([u_rand, grid_points.expand(4, m)], dim=1)

    ref = _batched_cell_index(grid_points, u)  # searchsorted path
    got = _batched_cell_index(
      grid_points, u, cell_lookup=cell_lookup, max_advance=max_advance
    )
    assert torch.equal(got, ref)


def test_cell_lookup_wired_grid_matches_searchsorted() -> None:
  """A grid with the bucket table attached (the ``CELL_LOOKUP_ENABLED``
  production path) produces byte-identical eval outputs to the default
  searchsorted grid — guards the opt-in wiring even though it is off by
  default."""
  cop = _fit_tll_bicop(seed=7)
  bc = TorchBicop.from_bicop(cop, cache_integrals=False)
  grid = bc.interp_grid  # searchsorted (cell_lookup is None by default)
  assert grid.cell_lookup is None

  # Same grid with the table attached (what CELL_LOOKUP_ENABLED builds).
  wired = TorchBicop.from_bicop(cop, cache_integrals=False).interp_grid
  wired.cell_lookup, wired.max_advance = _build_cell_lookup(wired.grid_points)

  u = torch.from_numpy(_eval_grid(400, seed=8))
  np.testing.assert_array_equal(
    grid.interpolate(u).numpy(), wired.interpolate(u).numpy()
  )
  for cond_var in (1, 2):
    np.testing.assert_array_equal(
      grid.integrate_1d(u, cond_var).numpy(),
      wired.integrate_1d(u, cond_var).numpy(),
    )
    np.testing.assert_array_equal(
      grid.inverse_integrate_1d(u, cond_var).numpy(),
      wired.inverse_integrate_1d(u, cond_var).numpy(),
    )
  np.testing.assert_array_equal(
    grid.integrate_2d(u).numpy(), wired.integrate_2d(u).numpy()
  )
