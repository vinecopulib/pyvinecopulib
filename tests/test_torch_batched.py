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
  int_on_grid_batched,
  integrate_1d_batched,
  integrate_2d_batched,
  interpolate_batched,
  interpolate_batched_many,
)


@pytest.mark.parametrize("is_linear", [False, True])
def test_interpolate_batched_many_matches_separate(is_linear: bool) -> None:
  """The fused K-grid lookup is bit-identical to K separate calls.

  ``interpolate_batched_many`` shares the cell search / weights across the K
  stacked grids; ``[:, k, :]`` must equal ``interpolate_batched(caches[:, k])``
  exactly (pure gathers + the same weighted-sum arithmetic), which is what lets
  the fused TLL cascade keep the 1e-15 batched-parity gate.
  """
  n_pairs, k_grids, m, n = 3, 3, 12, 40
  gen = torch.Generator().manual_seed(0)
  grid = (
    torch.linspace(0.0, 1.0, m, dtype=torch.float64)
    if is_linear
    else torch.sort(torch.rand(m, generator=gen, dtype=torch.float64)).values
  )
  grid[0], grid[-1] = 0.0, 1.0
  caches = torch.rand(
    n_pairs, k_grids, m, m, generator=gen, dtype=torch.float64
  )
  u = torch.rand(n_pairs, n, 2, generator=gen, dtype=torch.float64)

  many = interpolate_batched_many(grid, caches, u, is_linear)
  for k in range(k_grids):
    sep = interpolate_batched(grid, caches[:, k], u, is_linear)
    torch.testing.assert_close(many[:, k, :], sep, atol=0.0, rtol=0.0)


def _fit_tll_bicop(seed: int) -> pv.Bicop:
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  u_fit = cop.sample(2000, seeds=[seed, seed + 1, seed + 2])
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
