"""Bake-+-batched-eval equivalence with :class:`TorchBicop` across all four rotations.

The batched cascade in :class:`TorchVinecop` (see :mod:`pyvinecopulib.torch._batched`)
collapses the rotation-dispatch branches inside :class:`TorchBicop` by baking
each rotation into the stacked grid tensors at construction time. This test
file verifies that for each ``rotation in {0, 90, 180, 270}`` and each method
in ``{pdf, cdf, hfunc1, hfunc2}``:

  bake(grid, rotation) + batched_eval(u)  ==  TorchBicop(rotation).method(u)

at machine precision. The bake-then-evaluate path is used by the
``cache_integrals=True`` cascade (interp_at on baked caches + the linear
sign/offset/cdf_correction post-fold). The ``cache_integrals=False``
cascade uses a different rotate-at-eval-time path because trapezoidal
integration is not invariant under coordinate-permuting bakes (visits
different partial cells); that path is exercised by the higher-level
TorchVinecop tests in :mod:`test_torch_vinecop` against C++ pv.Vinecop.

The batched primitives themselves are also compared against their per-pair
counterparts in :mod:`_interp`.
"""

from __future__ import annotations

import numpy as np
import pytest

import pyvinecopulib as pv

torch = pytest.importorskip("torch")

from pyvinecopulib.torch import TorchBicop  # noqa: E402
from pyvinecopulib.torch._batched import (  # noqa: E402
  bake_pair_grid,
  bake_pair_metadata,
  int_on_grid_batched,
  integrate_1d_batched,
  integrate_2d_batched,
  interp_at_batched,
  interpolate_batched,
)

ROTATIONS = [0, 90, 180, 270]


def _fit_tll_bicop(seed: int) -> pv.Bicop:
  """Fit a TLL Bicop on Gaussian-sampled data — yields a non-trivial grid."""
  cop = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.6]]))
  u_fit = cop.simulate(2000, seeds=[seed, seed + 1, seed + 2])
  return pv.Bicop.from_data(
    u_fit,
    controls=pv.FitControlsBicop(family_set=[pv.families.tll], num_threads=1),
  )


def _eval_grid(n: int, seed: int) -> np.ndarray:
  rng = np.random.default_rng(seed)
  return rng.uniform(0.02, 0.98, size=(n, 2))


# --------------------------------------------------------------------------- #
# Batched-vs-unbatched primitive equivalence                                   #
# --------------------------------------------------------------------------- #


def test_interpolate_batched_matches_unbatched() -> None:
  """N=3 stacked grids vs three independent ``InterpolationGrid2D.interpolate``."""
  cops = [_fit_tll_bicop(seed) for seed in (10, 20, 30)]
  bcs = [TorchBicop.from_bicop(c) for c in cops]
  grid_points = bcs[0].interp_grid.grid_points
  values = torch.stack(
    [bc.interp_grid.values for bc in bcs], dim=0
  )  # (3, m, m)

  u_np = _eval_grid(500, seed=11)
  u_t = torch.from_numpy(u_np)
  u_batch = u_t.unsqueeze(0).expand(3, -1, -1).contiguous()  # (3, 500, 2)

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


# --------------------------------------------------------------------------- #
# Bake + batched-eval equivalence with TorchBicop, across all rotations         #
# --------------------------------------------------------------------------- #


def _rotated_bicop(base: pv.Bicop, rotation: int, cache: bool) -> TorchBicop:
  """Build a TorchBicop with the given rotation, sharing the underlying grid."""
  bc_zero = TorchBicop.from_bicop(base, cache_integrals=cache)
  return TorchBicop(
    grid_points=bc_zero.interp_grid.grid_points,
    values=bc_zero.interp_grid.values,
    rotation=rotation,
    cache_integrals=cache,
    norm_times=0,
  )


def _bake_single(
  bc: TorchBicop, rotation: int, cache: bool
) -> tuple[dict, dict]:
  """Bake the source TorchBicop's grids for a rotation; return (grids, meta)."""
  grids = bake_pair_grid(
    values=bc.interp_grid.values,
    rotation=rotation,
    h1_cache=bc._hfunc1_cache if cache else None,
    h2_cache=bc._hfunc2_cache if cache else None,
    cdf_cache=bc._cdf_cache if cache else None,
  )
  meta = bake_pair_metadata(rotation, dtype=bc.interp_grid.values.dtype)
  return grids, meta


@pytest.mark.parametrize("rotation", ROTATIONS)
def test_baked_pdf_matches_torchbicop(rotation: int) -> None:
  cache = True  # pdf bake is identical for both cache modes (no integration).
  base = _fit_tll_bicop(50 + rotation)
  bc_zero = TorchBicop.from_bicop(base, cache_integrals=cache)
  bc_rot = _rotated_bicop(base, rotation, cache)
  grid_points = bc_zero.interp_grid.grid_points

  grids, _ = _bake_single(bc_zero, rotation, cache)
  values_baked = grids["values_baked"].unsqueeze(0)  # (1, m, m)

  u_np = _eval_grid(400, seed=110 + rotation)
  u_t = torch.from_numpy(u_np)
  u_batch = u_t.unsqueeze(0)  # (1, 400, 2)

  expected = bc_rot.pdf(u_t).numpy()
  got = interpolate_batched(grid_points, values_baked, u_batch).clamp_min(1e-20)
  np.testing.assert_allclose(got.numpy()[0], expected, atol=1e-13, rtol=1e-13)


@pytest.mark.parametrize("rotation", ROTATIONS)
def test_baked_cdf_cache_true_matches_torchbicop(rotation: int) -> None:
  base = _fit_tll_bicop(60 + rotation)
  bc_zero = TorchBicop.from_bicop(base, cache_integrals=True)
  bc_rot = _rotated_bicop(base, rotation, cache=True)
  grid_points = bc_zero.interp_grid.grid_points

  grids, meta = _bake_single(bc_zero, rotation, cache=True)
  cdf_baked = grids["cdf_baked"].unsqueeze(0)
  corr = meta["cdf_correction"]  # (4,)

  u_np = _eval_grid(400, seed=120 + rotation)
  u_t = torch.from_numpy(u_np)
  u_batch = u_t.unsqueeze(0)

  p = interp_at_batched(grid_points, cdf_baked, u_batch)  # (1, n)
  got = corr[0] * p[0] + corr[1] * u_t[:, 0] + corr[2] * u_t[:, 1] + corr[3]
  expected = bc_rot.cdf(u_t).numpy()
  np.testing.assert_allclose(
    got.clamp(1e-10, 1.0 - 1e-10).numpy(), expected, atol=1e-13, rtol=1e-13
  )


@pytest.mark.parametrize("rotation", ROTATIONS)
@pytest.mark.parametrize("which", [1, 2])
def test_baked_hfunc_cache_true_matches_torchbicop(
  rotation: int, which: int
) -> None:
  base = _fit_tll_bicop(80 + rotation)
  bc_zero = TorchBicop.from_bicop(base, cache_integrals=True)
  bc_rot = _rotated_bicop(base, rotation, cache=True)
  grid_points = bc_zero.interp_grid.grid_points

  grids, meta = _bake_single(bc_zero, rotation, cache=True)
  cache_b = grids[f"h{which}_baked"].unsqueeze(0)
  sign = meta[f"h{which}_sign"]
  offset = meta[f"h{which}_offset"]

  u_np = _eval_grid(400, seed=140 + rotation + 17 * which)
  u_t = torch.from_numpy(u_np)
  u_batch = u_t.unsqueeze(0)

  raw = interp_at_batched(grid_points, cache_b, u_batch)  # (1, n)
  got = (sign * raw[0] + offset).clamp(0.0, 1.0)
  expected = (
    bc_rot.hfunc1(u_t).numpy() if which == 1 else bc_rot.hfunc2(u_t).numpy()
  )
  np.testing.assert_allclose(got.numpy(), expected, atol=1e-13, rtol=1e-13)
