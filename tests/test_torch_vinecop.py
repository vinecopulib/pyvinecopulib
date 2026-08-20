"""Tests for the optional ``pyvinecopulib.torch.TorchVinecop`` class.

Skipped when PyTorch isn't installed. Compares the torch ``TorchVinecop``
against the C++ ``pv.Vinecop`` (with TLL pair copulas) for ``pdf`` /
``rosenblatt`` / ``inverse_rosenblatt``; verifies the round-trip
``inverse_rosenblatt(rosenblatt(u)) ≈ u``; spot-checks an
independent vine; and confirms the truncation and discrete paths.
"""

from __future__ import annotations

import pickle

import numpy as np
import pytest

import pyvinecopulib as pv

torch = pytest.importorskip("torch")

from pyvinecopulib.torch import FitControlsTorchVinecop, TorchVinecop  # noqa: E402

_TLL_CONTROLS = pv.FitControlsVinecop(
  family_set=[pv.families.tll], num_threads=1
)


def _fit_tll_vine(u: np.ndarray, **extra) -> pv.Vinecop:
  ctl = pv.FitControlsVinecop(
    family_set=[pv.families.tll], num_threads=1, **extra
  )
  return pv.Vinecop.from_data(u, controls=ctl)


def _eval_grid(n: int, d: int, seed: int = 0) -> np.ndarray:
  rng = np.random.default_rng(seed)
  return rng.uniform(0.02, 0.98, size=(n, d))


def _simulate(d: int, n: int, seed: int = 0) -> np.ndarray:
  # Use a non-trivial structure: Gaussian pair-copula vine via Clayton
  # samples per pair would require building a parametric Vinecop. Easier:
  # take pseudo-observations from a multivariate normal with banded
  # correlation, which gives a smoothly varying density per pair.
  rng = np.random.default_rng(seed)
  base = rng.standard_normal(size=(n, 1))
  noise = rng.standard_normal(size=(n, d))
  x = 0.6 * base + 0.4 * noise
  return pv.to_pseudo_obs(x)


def test_pdf_matches_pvvinecop() -> None:
  u_fit = _simulate(d=5, n=2000, seed=1)
  cop_tll = _fit_tll_vine(u_fit)

  # Pin cache=False — this test verifies vine-level parity with the C++
  # on-the-fly cascade at 1e-10. The default cache_integrals=True uses
  # bilinear interp; its precision floor is pinned separately in
  # test_cache_integrals_precision_floor.
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=False)
  u_eval = _eval_grid(500, d=5, seed=11)

  out_torch = bc.pdf(torch.from_numpy(u_eval)).numpy()
  out_cpp = cop_tll.pdf(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


def test_rosenblatt_matches_pvvinecop() -> None:
  u_fit = _simulate(d=5, n=2000, seed=2)
  cop_tll = _fit_tll_vine(u_fit)

  # Pin cache=False — see test_pdf_matches_pvvinecop.
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=False)
  u_eval = _eval_grid(500, d=5, seed=12)

  out_torch = bc.rosenblatt(torch.from_numpy(u_eval)).numpy()
  out_cpp = cop_tll.rosenblatt(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


@pytest.mark.parametrize(
  "method,cache,max_abs,mean_abs",
  [
    # cache=False is bisection-precise at every step (forward cascades are
    # exact trapezoidal). Vine-level error stays at machine precision for
    # pdf/rosenblatt and at bisection convergence for inverse_rosenblatt.
    ("pdf", False, 1e-11, 1e-13),
    ("rosenblatt", False, 1e-10, 1e-13),
    ("inverse_rosenblatt", False, 1e-8, 1e-9),
    # cache=True replaces every integration / inversion with a bilinear
    # interp on a 30x30 grid — the per-pair ~1e-3 mean error compounds
    # through the d-1 cascade levels. pdf is unbounded so the *max* error
    # in absolute terms can be O(1); rosenblatt / inverse_rosenblatt are
    # u-space and bounded so their max stays around the bicop-level cap.
    # These bounds are loose — they pin the floor, not optimize against it.
    ("pdf", True, 1e1, 1e-1),
    ("rosenblatt", True, 1.0, 5e-2),
    ("inverse_rosenblatt", True, 5e-2, 5e-3),
  ],
)
def test_cache_integrals_precision_floor(
  method: str, cache: bool, max_abs: float, mean_abs: float
) -> None:
  """Pin the precision floor for each cache mode at a moderately-deep vine.

  ``cache_integrals=False`` is the right choice for likelihood / sampling
  applications: every per-pair integration and inversion is exact (modulo
  bisection convergence), so the cascade preserves precision. Setting
  ``cache=True`` is the right choice when ~1e-3 u-space error is
  acceptable and per-call speed matters — typically pdf-only workloads
  where the bilinear-interp gap is small compared to downstream noise.
  """
  d = 10
  u_fit = _simulate(d=d, n=2000, seed=1)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=cache)
  u_eval = _eval_grid(500, d=d, seed=2)
  out_cpp = getattr(cop_tll, method)(u_eval)
  out_torch = getattr(bc, method)(torch.from_numpy(u_eval)).numpy()
  diff = np.abs(out_torch - out_cpp)
  assert diff.max() < max_abs, f"{method}/{cache}: max {diff.max():.2e}"
  assert diff.mean() < mean_abs, f"{method}/{cache}: mean {diff.mean():.2e}"


def test_inverse_rosenblatt_matches_pvvinecop() -> None:
  u_fit = _simulate(d=5, n=2000, seed=3)
  cop_tll = _fit_tll_vine(u_fit)

  # Pin cache=False — C++ parity at 1e-9 via the bisection-precise
  # on-the-fly cascade. The cached path's precision floor is checked
  # by test_cache_integrals_precision_floor.
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=False)
  # Inverse-rosenblatt takes independent uniforms; sample fresh.
  w = _eval_grid(400, d=5, seed=13)

  out_torch = bc.inverse_rosenblatt(torch.from_numpy(w)).numpy()
  out_cpp = cop_tll.inverse_rosenblatt(w)
  # Bisection at n_iter=35 gives bracket-width accuracy 0.5**35 ≈ 6e-11.
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-9, rtol=1e-9)


def test_inverse_rosenblatt_roundtrip() -> None:
  u_fit = _simulate(d=5, n=2000, seed=4)
  cop_tll = _fit_tll_vine(u_fit)

  # Pin cache=False so the per-pair forward/inverse cascade is
  # bisection-precise; the round-trip identity requires the on-the-fly
  # path to hit the 99-percentile tolerance below.
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=False)
  rng = np.random.default_rng(14)
  u_eval = rng.uniform(0.1, 0.9, size=(400, 5))
  u_t = torch.from_numpy(u_eval)
  back = bc.inverse_rosenblatt(bc.rosenblatt(u_t)).numpy()
  diff = np.abs(back - u_eval)
  # ITP makes per-call precision ~1e-15. The cascade's irreducible noise
  # is the [TRIM_LO, TRIM_HI] = [1e-10, 1-1e-10] clamp inside hfunc /
  # hinv, which creates plateaus the solver cannot resolve at the
  # boundaries. The bulk of u_eval roundtrips at machine precision; a
  # tiny fraction (~0.1%) of points whose forward pass saturates against
  # the clamp don't, so we measure the 99th percentile rather than the
  # max.
  assert np.percentile(diff, 99) < 1e-10, (
    f"99th percentile diff = {np.percentile(diff, 99):.2e}; "
    f"max = {diff.max():.2e}"
  )


@pytest.mark.parametrize("batched", [False, True])
def test_independent_vine_short_circuits(batched: bool) -> None:
  d = 4
  structure = pv.RVineStructure.from_order([1, 2, 3, 4])
  # Build a vine with all-independent pair copulas.
  pair_copulas = [
    [pv.Bicop(family=pv.families.indep) for _ in range(d - 1 - t)]
    for t in range(d - 1)
  ]
  cop = pv.Vinecop.from_structure(
    structure=structure, pair_copulas=pair_copulas
  )

  # cache_integrals doesn't matter for the all-independence short-circuit
  # (every pair short-circuits to pdf=1, hfunc=u, etc.), but we pin
  # cache=False so the 1e-10 / 1e-12 tolerances below hold under either
  # default.
  bc = TorchVinecop.from_vinecop(cop, cache_integrals=False)
  u = _eval_grid(200, d=d, seed=21)
  u_t = torch.from_numpy(u)

  np.testing.assert_allclose(
    bc.pdf(u_t, batched=batched).numpy(), np.ones(200), atol=1e-12, rtol=0
  )
  np.testing.assert_allclose(
    bc.rosenblatt(u_t, batched=batched).numpy(), u, atol=1e-10, rtol=1e-10
  )
  # inverse_rosenblatt with batched=True now raises; only check batched=False.
  if not batched:
    np.testing.assert_allclose(
      bc.inverse_rosenblatt(u_t, batched=batched).numpy(),
      u,
      atol=1e-10,
      rtol=1e-10,
    )


@pytest.mark.parametrize("batched", [False, True])
def test_truncated_vine_pdf(batched: bool) -> None:
  u_fit = _simulate(d=6, n=2000, seed=5)
  cop_tll = _fit_tll_vine(u_fit, trunc_lvl=2)
  assert cop_tll.structure.trunc_lvl == 2

  # cache=False for 1e-10 C++ parity.
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=False)
  assert bc.trunc_lvl == 2

  u_eval = _eval_grid(300, d=6, seed=15)
  np.testing.assert_allclose(
    bc.pdf(torch.from_numpy(u_eval), batched=batched).numpy(),
    cop_tll.pdf(u_eval),
    atol=1e-10,
    rtol=1e-10,
  )


def test_from_vinecop_rejects_unsupported_family() -> None:
  d = 3
  structure = pv.RVineStructure.from_order([1, 2, 3])
  pair_copulas = [
    [
      pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.5]]))
      for _ in range(d - 1 - t)
    ]
    for t in range(d - 1)
  ]
  cop = pv.Vinecop.from_structure(
    structure=structure, pair_copulas=pair_copulas
  )
  with pytest.raises(ValueError, match="tll and indep"):
    TorchVinecop.from_vinecop(cop)


def test_pdf_shape_validation() -> None:
  u_fit = _simulate(d=4, n=1000, seed=6)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll)

  with pytest.raises(ValueError, match=r"shape \(n, 4\)"):
    bc.pdf(torch.zeros(10, 5, dtype=torch.float64))


# --------------------------------------------------------------------- #
# Batched cascade                                                        #
# --------------------------------------------------------------------- #


@pytest.mark.parametrize("cache", [False, True])
def test_pdf_batched_matches_unbatched(cache: bool) -> None:
  u_fit = _simulate(d=10, n=800, seed=300)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=cache)
  u_t = torch.from_numpy(_eval_grid(500, d=10, seed=310))
  ref = bc.pdf(u_t, batched=False).numpy()
  got = bc.pdf(u_t, batched=True).numpy()
  np.testing.assert_allclose(got, ref, atol=1e-12, rtol=1e-12)


@pytest.mark.parametrize("cache", [False, True])
def test_rosenblatt_batched_matches_unbatched(cache: bool) -> None:
  u_fit = _simulate(d=10, n=800, seed=301)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=cache)
  u_t = torch.from_numpy(_eval_grid(500, d=10, seed=311))
  ref = bc.rosenblatt(u_t, batched=False).numpy()
  got = bc.rosenblatt(u_t, batched=True).numpy()
  np.testing.assert_allclose(got, ref, atol=1e-13, rtol=1e-13)


def test_inverse_rosenblatt_batched_raises() -> None:
  """``batched=True`` on inverse_rosenblatt raises NotImplementedError.

  The inverse cascade's dependency graph is genuinely 2-D (some
  iterations depend on values at a *different* tree level), so the
  per-tree-level wavefront that works for pdf / rosenblatt doesn't
  apply. Raising surfaces the limitation explicitly instead of silently
  routing to the slower non-batched path.
  """
  u_fit = _simulate(d=6, n=600, seed=302)
  cop_tll = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop_tll)
  w_t = torch.from_numpy(_eval_grid(300, d=6, seed=312))
  with pytest.raises(NotImplementedError, match="batched=True"):
    bc.inverse_rosenblatt(w_t, batched=True)


def test_batched_matches_cpp_pdf() -> None:
  """End-to-end: torch batched pdf vs pv.Vinecop on the same fit."""
  u_fit = _simulate(d=8, n=1000, seed=320)
  cop_tll = _fit_tll_vine(u_fit)
  # cache=False for 1e-10 C++ parity.
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=False)
  u_eval = _eval_grid(400, d=8, seed=330)
  out_torch = bc.pdf(torch.from_numpy(u_eval), batched=True).numpy()
  out_cpp = cop_tll.pdf(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


def test_batched_matches_cpp_rosenblatt() -> None:
  u_fit = _simulate(d=8, n=1000, seed=321)
  cop_tll = _fit_tll_vine(u_fit)
  # cache=False for 1e-10 C++ parity.
  bc = TorchVinecop.from_vinecop(cop_tll, cache_integrals=False)
  u_eval = _eval_grid(400, d=8, seed=331)
  out_torch = bc.rosenblatt(torch.from_numpy(u_eval), batched=True).numpy()
  out_cpp = cop_tll.rosenblatt(u_eval)
  np.testing.assert_allclose(out_torch, out_cpp, atol=1e-10, rtol=1e-10)


@pytest.mark.parametrize("d", [5, 10])
def test_from_data_matches_from_vinecop(d: int) -> None:
  """``TorchVinecop.from_data(u, structure)`` should agree with
  ``pv.Vinecop.from_data(u, controls=..., structure=structure)`` to machine
  precision. Note: this is the *fixed-structure* C++ path. Freely-selected
  C++ Vinecop.from_data produces a different pair-copula orientation per
  edge (the graph-traversal order during selection differs from the order
  used when the structure is passed in), so the comparison target is the
  fixed-structure refit.
  """
  u_fit = _simulate(d=d, n=1500, seed=400 + d)
  # First, get a structure (any TLL fit will do).
  structure = _fit_tll_vine(u_fit).structure
  # Now refit C++ with that structure passed in.
  cop_cpp_fixed = pv.Vinecop.from_data(
    u_fit, controls=_TLL_CONTROLS, structure=structure
  )

  # Pin cache=False so the per-pair cascade is bisection-precise — needed
  # to match at 1e-9 (the cached path round-trips only to ~1e-3).
  bc_cpp = TorchVinecop.from_vinecop(cop_cpp_fixed, cache_integrals=False)
  bc_torch = TorchVinecop.from_data(
    torch.from_numpy(u_fit),
    structure,
    controls=FitControlsTorchVinecop(cache_integrals=False),
  )

  u_eval = torch.from_numpy(_eval_grid(400, d=d, seed=410 + d))
  np.testing.assert_allclose(
    bc_torch.pdf(u_eval).numpy(),
    bc_cpp.pdf(u_eval).numpy(),
    atol=1e-9,
    rtol=1e-9,
  )
  np.testing.assert_allclose(
    bc_torch.rosenblatt(u_eval).numpy(),
    bc_cpp.rosenblatt(u_eval).numpy(),
    atol=1e-9,
    rtol=1e-9,
  )
  np.testing.assert_allclose(
    bc_torch.inverse_rosenblatt(u_eval).numpy(),
    bc_cpp.inverse_rosenblatt(u_eval).numpy(),
    atol=1e-9,
    rtol=1e-9,
  )


def test_from_data_auto_selects_structure() -> None:
  """``from_data(structure=None)`` selects the R-vine natively in torch
  (``VinecopBase.select``, TLL, ``trunc_lvl=20``) and reuses the
  selection-time pairs — an exact port of ``pv.Vinecop.from_data``, so the
  selected structure (identical matrix encoding), the reoriented pairs, and
  the resulting density must all match."""
  u_fit = _simulate(d=5, n=1500, seed=777)
  cop_cpp = pv.Vinecop.from_data(
    u_fit,
    structure=None,
    var_types=[],
    controls=pv.FitControlsVinecop(
      family_set=[pv.families.tll], trunc_lvl=20, num_threads=1
    ),
  )
  got = TorchVinecop.from_data(
    torch.from_numpy(u_fit),
    controls=FitControlsTorchVinecop(cache_integrals=False),
  )
  # Identical structure, down to the matrix encoding.
  np.testing.assert_array_equal(
    np.asarray(got.structure.matrix), np.asarray(cop_cpp.structure.matrix)
  )
  # Identical density: same selection, same reused pairs (torch TLL matches
  # the ``Bicop`` TLL fit to machine precision), same flips.
  reference = TorchVinecop.from_vinecop(cop_cpp, cache_integrals=False)
  u_eval = torch.from_numpy(_eval_grid(400, d=5, seed=778))
  np.testing.assert_allclose(
    got.pdf(u_eval).numpy(),
    reference.pdf(u_eval).numpy(),
    atol=1e-9,
    rtol=1e-9,
  )


def test_batched_to_device_invalidates() -> None:
  """``.to()`` should drop the lazily-built BatchedVine so the next
  batched call rebuilds it on the new device."""
  u_fit = _simulate(d=5, n=500, seed=322)
  bc = TorchVinecop.from_vinecop(_fit_tll_vine(u_fit))
  u_t = torch.from_numpy(_eval_grid(50, d=5, seed=332))
  bc.pdf(u_t, batched=True)  # build the BatchedVine
  assert bc._batched is not None
  bc.to(torch.float32)
  assert bc._batched is None
  # And the re-build on the new dtype produces a finite output.
  out = bc.pdf(u_t.to(torch.float32), batched=True)
  assert torch.isfinite(out).all()


# --------------------------------------------------------------------------- #
# User-facing input validation for TorchVinecop                                #
# --------------------------------------------------------------------------- #


def test_from_data_rejects_non_2d_u() -> None:
  """``TorchVinecop.from_data`` rejects 1-D / 3-D inputs."""
  structure = pv.RVineStructure.sample(5, seeds=[1])
  with pytest.raises(ValueError, match="must be 2-D"):
    TorchVinecop.from_data(torch.zeros(100, dtype=torch.float64), structure)
  with pytest.raises(ValueError, match="must be 2-D"):
    TorchVinecop.from_data(
      torch.zeros(100, 5, 2, dtype=torch.float64), structure
    )


def test_from_data_rejects_structure_dim_mismatch() -> None:
  """``TorchVinecop.from_data`` rejects ``structure.dim != u.shape[1]``."""
  structure = pv.RVineStructure.sample(5, seeds=[1])  # d=5
  with pytest.raises(ValueError, match="does not match"):
    TorchVinecop.from_data(torch.zeros(100, 7, dtype=torch.float64), structure)


@pytest.mark.parametrize("op", ["pdf", "rosenblatt", "inverse_rosenblatt"])
def test_eval_rejects_wrong_input_shape(op: str) -> None:
  """Every eval entry point requires ``(n, d)`` and raises otherwise."""
  d = 5
  u_fit = _simulate(d=d, n=400, seed=520)
  bc = TorchVinecop.from_vinecop(_fit_tll_vine(u_fit))
  fn = getattr(bc, op)
  with pytest.raises(ValueError):
    fn(torch.zeros(100, dtype=torch.float64))  # 1-D
  with pytest.raises(ValueError):
    fn(torch.zeros(100, d + 1, dtype=torch.float64))  # wrong second dim


# ---------------------------------------------------------------------------
# API alignment with pv.Vinecop (from_structure, sample, cdf, num_threads)
# ---------------------------------------------------------------------------


def test_from_structure_pair_copulas_matches_init():
  """from_structure(structure, pair_copulas) builds the same vine as
  TorchVinecop(pair_copulas, structure)."""
  rng = np.random.default_rng(0)
  U = rng.uniform(0.001, 0.999, (300, 3))
  cop = _fit_tll_vine(U)
  from pyvinecopulib.torch import TorchBicop

  pcs = [[TorchBicop.from_bicop(b) for b in row] for row in cop.pair_copulas]
  via_init = TorchVinecop(pair_copulas=pcs, structure=cop.structure)
  via_factory = TorchVinecop.from_structure(
    structure=cop.structure, pair_copulas=pcs
  )

  u_t = torch.from_numpy(U[:10])
  assert torch.allclose(via_init.pdf(u_t), via_factory.pdf(u_t))


def test_from_structure_empty_is_independence_vine():
  """No pair_copulas -> every edge is an independence TorchBicop;
  pdf is identically 1."""
  s = pv.RVineStructure.sample(4, seeds=[1, 2, 3, 4, 5])
  tv = TorchVinecop.from_structure(structure=s)
  u = torch.rand(50, 4, dtype=torch.float64)
  assert torch.allclose(
    tv.pdf(u), torch.ones(50, dtype=torch.float64), atol=1e-12
  )


def test_from_structure_matrix_path():
  """Matrix kwarg routes through RVineStructure.from_matrix."""
  s = pv.RVineStructure.sample(3, seeds=[1, 2, 3, 4, 5])
  mat = np.asarray(s.matrix, dtype=np.uint64)
  tv = TorchVinecop.from_structure(matrix=mat)
  assert tv.d == 3
  assert tv.trunc_lvl == 2


def test_from_structure_rejects_both_or_neither():
  s = pv.RVineStructure.sample(3, seeds=[1, 2, 3, 4, 5])
  with pytest.raises(ValueError, match="exactly one"):
    TorchVinecop.from_structure()
  with pytest.raises(ValueError, match="exactly one"):
    TorchVinecop.from_structure(
      structure=s, matrix=np.asarray(s.matrix, dtype=np.uint64)
    )


def test_simulate_seeded_reproducible():
  rng = np.random.default_rng(0)
  U = rng.uniform(0.001, 0.999, (300, 3))
  tv = TorchVinecop.from_vinecop(_fit_tll_vine(U))

  a = tv.sample(n=200, seeds=[42])
  b = tv.sample(n=200, seeds=[42])
  assert torch.allclose(a, b)
  c = tv.sample(n=200, seeds=[43])
  assert not torch.allclose(a, c)


def test_simulate_qrng_runs_and_returns_uniforms():
  """qrng=True draws Halton/Sobol via pv.utils.sample_uniform, then
  pushes through inverse_rosenblatt. We check basic shape + range."""
  rng = np.random.default_rng(0)
  U = rng.uniform(0.001, 0.999, (300, 3))
  tv = TorchVinecop.from_vinecop(_fit_tll_vine(U))

  s = tv.sample(n=500, qrng=True, seeds=[1, 2, 3])
  assert s.shape == (500, 3)
  assert (s > 0).all() and (s < 1).all()


def test_cdf_matches_cpp_within_mc_error():
  rng = np.random.default_rng(0)
  U = rng.uniform(0.001, 0.999, (300, 3))
  cop = _fit_tll_vine(U)
  tv = TorchVinecop.from_vinecop(cop)
  query = torch.from_numpy(U[:10])

  cdf_torch = tv.cdf(query, N=20000, qrng=True, seeds=[1, 2, 3]).numpy()
  cdf_cpp = np.asarray(cop.cdf(U[:10], N=20000, seeds=[1, 2, 3]))
  # Both are MC estimates with N=20000; agree to ~3 sig figs.
  np.testing.assert_allclose(cdf_torch, cdf_cpp, atol=2e-2)


def test_pdf_accepts_num_threads_no_op():
  rng = np.random.default_rng(0)
  U = rng.uniform(0.001, 0.999, (200, 3))
  tv = TorchVinecop.from_vinecop(_fit_tll_vine(U))
  u_t = torch.from_numpy(U[:10])
  p1 = tv.pdf(u_t, num_threads=1)
  p4 = tv.pdf(u_t, num_threads=4)
  assert torch.allclose(p1, p4)


def test_rosenblatt_inverse_rosenblatt_accept_num_threads():
  rng = np.random.default_rng(0)
  U = rng.uniform(0.001, 0.999, (200, 3))
  tv = TorchVinecop.from_vinecop(_fit_tll_vine(U))
  u_t = torch.from_numpy(U[:5])
  # Just confirm the kwarg is accepted; behavior is identical.
  r1 = tv.rosenblatt(u_t, num_threads=1)
  r4 = tv.rosenblatt(u_t, num_threads=4)
  assert torch.allclose(r1, r4)
  ir1 = tv.inverse_rosenblatt(u_t, num_threads=1)
  ir4 = tv.inverse_rosenblatt(u_t, num_threads=4)
  assert torch.allclose(ir1, ir4)


def test_pdf_autograd_through_grid_param() -> None:
  """Grad flows through ``TorchVinecop.pdf`` to a ``TorchBicop`` grid value.

  Guards the Stage-2a cascade extraction onto ``VinecopBase``: ``pdf`` /
  ``rosenblatt`` must stay autograd-capable (only inverse / sample / cdf are
  wrapped in ``no_grad``), and the cascade must not detach pair outputs.
  """
  u_fit = _simulate(d=4, n=800, seed=700)
  cop = _fit_tll_vine(u_fit)
  bc = TorchVinecop.from_vinecop(cop, cache_integrals=False)
  # The stored module, not what `_get_pair_copula` hands the cascade: that is
  # the grid the parameters live on, and it is the same object here (this vine
  # is continuous, so nothing is wrapped).
  pair = bc._pair_module(0, 0)
  pair.interp_grid.values.requires_grad_(True)
  u = torch.from_numpy(_eval_grid(64, d=4, seed=701))
  out = bc.pdf(u, batched=False)
  assert out.requires_grad
  out.sum().backward()
  grad = pair.interp_grid.values.grad
  assert grad is not None
  assert torch.isfinite(grad).all()
  assert grad.abs().sum() > 0


def test_conditional_cdf_raises() -> None:
  """``cdf`` with external covariates ``x`` is unsupported (raises)."""
  u_fit = _simulate(d=4, n=400, seed=710)
  bc = TorchVinecop.from_vinecop(_fit_tll_vine(u_fit))
  u = torch.from_numpy(_eval_grid(5, d=4, seed=711))
  x = torch.zeros(5, 1, dtype=torch.float64)
  with pytest.raises(NotImplementedError, match="Conditional cdf"):
    bc.cdf(u, x=x)


def test_sample_conditional_requires_row_per_sample() -> None:
  """``sample(n, x=...)`` requires exactly one covariate row per sample."""
  u_fit = _simulate(d=3, n=400, seed=712)
  bc = TorchVinecop.from_vinecop(_fit_tll_vine(u_fit))
  x = torch.zeros(4, 1, dtype=torch.float64)
  with pytest.raises(ValueError, match="one covariate row per sample"):
    bc.sample(10, x=x)


def test_batched_cache_stays_out_of_the_state_dict() -> None:
  # The batched state is a memo derived from the pair copulas, not part of
  # the model. Registering it as a child module would put derived buffers in
  # every checkpoint taken after a batched call, and `load_state_dict` on a
  # fresh model would then reject them as unexpected keys.
  rng = np.random.default_rng(0)
  d = 4
  u = pv.to_pseudo_obs(
    rng.standard_normal((400, d)) @ rng.standard_normal((d, d))
  )
  cpp = pv.Vinecop.from_data(
    u, controls=pv.FitControlsVinecop(family_set=[pv.families.tll])
  )
  vine = TorchVinecop.from_vinecop(cpp)
  data = torch.as_tensor(u)

  keys_before = set(vine.state_dict())
  vine.pdf(data, batched=True)
  vine.rosenblatt(data, batched=True)
  assert set(vine.state_dict()) == keys_before

  fresh = TorchVinecop.from_vinecop(cpp)
  fresh.load_state_dict(vine.state_dict(), strict=True)
  torch.testing.assert_close(fresh.pdf(data), vine.pdf(data))


#: Cumulative atom masses of a Binomial(4, 0.5) count variable: its top atom
#: reaches ``F(x) = 1`` and its bottom ``F(x^-) = 0``, both places where the
#: cascade's clamp and a pair copula's own trimming have to agree.
_COUNT_CDF = np.cumsum([1.0, 4.0, 6.0, 4.0, 1.0]) / 16.0


def _discrete_data(
  var_types: list[str], n: int = 1500, seed: int = 3
) -> np.ndarray:
  """An ``(n, d + k)`` sample: counts at ``"d"``, continuous ranks at ``"c"``."""
  rng = np.random.default_rng(seed)
  d = len(var_types)
  base = rng.standard_normal((n, 1))
  z = 0.7 * base + 0.7 * rng.standard_normal((n, d))
  values, limits = [], []
  for j, t in enumerate(var_types):
    p = pv.to_pseudo_obs(z[:, [j]]).ravel()
    if t == "d":
      level = np.searchsorted(_COUNT_CDF[:-1], p)
      values.append(_COUNT_CDF[level])
      limits.append(np.where(level > 0, _COUNT_CDF[level - 1], 0.0))
    else:
      values.append(p)
  return np.column_stack(values + limits)


def _discrete_vinecop(var_types: list[str], u: np.ndarray) -> pv.Vinecop:
  return pv.Vinecop.from_data(u, var_types=var_types, controls=_TLL_CONTROLS)


@pytest.mark.parametrize(
  "var_types", [["d", "c", "c"], ["c", "d", "d"], ["d", "d", "d"]]
)
def test_from_vinecop_matches_discrete_vinecop(var_types: list[str]) -> None:
  # A discrete C++ vine lifted into torch. The stored grids are continuous and
  # `DiscretePair` supplies the mixed-discrete surface, so the difference is the
  # grid interpolation's, not the cascade's: `TorchBicop.cdf(cache=False)`
  # matches `Bicop.cdf` at 1e-10, and a difference quotient divides that by the
  # atom width (>= 0.06 for these counts), so ~1e-8 is the floor.
  u = _discrete_data(var_types)
  cop = _discrete_vinecop(var_types, u)
  bc = TorchVinecop.from_vinecop(cop, cache_integrals=False)
  assert bc.var_types == var_types
  u_t = torch.from_numpy(u)
  np.testing.assert_allclose(
    bc.pdf(u_t).numpy(), cop.pdf(u), atol=1e-7, rtol=1e-7
  )
  np.testing.assert_allclose(
    bc.rosenblatt(u_t, randomize_discrete=False).numpy(),
    cop.rosenblatt(u, randomize_discrete=False),
    atol=1e-7,
    rtol=1e-7,
  )


def test_from_structure_keeps_discrete_var_types() -> None:
  s = pv.RVineStructure.sample(3, seeds=[1, 2, 3, 4, 5])
  vine = TorchVinecop.from_structure(structure=s, var_types=["c", "d", "c"])
  assert vine.var_types == ["c", "d", "c"]
  # An independence vine has density 1 whatever the variable types are, and the
  # discrete edges still consume the four-column input.
  u = _discrete_data(["c", "d", "c"], n=64, seed=5)
  np.testing.assert_allclose(
    vine.pdf(torch.from_numpy(u)).numpy(), np.ones(64), atol=1e-12
  )


def test_from_data_matches_discrete_vinecop() -> None:
  # End to end: the torch TLL fit on data with atoms. `TllBicop::fit` ranks with
  # random tie-breaking from a fixed seed, which `TorchBicop.from_data` reuses,
  # so the fitted grids agree to machine precision and the vines to the
  # interpolation floor above.
  #
  # DELIBERATE TRIPWIRE. This holds because upstream computes a latent sample
  # for discrete data and then never fits on it (vinecopulib#736), which makes
  # tie-breaking the only thing an atom changes about the grid. If this test
  # starts failing right after a submodule bump, upstream fixed #736 and the
  # latent-sample step has to be ported into `_fit_tll` in the same PR --
  # regenerate the expected values from the new C++ output, and do not widen the
  # tolerances to absorb it.
  var_types = ["d", "c", "c"]
  u = _discrete_data(var_types)
  cop = _discrete_vinecop(var_types, u)
  bc = TorchVinecop.from_data(
    torch.from_numpy(u),
    var_types=var_types,
    controls=FitControlsTorchVinecop(cache_integrals=False),
  )
  assert np.array_equal(
    np.asarray(bc.structure.matrix), np.asarray(cop.structure.matrix)
  )
  np.testing.assert_allclose(
    bc.pdf(torch.from_numpy(u)).numpy(), cop.pdf(u), atol=1e-6, rtol=1e-6
  )


def test_discrete_vine_declines_the_batched_path() -> None:
  # The batched level carries no distribution-function grid, which a discrete
  # edge's h-functions are difference quotients of, so an explicit
  # `batched=True` must fall back rather than evaluate something else.
  var_types = ["d", "c", "c"]
  u = _discrete_data(var_types, n=400, seed=8)
  bc = TorchVinecop.from_vinecop(_discrete_vinecop(var_types, u))
  u_t = torch.from_numpy(u)
  np.testing.assert_array_equal(
    bc.pdf(u_t, batched=True).numpy(), bc.pdf(u_t, batched=False).numpy()
  )


def test_discrete_vine_round_trips_through_pickle() -> None:
  # `var_types` is plain Python state on the vine, and the per-edge type table
  # is derived from it at bind time — so both have to survive a round trip, or a
  # reloaded vine would evaluate a continuous density on the same data.
  var_types = ["d", "c", "c"]
  u = _discrete_data(var_types, n=400, seed=9)
  bc = TorchVinecop.from_vinecop(_discrete_vinecop(var_types, u))
  clone = pickle.loads(pickle.dumps(bc))
  assert clone.var_types == var_types
  assert clone.pair_var_types(0, 0) == bc.pair_var_types(0, 0)
  u_t = torch.from_numpy(u)
  np.testing.assert_array_equal(clone.pdf(u_t).numpy(), bc.pdf(u_t).numpy())
