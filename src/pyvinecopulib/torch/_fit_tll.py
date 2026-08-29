"""Pure-PyTorch port of the C++ ``TllBicop::fit`` for the ``constant`` method.

Mirrors the algorithm in ``lib/vinecopulib/.../bicop/implementation/tll.ipp``
end-to-end so :meth:`TorchBicop.from_data` produces an ``(m, m)`` density
grid that matches the C++ ``pv.Bicop.from_data`` output to machine
precision after :meth:`InterpolationGrid2D.normalize_margins(25)`.

Three pieces:

* :func:`_to_pseudo_obs_continuous` — empirical CDF on continuous data
  (rank/(n+1); no tie handling, mirrors C++ ``wdm`` ranks for
  jitter-free input).
* :func:`_ace` — alternating conditional expectations for the maximal-
  correlation coefficient. Outer/inner convergence loop matches the C++
  tolerances (``2e-15`` / ``1e-4``); the moving-average window smoother is
  a prefix-sum box mean rather than the C++ FFT, which is the same quantity
  in ``O(n)`` instead of ``O(n log n)``.
* :func:`fit_tll_constant` — the user-facing entry point. Wraps the
  bandwidth selection (``select_bandwidth_constant``) and the local-
  constant kernel density evaluation (``fit_local_likelihood_constant``)
  into a single call returning the unnormalized ``(m, m)`` values grid.

Only the ``constant`` method is supported here; the ``linear`` and
``quadratic`` variants need more elaborate per-grid-point bookkeeping
(see C++ ``TllBicop::fit_local_likelihood`` and ``calculate_infl``).
"""

from __future__ import annotations

import math
from typing import Callable, Optional

import torch
from torch import Tensor

from ._interp import InterpolationGrid2D

_SQRT_2PI_INV = 1.0 / math.sqrt(2.0 * math.pi)


# --------------------------------------------------------------------------- #
# Helpers                                                                      #
# --------------------------------------------------------------------------- #


def _qnorm(p: Tensor) -> Tensor:
  """Inverse standard normal CDF; matches C++ ``tools_stats::qnorm``."""
  return torch.special.ndtri(p)


def _to_pseudo_obs_continuous(x: Tensor) -> Tensor:
  """Empirical CDF on continuous data: ``rank/(n+1)``.

  Mirrors ``tools_stats::to_pseudo_obs`` with no ties (continuous input,
  no jittering needed). C++ uses ``wdm`` ranks; for unique data those are
  the same as :func:`torch.argsort` ranks.

  Both sorts are stable, so tied inputs get a defined order rather than
  whichever one the sort happens to produce for the shape it was handed --
  a leading pair axis otherwise makes the tie-break shape-dependent. Ties
  reach here even on continuous data: callers trim to
  ``(1e-10, 1 - 1e-10)`` first, which collapses everything beyond the clamp
  into one group.
  """
  n = x.shape[-2]
  ranks = x.argsort(dim=-2, stable=True).argsort(dim=-2, stable=True)
  return (ranks + 1).to(x.dtype) / (n + 1)


def _win_smoother(x: Tensor, wl: int) -> Tensor:
  """Centered moving average of half-window ``wl``, with edge clamps.

  Mirrors ``tools_stats::win``: the window is a box of length ``2*wl + 1``,
  the data are zero-padded outside their range, and the leading / trailing
  ``wl`` entries are then clamped to ``out[wl]`` and ``out[n - wl - 1]`` so
  the edges are flat. Batched over any leading dimensions; the window runs
  along the last one.

  A box mean is a difference of prefix sums, which is what this computes in
  ``O(n)``; ``tools_stats::win`` reaches the same quantity by FFT. Since
  ``wl`` is ``ceil(n / 5)`` the window grows with the data, so a direct
  convolution against the ``2*wl + 1``-tap kernel would be ``O(n**2)``.

  Differencing two prefix sums costs relative accuracy, because both
  operands carry the magnitude of the whole partial sum rather than of the
  window. That is immaterial for ACE, which smooths standardized scores so
  the partial sums stay ``O(sqrt(n))``, and it is the constraint a future
  caller has to check: on mean-shifted input the ratio grows.
  """
  n = x.shape[-1]
  zero = torch.zeros(x.shape[:-1] + (1,), dtype=x.dtype, device=x.device)
  prefix = torch.cat([zero, x.cumsum(-1)], dim=-1)
  idx = torch.arange(n, device=x.device)
  hi = (idx + wl + 1).clamp(max=n)
  lo = (idx - wl).clamp(min=0)
  # Dividing by the full window width, not by the number of terms in it, is
  # what makes the zero-padding of `tools_stats::win` show up: the tails are
  # averaged over fewer data points than the divisor counts. The clamps
  # below then overwrite exactly the region where that matters.
  out = (prefix[..., hi] - prefix[..., lo]) / (2 * wl + 1)
  if wl > 0:
    out[..., :wl] = out[..., wl : wl + 1]
    out[..., -wl:] = out[..., n - wl - 1 : n - wl]
  return out


def _cef(x: Tensor, ind: Tensor, ranks: Tensor, wl: int) -> Tensor:
  """``cef`` helper: ``win(x[ind], wl)[ranks]``.

  Smooths ``x`` in sorted order, then maps back to the original order.
  Mirrors ``tools_stats::cef``. Batched over any leading dimensions: the
  two permutations are per-lane, so they are applied with ``gather``
  rather than by fancy indexing.
  """
  sorted_x = torch.gather(x, -1, ind)
  return torch.gather(_win_smoother(sorted_x, wl), -1, ranks)


#: ``tools_stats::ace``'s outer tolerance, at float64. Scaled by the eps
#: ratio for narrower dtypes; see ``_ace``.
_ACE_OUTER_ABS_TOL: float = 2e-15


def _ace_step(
  source: Tensor,
  target: Tensor,
  live: Tensor,
  eps_prev: Tensor,
  abs_err: Tensor,
  trip: Tensor,
  gate: Tensor,
  n: int,
  ind: Tensor,
  ranks: Tensor,
  wl: int,
  iter_max: int,
  abs_tol: float,
) -> "tuple[Tensor, Tensor, Tensor, Tensor, Tensor]":
  """One ACE pass, including the per-lane loop state it advances.

  Both of :func:`_ace`'s loops take this same step, differing in which score
  they smooth, which permutation they smooth it under, and what gates their
  liveness. The state moves with the step rather than beside it so that the
  whole pass is one region -- the counter updates are four more kernels on
  ``(P,)`` tensors, which cost the same as the arithmetic they bookkeep.

  Every write is masked by ``live``, so a pass with a lane dead leaves that
  lane exactly where it was.

  Parameters
  ----------
  source : Tensor, shape (..., n)
      The score being conditioned on.
  target : Tensor, shape (..., n)
      The score being rewritten.
  live : Tensor, shape (...), dtype bool
      Lanes still iterating.
  eps_prev, abs_err, trip : Tensor, shape (...)
      The previous criterion, the last change in it, and the trip counter.
  gate : Tensor, shape (...), dtype bool
      An outer condition every lane must also satisfy to stay live. The
      inner loop passes the outer loop's liveness; the outer passes all-true.
  n : int
      Observations per lane.
  ind, ranks : Tensor, shape (..., n), dtype int
      The sort permutation of the conditioning argument, and its inverse.
  wl : int
      Smoothing half-window.
  iter_max : int
      Per-lane cap on ``trip``.
  abs_tol : float
      Convergence tolerance on ``abs_err``.

  Returns
  -------
  tuple of Tensor
      ``target``, ``eps_prev``, ``abs_err``, ``trip`` and the new ``live``.
  """
  v = _cef(source, ind, ranks, wl)
  v = v - v.sum(-1, keepdim=True) / n
  cand = v / ((v**2).sum(-1, keepdim=True) / (n - 1)).sqrt()
  eps = ((cand - source) ** 2).sum(-1) / n
  target = torch.where(live.unsqueeze(-1), cand, target)
  abs_err = torch.where(live, (eps_prev - eps).abs(), abs_err)
  eps_prev = torch.where(live, eps, eps_prev)
  trip = torch.where(live, trip + 1, trip)
  return (
    target,
    eps_prev,
    abs_err,
    trip,
    gate & (trip <= iter_max) & (abs_err > abs_tol),
  )


#: What :func:`_ace_step` and its compiled twin both are.
_AceStep = Callable[..., "tuple[Tensor, Tensor, Tensor, Tensor, Tensor]"]

#: The compiled :func:`_ace_step`, built on first use. Compiling costs seconds
#: of Inductor, so it is not paid by a process that never asks for it.
_COMPILED_ACE_STEP: Optional[_AceStep] = None


def _compiled_ace_step() -> _AceStep:
  """:func:`_ace_step` fused into one kernel, compiled once per process.

  ``dynamic=True`` because a vine hands the step one lane count per tree
  level -- nineteen of them at ``d = 20``. Specializing on each would exceed
  ``torch._dynamo.config.cache_size_limit`` (8 by default) and silently drop
  the widest levels back to eager, and raising that limit from inside a
  library would reach every other compiled function in the process. Dynamic
  shapes cover every width with about two graphs, for roughly a fifth of the
  speedup a specialized one gets.

  Returns
  -------
  callable
      The compiled step, with :func:`_ace_step`'s signature.
  """
  global _COMPILED_ACE_STEP
  if _COMPILED_ACE_STEP is None:
    _COMPILED_ACE_STEP = torch.compile(_ace_step, dynamic=True)
  return _COMPILED_ACE_STEP


def _ace(
  data: Tensor,
  *,
  outer_iter_max: int = 100,
  inner_iter_max: int = 10,
  # ``None`` resolves to upstream's tolerance scaled to the working
  # precision. The float64 value is ``tools_stats::ace``'s verbatim, so the
  # C++<->torch parity guard (tests/test_torch_bicop.py) is unmoved; a lower
  # precision cannot resolve it, and would iterate to `outer_iter_max`
  # against its own rounding noise instead of converging.
  outer_abs_tol: Optional[float] = None,
  inner_abs_tol: float = 1e-4,
  compile_step: bool = False,
) -> Tensor:
  """Alternating conditional expectations.

  Mirrors ``tools_stats::ace`` for the unweighted bivariate case (the
  only one that the TLL ``constant`` bandwidth path needs).

  Batched over any leading dimensions. Both convergence loops advance every
  lane together and *freeze* each as it converges -- its state, including
  its two trip counters, stops moving while the others carry on -- so
  *which* lanes a pair travelled with does not change its answer at all,
  and its iteration count is the one its own data earns.

  *How many* it travelled with is a separate matter and does move the last
  bits, here as everywhere in the batched fit: torch selects elementwise
  kernels by element count. The exact half of that promise is pinned by
  ``test_ace_freezes_each_lane_independently``, which varies the companions
  at a fixed shape.

  A batch runs until its slowest lane converges, so it trades the sum of the
  lanes' iterations for their maximum while doing the full batch's arithmetic
  at every step. How much that wins depends on how alike the lanes are, and
  iteration counts rise as dependence falls -- with no signal to find, the
  outer criterion grinds against its own rounding noise until it hits
  ``outer_iter_max``. One near-independent lane therefore paces a whole
  level, and the deeper trees, whose pairs sit nearer independence, are
  where that bites.

  Args:
    data: shape ``(..., n, 2)``.
    compile_step: fuse the per-pass body with ``torch.compile``. The pass is
      39 kernel launches over tensors small enough that the arithmetic is
      invisible, so fusing it is worth several times its cost at the sizes a
      vine is usually fitted at. It perturbs the last bits, and the outer
      criterion sits at the float64 noise floor, so a lane's iteration count
      can move; the C++ grid gate covers the result.

  Returns:
    ``(..., n, 2)`` tensor of the ACE-transformed scores ``phi``.
  """
  n = data.shape[-2]
  batch = data.shape[:-2]
  dtype, device = data.dtype, data.device
  if outer_abs_tol is None:
    outer_abs_tol = _ACE_OUTER_ABS_TOL * float(
      torch.finfo(dtype).eps / torch.finfo(torch.float64).eps
    )
  wl = int(math.ceil(n / 5.0))
  step = _compiled_ace_step() if compile_step else _ace_step
  # The outer loop has no condition beyond its own; the inner passes it.
  always = torch.ones(batch, dtype=torch.bool, device=device)

  # Two contiguous buffers rather than the columns of one `(..., n, 2)`
  # matrix: every operation below then runs on a contiguous buffer instead
  # of a stride-2 view, which is the same arithmetic in the same order but
  # fewer and cheaper kernels. The pair is stacked back on return.
  ind0 = data[..., 0].argsort(dim=-1, stable=True)
  ind1 = data[..., 1].argsort(dim=-1, stable=True)
  positions = torch.arange(n, device=device).expand(batch + (n,))
  ranks0 = torch.empty_like(ind0).scatter_(-1, ind0, positions)
  ranks1 = torch.empty_like(ind1).scatter_(-1, ind1, positions)

  shift = (n - 1) / 2.0 - 1.0
  scale = math.sqrt(n * (n - 1) / 12.0)
  phi0 = (ranks0.to(dtype) - shift) / scale
  phi1 = (ranks1.to(dtype) - shift) / scale

  # Per-lane loop state. A lane that has converged is *frozen*: its state
  # stops moving while the others keep iterating, which is what reproduces
  # the per-lane sequential answer exactly rather than approximately. The
  # trip counters are frozen too, because the caps are per-lane.
  ones = torch.ones(batch, dtype=dtype, device=device)
  longs = torch.ones(batch, dtype=torch.long, device=device)
  # Clones, not aliases: `ones` / `longs` are the templates every outer
  # iteration re-clones the inner state from, so an in-place write to either
  # of these would corrupt them.
  outer_eps, outer_abs_err, outer_iter = (
    ones.clone(),
    ones.clone(),
    longs.clone(),
  )
  outer_live = (outer_iter <= outer_iter_max) & (outer_abs_err > outer_abs_tol)

  while bool(outer_live.any()):
    inner_eps, inner_abs_err = ones.clone(), ones.clone()
    inner_iter = longs.clone()
    inner_live = (
      outer_live
      & (inner_iter <= inner_iter_max)
      & (inner_abs_err > inner_abs_tol)
    )
    while bool(inner_live.any()):
      # `phi1` is written by the inner loop, so it is gated on `inner_live`
      # -- which carries `outer_live`. A lane whose *outer* loop converged
      # must not have its `phi1` moved by other lanes' inner iterations.
      phi1, inner_eps, inner_abs_err, inner_iter, inner_live = step(
        phi0,
        phi1,
        inner_live,
        inner_eps,
        inner_abs_err,
        inner_iter,
        outer_live,
        n,
        ind1,
        ranks1,
        wl,
        inner_iter_max,
        inner_abs_tol,
      )
    phi0, outer_eps, outer_abs_err, outer_iter, outer_live = step(
      phi1,
      phi0,
      outer_live,
      outer_eps,
      outer_abs_err,
      outer_iter,
      always,
      n,
      ind0,
      ranks0,
      wl,
      outer_iter_max,
      outer_abs_tol,
    )

  return torch.stack([phi0, phi1], dim=-1)


def _pearson_cor(x: Tensor) -> Tensor:
  """Pearson correlation of the two columns of ``x: (..., n, 2)``.

  Returns one value per leading lane, so a 0-D tensor for a single pair.
  """
  x0 = x[..., 0] - x[..., 0].mean(-1, keepdim=True)
  x1 = x[..., 1] - x[..., 1].mean(-1, keepdim=True)
  return (x0 * x1).sum(-1) / ((x0**2).sum(-1).sqrt() * (x1**2).sum(-1).sqrt())


def _pairwise_mcor(x: Tensor, *, compile_step: bool = False) -> Tensor:
  """Maximal correlation via ACE + Pearson, one value per leading lane."""
  return _pearson_cor(_ace(x, compile_step=compile_step))


def _chol22(B: Tensor) -> Tensor:
  """Cholesky factor of a 2x2 SPD matrix; lower-triangular.

  Assembled with one ``stack`` rather than four writes into a zero matrix:
  each of those was a kernel launch for a single element, and on a device
  the launches cost more than the arithmetic.
  """
  a = B[..., 0, 0].sqrt()
  b = B[..., 1, 0] / a
  c = (B[..., 1, 1] - b**2).sqrt()
  zero = torch.zeros_like(a)
  return torch.stack(
    [torch.stack([a, zero], dim=-1), torch.stack([b, c], dim=-1)], dim=-2
  )


def _select_bandwidth_constant(
  z: Tensor, *, compile_step: bool = False
) -> Tensor:
  """Bandwidth matrix for the constant-method local-likelihood KDE.

  Mirrors ``TllBicop::select_bandwidth`` for ``method == "constant"``.
  Batched over any leading dimensions of ``z``, and entirely on device.

  Staying on device is what matches C++ at ``mcor == 0``, where
  ``std::pow(std::fabs(cor / mcor), 0.5 * mcor)`` is ``pow(inf, 0.0)`` --
  ``1.0`` by IEEE, where a host-float division raises instead. The division
  must stay unguarded: a ``clamp_min`` would diverge from upstream rather
  than agree with it.
  """
  n = z.shape[-2]
  cor = _pearson_cor(z).clamp(-0.95, 0.95)
  eye = torch.eye(2, dtype=z.dtype, device=z.device)
  cov = eye + cor[..., None, None] * (1.0 - eye)
  mult = n ** (-1.0 / 3.0)
  mcor = _pairwise_mcor(z, compile_step=compile_step)
  scale = (cor / mcor).abs() ** (0.5 * mcor)
  return mult * cov * scale[..., None, None]


#: Peak the grid blocking in :func:`_fit_local_likelihood_constant` aims to
#: stay under. The block's live ``(P, block, n)`` temporaries dominate the
#: fit's footprint, and six of them are in flight at once, so the peak is
#: ``6 * itemsize * block * n * P`` -- measured 48.5 bytes per element at
#: float64, against 48 predicted. Blocking is exact, grid points not
#: interacting, so trading block size for memory changes only the launch
#: count: a level 38 blocks deep instead of 4 costs launches that the
#: arithmetic hides, where a fixed block costs gigabytes that it does not.
_KDE_MEM_BUDGET_BYTES: int = 256 * 1024 * 1024


def _kde_grid_block(n: int, lanes: int, itemsize: int, grid: int) -> int:
  """Grid points per block that keep the live temporaries inside the budget.

  Parameters
  ----------
  n : int
      Observations per lane.
  lanes : int
      Pairs fitted together.
  itemsize : int
      Bytes per element of the working dtype.
  grid : int
      Total grid points; the block never exceeds it.

  Returns
  -------
  int
      Block size, at least one.
  """
  per_point = 6 * itemsize * n * max(lanes, 1)
  return max(1, min(grid, _KDE_MEM_BUDGET_BYTES // max(per_point, 1)))


def _fit_local_likelihood_constant(
  z: Tensor, z_data: Tensor, B: Tensor
) -> Tensor:
  """Local-constant kernel density estimate at ``z`` from ``z_data``.

  Mirrors ``TllBicop::fit_local_likelihood`` for ``method == "constant"``.
  ``z`` is ``(G, 2)`` grid points shared by every lane; ``z_data`` is
  ``(..., n, 2)`` and ``B`` is ``(..., 2, 2)``, so the result is
  ``(..., G)``.

  Evaluated in blocks of grid points rather than all at once. C++ loops one
  grid point at a time (``tll.ipp``); the two extremes bracket the same
  computation, because a grid point's density is a mean over the data and
  nothing couples one grid point to another. Blocking is therefore exact,
  and the block is sized from ``n`` and the lane count so the peak stays
  near :data:`_KDE_MEM_BUDGET_BYTES` however wide the level is -- a fixed
  block makes the footprint grow with both, which is what a level of a
  large vine on a large sample runs out of memory doing.
  """
  # Closed forms rather than `torch.linalg.inv` / `det`: the factor is 2x2
  # lower-triangular, and `inv` dispatches to cuSOLVER and reads its `info`
  # back to the host -- a device synchronization for a four-element matrix.
  rB = _chol22(B)
  ra, rb, rc = rB[..., 0, 0], rB[..., 1, 0], rB[..., 1, 1]
  ia = ra.reciprocal()
  ic = rc.reciprocal()
  ib = -rb * ia * ic
  det_irB = ia * ic
  # Applying the lower-triangular factor in closed form rather than as a
  # matmul: row 0 is `ia * z0 + 0 * z1` and row 1 sums `ib * z0 + ic * z1`
  # in that order, which is what the matmul did, without assembling the
  # matrix or materializing a trailing axis of 2.
  zd0 = ia[..., None] * z_data[..., 0]
  zd1 = ib[..., None] * z_data[..., 0] + ic[..., None] * z_data[..., 1]
  const = (_SQRT_2PI_INV * _SQRT_2PI_INV) * det_irB
  lanes = 1
  for extent in z_data.shape[:-2]:
    lanes *= int(extent)
  block = _kde_grid_block(
    int(z_data.shape[-2]), lanes, z_data.element_size(), int(z.shape[0])
  )
  out = []
  for lo in range(0, z.shape[0], block):
    zb = z[lo : lo + block]
    zg0 = ia[..., None] * zb[:, 0]
    zg1 = ib[..., None] * zb[:, 0] + ic[..., None] * zb[:, 1]
    d0 = zd0[..., None, :] - zg0[..., :, None]
    d1 = zd1[..., None, :] - zg1[..., :, None]
    kernels = torch.exp(-0.5 * (d0**2 + d1**2)) * const[..., None, None]
    out.append(kernels.mean(dim=-1))
  return torch.cat(out, dim=-1) if len(out) > 1 else out[0]


# --------------------------------------------------------------------------- #
# Public entry point                                                           #
# --------------------------------------------------------------------------- #


def fit_tll_constant(
  u: Tensor,
  grid_size: int = 30,
  mult: float = 1.0,
  grid_type: str = "normal",
  pseudo_obs: Optional[Tensor] = None,
  discrete_data: Optional[Tensor] = None,
  compile_fit: bool = False,
) -> tuple[Tensor, Tensor]:
  """Fit a TLL pair-copula via local-constant kernel density estimation.

  Pure-PyTorch port of ``TllBicop::fit`` for ``method == "constant"``.

  Args:
    u: ``(n, 2)`` pseudo-observations in ``(0, 1)``, or ``(P, n, 2)`` to fit
      ``P`` pairs together -- they share the grid, the window length and the
      eval points. A discrete edge is refused a pair axis.
    grid_size: number of grid points per axis (default 30; matches C++).
    mult: bandwidth multiplier passed through to ``select_bandwidth``;
      the C++ default is 1.
    pseudo_obs: ranks to fit on, overriding the ones derived from ``u``. C++
      ranks with ``ties_method="random"``, which differs from the argsort ranks
      only when the data has ties — i.e. when a margin has atoms, where the
      caller supplies them. On a discrete edge these ranks only seed the
      bandwidth; see ``discrete_data``.
    discrete_data: the ``(n, 4)`` layout ``[u1, u2, u1^-, u2^-]`` of a discrete
      or mixed edge. When given, the fit runs on the *latent* sample recovered
      from it rather than on the ranks, which is what ``TllBicop::fit`` does:
      the bandwidth is selected from the jittered ranks and only then is the
      latent sample drawn, with ``(B00 * B11) ** 0.25`` as its own bandwidth and
      ``B`` left as selected.
    compile_fit: fuse the bandwidth search's per-pass body with
      ``torch.compile``. Off by default: the first call compiles for seconds,
      which only a caller fitting many pairs in one process earns back. See
      ``FitControlsTorchBicop.compile_fit``.
    grid_type: either ``"normal"`` (default, matches C++ — grid_points are
      ``pnorm(linspace(-3.25, 3.25, m))``) or ``"linear"`` (grid_points
      are uniformly spaced on ``[Phi(-3.25), 1-Phi(-3.25)]``, with the
      InterpolationGrid2D constructor forcing the endpoints to exactly
      0 / 1). The linear option matches the z-range of the normal grid so
      bandwidth selection is unaffected; only the *storage* grid changes,
      enabling O(1) cell-finding at eval time.

  Returns:
    A ``(grid_points, values)`` pair. ``values`` is the unnormalized
    ``(m, m)`` density; callers should pass it through
    ``InterpolationGrid2D(grid_points, values, norm_maxiter=25,
    is_linear=(grid_type == "linear"))`` to match the C++
    ``Bicop.parameters`` output to machine precision for ``"normal"``.

    ``u`` may carry a leading pair axis, ``(P, n, 2)``, in which case
    ``values`` is ``(P, m, m)`` and the P pairs are fitted together: they
    share the grid, the window length and the eval points, and the two
    convergence loops in :func:`_ace` advance every lane at once, freezing
    each as it converges. ``grid_points`` is shared, so it is returned
    once. A discrete edge is refused with a pair axis -- see below.
  """
  if u.ndim not in (2, 3) or u.shape[-1] != 2:
    raise ValueError(
      f"u must have shape (n, 2) or (P, n, 2); got {tuple(u.shape)}"
    )
  if grid_size < 2:
    raise ValueError(f"grid_size must be >= 2; got {grid_size}")
  if mult <= 0:
    raise ValueError(f"mult must be > 0; got {mult}")
  dtype, device = u.dtype, u.device

  # Pseudo-observations + qnorm to z-space.
  psobs = _to_pseudo_obs_continuous(u) if pseudo_obs is None else pseudo_obs
  z_data = _qnorm(psobs)

  # Bandwidth selection.
  B = _select_bandwidth_constant(z_data, compile_step=compile_fit) * mult

  if discrete_data is not None:
    # Reuse the compiled draw, as structure selection reuses `wdm`: it is a
    # stochastic iterative reconstruction over a spatial index, and reproducing
    # it in torch would put the torch<->C++ grid parity at the mercy of a
    # reimplementation rather than of the same code.
    from ..pyvinecopulib_ext import find_latent_sample

    if u.ndim != 2:
      raise ValueError(
        "a discrete edge cannot be fitted with a leading pair axis: "
        "`find_latent_sample` is a compiled per-pair draw over a spatial "
        "index with a fixed-seed generator, so it has no batch axis to "
        "give. Fit discrete edges one at a time."
      )
    latent = find_latent_sample(
      discrete_data.detach().cpu().numpy(),
      float((B[0, 0] * B[1, 1]).item() ** 0.25),
    )
    z_data = _qnorm(torch.as_tensor(latent, dtype=dtype, device=device))

  # Storage and KDE-eval grids come from the centralized factory on
  # InterpolationGrid2D so any future grid type lives in one place.
  # For "normal" the two coincide and match C++; for "linear" the KDE
  # eval grid is uniform-then-clamped so qnorm stays finite while the
  # storage grid is the unclipped linspace(0, 1, m) — the resulting
  # densities at u=0/1 are stored at z=±3.25 (same trick as C++).
  grid_points = InterpolationGrid2D.make_grid_points(
    grid_type, grid_size, dtype=dtype, device=device
  )
  kde_points = InterpolationGrid2D.make_kde_eval_points(
    grid_type, grid_size, dtype=dtype, device=device
  )

  # Expand grid (row-major matches the C++ Eigen::Map<...>().transpose()
  # reshape in TllBicop::fit at the values-assembly step).
  g_i = kde_points.repeat_interleave(grid_size)
  g_j = kde_points.repeat(grid_size)
  grid_2d = torch.stack([g_i, g_j], dim=-1)
  z = _qnorm(grid_2d)

  # Local-constant KDE on the z-scale grid.
  f0 = _fit_local_likelihood_constant(z, z_data, B)

  # Transform z-space density to copula-scale density.
  phi_z = (
    torch.exp(-0.5 * (z[:, 0] ** 2 + z[:, 1] ** 2))
    * _SQRT_2PI_INV
    * _SQRT_2PI_INV
  )
  c = f0 / phi_z
  values = c.reshape(c.shape[:-1] + (grid_size, grid_size))

  # The canonical TorchBicop builds an InterpolationGrid2D from
  # (grid_points, values) with norm_maxiter=25 — that's what matches C++
  # to machine precision. The grid_points returned here are the
  # un-forced ones used for the fit positions; InterpolationGrid2D's
  # constructor will clamp the endpoints to 0/1 internally.
  return grid_points, values
