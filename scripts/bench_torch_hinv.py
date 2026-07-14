#!/usr/bin/env python
"""Speed / accuracy benchmark for the torch inverse h-functions.

Compares three ways of evaluating ``TorchBicop.hinv1`` / ``hinv2`` on a
fitted TLL bicop:

- ``closed_form`` — the current no-cache path: exact inversion of the
  piecewise-quadratic conditional cdf
  (``InterpolationGrid2D.inverse_integrate_1d``, port of vinecopulib#691).
- ``cached`` — ``cache_integrals=True``: one bilinear interpolation on the
  precomputed inverse grid.
- ``itp`` — the pre-#691-parity reference: a fixed-iter vectorized ITP
  root-finder over the on-the-fly h-function (self-contained copy of the
  removed ``pyvinecopulib.torch._util.solve_itp``), kept here so the
  old-vs-new comparison stays reproducible.

Accuracy is the max absolute difference against the C++ ``Bicop`` fitted on
the same data (which itself uses the #691 closed form).

Outputs a table to stdout with columns:
    device, n, method, time_ms, max_abs_diff_vs_cpp
"""

from __future__ import annotations

import argparse
import time
from typing import Callable

import numpy as np
import torch

import pyvinecopulib as pv
from pyvinecopulib.torch import TorchBicop

# --- reference ITP root-finder (verbatim copy of the removed
# pyvinecopulib.torch._util.solve_itp; Oliveira & Takahashi, 2020) -------- #


@torch.no_grad()
def _solve_itp(
  fun: Callable[[torch.Tensor], torch.Tensor],
  x_a: torch.Tensor,
  x_b: torch.Tensor,
  n_iter: int = 30,
  k1: float = 0.2,
  k2: float = 2.0,
  eps: float = 1e-12,
) -> torch.Tensor:
  x_a = x_a.clone()
  x_b = x_b.clone()
  fa = fun(x_a)
  fb = fun(x_b)
  swap = fa > 0
  x_a, x_b = torch.where(swap, x_b, x_a), torch.where(swap, x_a, x_b)
  fa, fb = torch.where(swap, fb, fa), torch.where(swap, fa, fb)
  for _ in range(n_iter):
    width = x_b - x_a
    mid = 0.5 * (x_a + x_b)
    denom = (fb - fa).clamp_min(1e-300)
    x_f = (fb * x_a - fa * x_b) / denom
    diff = mid - x_f
    sigma = torch.sign(diff)
    delta = k1 * width.abs().pow(k2)
    x_t = torch.where(delta <= diff.abs(), x_f + sigma * delta, mid)
    r = (0.5 * width.abs() - eps).clamp_min(0.0)
    x_itp = torch.where((x_t - mid).abs() <= r, x_t, mid - sigma * r)
    fx = fun(x_itp)
    below = fx < 0.0
    x_a = torch.where(below, x_itp, x_a)
    x_b = torch.where(below, x_b, x_itp)
    fa = torch.where(below, fx, fa)
    fb = torch.where(below, fb, fx)
  return 0.5 * (x_a + x_b)


def _hinv1_itp(bc: TorchBicop, u: torch.Tensor) -> torch.Tensor:
  fixed, p = u[:, 0:1], u[:, 1:2]

  def fun(x: torch.Tensor) -> torch.Tensor:
    u_eval = torch.cat([fixed, x], dim=-1)
    return bc.interp_grid.integrate_1d(u_eval, cond_var=1).unsqueeze(-1) - p

  x = _solve_itp(fun, torch.zeros_like(p), torch.ones_like(p))
  return x.squeeze(-1).clamp(0.0, 1.0)


# --- benchmark ------------------------------------------------------------ #


def _time(fn: Callable[[], torch.Tensor], device: str, reps: int) -> float:
  fn()  # warm-up
  if device == "cuda":
    torch.cuda.synchronize()
  t0 = time.perf_counter()
  for _ in range(reps):
    fn()
  if device == "cuda":
    torch.cuda.synchronize()
  return (time.perf_counter() - t0) / reps * 1e3


def main() -> None:
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument(
    "--devices", nargs="+", default=["cpu"], choices=["cpu", "cuda"]
  )
  parser.add_argument("--sizes", nargs="+", type=int, default=[10_000, 100_000])
  parser.add_argument("--reps", type=int, default=5)
  parser.add_argument("--seed", type=int, default=42)
  args = parser.parse_args()

  rng = np.random.default_rng(args.seed)
  cop_true = pv.Bicop(family=pv.families.gaussian, parameters=np.array([[0.7]]))
  data = cop_true.simulate(2_000, seeds=[args.seed])
  cpp = pv.Bicop.from_data(
    data, controls=pv.FitControlsBicop(family_set=[pv.families.tll])
  )

  print(
    f"{'device':<8}{'n':>10}{'method':>14}{'time_ms':>12}"
    f"{'max_abs_diff_vs_cpp':>22}"
  )
  for device in args.devices:
    bc_plain = TorchBicop.from_bicop(cpp, cache_integrals=False).to(device)
    bc_cached = TorchBicop.from_bicop(cpp, cache_integrals=True).to(device)
    for n in args.sizes:
      u = rng.uniform(0.02, 0.98, size=(n, 2))
      ref = cpp.hinv1(u)
      u_t = torch.as_tensor(u, dtype=torch.float64, device=device)
      methods: dict[str, Callable[[], torch.Tensor]] = {
        "closed_form": lambda: bc_plain.hinv1(u_t),
        "cached": lambda: bc_cached.hinv1(u_t),
        "itp": lambda: _hinv1_itp(bc_plain, u_t),
      }
      for name, fn in methods.items():
        ms = _time(fn, device, args.reps)
        diff = float(np.max(np.abs(fn().cpu().numpy() - ref)))
        print(f"{device:<8}{n:>10}{name:>14}{ms:>12.2f}{diff:>22.2e}")


if __name__ == "__main__":
  main()
