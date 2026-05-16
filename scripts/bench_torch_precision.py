#!/usr/bin/env python
"""Precision-vs-truth benchmark for the torch backend.

Section 1 — Bicop precision. For each (family, params, n):
  - sample u_true ~ cop_true.simulate(n) from the parametric family,
  - fit two TorchBicops on u_true: one with cache_integrals=False, one
    with cache_integrals=True,
  - on M=10k iid eval points in [0.02, 0.98]^2 compute the integrated
    absolute error (IAE = mean |fit - true|) for pdf, cdf, hfunc1,
    hfunc2, hinv1, hinv2, where "true" comes from the parametric
    cop_true.

Section 2 — Vine precision. For each d:
  - simulate a random R-vine structure,
  - build a parametric Gaussian vine with ρ_t = 0.6 · 0.7^t per tree,
  - sample u_true ~ cop_true.simulate(n) and fit
    TorchVinecop.from_data(u_true, structure=cop_true.structure) twice
    (cache=False, cache=True),
  - compute IAE = mean |fit.pdf - true.pdf| on M=10k iid points.

Outputs a long-format CSV to --output (default stdout) with columns:
    section, family, params, d, n, quantity, cache, grid_type, IAE
"""

from __future__ import annotations

import argparse
import csv
import sys

import numpy as np
import torch

import pyvinecopulib as pv
from pyvinecopulib.torch import TorchBicop, TorchVinecop


# --- shared helpers ------------------------------------------------------ #


def _u_eval(rng: np.random.Generator, m: int, d: int) -> np.ndarray:
  return rng.uniform(0.02, 0.98, size=(m, d))


def _iae(fit: np.ndarray, truth: np.ndarray) -> float:
  return float(np.mean(np.abs(fit - truth)))


# --- Section 1: bicop ---------------------------------------------------- #


_BICOP_SPECS = [
  ("gaussian", pv.families.gaussian, [0.3]),
  ("gaussian", pv.families.gaussian, [0.7]),
  ("clayton", pv.families.clayton, [2.0]),
  ("clayton", pv.families.clayton, [5.0]),
]


def _bicop_truth(cop: pv.Bicop, u: np.ndarray) -> dict[str, np.ndarray]:
  return {
    "pdf": cop.pdf(u),
    "cdf": cop.cdf(u),
    "hfunc1": cop.hfunc1(u),
    "hfunc2": cop.hfunc2(u),
    "hinv1": cop.hinv1(u),
    "hinv2": cop.hinv2(u),
  }


def _bicop_fit(bc: TorchBicop, u_t: torch.Tensor) -> dict[str, np.ndarray]:
  return {
    "pdf": bc.pdf(u_t).numpy(),
    "cdf": bc.cdf(u_t).numpy(),
    "hfunc1": bc.hfunc1(u_t).numpy(),
    "hfunc2": bc.hfunc2(u_t).numpy(),
    "hinv1": bc.hinv1(u_t).numpy(),
    "hinv2": bc.hinv2(u_t).numpy(),
  }


def _run_bicop_section(
  n_list: list[int], grid_types: list[str], m_eval: int, seed: int
) -> list[dict]:
  rows: list[dict] = []
  rng_master = np.random.default_rng(seed)
  for fam_label, fam, params in _BICOP_SPECS:
    cop_true = pv.Bicop(family=fam, parameters=np.array([params]))
    for n in n_list:
      sim_seeds = list(
        int(x)
        for x in rng_master.integers(1, 2**31 - 1, size=3, endpoint=False)
      )
      u_true = cop_true.simulate(n, seeds=sim_seeds)
      eval_seed = int(rng_master.integers(1, 2**31 - 1))
      u_eval = _u_eval(np.random.default_rng(eval_seed), m_eval, d=2)
      u_eval_t = torch.from_numpy(u_eval)

      truths = _bicop_truth(cop_true, u_eval)
      for grid_type in grid_types:
        for cache in (False, True):
          bc_fit = TorchBicop.from_data(
            torch.from_numpy(u_true),
            cache_integrals=cache,
            grid_type=grid_type,
          )
          fits = _bicop_fit(bc_fit, u_eval_t)
          for q in truths:
            rows.append(
              {
                "section": "bicop",
                "family": fam_label,
                "params": ";".join(f"{p:g}" for p in params),
                "d": 2,
                "n": n,
                "quantity": q,
                "cache": int(cache),
                "grid_type": grid_type,
                "IAE": _iae(fits[q], truths[q]),
              }
            )
      print(
        f"# bicop {fam_label}({params}) n={n} done", file=sys.stderr, flush=True
      )
  return rows


# --- Section 2: vine ----------------------------------------------------- #


def _build_gaussian_vine(d: int, structure_seed: int) -> pv.Vinecop:
  structure = pv.RVineStructure.simulate(d, seeds=[structure_seed])
  pair_copulas = []
  for t in range(d - 1):
    rho_t = 0.6 * (0.7**t)
    pair_copulas.append(
      [
        pv.Bicop(family=pv.families.gaussian, parameters=np.array([[rho_t]]))
        for _ in range(d - t - 1)
      ]
    )
  return pv.Vinecop.from_structure(
    structure=structure, pair_copulas=pair_copulas
  )


def _run_vine_section(
  d_list: list[int], grid_types: list[str], n: int, m_eval: int, seed: int
) -> list[dict]:
  rows: list[dict] = []
  rng_master = np.random.default_rng(seed)
  for d in d_list:
    structure_seed = int(rng_master.integers(1, 2**31 - 1))
    cop_true = _build_gaussian_vine(d, structure_seed=structure_seed)
    sim_seeds = list(
      int(x) for x in rng_master.integers(1, 2**31 - 1, size=3, endpoint=False)
    )
    u_true = cop_true.simulate(n, seeds=sim_seeds)
    eval_seed = int(rng_master.integers(1, 2**31 - 1))
    u_eval = _u_eval(np.random.default_rng(eval_seed), m_eval, d=d)
    u_eval_t = torch.from_numpy(u_eval)
    truth_pdf = cop_true.pdf(u_eval)

    for grid_type in grid_types:
      for cache in (False, True):
        vc_fit = TorchVinecop.from_data(
          torch.from_numpy(u_true),
          cop_true.structure,
          cache_integrals=cache,
          grid_type=grid_type,
        )
        fit_pdf = vc_fit.pdf(u_eval_t).numpy()
        rows.append(
          {
            "section": "vine",
            "family": "gaussian",
            "params": "rho_t=0.6*0.7^t",
            "d": d,
            "n": n,
            "quantity": "pdf",
            "cache": int(cache),
            "grid_type": grid_type,
            "IAE": _iae(fit_pdf, truth_pdf),
          }
        )
    print(f"# vine d={d} done", file=sys.stderr, flush=True)
  return rows


# --- driver -------------------------------------------------------------- #


def _parse_int_list(s: str) -> list[int]:
  return [int(x) for x in s.split(",") if x.strip()]


def main() -> None:
  ap = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
  )
  ap.add_argument(
    "--n-bicop",
    default="500,2000,10000",
    type=_parse_int_list,
    help="Sample sizes for the bicop section (default: 500,2000,10000).",
  )
  ap.add_argument(
    "--d-vine",
    default="5,10,20",
    type=_parse_int_list,
    help="Vine dimensions for the vine section (default: 5,10,20).",
  )
  ap.add_argument(
    "--n-vine",
    default=10000,
    type=int,
    help="Sample size for vine fits (default: 10000).",
  )
  ap.add_argument(
    "--m-eval",
    default=10000,
    type=int,
    help="Number of iid eval points for IAE Monte-Carlo (default: 10000).",
  )
  ap.add_argument("--seed", default=42, type=int)
  ap.add_argument(
    "--sections",
    default="bicop,vine",
    help="Which sections to run (default: bicop,vine).",
  )
  ap.add_argument(
    "--grid-types",
    default="normal,linear",
    help="Storage-grid types to sweep (default: normal,linear).",
  )
  ap.add_argument("--output", default="-")
  args = ap.parse_args()

  torch.set_num_threads(1)
  sections = {s.strip() for s in args.sections.split(",") if s.strip()}
  grid_types = [g.strip() for g in args.grid_types.split(",") if g.strip()]

  rows: list[dict] = []
  if "bicop" in sections:
    rows += _run_bicop_section(args.n_bicop, grid_types, args.m_eval, args.seed)
  if "vine" in sections:
    rows += _run_vine_section(
      args.d_vine, grid_types, args.n_vine, args.m_eval, args.seed + 1
    )

  out = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
  writer = csv.DictWriter(
    out,
    fieldnames=[
      "section",
      "family",
      "params",
      "d",
      "n",
      "quantity",
      "cache",
      "grid_type",
      "IAE",
    ],
  )
  writer.writeheader()
  for r in rows:
    r["IAE"] = f"{r['IAE']:.6e}"
    writer.writerow(r)
  if args.output != "-":
    out.close()


if __name__ == "__main__":
  main()
