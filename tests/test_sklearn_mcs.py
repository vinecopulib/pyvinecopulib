"""Smoke tests for the (private) MCS selectors used by the forest classes.

The forest tests cover the MCS path transitively. These cases hit the
two selectors with hand-crafted loss matrices where the survivor set
is obvious by construction, so a future refactor of `_mcs.py` is
caught at a small unit scale instead of via integration drift.
"""

import numpy as np
import pytest

pytest.importorskip("sklearn")
pytest.importorskip("scipy")

from pyvinecopulib.sklearn._mcs import (  # noqa: E402
  argmin_except_self,
  da_mcs_marg,
  da_mcs_unif,
  pairwise_sd_for_pairs,
  pairwise_sd_for_pairs_einsum,
)


def test_argmin_except_self_unique_min():
  mu = np.array([0.5, 0.1, 0.7, 0.3])
  k = argmin_except_self(mu)
  # Index 1 holds the global min; every other i picks 1 as partner.
  assert k[0] == 1
  assert k[2] == 1
  assert k[3] == 1
  # Index 1 picks the second min (index 3).
  assert k[1] == 3


def test_argmin_except_self_tied_min():
  mu = np.array([0.1, 0.5, 0.1])
  k = argmin_except_self(mu)
  # Two rows attain the min; they partner with each other.
  assert k[0] == 2
  assert k[2] == 0
  assert k[1] == 0


def _losses_with_two_dominant_models(n=200, p=6, rng_seed=0):
  """Construct a loss matrix where models 0–1 dominate; 2–5 lose."""
  rng = np.random.default_rng(rng_seed)
  losses = rng.normal(loc=0.0, scale=1.0, size=(n, p))
  losses[:, 2:] += 1.5  # models 2..5 have a clearly higher mean loss
  return losses


def test_da_mcs_marg_selects_top_models():
  losses = _losses_with_two_dominant_models(n=400, p=6, rng_seed=0)
  result = da_mcs_marg(losses, alpha=0.05, rng=np.random.default_rng(123))

  assert set(result.keys()) == {"decision", "stats", "threshold"}
  assert result["decision"].shape == (6,)
  # Among the two near-equal best models, at least one must survive
  # (the DA test can flip between tied competitors when paired against
  # each other). All clearly-worse models must be excluded.
  assert result["decision"][:2].any()
  assert not result["decision"][2:].any()


def test_da_mcs_unif_selects_top_models():
  losses = _losses_with_two_dominant_models(n=400, p=6, rng_seed=1)
  result = da_mcs_unif(losses, alpha=0.05, rng=np.random.default_rng(123))

  assert set(result.keys()) == {"decision", "stats", "stats_pre"}
  assert result["decision"].shape == (6,)
  assert result["decision"][:2].any()
  assert not result["decision"][2:].any()


def test_da_mcs_marg_deterministic_split():
  losses = _losses_with_two_dominant_models(n=400, p=6, rng_seed=2)
  r1 = da_mcs_marg(losses, alpha=0.05, randomize=False)
  r2 = da_mcs_marg(losses, alpha=0.05, randomize=False)
  np.testing.assert_array_equal(r1["decision"], r2["decision"])
  np.testing.assert_allclose(r1["stats"], r2["stats"])


def test_da_mcs_unif_deterministic_split():
  losses = _losses_with_two_dominant_models(n=400, p=6, rng_seed=3)
  r1 = da_mcs_unif(losses, alpha=0.05, randomize=False)
  r2 = da_mcs_unif(losses, alpha=0.05, randomize=False)
  np.testing.assert_array_equal(r1["decision"], r2["decision"])
  np.testing.assert_allclose(r1["stats"], r2["stats"])


@pytest.mark.parametrize("n,p", [(50, 4), (400, 6), (5, 3)])
def test_pairwise_sd_helpers_agree(n, p):
  """The plain and einsum SD helpers must compute the same thing.

  `da_mcs_unif` defaults to the einsum variant for speed; the plain
  variant is the readable reference. Nothing else asserts they match,
  so this guards the intentional duplication in `_mcs.py`.
  """
  rng = np.random.default_rng(7 * n + p)
  X = rng.standard_normal((n, p))
  mu = X.mean(axis=0)
  k = argmin_except_self(mu)
  sd_plain = pairwise_sd_for_pairs(X, mu, k)
  sd_einsum = pairwise_sd_for_pairs_einsum(X, mu, k)
  np.testing.assert_allclose(sd_einsum, sd_plain, rtol=1e-10, atol=1e-12)
