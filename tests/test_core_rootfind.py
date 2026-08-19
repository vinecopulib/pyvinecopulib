"""Tests for the backend-agnostic monotone root-finder.

`solve_increasing` backs two different inverses: the pair-copula h-inverses,
which always search the unit interval, and the marginal `icdf`, which may have
to search an unbounded support. The unit-interval path is load-bearing for
torch<->C++ parity, so the first test pins it against a literal transcription of
the bisection loop; the rest cover the bracket search that unbounded supports
need.
"""

from __future__ import annotations

from typing import Any, Callable

import numpy as np
import pytest
from array_api_compat import array_namespace

from pyvinecopulib.core._rootfind import solve_increasing


def _reference_unit_interval(
  f: Callable[[Any], Any], p: Any, n_iter: int = 50
) -> Any:
  """The unit-interval bisection, transcribed literally.

  Guards the h-inverse hot path: generalizing the bracket must not perturb it.
  """
  xp = array_namespace(p)
  a = xp.full_like(p, 0.0)
  b = xp.full_like(p, 1.0)
  for _ in range(n_iter):
    mid = 0.5 * (a + b)
    lower = f(mid) < p
    a = xp.where(lower, mid, a)
    b = xp.where(lower, b, mid)
  return xp.clip(0.5 * (a + b), 0.0, 1.0)


def _norm_cdf(x: Any) -> Any:
  from scipy.special import ndtr

  return ndtr(x)


def test_unit_interval_is_bit_for_bit_unchanged() -> None:
  """The default bracket reproduces the previous implementation exactly."""
  p = np.linspace(1e-10, 1 - 1e-10, 97)
  for f in (lambda x: x, lambda x: x**3, lambda x: np.sqrt(x)):
    got = solve_increasing(f, p)
    want = _reference_unit_interval(f, p)
    assert np.array_equal(got, want)


def test_identity_recovers_the_target() -> None:
  """Inverting the identity returns the target to bisection precision."""
  p = np.linspace(0.01, 0.99, 50)
  np.testing.assert_allclose(solve_increasing(lambda x: x, p), p, atol=1e-14)


def test_unbounded_support_inverts_the_normal_cdf() -> None:
  """An infinite bracket on both sides is widened until it brackets."""
  pytest.importorskip("scipy")
  from scipy.special import ndtri

  p = np.array([1e-8, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-3, 1 - 1e-8])
  got = solve_increasing(_norm_cdf, p, lo=-np.inf, hi=np.inf)
  np.testing.assert_allclose(got, ndtri(p), atol=1e-8)


def test_half_bounded_support_inverts_the_exponential_cdf() -> None:
  """A finite lower bound is kept while the upper one is searched."""
  p = np.array([0.01, 0.25, 0.5, 0.9, 0.999])
  got = solve_increasing(
    lambda x: 1.0 - np.exp(-x), p, lo=0.0, hi=np.inf, n_iter=80
  )
  np.testing.assert_allclose(got, -np.log1p(-p), atol=1e-10)
  assert np.all(got >= 0.0)


def test_heavy_tail_needs_many_doublings() -> None:
  """A Cauchy tail forces a wide search and still converges."""
  p = np.array([1e-6, 1e-4, 0.5, 1 - 1e-4])
  got = solve_increasing(
    lambda x: 0.5 + np.arctan(x) / np.pi, p, lo=-np.inf, hi=np.inf, n_iter=200
  )
  # atol covers p = 0.5, where the exact answer is 0 and rtol is meaningless.
  np.testing.assert_allclose(
    got, np.tan(np.pi * (p - 0.5)), rtol=1e-6, atol=1e-9
  )


def test_array_valued_brackets_broadcast() -> None:
  """Per-element supports are honored independently."""
  lo = np.array([0.0, 10.0, -5.0])
  hi = np.array([1.0, 11.0, 5.0])
  p = np.array([0.5, 10.5, 0.0])
  # f is the identity, so the solution is the target itself.
  got = solve_increasing(lambda x: x, p, lo=lo, hi=hi)
  np.testing.assert_allclose(got, p, atol=1e-12)


def test_result_is_clamped_to_the_bracket() -> None:
  """An unreachable target saturates at the endpoint rather than escaping."""
  p = np.array([-1.0, 2.0])
  got = solve_increasing(lambda x: x, p, lo=0.0, hi=1.0)
  assert np.all(got >= 0.0) and np.all(got <= 1.0)


def test_infinite_bracket_without_expansion_raises() -> None:
  """Refused rather than silently returning nan."""
  with pytest.raises(ValueError, match="max_expand"):
    solve_increasing(_norm_cdf, np.array([0.5]), lo=-np.inf, max_expand=0)


def test_torch_backend_matches_numpy() -> None:
  """The same code path runs on torch tensors and agrees with numpy."""
  torch = pytest.importorskip("torch")
  p_np = np.array([0.05, 0.5, 0.95])
  p_t = torch.tensor(p_np, dtype=torch.float64)

  got_np = solve_increasing(lambda x: 1.0 - np.exp(-x), p_np, lo=0.0, hi=np.inf)
  got_t = solve_increasing(
    lambda x: 1.0 - torch.exp(-x),
    p_t,
    lo=0.0,
    hi=float("inf"),
  )
  np.testing.assert_allclose(got_t.numpy(), got_np, atol=1e-12)
