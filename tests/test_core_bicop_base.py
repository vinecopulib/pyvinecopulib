"""Tests for the backend-neutral ``pyvinecopulib.core`` pair-copula base.

Exercises :class:`pyvinecopulib.core.BicopBase` — the canonical,
array-backend-agnostic partial implementation of the ``BicopLike`` contract —
purely on NumPy, so it also confirms that the neutral ``core`` layer runs
without PyTorch. A separate subprocess test pins the torch-free import
guarantee (downstream packages can build custom pairs on ``BicopBase`` in a
torch-less environment); a conformance test pins that the nanobind ``Bicop`` /
``Vinecop`` satisfy the neutral ``BicopLike`` / ``VinecopLike`` contracts.
"""

from __future__ import annotations

import subprocess
import sys
from typing import Any, Optional

import numpy as np
import pytest

import pyvinecopulib as pv
from pyvinecopulib.core import (
  BicopBase,
  BicopLike,
  IndependencePair,
  VinecopLike,
)


class _IndepPair(BicopBase[np.ndarray]):
  """Independence pair copula (``c == 1``); inherits the BicopBase defaults.

  Implements only the abstract surface (``pdf`` / ``hfunc1`` / ``hfunc2``), so
  ``hinv1`` / ``hinv2`` / ``cdf`` come from :class:`BicopBase` and are what these
  tests exercise.
  """

  def pdf(self, u: np.ndarray, x: Optional[np.ndarray] = None) -> np.ndarray:
    return np.ones(u.shape[0], dtype=u.dtype)

  def hfunc1(self, u: np.ndarray, x: Optional[np.ndarray] = None) -> np.ndarray:
    return u[:, 1]

  def hfunc2(self, u: np.ndarray, x: Optional[np.ndarray] = None) -> np.ndarray:
    return u[:, 0]

  def _sample_uniform(self, n: int, qrng: bool, seeds: list[int]) -> np.ndarray:
    rng = np.random.default_rng(seeds[0] if seeds else 0)
    return rng.uniform(size=(n, 2))


class _SqrtPair(BicopBase[np.ndarray]):
  """Toy pair with monotone ``hfunc == (free arg)**2``.

  Non-identity h-functions, so :meth:`BicopBase.hinv1` / :meth:`BicopBase.hinv2`
  genuinely exercise the numerical (bisection) inverse.
  """

  def pdf(self, u: np.ndarray, x: Optional[np.ndarray] = None) -> np.ndarray:
    return np.ones(u.shape[0], dtype=u.dtype)

  def hfunc1(self, u: np.ndarray, x: Optional[np.ndarray] = None) -> np.ndarray:
    return u[:, 1] ** 2

  def hfunc2(self, u: np.ndarray, x: Optional[np.ndarray] = None) -> np.ndarray:
    return u[:, 0] ** 2


def test_bicopbase_numerical_hinv_identity() -> None:
  """Numerical ``hinv`` inverts the identity h-functions to the target."""
  cop = _IndepPair()
  u = np.array([[0.3, 0.7], [0.5, 0.2], [0.1, 0.9]])
  np.testing.assert_allclose(cop.hinv1(u), u[:, 1], atol=1e-9)
  np.testing.assert_allclose(cop.hinv2(u), u[:, 0], atol=1e-9)


def test_bicopbase_numerical_hinv_nontrivial() -> None:
  """Numerical ``hinv`` solves ``x**2 = p`` -> ``sqrt(p)`` via bisection."""
  cop = _SqrtPair()
  p = np.array([0.04, 0.25, 0.81])
  u1 = np.full_like(p, 0.5)
  np.testing.assert_allclose(
    cop.hinv1(np.stack([u1, p], axis=-1)), np.sqrt(p), atol=1e-9
  )
  np.testing.assert_allclose(
    cop.hinv2(np.stack([p, u1], axis=-1)), np.sqrt(p), atol=1e-9
  )


def test_bicopbase_cdf_raises() -> None:
  """The base ``cdf`` raises (the vine cdf uses Monte-Carlo, not per-pair)."""
  cop = _IndepPair()
  with pytest.raises(NotImplementedError):
    cop.cdf(np.array([[0.5, 0.5]]))


def test_bicopbase_loglik() -> None:
  """``loglik`` sums the log-density; the independence pair gives 0."""
  cop = _IndepPair()
  u = np.array([[0.3, 0.7], [0.5, 0.5], [0.9, 0.1]])
  assert float(cop.loglik(u)) == pytest.approx(0.0, abs=1e-12)


def test_bicopbase_loglik_preserves_extreme_tail_density() -> None:
  """Valid densities below 1e-20 remain part of the likelihood."""
  ref = pv.Bicop.from_family(
    family=pv.families.gaussian, parameters=np.array([[0.99]])
  )

  class _Hosted(BicopBase[np.ndarray]):
    def pdf(self, u: np.ndarray, *, x: Any = None) -> np.ndarray:
      return np.asarray(ref.pdf(u))

    def hfunc1(self, u: np.ndarray, *, x: Any = None) -> np.ndarray:
      return np.asarray(ref.hfunc1(u))

    def hfunc2(self, u: np.ndarray, *, x: Any = None) -> np.ndarray:
      return np.asarray(ref.hfunc2(u))

  u = np.array([[1e-6, 1 - 1e-6], [1e-5, 1 - 1e-5]])
  assert np.all(ref.pdf(u) < 1e-20)
  np.testing.assert_allclose(_Hosted().loglik(u), ref.loglik(u), rtol=1e-14)


def test_bicopbase_simulate_default() -> None:
  """Default ``sample`` (inverse Rosenblatt) returns (n, 2) samples in (0, 1)."""
  cop = _IndepPair()
  s = cop.sample(50, seeds=[7])
  assert s.shape == (50, 2)
  assert bool((s > 0).all()) and bool((s < 1).all())
  # independence hfunc1 is the identity -> the sample is the base uniforms.
  base = _IndepPair()._sample_uniform(50, False, [7])
  np.testing.assert_allclose(s, base, atol=1e-9)


def test_bicopbase_simulate_requires_draw_hook() -> None:
  """``sample`` raises when the backend has not provided ``_sample_uniform``."""
  cop = _SqrtPair()
  with pytest.raises(NotImplementedError):
    cop.sample(5)


def test_bicopbase_requires_row_aligned_covariates() -> None:
  """Inherited pair operations do not broadcast a one-row conditioning design."""
  cop = _IndepPair()
  u = np.full((3, 2), 0.5)
  for x in (np.zeros(3), np.zeros((1, 1))):
    for call in (
      lambda: cop.loglik(u, x=x),
      lambda: cop.hinv1(u, x=x),
      lambda: cop.hinv2(u, x=x),
      lambda: cop.sample(3, x=x),
    ):
      with pytest.raises(ValueError, match="one row per observation|shape"):
        call()


def test_bicopbase_plot_runs() -> None:
  """The inherited ``plot`` delegates to the shared helper without error (Agg)."""
  import matplotlib.pyplot as plt

  _IndepPair().plot(plot_type="contour")
  plt.close("all")


def test_cpp_classes_satisfy_neutral_protocols() -> None:
  """The nanobind ``Bicop`` / ``Vinecop`` satisfy ``BicopLike`` / ``VinecopLike``.

  ``BicopLike`` mirrors the C++ ``Bicop`` evaluation surface (``pdf`` / ``cdf`` /
  ``hfunc1`` / ``hfunc2`` / ``hinv1`` / ``hinv2``, no ``dtype`` / ``device``), so
  a fitted C++ pair / vine conforms structurally; this guards against future
  contract drift.
  """
  bicop = pv.Bicop(family=pv.families.indep)
  assert isinstance(bicop, BicopLike)

  structure = pv.RVineStructure.from_order([1, 2])
  vine = pv.Vinecop.from_structure(structure=structure, pair_copulas=[[bicop]])
  assert isinstance(vine, VinecopLike)


def test_core_import_is_torch_free() -> None:
  """Importing ``pyvinecopulib.core`` must not pull in PyTorch."""
  code = (
    "import sys; import pyvinecopulib.core; "
    "sys.exit(0 if 'torch' not in sys.modules else 1)"
  )
  result = subprocess.run(  # noqa: S603
    [sys.executable, "-c", code], capture_output=True, text=True
  )
  assert result.returncode == 0, result.stderr


def test_conditioning_matrix_is_keyword_only() -> None:
  """``x`` must not be passable where ``Bicop`` expects ``parameters``.

  ``BicopLike`` is ``runtime_checkable``, so ``pv.Bicop`` satisfies it on
  method names alone -- while its second positional argument is per-row
  ``parameters``, not a conditioning matrix. If the cascade passed ``x``
  positionally, hosting a ``pv.Bicop`` in a non-simplified vine would feed the
  conditioning values in as parameters and return a wrong density instead of
  raising.
  """
  import inspect

  from pyvinecopulib.core import BicopBase, BicopLike

  for owner in (BicopLike, BicopBase):
    for name in ("pdf", "cdf", "hfunc1", "hfunc2", "hinv1", "hinv2"):
      kind = inspect.signature(getattr(owner, name)).parameters["x"].kind
      assert kind is inspect.Parameter.KEYWORD_ONLY, f"{owner.__name__}.{name}"

  u = np.full((4, 2), 0.5)
  x = np.ones((4, 1))
  compiled = pv.Bicop(family=pv.families.clayton, parameters=np.array([[2.0]]))
  assert isinstance(compiled, BicopLike)
  # Dispatched through getattr so `ty` does not reject the call it is meant to
  # reject -- a static error here is the same guarantee, one step earlier.
  with pytest.raises(TypeError):
    getattr(compiled, "pdf")(u, x=x)


def test_independence_pair_is_the_independence_copula() -> None:
  """Every member of `IndependencePair`, against the compiled `indep` `Bicop`.

  The class is public and `VinecopBase.select` hands it out, but the torch
  vine substitutes its own grid for storage, so nothing else here calls its
  `pdf`, `cdf` or inverses. This does, on NumPy, where it is used as written.
  """
  rng = np.random.default_rng(11)
  u = rng.uniform(0.01, 0.99, size=(500, 2))
  pair: BicopLike[np.ndarray] = IndependencePair()
  ref = pv.Bicop(family=pv.families.indep)

  np.testing.assert_array_equal(pair.pdf(u), np.ones(len(u)))
  np.testing.assert_array_equal(pair.hfunc1(u), u[:, 1])
  np.testing.assert_array_equal(pair.hfunc2(u), u[:, 0])
  np.testing.assert_array_equal(pair.hinv1(u), u[:, 1])
  np.testing.assert_array_equal(pair.hinv2(u), u[:, 0])
  # `cdf` is the one product, so it rounds where the others cannot.
  np.testing.assert_allclose(pair.cdf(u), ref.cdf(u), rtol=0.0, atol=2.3e-16)
  for name in ("pdf", "hfunc1", "hfunc2", "hinv1", "hinv2"):
    np.testing.assert_array_equal(
      getattr(pair, name)(u), getattr(ref, name)(u), err_msg=name
    )

  # Symmetric, so flipping is a no-op rather than a new object's behavior.
  assert pair.flip() is pair
  assert repr(pair) == "IndependencePair()"

  # The public concrete pair also fulfills the sampling member of BicopLike.
  np.testing.assert_array_equal(
    pair.sample(20, seeds=[7]), ref.sample(20, seeds=[7])
  )

  # A wider layout is accepted: the extra left-limit columns are ignored,
  # which is what a discrete edge below the threshold would hand it.
  wide = np.hstack([u, u - 1e-3])
  np.testing.assert_array_equal(pair.pdf(wide), np.ones(len(u)))
  np.testing.assert_array_equal(pair.hfunc1(wide), u[:, 1])
