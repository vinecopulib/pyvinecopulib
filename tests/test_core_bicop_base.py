"""Tests for the backend-neutral ``pyvinecopulib.core`` pair-copula base.

Exercises :class:`pyvinecopulib.core.BicopBase` — the canonical,
array-backend-agnostic partial implementation of the ``BicopLike`` contract —
purely on NumPy, so it also confirms that the neutral ``core`` layer runs
without PyTorch. A separate subprocess test pins the torch-free import
guarantee (downstream packages can build custom pairs on ``BicopBase`` in a
torch-less environment).
"""

from __future__ import annotations

import subprocess
import sys
from typing import Any, Optional

import numpy as np
import pytest

from pyvinecopulib.core import BicopBase


class _IndepPair(BicopBase[np.ndarray]):
  """Independence pair copula (``c == 1``); inherits every BicopBase default.

  Implements only the abstract surface (``log_pdf`` / ``hfunc1`` / ``hfunc2`` +
  ``dtype`` / ``device``), so ``pdf`` / ``hinv1`` / ``hinv2`` / ``cdf`` come
  from :class:`BicopBase` and are what these tests exercise.
  """

  def __init__(self, ref: np.ndarray) -> None:
    self._ref: Any = ref

  @property
  def dtype(self) -> object:
    return self._ref.dtype

  @property
  def device(self) -> object:
    return self._ref.device

  def log_pdf(
    self, u: np.ndarray, cond: Optional[np.ndarray] = None
  ) -> np.ndarray:
    return np.zeros(u.shape[0], dtype=u.dtype)

  def hfunc1(
    self, u: np.ndarray, cond: Optional[np.ndarray] = None
  ) -> np.ndarray:
    return u[:, 1]

  def hfunc2(
    self, u: np.ndarray, cond: Optional[np.ndarray] = None
  ) -> np.ndarray:
    return u[:, 0]


class _SqrtPair(BicopBase[np.ndarray]):
  """Toy pair with monotone ``hfunc == (free arg)**2``.

  Non-identity h-functions, so :meth:`BicopBase.hinv1` / :meth:`BicopBase.hinv2`
  genuinely exercise the numerical (bisection) inverse.
  """

  def __init__(self, ref: np.ndarray) -> None:
    self._ref: Any = ref

  @property
  def dtype(self) -> object:
    return self._ref.dtype

  @property
  def device(self) -> object:
    return self._ref.device

  def log_pdf(
    self, u: np.ndarray, cond: Optional[np.ndarray] = None
  ) -> np.ndarray:
    return np.zeros(u.shape[0], dtype=u.dtype)

  def hfunc1(
    self, u: np.ndarray, cond: Optional[np.ndarray] = None
  ) -> np.ndarray:
    return u[:, 1] ** 2

  def hfunc2(
    self, u: np.ndarray, cond: Optional[np.ndarray] = None
  ) -> np.ndarray:
    return u[:, 0] ** 2


def test_bicopbase_pdf_default() -> None:
  """``pdf`` defaults to ``exp(log_pdf)`` (== 1 for the independence pair)."""
  cop = _IndepPair(np.empty(0, dtype=np.float64))
  u = np.array([[0.3, 0.7], [0.5, 0.5], [0.9, 0.1]])
  np.testing.assert_allclose(cop.pdf(u), np.ones(3))


def test_bicopbase_numerical_hinv_identity() -> None:
  """Numerical ``hinv`` inverts the identity h-functions to the target."""
  cop = _IndepPair(np.empty(0, dtype=np.float64))
  u = np.array([[0.3, 0.7], [0.5, 0.2], [0.1, 0.9]])
  np.testing.assert_allclose(cop.hinv1(u), u[:, 1], atol=1e-9)
  np.testing.assert_allclose(cop.hinv2(u), u[:, 0], atol=1e-9)


def test_bicopbase_numerical_hinv_nontrivial() -> None:
  """Numerical ``hinv`` solves ``x**2 = p`` -> ``sqrt(p)`` via bisection."""
  cop = _SqrtPair(np.empty(0, dtype=np.float64))
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
  cop = _IndepPair(np.empty(0, dtype=np.float64))
  with pytest.raises(NotImplementedError):
    cop.cdf(np.array([[0.5, 0.5]]))


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
