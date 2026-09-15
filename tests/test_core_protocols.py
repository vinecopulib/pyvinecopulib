"""The four contracts, and that their reference implementations satisfy them.

``runtime_checkable`` compares member *names* only, so ``isinstance`` says
nothing about signatures. What this file pins is the other half: that ``Bicop``,
``Vinecop`` and ``Kde1d`` satisfy their contracts **structurally**, which is
what lets a consumer type against the contract instead of a concrete class.

The static half is checked by ``ty`` rather than at runtime -- the annotated
assignments below are the assertion, and a signature drifting apart from its
contract fails ``make check``. The runtime half is asserted here as well, since
the two can disagree.
"""

from typing import Any, Optional

import numpy as np
import pytest

import pyvinecopulib as pv
from pyvinecopulib.core import (
  BicopBase,
  BicopLike,
  MarginBase,
  MarginLike,
  VinecopLike,
  VinedistLike,
)


def test_the_compiled_classes_satisfy_their_contracts_structurally() -> None:
  """``ty`` checks these assignments; the ``isinstance`` calls are the twin."""
  bicop: BicopLike[np.ndarray] = pv.Bicop(
    family=pv.families.clayton, parameters=np.array([[2.0]])
  )
  vinecop: VinecopLike[np.ndarray] = pv.Vinecop.from_structure(
    structure=pv.RVineStructure.from_order([1, 2, 3])
  )
  margin: MarginLike[np.ndarray] = pv.core.Kde1d()
  dist: VinedistLike[np.ndarray] = pv.Vinedist(
    vinecop, margins=[pv.core.Kde1d() for _ in range(3)]
  )
  for obj, proto in (
    (bicop, BicopLike),
    (vinecop, VinecopLike),
    (margin, MarginLike),
    (dist, VinedistLike),
  ):
    assert isinstance(obj, proto), type(obj).__name__


def test_the_contracts_carry_no_covariate_argument() -> None:
  """``x`` lives on the bases, which is what lets the compiled classes conform.

  A protocol may ask for *fewer* parameters than an implementation provides,
  never more -- so declaring ``x`` here would put ``Bicop``, ``Vinecop`` and
  ``Kde1d``, none of which model covariates, outside their own contracts.
  """
  import inspect

  for proto in (BicopLike, VinecopLike, MarginLike, VinedistLike):
    for name in dir(proto):
      member = getattr(proto, name, None)
      if name.startswith("_") or not callable(member):
        continue
      assert "x" not in inspect.signature(member).parameters, f"{proto}.{name}"


def test_an_unconditional_margin_needs_no_covariate_parameter() -> None:
  """The documented two-primitive margin, written the way the docs write it."""

  class ShiftedExp(MarginBase[np.ndarray]):
    def __init__(self, rate: float = 1.0, shift: float = 0.0) -> None:
      self.rate, self.shift = rate, shift

    @property
    def support(self) -> tuple[float, float]:
      return (self.shift, float("inf"))

    def pdf(self, y: Any) -> Any:
      return self.rate * np.exp(-self.rate * (y - self.shift))

    def cdf(self, y: Any) -> Any:
      return 1.0 - np.exp(-self.rate * (y - self.shift))

  m = ShiftedExp(rate=2.0, shift=1.0)
  assert isinstance(m, MarginLike)
  np.testing.assert_allclose(m.icdf(np.array([0.5])), [1.34657359], rtol=1e-8)


def test_an_unconditional_pair_needs_no_covariate_parameter() -> None:
  """The same, one level up: three primitives, none of them declaring ``x``."""

  class Independence(BicopBase[np.ndarray]):
    def _pdf_raw(self, u: Any) -> Any:
      return np.ones(u.shape[0])

    def _hfunc1_raw(self, u: Any) -> Any:
      return u[:, 1]

    def _hfunc2_raw(self, u: Any) -> Any:
      return u[:, 0]

  pair = Independence()
  assert isinstance(pair, BicopLike)
  u = np.array([[0.3, 0.7], [0.5, 0.5]])
  np.testing.assert_array_equal(pair.pdf(u), np.ones(2))
  np.testing.assert_allclose(pair.hinv1(u), u[:, 1], atol=1e-10)

  # Undeclared means "cannot be hosted in a conditional vine", and the refusal
  # is immediate rather than an unconditional answer under a conditional call.
  with pytest.raises(TypeError):
    pair.pdf(u, x=np.zeros((2, 1)))


def test_a_conditional_pair_declares_the_covariates_it_reads() -> None:
  """Widening the primitive is how a pair says it models covariates."""

  class Conditional(BicopBase[np.ndarray]):
    def _pdf_raw(self, u: Any, *, x: Optional[Any] = None) -> Any:
      scale = 1.0 if x is None else 1.0 + float(x[0, 0])
      return np.full(u.shape[0], scale)

    def _hfunc1_raw(self, u: Any, *, x: Optional[Any] = None) -> Any:
      del x
      return u[:, 1]

    def _hfunc2_raw(self, u: Any, *, x: Optional[Any] = None) -> Any:
      del x
      return u[:, 0]

  pair = Conditional()
  u = np.array([[0.3, 0.7], [0.5, 0.5]])
  np.testing.assert_array_equal(pair.pdf(u), np.ones(2))
  np.testing.assert_array_equal(
    pair.pdf(u, x=np.full((2, 1), 3.0)), np.full(2, 4.0)
  )
