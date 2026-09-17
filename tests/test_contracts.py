"""One conformance suite per level, over every implementation that ships.

The library has four levels -- margin, pair copula, vine, distribution -- and
each states its contract as a ``runtime_checkable`` protocol. What a consumer
may ask of a part is therefore the same question whatever the part is, and
asking it once per implementation is what this file does.

It replaces the same clauses asserted per class: `isinstance(_, <Protocol>)`
appeared at 23 sites across 9 files, and the margin identities
(``cdf(icdf(p)) == p``, ``logpdf == log(pdf)``, a 0-d ``loglik``) at 5 more.
`tests/test_margins_contract.py` already proves the shape works; this widens it
from the two parametric adapters to every level.

What stays in the per-class files is what is *not* contract: a family's
arithmetic, a lane's `nn.Module` behavior, an implementation's own errors.
"""

from __future__ import annotations

from typing import Any, NamedTuple

import numpy as np
import pytest

import pyvinecopulib as pv
from pyvinecopulib.core import (
  BicopLike,
  IndependenceBicop,
  Kde1d,
  MarginLike,
  VinecopLike,
  Vinedist,
  VinedistLike,
)

_RNG = np.random.default_rng(0)
_N = 400
_D = 3
_Y = _RNG.normal(size=(_N, _D))
_U = pv.to_pseudo_obs(_Y)
_P = np.array([0.1, 0.25, 0.5, 0.75, 0.9])


class Impl(NamedTuple):
  """One implementation of one contract, and how to build a fitted instance."""

  name: str
  build: Any
  extra: str | None = None


def _built(impl: Impl) -> Any:
  """Build it, skipping where its optional extra is absent."""
  if impl.extra is not None:
    pytest.importorskip(impl.extra)
  return impl.build()


def _torch_u() -> Any:
  import torch

  return torch.as_tensor(_U, dtype=torch.float64)


# --- the four tables --------------------------------------------------------- #

MARGINS = [
  Impl("Kde1d", lambda: Kde1d().fit(_Y[:, 0])),
  Impl(
    "SciPyMargin",
    lambda: (
      __import__("pyvinecopulib.margins", fromlist=["SciPyMargin"])
      .SciPyMargin("norm")
      .fit(_Y[:, 0])
    ),
    extra="scipy",
  ),
  Impl(
    "TorchKde1d",
    lambda: (
      __import__("pyvinecopulib.torch", fromlist=["TorchKde1d"])
      .TorchKde1d()
      .fit(__import__("torch").as_tensor(_Y[:, 0]))
    ),
    extra="torch",
  ),
]

BICOPS = [
  Impl("Bicop", lambda: pv.Bicop.from_data(_U[:, :2])),
  Impl("IndependenceBicop", IndependenceBicop),
  Impl(
    "TorchTllBicop",
    lambda: __import__(
      "pyvinecopulib.torch", fromlist=["TorchTllBicop"]
    ).TorchTllBicop.from_data(_torch_u()[:, :2]),
    extra="torch",
  ),
]

VINECOPS = [
  Impl("Vinecop", lambda: pv.Vinecop.from_data(_U)),
  Impl(
    "TorchVinecop",
    lambda: __import__(
      "pyvinecopulib.torch", fromlist=["TorchVinecop"]
    ).TorchVinecop.from_data(_torch_u()),
    extra="torch",
  ),
]

DISTS = [
  Impl("Vinedist", lambda: Vinedist.from_data(_Y)),
  Impl(
    "TorchVinedist",
    lambda: __import__(
      "pyvinecopulib.torch", fromlist=["TorchVinedist"]
    ).TorchVinedist.from_data(__import__("torch").as_tensor(_Y)),
    extra="torch",
  ),
]

_LEVELS = [
  ("margin", MARGINS, MarginLike),
  ("bicop", BICOPS, BicopLike),
  ("vinecop", VINECOPS, VinecopLike),
  ("vinedist", DISTS, VinedistLike),
]
_ALL = [
  pytest.param(impl, proto, id=f"{level}-{impl.name}")
  for level, impls, proto in _LEVELS
  for impl in impls
]


def _numpy(value: Any) -> np.ndarray:
  """Whatever array a lane answered in, as NumPy."""
  from pyvinecopulib.core.extend import to_numpy

  return to_numpy(value, dtype=float)


# --- what every level owes --------------------------------------------------- #


@pytest.mark.parametrize(("impl", "proto"), _ALL)
def test_it_satisfies_its_contract(impl: Impl, proto: Any) -> None:
  """`isinstance` against the protocol, for every shipped implementation.

  Structural, not nominal: the compiled classes inherit nothing from these
  protocols, and satisfying them anyway is what lets a consumer type against
  the contract rather than a concrete class.
  """
  assert isinstance(_built(impl), proto)


@pytest.mark.parametrize(("impl", "proto"), _ALL)
def test_npars_is_a_number_or_says_it_has_none(impl: Impl, proto: Any) -> None:
  """`npars` always answers; `nan` is how it reports having no count.

  It may not raise: a `runtime_checkable` protocol evaluates its data members
  during `isinstance` on Python 3.11, so a member that raises makes the
  conformance check itself raise.
  """
  del proto
  value = float(_built(impl).npars)
  assert value >= 0.0 or np.isnan(value)


# --- per-level identities ---------------------------------------------------- #


@pytest.mark.parametrize("impl", MARGINS, ids=lambda i: i.name)
def test_a_margin_inverts_its_own_distribution(impl: Impl) -> None:
  """`cdf(icdf(p)) == p` on the levels every margin covers."""
  margin = _built(impl)
  p = _P if impl.extra != "torch" else __import__("torch").as_tensor(_P)
  np.testing.assert_allclose(_numpy(margin.cdf(margin.icdf(p))), _P, atol=1e-8)


@pytest.mark.parametrize("impl", MARGINS, ids=lambda i: i.name)
def test_a_margin_logs_its_own_density(impl: Impl) -> None:
  """`logpdf` is `log(pdf)` wherever the density is positive."""
  margin = _built(impl)
  y = (
    _Y[:20, 0]
    if impl.extra != "torch"
    else __import__("torch").as_tensor(_Y[:20, 0])
  )
  dens, logdens = _numpy(margin.pdf(y)), _numpy(margin.logpdf(y))
  positive = dens > 0
  np.testing.assert_allclose(
    logdens[positive], np.log(dens[positive]), rtol=1e-10
  )


@pytest.mark.parametrize("impl", BICOPS, ids=lambda i: i.name)
def test_a_pair_inverts_its_own_h_function(impl: Impl) -> None:
  """`hinv1(u1, hfunc1(u))` recovers `u2`, whatever the pair."""
  pair = _built(impl)
  u = _U[:20, :2] if impl.extra != "torch" else _torch_u()[:20, :2]
  h = pair.hfunc1(u)
  back = pair.hinv1(
    pair._prep(np.column_stack([_numpy(u)[:, 0], _numpy(h)]))
    if impl.extra == "torch"
    else np.column_stack([u[:, 0], h])
  )
  np.testing.assert_allclose(_numpy(back), _numpy(u)[:, 1], atol=1e-6)


@pytest.mark.parametrize("impl", VINECOPS, ids=lambda i: i.name)
def test_a_vine_inverts_its_own_rosenblatt(impl: Impl) -> None:
  """`inverse_rosenblatt(rosenblatt(u))` recovers `u`, whatever the lane."""
  vine = _built(impl)
  u = _U[:50] if impl.extra != "torch" else _torch_u()[:50]
  back = vine.inverse_rosenblatt(vine.rosenblatt(u))
  np.testing.assert_allclose(_numpy(back), _numpy(u), atol=1e-6)


@pytest.mark.parametrize("impl", DISTS, ids=lambda i: i.name)
def test_a_distribution_exponentiates_its_own_log_density(impl: Impl) -> None:
  """`pdf` is `exp(logpdf)`, and `loglik` is their 0-d total."""
  dist = _built(impl)
  y = (
    _Y[:20] if impl.extra != "torch" else __import__("torch").as_tensor(_Y[:20])
  )
  np.testing.assert_allclose(
    _numpy(dist.pdf(y)), np.exp(_numpy(dist.logpdf(y))), rtol=1e-10
  )
  total = dist.loglik(y)
  assert np.ndim(_numpy(total)) == 0
  np.testing.assert_allclose(
    _numpy(total), _numpy(dist.logpdf(y)).sum(), rtol=1e-10
  )
