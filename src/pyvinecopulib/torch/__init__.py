"""PyTorch backend for the TLL / kernel bivariate and vine copulas.

This subpackage mirrors :mod:`pyvinecopulib.core`'s evaluation chain in
pure PyTorch. It exposes three classes:

* :class:`TorchBicop` — a ``torch.nn.Module`` evaluator for a single TLL
  pair copula. Built either from a fitted :class:`pyvinecopulib.Bicop`
  (:meth:`TorchBicop.from_bicop`) or directly from data via a pure-torch
  TLL fit (:meth:`TorchBicop.from_data`). Exposes ``pdf`` / ``cdf`` /
  ``hfunc1`` / ``hfunc2`` / ``hinv1`` / ``hinv2`` / ``sample``.

* :class:`TorchVinecop` — a ``torch.nn.Module`` evaluator for an R-vine
  built on top of :class:`TorchBicop` pair copulas. Provides ``pdf`` /
  ``rosenblatt`` / ``inverse_rosenblatt``. Two implementations are
  available: a direct port of the C++ ``Vinecop`` cascade
  (``impl="legacy"``, byte-for-byte agreement with
  :class:`pyvinecopulib.Vinecop`) and a dict-based reformulation
  (``impl="lazy"``) following Cheng, Vatter, Nagler & Chen (2025).

* :class:`InterpolationGrid2D` — the underlying bilinear-interpolation
  grid (re-exported for advanced users). Provides static factories for
  the supported storage grids (``"normal"``, Phi-spaced like C++; or
  ``"linear"``, uniform on [0, 1] with O(1) cell-finding).

The whole evaluation chain stays in PyTorch, so models can move to GPU
with ``.to(device)`` and compose with autograd-aware downstream code.

References
----------
* The ``lazy`` and ``batched`` evaluation strategies are inspired by
  Cheng, Vatter, Nagler & Chen, *Vine Copulas as Differentiable
  Computational Graphs*, arXiv:2506.13318, 2025.

Installation
------------
Requires PyTorch. Install with ``pip install pyvinecopulib[torch]``.
"""

try:
  import torch  # noqa: F401
except ImportError as e:
  raise ImportError(
    "pyvinecopulib.torch requires PyTorch. "
    "Install it with `pip install pyvinecopulib[torch]`."
  ) from e

from .bicop import TorchBicop
from .vinecop import TorchVinecop
from ._interp import InterpolationGrid2D

__all__ = ["TorchBicop", "TorchVinecop", "InterpolationGrid2D"]
