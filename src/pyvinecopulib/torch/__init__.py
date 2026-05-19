"""PyTorch backend for the TLL / kernel bivariate and vine copulas.

This subpackage mirrors :mod:`pyvinecopulib.core`'s evaluation chain in
pure PyTorch. It exposes:

* :class:`TorchBicop` — a ``torch.nn.Module`` evaluator for a single
  bivariate pair copula on a density grid. Built either from a fitted
  :class:`pyvinecopulib.Bicop` (:meth:`TorchBicop.from_bicop`) or
  directly from data via :meth:`TorchBicop.from_data` (which dispatches
  on :class:`FitControlsTorchBicop` — pure-torch TLL by default;
  ``method="vdc"`` plugs in the optional `vine-denoising-copula
  <https://github.com/KempnerInstitute/vine-denoising-copula>`_
  pretrained estimator). Exposes ``pdf`` / ``cdf`` / ``hfunc1`` /
  ``hfunc2`` / ``hinv1`` / ``hinv2`` / ``simulate``.

* :class:`TorchVinecop` — a ``torch.nn.Module`` evaluator for an R-vine
  built on top of :class:`TorchBicop` pair copulas. Provides ``pdf`` /
  ``rosenblatt`` / ``inverse_rosenblatt``. Two implementations are
  available: a direct port of the C++ ``Vinecop`` cascade
  (``impl="legacy"``, byte-for-byte agreement with
  :class:`pyvinecopulib.Vinecop`) and a dict-based reformulation
  (``impl="lazy"``) following Cheng, Vatter, Nagler & Chen (2025).

* :class:`FitControlsTorchBicop` — fit-time controls for
  :meth:`TorchBicop.from_data`. Mirrors
  :class:`pyvinecopulib.FitControlsBicop`; carries the ``method``
  selector plus method-specific parameters.

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

The optional ``method="vdc"`` path additionally requires the
`vine-denoising-copula
<https://github.com/KempnerInstitute/vine-denoising-copula>`_ package,
which is not yet on PyPI. The ``[vdc]`` extra resolves it from GitHub
via ``[tool.uv.sources]`` for uv users (``uv sync --extra vdc``); for
plain pip, install directly::

    pip install "vine-denoising-copula @ git+https://github.com/KempnerInstitute/vine-denoising-copula"
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
from ._controls import FitControlsTorchBicop

__all__ = [
  "TorchBicop",
  "TorchVinecop",
  "InterpolationGrid2D",
  "FitControlsTorchBicop",
]
