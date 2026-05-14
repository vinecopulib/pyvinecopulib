"""PyTorch port of the TLL / kernel bivariate copula.

This subpackage exposes :class:`TorchBicop`, a ``torch.nn.Module`` that
takes a fitted :class:`pyvinecopulib.Bicop` (with the ``tll`` family) and
provides the standard copula evaluation API (``pdf`` / ``cdf`` / ``hfunc`` /
``hinv`` / ``sample``) on top of a PyTorch interpolation grid. The
evaluation chain stays in PyTorch end-to-end, so the bicop can be moved
to GPU and combined with autograd-aware downstream code.

Fitting is delegated to the C++ library; this subpackage is a thin
PyTorch evaluator over the resulting grid.

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
from ._interp import InterpolationGrid2D

__all__ = ["TorchBicop", "InterpolationGrid2D"]
