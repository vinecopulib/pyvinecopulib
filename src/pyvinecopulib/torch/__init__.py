"""PyTorch backend for vine copula evaluation.

This subpackage is a pure-PyTorch port of the evaluation chain in
:mod:`pyvinecopulib.core`. Pick it when you need any of:

* **GPU placement** — every class is a :class:`torch.nn.Module`, so
  ``.to("cuda")`` moves an entire vine to the GPU in one line.
* **Autograd compatibility** — the cascade is built from differentiable
  ops, so the joint density / Rosenblatt outputs flow gradients back to
  any upstream parameters (e.g. learned marginals or transforms).
* **Composition with PyTorch pipelines** — the modules drop in next to
  any other ``torch.nn.Module`` and respect the standard ``train`` /
  ``eval`` / ``state_dict`` protocol.

The default fits use the **TLL** family — *Transformed Local
Likelihood* (Geenens 2014 [2]_; Nagler 2018 [3]_) — a non-parametric
pair-copula estimator that fits a kernel density on a grid in the
inverse-normal-transformed copula space. This is the same family the
C++ library exposes as :data:`pyvinecopulib.families.tll` and is the
default everywhere in pyvinecopulib because it captures arbitrary
non-Gaussian-like dependence without picking a parametric form.

If you have not used vine copulas before, the
:doc:`concepts page </concepts>` introduces pair copulas and R-vines
in ~5 minutes.

Requires PyTorch. Install with ``pip install pyvinecopulib[torch]``.
The optional ``method="vdc"`` path additionally requires
`vine-denoising-copula <https://github.com/KempnerInstitute/vine-denoising-copula>`_,
which is not yet on PyPI. The ``[vdc]`` extra resolves it from GitHub
via ``[tool.uv.sources]`` for uv users (``uv sync --extra vdc``); for
plain pip, install directly::

    pip install "vine-denoising-copula @ git+https://github.com/KempnerInstitute/vine-denoising-copula"

See Also
--------
pyvinecopulib.core : C++/nanobind backend (default everywhere).
pyvinecopulib.sklearn : sklearn-compatible vine-copula estimators that route through either backend.

Notes
-----
**What's exposed.**

- :class:`TorchBicop` — evaluator for a single bivariate pair copula
  on a density grid. Build from a fitted
  :class:`pyvinecopulib.Bicop` via :meth:`TorchBicop.from_bicop`, or
  fit directly from data via :meth:`TorchBicop.from_data` (which
  dispatches on :class:`FitControlsTorchBicop`: pure-torch TLL by
  default; ``method="vdc"`` plugs in the optional amortized vine
  copula estimator — see *VDC* below). Exposes the standard surface:
  ``pdf`` / ``cdf`` / ``hfunc1`` / ``hfunc2`` / ``hinv1`` / ``hinv2``
  / ``simulate``.

- :class:`TorchVinecop` — evaluator for a full R-vine built on top of
  :class:`TorchBicop` pair copulas. Provides ``pdf`` / ``cdf`` /
  ``rosenblatt`` / ``inverse_rosenblatt`` / ``simulate`` with the same
  signatures as :class:`pyvinecopulib.Vinecop`. Two equivalent
  cascade implementations are available (``impl="legacy"`` mirrors the
  C++ cascade byte-for-byte; ``impl="lazy"`` is the dict-based
  reformulation of Cheng, Vatter, Nagler & Chen, 2025 [1]_).

- :class:`FitControlsTorchBicop`, :class:`FitControlsTorchVinecop` —
  fit-time control dataclasses mirroring
  :class:`pyvinecopulib.FitControlsBicop` /
  :class:`pyvinecopulib.FitControlsVinecop`. Bundle method selection
  and method-specific knobs.

- :class:`InterpolationGrid2D` — the bilinear-interpolation grid that
  backs :class:`TorchBicop` (re-exported for advanced users).

**VDC — the optional amortized estimator.** ``method="vdc"`` swaps
the per-fit TLL kernel for the pretrained amortized estimator of
Safaai (2026), reference [4]_ below. Reference implementation:
`KempnerInstitute/vine-denoising-copula <https://github.com/KempnerInstitute/vine-denoising-copula>`_.
VDC is useful when you want to skip per-pair bandwidth selection
entirely (the network is pretrained on a wide family of copulas).
See the dataclass field docs on :class:`FitControlsTorchBicop` for
checkpoint / DDIM / CFG knobs.

References
----------
.. [1] Cheng, B., Vatter, T., Nagler, T. & Chen, V. (2025).
       *Vine Copulas as Differentiable Computational Graphs.*
       arXiv:2506.13318 — basis for ``impl="lazy"`` / ``batched=True``.
.. [2] Geenens, G. (2014). *Probit Transformation for Kernel Density
       Estimation on the Unit Interval.* JASA 109(505), 346–358 —
       original TLL motivation.
.. [3] Nagler, T. (2018). *A Generic Approach to Nonparametric
       Function Estimation with Mixed Data.* Statistics & Probability
       Letters 137, 326–330 — multivariate TLL.
.. [4] Safaai, H. (2026). *Amortized Vine Copulas for High-Dimensional
       Density and Information Estimation.* arXiv:2604.20568 — VDC.
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
from ._controls import FitControlsTorchBicop, FitControlsTorchVinecop

__all__ = [
  "TorchBicop",
  "TorchVinecop",
  "InterpolationGrid2D",
  "FitControlsTorchBicop",
  "FitControlsTorchVinecop",
]
