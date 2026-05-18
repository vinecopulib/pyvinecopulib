"""Fit-time controls for :class:`TorchBicop`.

Mirrors :class:`pyvinecopulib.FitControlsBicop`: method-specific args live
on this dataclass; cross-cutting args (``cache_integrals`` / ``device`` /
``dtype``) stay on :meth:`TorchBicop.from_data`. Adding a new fitter to
the torch backend only requires extending this dataclass and the dispatch
in :meth:`TorchBicop.from_data` — the public ``from_data`` signature is
forward-stable.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

METHODS: tuple[str, ...] = ("tll", "vdc")


@dataclass
class FitControlsTorchBicop:
  """Controls for :meth:`TorchBicop.from_data`.

  Parameters
  ----------
  method:
    Which fitter to use. ``"tll"`` (default) is the pure-torch TLL KDE
    that matches the C++ ``pv.Bicop.from_data(u, family_set=[tll])`` fit
    to machine precision. ``"vdc"`` is the Kempner Institute's
    `vine-denoising-copula
    <https://github.com/KempnerInstitute/vine-denoising-copula>`_
    pretrained denoiser / diffusion estimator (arXiv 2604.20568). vdc
    is not on PyPI yet; install via the ``[vdc]`` extra under uv
    (``uv sync --extra vdc``) or with pip directly:
    ``pip install "vine-denoising-copula @ git+https://github.com/KempnerInstitute/vine-denoising-copula"``.
  grid_size:
    *TLL only.* Density grid size per axis (default ``30``; matches C++).
  mult:
    *TLL only.* Bandwidth multiplier (default ``1.0``; matches C++).
  grid_type:
    *TLL only.* Storage grid type — ``"normal"`` (default, Phi-spaced,
    byte-for-byte parity with C++) or ``"linear"`` (uniform on ``[0, 1]``
    with the O(1) cell-finding fast-path in
    :meth:`InterpolationGrid2D._cell_index`).
  vdc_bundle:
    *VDC only.* A pre-loaded ``vdc.LoadedPretrainedModel``. If ``None``,
    the loader is invoked with ``vdc_model_id`` and the result is cached
    at module scope (keyed on model_id × device) so repeated fits
    amortize the HuggingFace download.
  vdc_model_id:
    *VDC only.* Pretrained checkpoint identifier (default
    ``"vdc-denoiser-m64-v1"``; ``m=64`` density grid).
  vdc_diffusion_steps:
    *VDC only.* DDIM step count for diffusion checkpoints. ``None``
    (default) lets vdc pick its own; ignored for the denoiser
    checkpoint, which is a single forward pass.
  vdc_cfg_scale:
    *VDC only.* Classifier-free guidance scale for diffusion checkpoints
    (default ``1.0``); ignored for the denoiser checkpoint.
  vdc_projection_iters:
    *VDC only.* IPFP / Sinkhorn iterations enforcing uniform marginals
    on the cell-center grid (default ``50``). With this many iterations
    vdc's output is already at machine-precision uniform marginals, so
    the resulting :class:`TorchBicop` is built with ``norm_times=0``.
  """

  method: str = "tll"

  # TLL-only
  grid_size: int = 30
  mult: float = 1.0
  grid_type: str = "normal"

  # VDC-only
  vdc_bundle: Optional[Any] = None
  vdc_model_id: str = "vdc-denoiser-m64-v1"
  vdc_diffusion_steps: Optional[int] = None
  vdc_cfg_scale: float = 1.0
  vdc_projection_iters: int = 50

  def __post_init__(self) -> None:
    if self.method not in METHODS:
      raise ValueError(
        f"unknown method={self.method!r}; expected one of {METHODS}"
      )
