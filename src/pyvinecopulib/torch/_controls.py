"""Fit-time controls for :class:`TorchBicop` and :class:`TorchVinecop`.

Mirrors :class:`pyvinecopulib.FitControlsBicop` /
:class:`pyvinecopulib.FitControlsVinecop`: method-specific args live on
the dataclass; cross-cutting args stay on the relevant ``from_data``
signature only where they don't fit naturally on the controls.

Adding a new fitter to the torch backend only requires extending the
relevant dataclass and the dispatch in the corresponding ``from_data``
— the public ``from_data`` signatures are forward-stable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

METHODS: tuple[str, ...] = ("tll", "vdc")
_VINE_IMPLS: tuple[str, ...] = ("legacy", "lazy")


@dataclass
class FitControlsTorchBicop:
  """Controls for :meth:`TorchBicop.from_data`.

  Mirrors :class:`pyvinecopulib.FitControlsBicop` on the torch side:
  ``method`` picks the pair-copula fitter and the remaining fields
  carry method-specific hyperparameters.

  Parameters
  ----------
  method:
    Which fitter to use.

    * ``"tll"`` (default) — pure-torch *Transformed Local Likelihood*
      kernel density estimator on a 2-D grid in the inverse-normal
      transformed copula space (Geenens 2014; Nagler 2018). Matches
      the C++ ``pv.Bicop.from_data(u, family_set=[pv.families.tll])``
      fit to machine precision; this is the same family the C++
      library exposes as :data:`pyvinecopulib.families.tll` and the
      pyvinecopulib default for non-parametric copula estimation.
    * ``"vdc"`` — the pretrained amortized vine-copula estimator of
      Safaai (2026), *Amortized Vine Copulas for High-Dimensional
      Density and Information Estimation*, arXiv:2604.20568
      (<https://arxiv.org/abs/2604.20568>). Reference implementation:
      `KempnerInstitute/vine-denoising-copula
      <https://github.com/KempnerInstitute/vine-denoising-copula>`_.
      Not on PyPI yet; install via the ``[vdc]`` extra under uv
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


@dataclass
class FitControlsTorchVinecop:
  """Controls for :meth:`TorchVinecop.from_data` and for the runtime
  pdf / rosenblatt cascade.

  Mirrors :class:`pyvinecopulib.FitControlsVinecop` in the sense that it
  bundles all vine-fit knobs into one object: a nested
  :class:`FitControlsTorchBicop` controls how each pair-copula is fit,
  while ``cache_integrals`` / ``device`` / ``dtype`` are vine-level
  placement / precision knobs, and ``impl`` / ``batched`` select the
  cascade variant used by the post-fit evaluation entry points.

  Parameters
  ----------
  bicop_controls:
    Controls applied to every pair-copula fit. Defaults to
    ``FitControlsTorchBicop()`` (TLL on a 30×30 normal-spaced grid).
  cache_integrals:
    If ``True``, precompute the cdf/hfunc/hinv caches on every pair
    copula's interpolation grid; see :class:`TorchBicop` for the
    accuracy/speed trade-off.
  device:
    Target torch device for the fitted pair copulas; ``None`` keeps the
    input's device.
  dtype:
    Target torch dtype; ``None`` defaults to ``torch.float64``
    (parity with the C++ implementation).
  impl:
    Cascade implementation used by :meth:`TorchVinecop.pdf` /
    :meth:`rosenblatt` / :meth:`inverse_rosenblatt`. ``"legacy"`` is
    the dense-scratch port of the C++ cascade; ``"lazy"`` is the
    dict-based variant with ref-counted GC. Both produce numerically
    identical outputs.
  batched:
    If ``True``, the cascade fires a single batched bicop call per
    tree level. Orthogonal to ``impl``; available for ``pdf`` /
    ``rosenblatt`` only (``inverse_rosenblatt(batched=True)`` raises).
  """

  bicop_controls: FitControlsTorchBicop = field(
    default_factory=FitControlsTorchBicop
  )
  cache_integrals: bool = False
  device: Optional[Any] = None
  dtype: Optional[Any] = None
  impl: str = "legacy"
  batched: bool = False

  def __post_init__(self) -> None:
    if self.impl not in _VINE_IMPLS:
      raise ValueError(
        f"unknown impl={self.impl!r}; expected one of {_VINE_IMPLS}"
      )
