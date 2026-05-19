"""vdc-based fitter for :class:`TorchBicop`.

Wraps the Kempner Institute's `vine-denoising-copula
<https://github.com/KempnerInstitute/vine-denoising-copula>`_ pretrained
estimator (arXiv 2604.20568): vdc takes a ``(n, 2)`` pseudo-observation
sample and produces an ``(m, m)`` density grid on the cell-centered
grid ``(0.5/m, 1.5/m, …, (m-0.5)/m)``. We replicate-pad that grid to
``(m+2, m+2)`` on ``[0, cell-centers..., 1]`` and hand it to the existing
:class:`InterpolationGrid2D`, which then provides ``pdf`` / ``cdf`` /
``hfunc`` / ``hinv`` / sampling exactly as for the TLL fit.

Why replicate-pad works: vdc's ``copula_project`` IPFP enforces
midpoint-rule uniform marginals on the cell-center grid (``sum_i
density[i, j] * (1/m) = 1``). The trapezoidal integral of the
replicate-padded density on ``[0, 1]`` coincides with that midpoint sum
at every cell center, so marginals stay uniform without re-running
``normalize_margins`` and the h-function values at cell centers match
vdc's ``(cumsum - 0.5*density) * (1/m)`` formula exactly.

The ``vdc`` package is imported lazily so ``pyvinecopulib.torch``
imports cleanly without it. With ``uv``, install via the ``[vdc]``
extra (pyvinecopulib's ``pyproject.toml`` maps it to the GitHub repo
through ``[tool.uv.sources]``)::

    uv sync --extra vdc                  # in a pyvinecopulib checkout
    uv pip install "pyvinecopulib[vdc]"  # from an arbitrary uv project

Plain pip doesn't read ``[tool.uv.sources]``, so install vdc directly::

    pip install "vine-denoising-copula @ git+https://github.com/KempnerInstitute/vine-denoising-copula"
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional

import torch
from torch import Tensor

if TYPE_CHECKING:  # pragma: no cover
  from ._controls import FitControlsTorchBicop

# (model_id, device-str) -> LoadedPretrainedModel.
_BUNDLE_CACHE: dict[tuple[str, str], Any] = {}


def _try_import_vdc():
  """Lazy-import vine-denoising-copula with a friendly error.

  .. note::

     As of vdc 0.1.0 on GitHub main, the published wheel references
     several submodules that aren't shipped (``vdc.inference``,
     ``vdc.vine``). We work around this with a ``sys.modules`` shim
     installed by :func:`_install_upstream_shims` before the first
     :func:`vdc.load_pretrained_model` call — see that function for the
     details. Once upstream restores the missing subpackages the shim
     detects them and becomes a no-op.
  """
  try:
    import vdc as _vdc
  except (
    ImportError
  ) as e:  # pragma: no cover - exercised in tests via monkeypatch
    raise ImportError(
      "method='vdc' requires the vine-denoising-copula package, which "
      "is not yet on PyPI. Install with `uv sync --extra vdc` (uses the "
      "GitHub source declared in pyproject.toml's [tool.uv.sources]) or "
      "with pip directly:\n"
      '    pip install "vine-denoising-copula @ '
      'git+https://github.com/KempnerInstitute/vine-denoising-copula"'
    ) from e
  return _vdc


_SHIMMED = False


def _install_upstream_shims() -> None:
  """Work around the vdc 0.1.0 upstream packaging bug.

  As shipped on GitHub ``main``, the vdc wheel references two submodules
  that aren't actually packaged:

  - ``vdc.inference.density`` (imported by ``vdc.data.onthefly`` for
    ``scatter_to_hist`` — the real function lives in ``vdc.data.hist`` —
    and also referenced inside ``vdc.pretrained.estimate_pair_density_from_samples``
    for ``sample_density_grid``, which is only used in the diffusion
    branch we don't touch).
  - ``vdc.vine.copula_diffusion`` (imported by ``vdc.data.onthefly`` for
    ``DiffusionCopulaModel`` — only used at training time, never
    instantiated during inference).

  Both are referenced at *import time*, so ``vdc.load_pretrained_model``
  raises ``ModuleNotFoundError`` before doing anything useful. We inject
  ``sys.modules`` stubs that satisfy the imports without changing
  inference behaviour. Once upstream restores the missing subpackages
  this function becomes a no-op (we detect that and skip).
  """
  global _SHIMMED
  if _SHIMMED:
    return
  import importlib
  import importlib.util
  import sys
  import types

  # If the real modules are present (post-upstream-fix), do nothing.
  try:
    importlib.import_module("vdc.inference.density")
    importlib.import_module("vdc.vine.copula_diffusion")
    _SHIMMED = True
    return
  except ModuleNotFoundError:
    pass

  vdc_spec = importlib.util.find_spec("vdc")
  if vdc_spec is None or vdc_spec.origin is None:
    _SHIMMED = True
    return
  vdc_dir = vdc_spec.origin.rsplit("/", 1)[0]

  # Load ``scatter_to_hist`` from ``vdc/data/hist.py`` directly via the
  # file-location loader — going through ``import vdc.data.hist`` would
  # trigger ``vdc/data/__init__.py``'s transitive broken imports.
  hist_spec = importlib.util.spec_from_file_location(
    "_vdc_hist_internal", f"{vdc_dir}/data/hist.py"
  )
  if hist_spec is None or hist_spec.loader is None:
    _SHIMMED = True
    return
  hist_mod = importlib.util.module_from_spec(hist_spec)
  hist_spec.loader.exec_module(hist_mod)

  def _missing_sample_density_grid(*args, **kwargs):
    raise NotImplementedError(
      "vdc.inference.density.sample_density_grid is missing from the "
      "installed vdc wheel (upstream packaging bug); diffusion-checkpoint "
      "inference is unavailable. Use a denoiser checkpoint such as "
      "'vdc-denoiser-m64-v1' instead."
    )

  # Use setattr for the dynamic module-attribute writes so the type
  # checker doesn't flag ModuleType as missing the attributes we're
  # creating on the fly.
  inference_pkg = types.ModuleType("vdc.inference")
  setattr(inference_pkg, "__path__", [])
  density_mod = types.ModuleType("vdc.inference.density")
  setattr(density_mod, "scatter_to_hist", getattr(hist_mod, "scatter_to_hist"))
  setattr(density_mod, "sample_density_grid", _missing_sample_density_grid)
  sys.modules["vdc.inference"] = inference_pkg
  sys.modules["vdc.inference.density"] = density_mod

  vine_pkg = types.ModuleType("vdc.vine")
  setattr(vine_pkg, "__path__", [])
  cdiff_mod = types.ModuleType("vdc.vine.copula_diffusion")

  class _StubDiffusionCopulaModel:
    """Stub for an upstream-missing class that's import-referenced but
    never instantiated during inference."""

    def __init__(self, *args, **kwargs) -> None:
      pass

  setattr(cdiff_mod, "DiffusionCopulaModel", _StubDiffusionCopulaModel)
  sys.modules["vdc.vine"] = vine_pkg
  sys.modules["vdc.vine.copula_diffusion"] = cdiff_mod

  _SHIMMED = True


def _load_bundle(model_id: str, device: Optional[torch.device | str]):
  """Return a cached ``vdc.LoadedPretrainedModel`` for ``model_id`` on ``device``.

  The cache is keyed on ``(model_id, str(device))`` so the same checkpoint
  loaded on cpu and cuda gives two distinct entries. The HuggingFace
  download is only paid on the first cache miss per process.
  """
  vdc = _try_import_vdc()
  _install_upstream_shims()
  key = (model_id, str(device))
  if key not in _BUNDLE_CACHE:
    _BUNDLE_CACHE[key] = vdc.load_pretrained_model(
      model_id, device=device or "cpu"
    )
  return _BUNDLE_CACHE[key]


def _replicate_pad(density: Tensor) -> tuple[Tensor, Tensor]:
  """Pad a cell-centered ``(m, m)`` density to a ``(m+2, m+2)`` grid that
  spans ``[0, 1]``.

  Returns ``(grid_points, padded_values)`` ready for the
  :class:`TorchBicop` constructor. The padded grid points are
  ``[0, 0.5/m, 1.5/m, …, (m-0.5)/m, 1]`` (length ``m+2``); padded values
  replicate the first/last row and column. With this padding the
  trapezoidal integral on the padded grid equals vdc's midpoint cumsum
  at every cell center — see this module's docstring for the algebra.
  """
  if density.ndim != 2 or density.shape[0] != density.shape[1]:
    raise ValueError(
      f"vdc density must be a square 2D array; got shape {tuple(density.shape)}"
    )
  m = density.shape[0]
  cc = (torch.arange(m, dtype=density.dtype, device=density.device) + 0.5) / m
  grid_points = torch.cat([cc.new_zeros(1), cc, cc.new_ones(1)])
  padded = (
    torch.nn.functional.pad(
      density.unsqueeze(0).unsqueeze(0), (1, 1, 1, 1), mode="replicate"
    )
    .squeeze(0)
    .squeeze(0)
  )
  return grid_points, padded


def fit_vdc(
  u: Tensor,
  controls: "FitControlsTorchBicop",
  *,
  device: Optional[torch.device | str],
  dtype: torch.dtype,
) -> tuple[Tensor, Tensor]:
  """Estimate a bivariate density via vdc and return ``(grid_points, values)``.

  ``u`` is an ``(n, 2)`` tensor of pseudo-observations in ``[0, 1]^2``.
  The returned ``(grid_points, values)`` are sized ``(m+2,)`` and
  ``(m+2, m+2)`` respectively, where ``m`` is the pretrained model's
  density-grid resolution (``64`` for ``vdc-denoiser-m64-v1``).
  """
  if u.ndim != 2 or u.shape[1] != 2:
    raise ValueError(f"u must have shape (n, 2); got {tuple(u.shape)}")
  vdc = _try_import_vdc()
  bundle = controls.vdc_bundle
  if bundle is None:
    bundle = _load_bundle(controls.vdc_model_id, device)
  u_np = u.detach().cpu().numpy()
  density_np = vdc.estimate_pair_density_from_samples(
    bundle,
    u_np,
    diffusion_steps=controls.vdc_diffusion_steps,
    cfg_scale=controls.vdc_cfg_scale,
    projection_iters=controls.vdc_projection_iters,
  )
  density = torch.as_tensor(density_np, dtype=dtype, device=device)
  return _replicate_pad(density)
