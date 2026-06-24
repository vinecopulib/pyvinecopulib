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

METHODS: tuple[str, ...] = ("tll",)


@dataclass
class FitControlsTorchBicop:
  """Controls for `TorchBicop.from_data`.

  Mirrors `pyvinecopulib.FitControlsBicop`: ``method`` picks the
  pair-copula fitter and the remaining fields carry method-specific
  hyperparameters.

  Attributes
  ----------
  method : {"tll"}, default="tll"
      Fitter to use. ``"tll"`` is a pure-torch *Transformed Local
      Likelihood* kernel density estimator on a 2-D grid in the
      inverse-normal-transformed copula space (Geenens, 2014;
      Nagler, 2018), matching the C++ TLL fit to machine precision.
      It is the only fitter currently shipped; ``method`` is kept as
      the extension seam for future torch fitters.
  grid_size : int, default=30
      *TLL only.* Density grid size per axis.
  mult : float, default=1.0
      *TLL only.* Bandwidth multiplier.
  grid_type : {"normal", "linear"}, default="normal"
      *TLL only.* Storage grid type — ``"normal"`` (Phi-spaced,
      the C++-parity default) or ``"linear"`` (uniform on
      ``[0, 1]`` with the O(1) cell-finding fast-path).
  """

  method: str = "tll"

  # TLL-only
  grid_size: int = 30
  mult: float = 1.0
  grid_type: str = "normal"

  def __post_init__(self) -> None:
    if self.method not in METHODS:
      raise ValueError(
        f"unknown method={self.method!r}; expected one of {METHODS}"
      )


@dataclass
class FitControlsTorchVinecop:
  """Controls for `TorchVinecop.from_data` and the runtime cascade.

  Mirrors `pyvinecopulib.FitControlsVinecop`: bundles all vine-fit
  knobs into one object. A nested `FitControlsTorchBicop` controls
  how each pair-copula is fit; vine-level fields below carry
  placement / precision / cascade-variant settings.

  Attributes
  ----------
  bicop_controls : FitControlsTorchBicop
      Controls applied to every pair-copula fit.
  cache_integrals : bool, default=True
      If ``True``, precompute the cdf / hfunc / hinv caches on
      every pair copula's interpolation grid. Cached lookups are
      1–2 orders of magnitude faster than the on-the-fly path
      with a ~1e-3 IAE cost.
  device : torch.device or None, default=None
      Target torch device for the fitted pair copulas. `None`
      keeps the input's device.
  dtype : torch.dtype or None, default=None
      Target torch dtype. `None` defaults to ``torch.float64``
      (parity with C++).
  batched : bool, default=False
      If ``True``, fires a single batched bicop call per tree
      level. Available for ``pdf`` / ``rosenblatt`` only
      (``inverse_rosenblatt(batched=True)`` raises).
  """

  bicop_controls: FitControlsTorchBicop = field(
    default_factory=FitControlsTorchBicop
  )
  cache_integrals: bool = True
  device: Optional[Any] = None
  dtype: Optional[Any] = None
  batched: bool = False
