"""Backends for the sklearn vine-copula estimators.

Each backend is a configured adapter that knows how to fit a vine on
pseudo-observations and how to evaluate ``pdf`` / ``cdf`` / ``simulate``
on that vine. Two concrete backends ship:

- :class:`VinecopBackend` (default) — wraps :class:`pyvinecopulib.Vinecop`.
  Holds an optional :class:`pyvinecopulib.FitControlsVinecop` and an
  optional :class:`pyvinecopulib.RVineStructure`. **No PyTorch
  dependency**; this is what the sklearn module uses out of the box.
- :class:`TorchVinecopBackend` — wraps
  :class:`pyvinecopulib.torch.TorchVinecop`. Holds an optional
  :class:`pyvinecopulib.torch.FitControlsTorchVinecop` and an
  optional :class:`pyvinecopulib.RVineStructure`. Constructing this
  class triggers the torch import (it's the explicit opt-in signal).

The :class:`VinecopLike` Protocol describes the shared post-fit
surface that both :class:`pyvinecopulib.Vinecop` and
:class:`pyvinecopulib.torch.TorchVinecop` satisfy structurally;
downstream code can type against it whenever it only needs runtime
evaluation, not fitting.

Notes
-----
**Which backend should I pick?**

Stay on the default (:class:`VinecopBackend`) when you want the
fastest CPU-bound vine fits, multi-threaded evaluation
(``controls.num_threads``), the full parametric pair-copula family
set (Gaussian, Student, Clayton, Gumbel, Frank, Joe, BB families, …)
**and** the non-parametric *Transformed Local Likelihood* (TLL)
family — all backed by the C++/nanobind core.

Switch to :class:`TorchVinecopBackend` when you need any of:

- *GPU placement* — drop a fitted vine on the GPU with
  ``.to("cuda")`` and evaluate batched ``pdf`` / ``cdf`` /
  ``simulate`` calls there.
- *Autograd* — the entire cascade is built from differentiable
  PyTorch ops, so gradients flow back through ``pdf`` / Rosenblatt
  outputs to any upstream parameters (e.g. learned marginals or
  feature transforms).
- *Composition with PyTorch pipelines* — the vine is a
  ``torch.nn.Module`` that drops into any other model.

The torch backend currently supports the TLL family only (which is
both the default and what the GPU path is built around) and the
optional ``"vdc"`` amortized estimator of Safaai (2026); other
parametric families require :class:`VinecopBackend`.

Examples
--------
Pass a configured backend instance to any sklearn estimator via
``backend=``::

    from pyvinecopulib.sklearn import VineDensity
    from pyvinecopulib.sklearn.backends import (
        VinecopBackend, TorchVinecopBackend,
    )
    from pyvinecopulib.torch import FitControlsTorchVinecop
    import pyvinecopulib as pv
    import torch

    VineDensity()                                          # default cpp
    VineDensity(backend=VinecopBackend(
        controls=pv.FitControlsVinecop(num_threads=4),
    ))
    VineDensity(backend=TorchVinecopBackend(
        controls=FitControlsTorchVinecop(
            device="cuda", dtype=torch.float32,
        ),
    ))
"""

from __future__ import annotations

import copy as _copy
from typing import Any, ClassVar, Optional, Protocol, runtime_checkable

import numpy as np

import pyvinecopulib as pv


@runtime_checkable
class VinecopLike(Protocol):
  """Post-fit duck-type for vine-copula evaluators.

  Satisfied structurally by both `pyvinecopulib.Vinecop` and
  `pyvinecopulib.torch.TorchVinecop`. Covers only the runtime
  methods sklearn estimators consume; backend strategy objects
  (`VinecopBackend`, `TorchVinecopBackend`) absorb signature
  differences for fit-time configuration.
  """

  structure: pv.RVineStructure

  def pdf(self, u, num_threads: int = 1, /, **kw) -> Any: ...
  def cdf(
    self,
    u,
    N: int = 10000,
    num_threads: int = 1,
    seeds: list[int] = [],
    /,
    **kw,
  ) -> Any: ...
  def simulate(
    self,
    n: int,
    qrng: bool = False,
    num_threads: int = 1,
    seeds: list[int] = [],
    /,
    **kw,
  ) -> Any: ...


def _default_cpp_controls() -> pv.FitControlsVinecop:
  return pv.FitControlsVinecop(
    family_set=[pv.families.tll], trunc_lvl=20, num_threads=1
  )


class VinecopBackend:
  """C++/nanobind backend (default). Wraps `pyvinecopulib.Vinecop`.

  Stores constructor arguments verbatim per the scikit-learn
  developer guide; the default `FitControlsVinecop` is materialised
  lazily in `fit_vine` when ``controls is None``.

  Parameters
  ----------
  controls :
      A `FitControlsVinecop` instance bundling pair-family / threading
      / structure-selection knobs. `None` defaults to TLL with
      `trunc_lvl=20`.
  structure :
      An optional pre-specified `RVineStructure`. When provided,
      `fit_vine` skips structure selection.
  """

  name: ClassVar[str] = "vinecop"
  supports_discrete: ClassVar[bool] = True
  supports_cdf: ClassVar[bool] = True
  supports_simulate: ClassVar[bool] = True

  def __init__(
    self,
    *,
    controls: Optional[pv.FitControlsVinecop] = None,
    structure: Optional[pv.RVineStructure] = None,
  ) -> None:
    self.controls = controls
    self.structure = structure

  # -- fit + evaluate ----------------------------------------------------- #

  def _effective_controls(self) -> pv.FitControlsVinecop:
    return (
      self.controls if self.controls is not None else _default_cpp_controls()
    )

  def fit_vine(self, U: np.ndarray, *, var_types: list[str]) -> Any:
    return pv.Vinecop.from_data(
      data=U,
      structure=self.structure,
      var_types=var_types,
      controls=self._effective_controls(),
    )

  def pdf(self, vine: Any, U: np.ndarray) -> np.ndarray:
    return np.asarray(
      vine.pdf(U, num_threads=self._effective_controls().num_threads)
    )

  def cdf(
    self, vine: Any, U: np.ndarray, *, N: int, seeds: list[int]
  ) -> np.ndarray:
    return np.asarray(
      vine.cdf(
        U,
        N=N,
        num_threads=self._effective_controls().num_threads,
        seeds=seeds,
      )
    )

  def simulate(
    self, vine: Any, n_samples: int, *, seeds: list[int]
  ) -> np.ndarray:
    return np.asarray(
      vine.simulate(
        n_samples,
        num_threads=self._effective_controls().num_threads,
        seeds=seeds,
      )
    )

  def structure_of(self, vine: Any) -> pv.RVineStructure:
    return vine.structure

  # -- forest plumbing (return new backend; copy-on-write) ---------------- #

  def with_random_structure(self, d: int, seeds: list[int]) -> "VinecopBackend":
    new = _copy.copy(self)
    new.structure = pv.RVineStructure.simulate(d, seeds=seeds)
    return new

  def with_local_random(self, seeds: list[int]) -> "VinecopBackend":
    new_controls = _copy.copy(self._effective_controls())
    new_controls.tree_algorithm = "random_weighted"
    new_controls.seeds = seeds
    new = _copy.copy(self)
    new.controls = new_controls
    new.structure = None
    return new

  def with_num_threads(self, num_threads: int) -> "VinecopBackend":
    new_controls = _copy.copy(self._effective_controls())
    new_controls.num_threads = num_threads
    new = _copy.copy(self)
    new.controls = new_controls
    return new


class TorchVinecopBackend:
  """PyTorch backend. Wraps `pyvinecopulib.torch.TorchVinecop`.

  Pick this backend for GPU placement (``.to("cuda")``), autograd
  through the vine cascade, or composition with other
  ``torch.nn.Module`` code. The default `VinecopBackend` is generally
  faster on CPU for the same problem.

  Constructing this class imports torch — the explicit opt-in signal
  that PyTorch is required. The default `FitControlsTorchVinecop` is
  materialised lazily in `fit_vine`.

  Parameters
  ----------
  controls :
      A `FitControlsTorchVinecop` instance bundling cascade /
      placement / precision knobs. `None` resolves to defaults at
      fit time.
  structure :
      An optional pre-specified `RVineStructure`. When provided,
      `fit_vine` skips structure selection.

  Notes
  -----
  `with_num_threads` is a no-op on this backend; for CPU intraop
  parallelism call ``torch.set_num_threads(N)`` globally before
  evaluating. The sklearn forest's outer parallelism (via joblib)
  still applies independently.
  """

  name: ClassVar[str] = "torch_vinecop"
  supports_discrete: ClassVar[bool] = False
  supports_cdf: ClassVar[bool] = True
  supports_simulate: ClassVar[bool] = True

  def __init__(
    self,
    *,
    controls: Optional[Any] = None,
    structure: Optional[pv.RVineStructure] = None,
  ) -> None:
    from pyvinecopulib.torch import TorchVinecop  # noqa: F401

    self.controls = controls
    self.structure = structure
    # When the user (or the forest) selects local random structure
    # sampling, structure selection still happens on the C++ side and
    # is lifted into torch via `TorchVinecop.from_vinecop`. This slot
    # holds the C++ controls used for that step; ``None`` means
    # ``_default_cpp_controls()`` (mst_prim / Dissmann).
    self._cpp_structure_controls: Optional[pv.FitControlsVinecop] = None

  # -- fit + evaluate ----------------------------------------------------- #

  def _effective_controls(self) -> Any:
    if self.controls is not None:
      return self.controls
    from pyvinecopulib.torch import FitControlsTorchVinecop

    return FitControlsTorchVinecop()

  def _effective_cpp_structure_controls(self) -> pv.FitControlsVinecop:
    return (
      self._cpp_structure_controls
      if self._cpp_structure_controls is not None
      else _default_cpp_controls()
    )

  def fit_vine(self, U: np.ndarray, *, var_types: list[str]) -> Any:
    if any(v != "c" for v in var_types):
      raise NotImplementedError(
        "TorchVinecopBackend is continuous-only; got var_types="
        f"{var_types!r}. Use VinecopBackend for discrete margins."
      )
    from pyvinecopulib.torch import TorchVinecop

    if self.structure is None:
      cpp_vine = pv.Vinecop.from_data(
        data=U,
        structure=None,
        var_types=var_types,
        controls=self._effective_cpp_structure_controls(),
      )
      ctrls = self._effective_controls()
      return TorchVinecop.from_vinecop(
        cpp_vine,
        cache_integrals=ctrls.cache_integrals,
        device=ctrls.device,
        dtype=ctrls.dtype,
      )
    return TorchVinecop.from_data(
      U, self.structure, controls=self._effective_controls()
    )

  def pdf(self, vine: Any, U: np.ndarray) -> np.ndarray:
    import torch

    ctrls = self._effective_controls()
    ref = vine._ref_tensor()
    u_t = torch.as_tensor(U, dtype=ref.dtype, device=ref.device)
    out = vine.pdf(u_t, impl=ctrls.impl, batched=ctrls.batched)
    return out.detach().cpu().numpy()

  def cdf(
    self, vine: Any, U: np.ndarray, *, N: int, seeds: list[int]
  ) -> np.ndarray:
    import torch

    ctrls = self._effective_controls()
    ref = vine._ref_tensor()
    u_t = torch.as_tensor(U, dtype=ref.dtype, device=ref.device)
    out = vine.cdf(u_t, N=N, qrng=True, seeds=seeds, impl=ctrls.impl)
    return out.detach().cpu().numpy()

  def simulate(
    self, vine: Any, n_samples: int, *, seeds: list[int]
  ) -> np.ndarray:
    ctrls = self._effective_controls()
    out = vine.simulate(n_samples, qrng=False, seeds=seeds, impl=ctrls.impl)
    return out.detach().cpu().numpy()

  def structure_of(self, vine: Any) -> pv.RVineStructure:
    return vine.structure

  # -- forest plumbing ---------------------------------------------------- #

  def with_random_structure(
    self, d: int, seeds: list[int]
  ) -> "TorchVinecopBackend":
    new = _copy.copy(self)
    new.structure = pv.RVineStructure.simulate(d, seeds=seeds)
    return new

  def with_local_random(self, seeds: list[int]) -> "TorchVinecopBackend":
    new_cpp = _copy.copy(self._effective_cpp_structure_controls())
    new_cpp.tree_algorithm = "random_weighted"
    new_cpp.seeds = seeds
    new = _copy.copy(self)
    new._cpp_structure_controls = new_cpp
    new.structure = None
    return new

  def with_num_threads(self, num_threads: int) -> "TorchVinecopBackend":
    # No-op: torch threading is global / device-bound.
    return self


def resolve_backend(backend: Any) -> Any:
  """Coerces a user-supplied ``backend=`` value to a concrete backend.

  Parameters
  ----------
  backend :
      `None` (returns a default-constructed `VinecopBackend`); a
      backend-like instance implementing ``fit_vine`` / ``pdf`` /
      ``cdf`` / ``simulate`` / ``structure_of`` (returned as-is);
      anything else raises `TypeError`.

  Returns
  -------
  A concrete backend instance.

  Raises
  ------
  TypeError
      If ``backend`` does not satisfy the duck-typed interface.
  """
  if backend is None:
    return VinecopBackend()
  for attr in (
    "fit_vine",
    "pdf",
    "cdf",
    "simulate",
    "structure_of",
    "with_random_structure",
    "with_local_random",
    "with_num_threads",
  ):
    if not callable(getattr(backend, attr, None)):
      raise TypeError(
        "backend must be a VinecopBackend, TorchVinecopBackend, or "
        f"compatible duck-typed object; got {type(backend).__name__} "
        f"missing method {attr!r}"
      )
  return backend


__all__ = [
  "VinecopLike",
  "VinecopBackend",
  "TorchVinecopBackend",
  "resolve_backend",
]
