"""Backends for the sklearn vine-copula estimators.

Each backend is a configured adapter that knows how to fit a vine on
pseudo-observations and how to evaluate ``pdf`` / ``cdf`` / ``sample``
on that vine. Two concrete backends ship:

- :class:`VinecopBackend` (default) — wraps ``Vinecop``.
  Holds an optional ``FitControlsVinecop`` and an
  optional ``RVineStructure``. **No PyTorch
  dependency**; this is what the sklearn module uses out of the box.
- :class:`TorchVinecopBackend` — wraps
  :class:`pyvinecopulib.torch.TorchVinecop`. Holds an optional
  :class:`pyvinecopulib.torch.FitControlsTorchVinecop` and an
  optional ``RVineStructure``. Constructing this
  class triggers the torch import (it's the explicit opt-in signal).

Notes
-----
**Which backend should I pick?**

Stay on the default (:class:`VinecopBackend`) when you want the
fastest CPU-bound vine fits, multi-threaded evaluation
(``controls.num_threads``), the full parametric pair-copula family
set (Gaussian, Student, Clayton, Gumbel, Frank, Joe, BB families, …)
**and** the non-parametric *Transformed Local Likelihood* (TLL)
family — all backed by ``Vinecop``.

Switch to :class:`TorchVinecopBackend` when you need any of:

- *GPU placement* — drop a fitted vine on the GPU with
  ``.to("cuda")`` and evaluate batched ``pdf`` / ``cdf`` /
  ``sample`` calls there.
- *Autograd* — the entire cascade is built from differentiable
  PyTorch ops, so gradients flow back through ``pdf`` / Rosenblatt
  outputs to any upstream parameters (e.g. learned marginals or
  feature transforms).
- *Composition with PyTorch pipelines* — the vine is a
  ``torch.nn.Module`` that drops into any other model.

The torch backend currently supports the TLL family only (which is
both the default and what the GPU path is built around); other
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

    VineDensity()                                          # default
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
from typing import Any, Optional

import numpy as np

import pyvinecopulib as pv


def _default_cpp_controls() -> pv.FitControlsVinecop:
  return pv.FitControlsVinecop(
    family_set=[pv.families.tll], trunc_lvl=20, num_threads=1
  )


class _VinecopBackendBase:
  """Shared adapter surface for the vine-copula backends.

  Holds the vine-fit configuration and the fitted vine's evaluation surface,
  factoring the parts identical across ``VinecopBackend`` and
  ``TorchVinecopBackend`` (``structure_of`` and the copy-on-write
  ``with_*`` derivations) here. Both controls types expose the same
  tree-selection fields (``tree_algorithm`` / ``seeds`` / ``trunc_lvl``),
  so the ``with_*`` helpers act on ``_effective_controls()`` uniformly.
  Concrete backends override only the genuinely divergent members: which vine
  class ``fit_vine`` builds, the per-op evaluation kwargs / output conversion,
  and the default controls.

  Parameters
  ----------
  controls : object, or None, optional
      Fit-time controls for the concrete backend (a ``FitControlsVinecop`` or a
      ``FitControlsTorchVinecop``). `None` resolves to the backend's default at
      fit time.
  structure : RVineStructure, or None, optional
      A pre-specified vine structure; when provided, `fit_vine` skips structure
      selection.
  """

  def __init__(
    self,
    *,
    controls: Optional[Any] = None,
    structure: Optional[pv.RVineStructure] = None,
  ) -> None:
    self.controls = controls
    self.structure = structure

  # -- hooks a concrete backend provides ---------------------------------- #
  def _default_controls(self) -> Any:
    raise NotImplementedError

  def _effective_controls(self) -> Any:
    return (
      self.controls if self.controls is not None else self._default_controls()
    )

  def fit_vine(self, U: np.ndarray, *, var_types: list[str]) -> Any:
    raise NotImplementedError

  def pdf(self, vine: Any, U: np.ndarray) -> np.ndarray:
    raise NotImplementedError

  def cdf(
    self, vine: Any, U: np.ndarray, *, N: int, seeds: list[int]
  ) -> np.ndarray:
    raise NotImplementedError

  def sample(
    self, vine: Any, n_samples: int, *, seeds: list[int]
  ) -> np.ndarray:
    raise NotImplementedError

  # -- shared surface (single source of truth) ---------------------------- #
  def default_margin(
    self, var_type: str, bounds: Optional[tuple[float, float]]
  ) -> Any:
    """The margin an estimator should fit when the caller named none.

    A seam rather than an `isinstance` check, so a backend whose vine lives on
    another array namespace can supply a matching margin.

    Note that ``TorchVinecopBackend`` deliberately does **not** override this.
    The sklearn layer wraps its vine in ``_BackendVinecop``, which converts to
    NumPy at the backend boundary because the estimators' public methods return
    arrays; a torch margin here would put the two halves of one ``Vinedist`` on
    different namespaces. End-to-end torch is reached through
    ``pyvinecopulib.torch.TorchVinedist.from_data`` instead.

    Parameters
    ----------
    var_type : str
        ``Kde1d``'s spelling of the variable type.
    bounds : tuple of float, or None
        Declared support, or ``None`` where the input states none.

    Returns
    -------
    MarginLike
        An unfitted kernel-density margin.
    """
    lo, hi = (None, None) if bounds is None else (bounds[0], bounds[1])
    return pv.core.Kde1d(type=var_type, xmin=lo, xmax=hi)

  def structure_of(self, vine: Any) -> pv.RVineStructure:
    return vine.structure

  def with_random_structure(
    self, d: int, seeds: list[int]
  ) -> "_VinecopBackendBase":
    new = _copy.copy(self)
    new.structure = pv.RVineStructure.sample(d, seeds=seeds)
    return new

  def with_local_random(self, seeds: list[int]) -> "_VinecopBackendBase":
    # Both controls types carry ``tree_algorithm`` / ``seeds``; set them on a
    # copy of the effective controls (copy-on-write) and clear the structure so
    # ``fit_vine`` selects a fresh Kendall-tau-weighted random tree.
    new_ctrls = _copy.copy(self._effective_controls())
    new_ctrls.tree_algorithm = "random_weighted"
    new_ctrls.seeds = seeds
    new = _copy.copy(self)
    new.controls = new_ctrls
    new.structure = None
    return new

  def with_num_threads(self, num_threads: int) -> "_VinecopBackendBase":
    new_controls = _copy.copy(self._effective_controls())
    new_controls.num_threads = num_threads
    new = _copy.copy(self)
    new.controls = new_controls
    return new


class VinecopBackend(_VinecopBackendBase):
  """Default backend. Wraps ``Vinecop``.

  Stores constructor arguments verbatim per the scikit-learn developer guide;
  the default ``FitControlsVinecop`` is materialized lazily in `fit_vine` when
  ``controls is None``.

  Parameters
  ----------
  controls : FitControlsVinecop, or None, optional
      Pair-family / threading / structure-selection knobs. `None` defaults to
      TLL with `trunc_lvl=20`.
  structure : RVineStructure, or None, optional
      A pre-specified vine structure; when provided, `fit_vine` skips structure
      selection.
  """

  def __init__(
    self,
    *,
    controls: Optional[pv.FitControlsVinecop] = None,
    structure: Optional[pv.RVineStructure] = None,
  ) -> None:
    super().__init__(controls=controls, structure=structure)

  def _default_controls(self) -> pv.FitControlsVinecop:
    return _default_cpp_controls()

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

  def sample(
    self, vine: Any, n_samples: int, *, seeds: list[int]
  ) -> np.ndarray:
    return np.asarray(
      vine.sample(
        n_samples,
        num_threads=self._effective_controls().num_threads,
        seeds=seeds,
      )
    )


class TorchVinecopBackend(_VinecopBackendBase):
  """PyTorch backend. Wraps :class:`pyvinecopulib.torch.TorchVinecop`.

  Pick this backend for GPU placement (``.to("cuda")``), autograd through the
  vine cascade, or composition with other ``torch.nn.Module`` code. The default
  :class:`VinecopBackend` is generally faster on CPU for the same problem.

  Constructing this class imports torch — the explicit opt-in signal that
  PyTorch is required. The default ``FitControlsTorchVinecop`` is materialized
  lazily at fit time.

  Parameters
  ----------
  controls : FitControlsTorchVinecop, or None, optional
      Cascade / placement / precision knobs plus the native (pure-torch)
      structure-selection knobs (``tree_algorithm`` / ``seeds`` / ``trunc_lvl``
      / ``tree_criterion`` / ``threshold``). `None` resolves to defaults at fit
      time.
  structure : RVineStructure, or None, optional
      A pre-specified vine structure; when provided, `fit_vine` skips structure
      selection.

  Notes
  -----
  `with_num_threads` is a no-op on this backend; for CPU intraop parallelism
  call ``torch.set_num_threads(N)`` globally before evaluating. Any outer
  parallelism a caller wraps around the estimator still applies
  independently.
  """

  def __init__(
    self,
    *,
    controls: Optional[Any] = None,
    structure: Optional[pv.RVineStructure] = None,
  ) -> None:
    from pyvinecopulib.torch import TorchVinecop  # noqa: F401  (torch opt-in)

    super().__init__(controls=controls, structure=structure)

  def _default_controls(self) -> Any:
    from pyvinecopulib.torch import FitControlsTorchVinecop

    return FitControlsTorchVinecop()

  def fit_vine(self, U: np.ndarray, *, var_types: list[str]) -> Any:
    from pyvinecopulib.torch import TorchVinecop

    return TorchVinecop.from_data(
      U,
      self.structure,
      controls=self._effective_controls(),
      var_types=var_types,
    )

  def pdf(self, vine: Any, U: np.ndarray) -> np.ndarray:
    out = vine.pdf(U, batched=self._effective_controls().batched)
    return out.detach().cpu().numpy()

  def cdf(
    self, vine: Any, U: np.ndarray, *, N: int, seeds: list[int]
  ) -> np.ndarray:
    out = vine.cdf(U, N=N, qrng=True, seeds=seeds)
    return out.detach().cpu().numpy()

  def sample(
    self, vine: Any, n_samples: int, *, seeds: list[int]
  ) -> np.ndarray:
    out = vine.sample(n_samples, qrng=False, seeds=seeds)
    return out.detach().cpu().numpy()

  def with_num_threads(self, num_threads: int) -> "TorchVinecopBackend":
    # No-op: torch threading is global / device-bound.
    return self


class _BackendVinecop:
  """A fitted vine that evaluates through the backend that fitted it.

  The estimators hand this to :class:`~pyvinecopulib.core.Vinedist` in place of
  the raw vine, so the distribution object evaluates exactly as the estimator
  does — with the backend's threading / batching arguments, and with results
  brought back to NumPy from wherever the vine computed them. Anything the
  backend does not adapt is read off the vine itself, so this still answers
  ``dim`` / ``var_types`` / ``structure`` and the Rosenblatt pair.

  Parameters
  ----------
  backend : object
      A resolved backend.
  vine : object
      The vine that backend fitted.
  """

  def __init__(self, backend: Any, vine: Any) -> None:
    self.backend = backend
    self.vine = vine

  def __getattr__(self, name: str) -> Any:
    # Underscored names are never forwarded: unpickling looks up `__setstate__`
    # before `vine` exists, and forwarding it would recurse.
    vine = self.__dict__.get("vine")
    if vine is None or name.startswith("_"):
      raise AttributeError(name)
    return getattr(vine, name)

  def pdf(self, u: np.ndarray) -> np.ndarray:
    """Copula density at ``u``.

    Parameters
    ----------
    u : ndarray, shape (n, d) or (n, d + k), dtype float
        Copula-scale data.

    Returns
    -------
    ndarray, shape (n,), dtype float
        Density values.
    """
    return self.backend.pdf(self.vine, u)

  def cdf(
    self,
    u: np.ndarray,
    *,
    N: int = 10000,
    seeds: Optional[list[int]] = None,
  ) -> np.ndarray:
    """Copula distribution function at ``u``.

    Parameters
    ----------
    u : ndarray, shape (n, d) or (n, d + k), dtype float
        Copula-scale data.
    N : int, optional
        Number of quasi-random points for the Monte-Carlo integration.
    seeds : list of int, or None, optional
        RNG seeds.

    Returns
    -------
    ndarray, shape (n,), dtype float
        Distribution values.
    """
    return self.backend.cdf(self.vine, u, N=N, seeds=list(seeds or []))

  def sample(self, n: int, *, seeds: Optional[list[int]] = None) -> np.ndarray:
    """Draw ``n`` samples on the copula scale.

    Parameters
    ----------
    n : int
        Number of samples.
    seeds : list of int, or None, optional
        RNG seeds.

    Returns
    -------
    ndarray, shape (n, d), dtype float
        Samples in ``[0, 1]^d``.
    """
    return self.backend.sample(self.vine, n, seeds=list(seeds or []))

  def __repr__(self) -> str:
    """Name the backend and the vine it wraps.

    Returns
    -------
    str
        The representation.
    """
    return (
      f"_BackendVinecop({type(self.backend).__name__}, "
      f"{type(self.vine).__name__})"
    )


def resolve_backend(backend: Any) -> Any:
  """Coerce a user-supplied ``backend=`` value to a concrete backend.

  Parameters
  ----------
  backend : object, or None
      `None` returns a default-constructed :class:`VinecopBackend`; any other
      value (a backend instance) is returned unchanged.

  Returns
  -------
  object
      A concrete backend instance.
  """
  return backend if backend is not None else VinecopBackend()


__all__ = [
  "VinecopBackend",
  "TorchVinecopBackend",
  "resolve_backend",
]
