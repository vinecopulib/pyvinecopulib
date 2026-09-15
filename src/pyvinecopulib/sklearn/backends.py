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
from typing import (
  TYPE_CHECKING,
  Any,
  Generic,
  Protocol,
  Optional,
  Self,
  Sequence,
  TypeVar,
  cast,
)

import numpy as np

import pyvinecopulib as pv

if TYPE_CHECKING:
  # Types only: the runtime import of `torch` stays inside
  # `TorchVinecopBackend`, where constructing the class is the opt-in signal.
  from pyvinecopulib.torch import (
    FitControlsTorchVinecop,
    TorchKde1d,
    TorchVinecop,
  )


def _default_cpp_controls() -> pv.FitControlsVinecop:
  return pv.FitControlsVinecop(
    family_set=[pv.families.tll], trunc_lvl=20, num_threads=1
  )


#: The vine a backend fits. Carried as a type parameter so a concrete backend
#: names its own class -- `Vinecop` or `TorchVinecop` -- in every signature
#: without narrowing an inherited parameter, which an override may not do.
#: A fitted vine, whichever lane fitted it. Bound to the public contract --
#: both lanes satisfy it -- while each concrete backend names its own vine
#: type, which is where the per-lane keywords (`num_threads`, batching) live.
_VineT = TypeVar("_VineT", bound=pv.core.VinecopLike[Any])


class _TreeControls(Protocol):
  """The tree-selection fields both lanes' controls carry.

  The bound on `_CtrlT`, and everything the shared `with_*` derivations
  read. Each concrete backend parameterizes the class with its own
  controls type, so a lane-specific field -- `num_threads` on the compiled
  lane, `device` / `dtype` on the torch one -- is read at the precise type
  rather than off an `Any`.
  """

  tree_algorithm: str
  seeds: list[int]

  def to_dict(self) -> dict[str, Any]: ...


#: The controls a backend holds. A type parameter for the same reason `_VineT`
#: is one: the two lanes share no field beyond those above, so the base can
#: only name what it reads and the subclass supplies the rest.
_CtrlT = TypeVar("_CtrlT", bound=_TreeControls)


class _VinecopBackendBase(Generic[_VineT, _CtrlT]):
  """Shared adapter surface for the vine-copula backends.

  Holds the vine-fit configuration and the fitted vine's evaluation surface,
  factoring the parts identical across ``VinecopBackend`` and
  ``TorchVinecopBackend`` (``structure_of`` and the copy-on-write
  ``with_*`` derivations) here. Both controls types expose the same
  tree-selection fields (``tree_algorithm`` / ``seeds`` / ``trunc_lvl``),
  so the ``with_*`` helpers act on ``_effective_controls()`` uniformly.
  Concrete backends override only the actually divergent members: which vine
  class ``fit_vine`` builds, the per-op evaluation kwargs / output conversion,
  and the default controls.

  Parameters
  ----------
  controls : ControlsLike, or None, optional
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
    controls: Optional[_CtrlT] = None,
    structure: Optional[pv.RVineStructure] = None,
  ) -> None:
    self.controls = controls
    self.structure = structure

  # -- hooks a concrete backend provides ---------------------------------- #
  def _default_controls(self) -> _CtrlT:
    raise NotImplementedError

  def _effective_controls(self) -> _CtrlT:
    return (
      self.controls if self.controls is not None else self._default_controls()
    )

  def fit_vine(self, U: np.ndarray, *, var_types: list[str]) -> _VineT:
    raise NotImplementedError

  def pdf(self, vine: _VineT, U: np.ndarray) -> np.ndarray:
    raise NotImplementedError

  def cdf(
    self, vine: _VineT, U: np.ndarray, *, N: int, seeds: list[int]
  ) -> np.ndarray:
    raise NotImplementedError

  def sample(
    self, vine: _VineT, n_samples: int, *, qrng: bool = False, seeds: list[int]
  ) -> np.ndarray:
    raise NotImplementedError

  # -- shared surface (single source of truth) ---------------------------- #
  # `MarginLike[Any]` rather than a namespace-specific one: a backend whose
  # vine lives on another array namespace supplies a margin on that namespace,
  # so no single parameterization describes the hook.
  def default_margin(
    self, var_type: str, bounds: Optional[tuple[float, float]]
  ) -> pv.core.MarginLike[Any]:
    """The margin an estimator should fit when the caller named none.

    A hook rather than an `isinstance` check, so a backend whose vine lives on
    another array namespace supplies a margin on that namespace and the two
    halves of one ``Vinedist`` stay together.

    Parameters
    ----------
    var_type : str
        ``Kde1d``'s spelling of the variable type.
    bounds : tuple of float, or None, optional
        Declared support, or ``None`` where the input states none.

    Returns
    -------
    MarginLike
        An unfitted kernel-density margin.
    """
    lo, hi = (None, None) if bounds is None else (bounds[0], bounds[1])
    return pv.core.Kde1d(type=var_type, xmin=lo, xmax=hi)

  def bind_distribution(
    self, vine: _VineT, margins: Sequence[pv.core.MarginLike[Any]]
  ) -> pv.core.VinedistLike[Any]:
    """Assemble the fitted vine and its margins into one distribution.

    The copula is wrapped in ``_BackendVinecop`` so the distribution evaluates
    with the backend's own threading and batching arguments, which is what keeps
    an estimator's ``distribution_`` numerically identical to its own methods.

    Parameters
    ----------
    vine : object
        The fitted vine.
    margins : sequence of MarginLike
        The fitted margins, in the vine's variable order.

    Returns
    -------
    Vinedist
        The joint distribution.
    """
    return pv.core.Vinedist(_BackendVinecop(self, vine), list(margins))

  def structure_of(self, vine: _VineT) -> pv.RVineStructure:
    # `VinecopLike[Any]` resolves its members to `Unknown`, so the read is
    # named rather than inferred.
    return cast("pv.RVineStructure", vine.structure)

  def with_random_structure(self, d: int, seeds: list[int]) -> Self:
    new = _copy.copy(self)
    new.structure = pv.RVineStructure.sample(d, seeds=seeds)
    return new

  def with_local_random(self, seeds: list[int]) -> Self:
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

  def with_fit_seeds(self, seeds: list[int]) -> Self:
    """Return a copy whose stochastic fit draws use ``seeds``.

    Parameters
    ----------
    seeds : list of int
        Seeds passed to the random-tree selector.

    Returns
    -------
    _VinecopBackendBase[Any, Any]
        Independent backend configuration for one fit, of this backend's own
        class.

    The estimator owns its ``random_state``. Copying the controls keeps a
    caller-owned backend reusable while making random-tree selection
    reproducible through the estimator's public parameter.
    """
    new = _copy.copy(self)
    controls = _copy.copy(self._effective_controls())
    controls.seeds = seeds
    new.controls = controls
    return new


class VinecopBackend(
  _VinecopBackendBase["pv.Vinecop", "pv.FitControlsVinecop"]
):
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

  #: This lane's controls, narrower than the base's ``ControlsLike``: the
  #: fields ``fit_vine`` and the evaluation methods read are this type's.
  controls: Optional[pv.FitControlsVinecop]

  def __init__(
    self,
    *,
    controls: Optional[pv.FitControlsVinecop] = None,
    structure: Optional[pv.RVineStructure] = None,
  ) -> None:
    super().__init__(controls=controls, structure=structure)

  def with_num_threads(self, num_threads: int) -> Self:
    """A copy of this backend fitting and evaluating on ``num_threads``.

    On this lane only: PyTorch threading is global and device-bound, so
    ``TorchVinecopBackend`` overrides it as a no-op.

    Parameters
    ----------
    num_threads : int
        Threads to use.

    Returns
    -------
    VinecopBackend
        A copy carrying the new thread count.
    """
    new_controls = _copy.copy(self._effective_controls())
    new_controls.num_threads = num_threads
    new = _copy.copy(self)
    new.controls = new_controls
    return new

  def _default_controls(self) -> pv.FitControlsVinecop:
    return _default_cpp_controls()

  def fit_vine(self, U: np.ndarray, *, var_types: list[str]) -> pv.Vinecop:
    return pv.Vinecop.from_data(
      data=U,
      structure=self.structure,
      var_types=var_types,
      controls=self._effective_controls(),
    )

  def pdf(self, vine: pv.Vinecop, U: np.ndarray) -> np.ndarray:
    return np.asarray(
      vine.pdf(U, num_threads=self._effective_controls().num_threads)
    )

  def cdf(
    self, vine: pv.Vinecop, U: np.ndarray, *, N: int, seeds: list[int]
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
    self,
    vine: pv.Vinecop,
    n_samples: int,
    *,
    qrng: bool = False,
    seeds: list[int],
  ) -> np.ndarray:
    return np.asarray(
      vine.sample(
        n_samples,
        qrng=qrng,
        num_threads=self._effective_controls().num_threads,
        seeds=seeds,
      )
    )


class TorchVinecopBackend(
  _VinecopBackendBase["TorchVinecop", "FitControlsTorchVinecop"]
):
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

  #: This lane's controls; see :class:`VinecopBackend`'s.
  controls: Optional["FitControlsTorchVinecop"]

  def __init__(
    self,
    *,
    controls: Optional["FitControlsTorchVinecop"] = None,
    structure: Optional[pv.RVineStructure] = None,
  ) -> None:
    from pyvinecopulib.torch import TorchVinecop  # noqa: F401  (torch opt-in)

    super().__init__(controls=controls, structure=structure)

  def _default_controls(self) -> "FitControlsTorchVinecop":
    from pyvinecopulib.torch import FitControlsTorchVinecop

    return FitControlsTorchVinecop()

  def fit_vine(self, U: np.ndarray, *, var_types: list[str]) -> "TorchVinecop":
    from pyvinecopulib.torch import TorchVinecop

    return TorchVinecop.from_data(
      U,
      self._effective_controls(),
      structure=self.structure,
      var_types=var_types,
    )

  # `vine` is a `TorchVinecop`, as it is on `sample` -- typed `Any` on the two
  # methods that hand it an array, because the cascade declares a `Tensor` and
  # what arrives is the estimator's NumPy array, which the vine's placement
  # hook brings across.
  def pdf(self, vine: "TorchVinecop", U: np.ndarray) -> np.ndarray:
    # No `batched=`: the vine resolves it per device, which is what every
    # other call on it does -- `cdf` and `sample` here, and `TorchVinedist`.
    out = vine.pdf(cast("Any", U))
    return out.detach().cpu().numpy()

  def cdf(
    self, vine: "TorchVinecop", U: np.ndarray, *, N: int, seeds: list[int]
  ) -> np.ndarray:
    out = vine.cdf(cast("Any", U), N=N, qrng=True, seeds=seeds)
    return out.detach().cpu().numpy()

  def sample(
    self,
    vine: "TorchVinecop",
    n_samples: int,
    *,
    qrng: bool = False,
    seeds: list[int],
  ) -> np.ndarray:
    out = vine.sample(n_samples, qrng=qrng, seeds=seeds)
    return out.detach().cpu().numpy()

  def default_margin(
    self, var_type: str, bounds: Optional[tuple[float, float]]
  ) -> "TorchKde1d":
    """A ``TorchKde1d`` placed and typed like the copula this backend fits.

    ``device`` and ``dtype`` come from the effective controls, so
    ``FitControlsTorchVinecop(dtype=torch.float32)`` does not leave float64
    margins on a float32 copula.

    Parameters
    ----------
    var_type : str
        ``Kde1d``'s spelling of the variable type.
    bounds : tuple of float, or None, optional
        Declared support, or ``None`` where the input states none.

    Returns
    -------
    TorchKde1d
        An unfitted torch kernel-density margin.
    """
    import torch

    from pyvinecopulib.torch import TorchKde1d

    controls = self._effective_controls()
    lo, hi = (None, None) if bounds is None else (bounds[0], bounds[1])
    return TorchKde1d(
      type=var_type,
      xmin=lo,
      xmax=hi,
      device=controls.device,
      dtype=controls.dtype if controls.dtype is not None else torch.float64,
    )

  def bind_distribution(
    self, vine: "TorchVinecop", margins: Sequence[pv.core.MarginLike[Any]]
  ) -> pv.core.VinedistLike[Any]:
    """Assemble a ``TorchVinedist`` from the fitted torch vine and its margins.

    The raw vine goes in rather than a ``_BackendVinecop`` wrapper: the point of
    publishing this object is that it is torch throughout, so ``.to(device)``
    moves it and a loss through ``logpdf`` reaches the margins' buffers -- which
    a wrapper converting to NumPy at every call would undo. A core ``Kde1d``
    margin (from ``margins="kde"``, or one the caller passed already fitted) is
    lifted with ``TorchKde1d.from_kde1d``, an exact transfer of the same
    grid.

    Parameters
    ----------
    vine : TorchVinecop
        The fitted vine.
    margins : sequence of MarginLike
        The fitted margins, in the vine's variable order.

    Returns
    -------
    TorchVinedist
        The joint distribution, entirely in torch.
    """
    import torch

    from pyvinecopulib.torch import TorchKde1d, TorchVinedist

    controls = self._effective_controls()
    dtype = controls.dtype if controls.dtype is not None else torch.float64
    lifted = [
      TorchKde1d.from_kde1d(m, device=controls.device, dtype=dtype)
      if isinstance(m, pv.core.Kde1d)
      else m
      for m in margins
    ]
    return TorchVinedist(vine, lifted)

  def with_num_threads(self, num_threads: int) -> Self:
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
  backend : _VinecopBackendBase[Any, Any]
      A resolved backend.
  vine : VinecopLike
      The vine that backend fitted.
  """

  def __init__(
    self, backend: _VinecopBackendBase[Any, Any], vine: pv.core.VinecopLike[Any]
  ) -> None:
    self.backend = backend
    self.vine = vine

  def __getattr__(self, name: str) -> Any:  # noqa: ANN401 - proxied attribute
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
    N : int, default=10000
        Number of quasi-random points for the Monte-Carlo integration.
    seeds : list of int, or None, optional
        RNG seeds.

    Returns
    -------
    ndarray, shape (n,), dtype float
        Distribution values.
    """
    return self.backend.cdf(self.vine, u, N=N, seeds=list(seeds or []))

  def sample(
    self,
    n: int,
    *,
    qrng: bool = False,
    seeds: Optional[list[int]] = None,
  ) -> np.ndarray:
    """Draw ``n`` samples on the copula scale.

    Parameters
    ----------
    n : int
        Number of samples.
    qrng : bool, default=False
        Draw quasi-random base uniforms instead of pseudo-random ones.
    seeds : list of int, or None, optional
        RNG seeds.

    Returns
    -------
    ndarray, shape (n, d), dtype float
        Samples in ``[0, 1]^d``.
    """
    return self.backend.sample(self.vine, n, qrng=qrng, seeds=list(seeds or []))

  # The three members below are read off the vine unchanged -- the backend has
  # nothing to adapt in them. They are written out rather than left to
  # `__getattr__` so the class satisfies `VinecopLike` statically, which is
  # what `Vinedist` asks for.
  @property
  def structure(self) -> pv.RVineStructure:
    """The R-vine structure the vine was fitted on.

    Returns
    -------
    RVineStructure
        The structure.
    """
    return cast("pv.RVineStructure", self.vine.structure)

  def rosenblatt(self, u: np.ndarray) -> np.ndarray:
    """Rosenblatt transform, as the vine computes it.

    Parameters
    ----------
    u : ndarray, shape (n, d), dtype float
        Copula-scale observations.

    Returns
    -------
    ndarray, shape (n, d), dtype float
        Independent uniforms.
    """
    return np.asarray(self.vine.rosenblatt(u))

  def inverse_rosenblatt(self, u: np.ndarray) -> np.ndarray:
    """Inverse Rosenblatt transform, as the vine computes it.

    Parameters
    ----------
    u : ndarray, shape (n, d), dtype float
        Independent uniforms.

    Returns
    -------
    ndarray, shape (n, d), dtype float
        Copula-scale observations.
    """
    return np.asarray(self.vine.inverse_rosenblatt(u))

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


# The return is the backend itself, and typed `Any` because it is what pins
# `VineBase.backend_`: `bind_distribution` can only promise a `VinedistLike`,
# which is narrower than the `Vinedist` the estimators publish as
# `distribution_`, so naming the type here would narrow every read of it.
def resolve_backend(
  backend: Optional[_VinecopBackendBase[Any, Any]],
) -> Any:  # noqa: ANN401 - see above
  """Coerce a user-supplied ``backend=`` value to a concrete backend.

  Parameters
  ----------
  backend : _VinecopBackendBase[Any, Any], or None, optional
      `None` returns a default-constructed :class:`VinecopBackend`; any other
      value (a backend instance) is returned unchanged.

  Returns
  -------
  _VinecopBackendBase[Any, Any]
      A concrete backend instance.
  """
  return backend if backend is not None else VinecopBackend()


__all__ = [
  "VinecopBackend",
  "TorchVinecopBackend",
  "resolve_backend",
]
