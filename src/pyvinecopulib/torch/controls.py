"""Fit-time controls for :class:`TorchBicop` and :class:`TorchVinecop`.

Mirrors ``FitControlsBicop`` /
``FitControlsVinecop``: method-specific args live on
the dataclass; cross-cutting args stay on the relevant ``from_data``
signature only where they don't fit naturally on the controls.

Adding a new fitter to ``TorchBicop`` only requires extending the
relevant dataclass and the dispatch in the corresponding ``from_data``
— the public ``from_data`` signatures are forward-stable.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields
from typing import Any, Optional

METHODS: tuple[str, ...] = ("tll",)

#: Structure-selection algorithms accepted by ``FitControlsTorchVinecop``,
#: mirroring ``FitControlsVinecop.tree_algorithm``.
TREE_ALGORITHMS: tuple[str, ...] = (
  "mst_prim",
  "mst_kruskal",
  "random_weighted",
  "random_unweighted",
)


@dataclass
class FitControlsTorchBicop:
  """Controls for ``TorchBicop.from_data()``.

  Mirrors ``FitControlsBicop``: ``method`` picks the
  pair-copula fitter and the remaining fields carry method-specific
  hyperparameters.

  Attributes
  ----------
  method : {"tll"}, default="tll"
      Fitter to use. ``"tll"`` is a pure-torch *Transformed Local
      Likelihood* kernel density estimator on a 2-D grid in the
      inverse-normal-transformed copula space (Geenens, 2014;
      Nagler, 2018), matching the ``Bicop`` TLL fit to
      machine precision. It is the only fitter currently shipped;
      ``method`` is kept as the extension seam for future torch fitters.
  grid_size : int, default=30
      *TLL only.* Density grid size per axis.
  mult : float, default=1.0
      *TLL only.* Bandwidth multiplier.
  grid_type : {"normal", "linear"}, default="normal"
      *TLL only.* Storage grid type — ``"normal"`` (Phi-spaced,
      the ``Bicop``-parity default) or ``"linear"``
      (uniform on ``[0, 1]`` with the O(1) cell-finding fast-path).
  compile_fit : bool, default=False
      *TLL only.* Fuse the bandwidth search's per-pass body with
      ``torch.compile``. The pass is 39 kernel launches over tensors small
      enough that the arithmetic is invisible even at ``n = 30000``, so
      fusing it roughly halves a pair fit: measured 2.0x end to end at
      ``n = 500``, 1.6x at 5000, 1.4x at 10000, tapering to 1.04x at 30000
      and never below 1.

      Off by default because the first call compiles for about 5 seconds and
      a second lane count costs another 8, which only a caller fitting many
      pairs in one process earns back -- around 500 pair fits at
      ``n = 2000``, so roughly 60 vines at ``d = 9``. Fitting one vine, it is
      a large loss. The same trade as
      ``FitControlsTorchVinecop.compile``, which fuses the evaluation
      cascades.

      Fusion reorders the arithmetic, so the fitted grid moves by about
      1e-15 and a lane's iteration count can change -- the outer criterion
      sits at the float64 noise floor. The ``Bicop`` parity gate is
      unaffected.
  """

  method: str = "tll"

  # TLL-only
  grid_size: int = 30
  mult: float = 1.0
  grid_type: str = "normal"
  compile_fit: bool = False

  def to_dict(self) -> dict:
    """Return the settings as a plain dictionary.

    Returns
    -------
    dict
        One entry per setting, keyed by the field name. Nested controls are
        kept as objects rather than flattened, so the entry round-trips.

    See Also
    --------
    pyvinecopulib.core.ControlsLike : The contract this satisfies.
    """
    return {f.name: getattr(self, f.name) for f in fields(self)}

  def __post_init__(self) -> None:
    if self.method not in METHODS:
      raise ValueError(
        f"unknown method={self.method!r}; expected one of {METHODS}"
      )


@dataclass
class FitControlsTorchVinecop(FitControlsTorchBicop):
  """Controls for :meth:`~pyvinecopulib.torch.TorchVinecop.from_data` and the cascade.

  Mirrors :class:`~pyvinecopulib.core.FitControlsVinecop`: bundles all vine-fit
  knobs into one object. A nested ``FitControlsTorchBicop`` controls
  how each pair-copula is fit; the vine-level fields below carry the
  structure-selection knobs plus placement / precision /
  cascade-variant settings.

  Attributes
  ----------
  trunc_lvl : int, default=20
      Maximum number of trees to select when
      :meth:`~pyvinecopulib.torch.TorchVinecop.from_data` is called with
      ``structure=None``.
  tree_criterion : {"tau", "rho", "hoeffd", "mcor", "cxi", "joe"}, \
default="tau"
      Dependence measure used to weight candidate edges during structure
      selection (passed to ``wdm``).
  threshold : float, default=0.0
      Dependence threshold, acting twice as it does in ``Vinecop``: a
      candidate edge below it is deprioritized during spanning-tree
      selection, and a surviving edge below it is left independent rather
      than fitted. It applies along a supplied structure too, where only the
      second half has anything to do. At the default nothing is below it.
  tree_algorithm : {"mst_prim", "mst_kruskal", "random_weighted", \
"random_unweighted"}, default="mst_prim"
      Spanning-tree algorithm for structure selection: Dissmann's maximum
      spanning tree (``mst_*``) or Wilson's (weighted / uniform) random
      spanning tree (``random_*``).
  seeds : list of int, default=[]
      RNG seeds for the random tree algorithms (ignored by the MST ones).
  conditioning_set : list of int, default=[]
      1-based variables to place at the tail of the selected order, so they can
      be conditioned on with
      :meth:`~pyvinecopulib.core.VinecopBase.sample_conditional`. Requires an MST
      ``tree_algorithm`` and no truncation; empty means unconditional selection.
  cache_integrals : bool or None, default=None
      If ``True``, precompute prefix tables for CDF and h-function evaluation.
      They reconstruct the uncached integrals up to summation order and rebuild
      in-graph when grid values require gradients. ``None`` resolves to
      ``True``, including for discrete vines.
  device : torch.device or None, default=None
      Target torch device for the fitted pair copulas. ``None``
      keeps the input's device.
  dtype : torch.dtype or None, default=None
      Target torch dtype. ``None`` defaults to ``torch.float64``
      (parity with :class:`~pyvinecopulib.core.Vinecop`).
  compile : bool, default=False
      If ``True``, wrap the batched cascades in
      :func:`torch.compile` with ``dynamic=False``. Inductor fuses the
      cascade's elementwise chains into a handful of kernels, which is worth
      a great deal on CUDA -- where the eager cascade is bound by
      kernel-launch count rather than arithmetic -- and costs a
      one-off compilation of tens of seconds on the first call at each
      input shape. Off by default for that reason: a single evaluation is
      slower compiled than not. Results agree with the eager path to
      floating point, not exactly.
  batched_fit : bool or None, default=None
      Fit a whole tree level in one call rather than one edge at a time.
      ``None`` resolves per device -- on for CUDA, off otherwise -- the same
      shape as the evaluation cascade's ``batched``. The per-level fitter
      advances every edge's bandwidth search together and freezes each lane
      as it converges, which trades a larger working set for far fewer
      kernel launches.

      That trade pays where launches cost something. On CUDA a whole vine
      fit is 1.3-4.0x faster over ``d`` in 5..20 and ``n`` in 2000..12000,
      on both the fixed-structure and the selecting path, and never slower
      than fitting edge at a time. It wins least at large ``n``, where the
      fit is arithmetic-bound and there was little overhead to remove.

      The working set does not grow with the level, so there is no size at
      which the trade inverts: the kernel evaluation blocks its grid axis
      against a memory budget, and the block shrinks as ``n`` and the level
      width rise. Peak stays near 280 MiB whether the vine is ``d = 9`` at
      ``n = 12000`` or ``d = 25`` at ``n = 20000``.

      Figures come from interleaving the two arms and taking medians: on a
      thermally throttling laptop a lone measurement of these cells moves by
      half again, which is enough to invent a reversal or hide one.

      On cpu it is 0.44-1.05x at one torch thread -- mostly slower, there
      being no launch overhead to amortize. Many threads favor it again,
      the larger kernels parallelizing better than many small ones, and
      torch uses every core by default, so the cpu default is the
      conservative reading of a measurement that moves with thread count
      rather than a claim that batching cannot pay there. Either way the
      torch fit is far from competitive with the core backend on cpu.

      A level carrying a discrete edge or a conditioning context is always
      fitted edge at a time: those cannot stack.

      The batched result agrees with the per-edge one to floating point on
      every device -- not bit for bit, on cpu included. A lane's iterations
      are independent of its batch-mates, each freezing as it converges,
      but the arithmetic passes through kernels torch selects by element
      count: the bandwidth search's `pow` takes a vectorized path past
      ``2 * Vectorized<double>::size()`` lanes, which is 8 on AVX2 and 4 on
      NEON. So the last bits depend on how many pairs travelled together,
      and a machine that vectorizes sooner diverges where another does not.
      Selected structures still match, the tree criterion being a function
      of ranks rather than of those bits.

  Notes
  -----
  Structure selection runs natively on the torch interpolation grids. It is
  TLL-only, and the criteria for automatic truncation / thresholding (``aic`` /
  ``bic`` / ``mbicv``) are not available here: ``trunc_lvl`` is a fixed cap.

  This *is* a :class:`FitControlsTorchBicop`, so the settings governing each
  pair-copula fit are its own attributes rather than a nested object, and one
  controls instance configures both halves of a vine fit. The same holds for
  :class:`~pyvinecopulib.core.FitControlsVinecop` and
  :class:`~pyvinecopulib.core.FitControlsBicop`.
  """

  trunc_lvl: int = 20
  tree_criterion: str = "tau"
  threshold: float = 0.0
  tree_algorithm: str = "mst_prim"
  seeds: list[int] = field(default_factory=list)
  conditioning_set: list[int] = field(default_factory=list)
  cache_integrals: Optional[bool] = None
  device: Optional[Any] = None
  dtype: Optional[Any] = None
  compile: bool = False
  batched_fit: Optional[bool] = None

  def to_dict(self) -> dict:
    """Return the settings as a plain dictionary.

    Returns
    -------
    dict
        One entry per setting, keyed by the field name. Nested controls are
        kept as objects rather than flattened, so the entry round-trips.

    See Also
    --------
    pyvinecopulib.core.ControlsLike : The contract this satisfies.
    """
    return {f.name: getattr(self, f.name) for f in fields(self)}

  def __post_init__(self) -> None:
    super().__post_init__()
    if self.tree_algorithm not in TREE_ALGORITHMS:
      raise ValueError(
        f"unknown tree_algorithm={self.tree_algorithm!r}; expected one of "
        f"{TREE_ALGORITHMS}"
      )
