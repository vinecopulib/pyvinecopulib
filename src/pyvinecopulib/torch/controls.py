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

from dataclasses import dataclass, field
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
  """Controls for :meth:`~pyvinecopulib.torch.TorchVinecop.from_data` and the cascade.

  Mirrors :class:`~pyvinecopulib.core.FitControlsVinecop`: bundles all vine-fit
  knobs into one object. A nested ``FitControlsTorchBicop`` controls
  how each pair-copula is fit; the vine-level fields below carry the
  structure-selection knobs plus placement / precision /
  cascade-variant settings.

  Attributes
  ----------
  bicop_controls : FitControlsTorchBicop
      Controls applied to every pair-copula fit.
  trunc_lvl : int, default=20
      Maximum number of trees to select when
      :meth:`~pyvinecopulib.torch.TorchVinecop.from_data` is called with
      ``structure=None``.
  tree_criterion : {"tau", "rho", "hoeffd", "mcor", "cxi", "joe"}, \
default="tau"
      Dependence measure used to weight candidate edges during structure
      selection (passed to ``wdm``).
  threshold : float, default=0.0
      Dependence threshold: candidate edges below it are deprioritized during
      spanning-tree selection.
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
      If ``True``, precompute the cdf / hfunc / hinv caches on
      every pair copula's interpolation grid. Cached lookups are
      1–2 orders of magnitude faster than the on-the-fly path
      with a ~1e-3 IAE cost. ``None`` resolves to ``True``, on a
      discrete vine too: the tables reconstruct the integral exactly,
      so a discrete edge can difference them.
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

      That trade pays where launches cost something, and it does not pay
      everywhere. On CUDA a whole vine fit is roughly 1.2-2.7x faster over
      ``d`` in 5..20 and ``n`` in 2000..12000, most at moderate size and
      least at large ``n``, where the fit is arithmetic-bound and there was
      little overhead to remove. **The gain is not monotone: at ``d = 20``
      and large ``n`` it reverses to about 0.95x**, the level's working set
      having grown past what batching buys, so that corner is a few percent
      slower unless the flag is switched off. Worth it against 2-3x at
      moderate size, but a real corner rather than a rounding error.

      Sizing that corner needs care, and the numbers here come from
      interleaving the two arms over five repetitions rather than from a
      single run: on a thermally throttling laptop a lone measurement of
      these cells moves by half again, enough to invent a reversal where
      there is none or hide one that is there.

      On cpu it is 0.44-1.05x at one torch thread -- mostly slower, there
      being no launch overhead to amortize. Many threads favor it again,
      the larger kernels parallelizing better than many small ones, and
      torch uses every core by default, so the cpu default is the
      conservative reading of a measurement that moves with thread count
      rather than a claim that batching cannot pay there. Either way the
      torch fit is far from competitive with the compiled backend on cpu.

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
  """

  bicop_controls: FitControlsTorchBicop = field(
    default_factory=FitControlsTorchBicop
  )
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

  def __post_init__(self) -> None:
    if self.tree_algorithm not in TREE_ALGORITHMS:
      raise ValueError(
        f"unknown tree_algorithm={self.tree_algorithm!r}; expected one of "
        f"{TREE_ALGORITHMS}"
      )
