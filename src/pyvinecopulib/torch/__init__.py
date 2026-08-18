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
Likelihood* (Geenens 2014 [1]_; Nagler 2018 [2]_) — a non-parametric
pair-copula estimator that fits a kernel density on a grid in the
inverse-normal-transformed copula space. This is the same family
exposed as ``tll`` and is the
default everywhere in pyvinecopulib because it captures arbitrary
non-Gaussian-like dependence without picking a parametric form.

If you have not used vine copulas before, the
:doc:`concepts page </concepts>` introduces pair copulas and R-vines
in ~5 minutes.

Requires PyTorch. Install with ``pip install pyvinecopulib[torch]``.

See Also
--------
pyvinecopulib.core : Reference vine-copula evaluators (default everywhere).
pyvinecopulib.sklearn : sklearn-compatible vine-copula estimators that route through either backend.

Notes
-----
**What's exposed.**

- :class:`TorchBicop` — evaluator for a single bivariate pair copula
  on a density grid. Build from a fitted
  ``Bicop`` via :meth:`TorchBicop.from_bicop`, or
  fit directly from data via :meth:`TorchBicop.from_data` (which
  dispatches on :class:`FitControlsTorchBicop`: pure-torch TLL).
  Exposes the standard surface:
  ``pdf`` / ``cdf`` / ``hfunc1`` / ``hfunc2`` / ``hinv1`` / ``hinv2``
  / ``simulate``.

- :class:`TorchVinecop` — evaluator for a full R-vine built on top of
  :class:`TorchBicop` pair copulas. Provides ``pdf`` / ``cdf`` /
  ``rosenblatt`` / ``inverse_rosenblatt`` / ``simulate`` with the same
  signatures as ``Vinecop``. The cascade mirrors
  ``Vinecop`` byte-for-byte; ``pdf`` / ``rosenblatt``
  also accept ``batched=True`` (one stacked bicop call per tree level).

- :class:`TorchMargin` — a univariate margin on a ``torch.distributions``
  family. Registers the family's parameters and rebuilds the distribution on
  each call, so ``.to(device)``, ``state_dict()`` and autograd all reach them.
  Continuous families only, and only those implementing ``cdf``.

- :class:`TorchVinedist` — a full multivariate distribution: a
  :class:`TorchVinecop` plus one margin per variable, with ``logpdf`` / ``pdf``
  / ``cdf`` / ``simulate`` / ``rosenblatt`` on the original data scale. Every
  part is a registered child module, so the joint log-density is differentiable
  with respect to the marginal parameters.

- :class:`FitControlsTorchBicop`, :class:`FitControlsTorchVinecop` —
  fit-time control dataclasses mirroring
  ``FitControlsBicop`` /
  ``FitControlsVinecop``. Bundle method selection
  and method-specific knobs.

References
----------
.. [1] Geenens, G. (2014). *Probit Transformation for Kernel Density
       Estimation on the Unit Interval.* JASA 109(505), 346–358 —
       original TLL motivation.
.. [2] Nagler, T. (2018). *A Generic Approach to Nonparametric
       Function Estimation with Mixed Data.* Statistics & Probability
       Letters 137, 326–330 — multivariate TLL.
"""

try:
  import torch  # noqa: F401
except ImportError as e:
  raise ImportError(
    "pyvinecopulib.torch requires PyTorch. "
    "Install it with `pip install pyvinecopulib[torch]`."
  ) from e

from .bicop import TorchBicop
from .margin import TorchMargin
from .vinecop import TorchVinecop
from .vinedist import TorchVinedist
from .controls import FitControlsTorchBicop, FitControlsTorchVinecop

__all__ = [
  "TorchBicop",
  "TorchMargin",
  "TorchVinecop",
  "TorchVinedist",
  "FitControlsTorchBicop",
  "FitControlsTorchVinecop",
]
