Concepts
========

A short tour of the ideas behind ``pyvinecopulib`` for readers
arriving from scikit-learn or PyTorch. Five minutes here will make
the rest of the API documentation much easier to read.

What is a copula?
-----------------

By **Sklar's theorem**, every joint distribution
:math:`F(x_1, \ldots, x_d)` factorises as

.. math::

   F(x_1, \ldots, x_d) \;=\; C\bigl(F_1(x_1), \ldots, F_d(x_d)\bigr),

where :math:`F_1, \ldots, F_d` are the **marginal CDFs** and
:math:`C: [0, 1]^d \to [0, 1]` is the **copula** — a multivariate CDF
on the unit cube with uniform marginals. The marginals carry the
shape of each variable on its own; the copula carries the dependence
between them.

Plugging the chain rule in gives an equally clean factorisation of
the joint density:

.. math::

   f(x_1, \ldots, x_d) \;=\; c\bigl(F_1(x_1), \ldots, F_d(x_d)\bigr)
   \,\prod_{j=1}^{d} f_j(x_j),

so estimating :math:`f` reduces to (a) fitting :math:`d` univariate
densities and (b) fitting one copula density :math:`c` on the
pseudo-observations :math:`u_j = \hat F_j(x_j)`.

The first half is well-trodden — pyvinecopulib uses
:class:`pyvinecopulib.utils.Kde1d` for the marginals. The interesting
part is :math:`c`.

Pair-copula construction and R-vines
------------------------------------

Estimating a fully-flexible :math:`d`-dimensional copula directly is
hard once :math:`d` is larger than a handful. The **pair-copula
construction** (PCC) of Bedford & Cooke (2002) sidesteps that by
decomposing :math:`c` into a product of **bivariate** building blocks
indexed by a tree structure called an **R-vine**:

* The bottom tree edges are unconditional pair copulas.
* Higher-tree edges are conditional pair copulas of the form
  :math:`c_{j, k \mid S}`, where :math:`S` is a small "conditioning"
  set.

An R-vine is a sequence of :math:`d - 1` trees encoding which
conditional pair copulas to use; each pair copula can be fit
independently of the others. Aas et al. (2009) popularised this
construction by demonstrating that pair-copula models scale far
better than direct multivariate alternatives in finance and beyond.

In ``pyvinecopulib`` the R-vine structure lives on
:class:`pyvinecopulib.RVineStructure`, and the fitted vine itself
lives on :class:`pyvinecopulib.Vinecop`. The :class:`Vinecop` exposes
the standard surface — ``pdf``, ``cdf``, ``simulate``,
``rosenblatt`` / ``inverse_rosenblatt`` — and the corresponding
``TorchVinecop`` from :mod:`pyvinecopulib.torch` mirrors that surface
in pure PyTorch for GPU / autograd workflows.

The TLL family
--------------

The default pair-copula family in ``pyvinecopulib`` is **TLL** —
*Transformed Local Likelihood*, a non-parametric kernel estimator
introduced by Geenens (2014) and extended by Nagler (2018). TLL
estimates the copula density on a grid in the inverse-normal
(:math:`\Phi^{-1}`) transformed space, where the local-likelihood
machinery is well-behaved at the boundary of the unit square. Each
grid cell's density is set to maximise a locally-weighted likelihood
with a bandwidth chosen automatically.

Why TLL is the default:

* It captures arbitrary non-Gaussian-like dependence shapes (heavy
  tails, asymmetric dependence) without committing to a parametric
  family.
* It composes cleanly with the vine: every pair copula on every edge
  uses the same evaluator, so fits scale predictably with :math:`d`.
* It is the family the PyTorch backend
  (:class:`pyvinecopulib.torch.TorchBicop`) supports natively — the
  density grid is the natural representation for GPU / autograd
  evaluation.

The C++ library exposes the family as
:data:`pyvinecopulib.families.tll`. If you have a clear parametric
prior — Gaussian, Clayton, Gumbel, etc. — pass it via
:class:`pyvinecopulib.FitControlsVinecop.family_set` instead.

Where to next
-------------

* :class:`pyvinecopulib.Vinecop` and :class:`pyvinecopulib.Bicop` —
  the C++ classes that do the actual fitting and evaluation.
* :mod:`pyvinecopulib.sklearn` — scikit-learn-compatible vine-copula
  estimators (:class:`~pyvinecopulib.sklearn.VineDensity`,
  :class:`~pyvinecopulib.sklearn.VineRegressor`, plus forest
  variants). The sklearn module wraps the marginal + copula
  pipeline behind a single ``.fit`` / ``.predict`` interface and
  can route either through the C++ default or the PyTorch backend.
* :mod:`pyvinecopulib.torch` — PyTorch evaluators for
  :class:`~pyvinecopulib.torch.TorchBicop` and
  :class:`~pyvinecopulib.torch.TorchVinecop`. Pick this when you
  need GPU placement, autograd, or composition with other
  ``torch.nn.Module`` code.
* :mod:`pyvinecopulib.utils` — :class:`Kde1d` for the marginals,
  plus quasi-random uniform generators (``sobol``, ``ghalton``,
  ``simulate_uniform``) used by quasi-MC integration.

References
----------

* Bedford, T. & Cooke, R. M. (2002). *Vines — a new graphical model
  for dependent random variables.* Annals of Statistics 30(4),
  1031–1068.
* Aas, K., Czado, C., Frigessi, A. & Bakken, H. (2009). *Pair-copula
  constructions of multiple dependence.* Insurance: Mathematics and
  Economics 44(2), 182–198.
* Geenens, G. (2014). *Probit Transformation for Kernel Density
  Estimation on the Unit Interval.* JASA 109(505), 346–358.
* Nagler, T. (2018). *A Generic Approach to Nonparametric Function
  Estimation with Mixed Data.* Statistics & Probability Letters 137,
  326–330.
* Cheng, B., Vatter, T., Nagler, T. & Chen, V. (2025). *Vine Copulas
  as Differentiable Computational Graphs.* arXiv:2506.13318 — basis
  for the lazy / batched torch cascades.
* Safaai, H. (2026). *Amortized Vine Copulas for High-Dimensional
  Density and Information Estimation.* arXiv:2604.20568 — the
  pretrained estimator behind ``method="vdc"`` on
  :class:`~pyvinecopulib.torch.TorchBicop`.
