Concepts
========

A self-contained tour of the ideas behind ``pyvinecopulib``. The
emphasis is on the formulas you'll see in docstrings and the API
classes that implement them. Skimming this page first should make
the rest of the documentation much easier to navigate; revisit it
whenever a docstring throws an unfamiliar symbol your way.

.. contents::
   :local:
   :depth: 2


.. _concepts-sklar:

Sklar's theorem and pseudo-observations
---------------------------------------

A *copula* is a multivariate distribution with uniform-:math:`[0, 1]`
marginals. By **Sklar's theorem** (Sklar, 1959), every joint
distribution :math:`F` on :math:`\mathbb{R}^d` factorizes as

.. math::

   F(x_1, \ldots, x_d) \;=\; C\bigl(F_1(x_1), \ldots, F_d(x_d)\bigr),

where :math:`F_1, \ldots, F_d` are the **marginal CDFs** and
:math:`C : [0, 1]^d \to [0, 1]` is the copula. When the marginals
are continuous :math:`C` is unique. Differentiating both sides gives
the corresponding density factorization

.. math::

   f(x_1, \ldots, x_d)
   \;=\;
   c\bigl(F_1(x_1), \ldots, F_d(x_d)\bigr) \,
   \prod_{j=1}^{d} f_j(x_j),

so the log-likelihood splits cleanly into a sum of marginal
log-likelihoods and a copula log-likelihood:

.. math::

   \log f(\mathbf x)
   \;=\;
   \log c\bigl(F_1(x_1), \ldots, F_d(x_d)\bigr)
   \;+\; \sum_{j=1}^{d} \log f_j(x_j).

This justifies a two-step **inference-functions-for-margins**
estimator (Joe & Xu, 1996): first fit the
:math:`d` univariate marginals; then fit the copula on the
**pseudo-observations**

.. math::

   \hat u_{ij} \;=\; \hat F_j(x_{ij}),
   \qquad i = 1, \ldots, n, \; j = 1, \ldots, d.

Both halves of that factorization are first-class here. The marginal
half is the :mod:`pyvinecopulib.margins` subpackage — nonparametric by
default (:class:`pyvinecopulib.utils.Kde1d`, a boundary-corrected 1-d
KDE that also handles ordered-discrete data), parametric with family
selection when you want it, or any distribution object you already
have. The copula half is :class:`pyvinecopulib.core.Bicop` and
:class:`pyvinecopulib.core.Vinecop`. The object that holds both, and so
gives you :math:`f` and :math:`F` on the scale of the data rather than
on the copula scale, is :class:`pyvinecopulib.core.Vinedist` — see
:ref:`concepts-distributions`.

Two convenience helpers convert raw data to pseudo-observations and
measure dependence on them:

* :func:`pyvinecopulib.utils.to_pseudo_obs` ranks each column to
  :math:`(0, 1)`, the canonical input shape for the copula classes.
* :func:`pyvinecopulib.utils.wdm` computes weighted versions of
  Kendall's :math:`\tau`, Spearman's :math:`\rho`, Blomqvist's
  :math:`\beta`, Hoeffding's :math:`D`, and Pearson correlation —
  margin-free dependence measures that drive both diagnostics and
  Dissmann's structure-selection heuristic.

Most of this page is about :math:`c` — how a vine represents and fits
it, which is where the modeling choices live. Keep in mind while
reading that the marginal factor never goes away: it is what makes the
model a distribution rather than a dependence structure, and it is the
factor that :ref:`concepts-distributions` puts back.


.. _concepts-bivariate:

Bivariate copulas
-----------------

The atomic object is the bivariate copula
:class:`pyvinecopulib.core.Bicop`. It exposes the standard surface
expected from a :math:`d = 2` distribution on
:math:`[0, 1]^2`:

* :math:`c(u_1, u_2)` — copula density —
  :meth:`pyvinecopulib.core.Bicop.pdf`;
* :math:`C(u_1, u_2)` — copula CDF —
  :meth:`pyvinecopulib.core.Bicop.cdf`;
* random sampling on :math:`[0, 1]^2` —
  :meth:`pyvinecopulib.core.Bicop.simulate`.

Pair copulas inside a vine also need **h-functions**, the partial
conditional CDFs

.. _concepts-h-functions:

.. math::

   h_1(u_1, u_2) \;&=\; \mathbb{P}(U_2 \le u_2 \mid U_1 = u_1)
   \;=\; \frac{\partial C(u_1, u_2)}{\partial u_1}, \\[4pt]
   h_2(u_1, u_2) \;&=\; \mathbb{P}(U_1 \le u_1 \mid U_2 = u_2)
   \;=\; \frac{\partial C(u_1, u_2)}{\partial u_2},

and their inverses :math:`h_1^{-1}`, :math:`h_2^{-1}`. These map to

* :meth:`pyvinecopulib.core.Bicop.hfunc1`,
  :meth:`pyvinecopulib.core.Bicop.hfunc2`,
* :meth:`pyvinecopulib.core.Bicop.hinv1`,
  :meth:`pyvinecopulib.core.Bicop.hinv2`.

H-functions are the workhorse of vine evaluation: they turn the
:math:`[0, 1]^2` outputs of one tree into the conditional
pseudo-observations consumed by the next, and their inverses drive
:meth:`pyvinecopulib.core.Vinecop.simulate` and
:meth:`pyvinecopulib.core.Vinecop.inverse_rosenblatt`.

Every :class:`~pyvinecopulib.core.Bicop` belongs to one of the families catalogued in
:ref:`concepts-families` below; switch families via
:class:`pyvinecopulib.core.FitControlsBicop`.


.. _concepts-vines:

Vine structures
---------------

Estimating a fully-flexible :math:`d`-dimensional copula directly is
hard once :math:`d` is larger than a handful. **Pair-copula
constructions** (Joe, 1996; Bedford & Cooke, 2001, 2002)
sidestep that by decomposing the joint copula density into a product
of :math:`d(d-1)/2` bivariate building blocks. The order of
conditioning is encoded by a graphical object called a *vine*.

.. admonition:: Definition (regular vine)

   A **regular vine** (R-vine) on :math:`d` variables is a sequence
   of trees :math:`(V_t, E_t)`, :math:`t = 1, \ldots, d - 1`, with

   * :math:`V_1 = \{1, \ldots, d\}`,
   * :math:`V_t = E_{t-1}` for :math:`t \ge 2`,
   * the **proximity condition**: two nodes in :math:`(V_{t+1}, E_{t+1})`
     are connected only if the corresponding edges in
     :math:`(V_t, E_t)` share a node.

A vine *copula* labels each edge :math:`e \in E_t` with a
conditioned set :math:`\{j_e, k_e\}` and a conditioning set
:math:`D_e \subset \{1, \ldots, d\} \setminus \{j_e, k_e\}`, plus a
bivariate **pair-copula** :math:`c_{j_e, k_e \mid D_e}` describing
the conditional dependence between :math:`U_{j_e}` and
:math:`U_{k_e}` given :math:`\mathbf U_{D_e}`. See Czado & Nagler
(2022) for a textbook treatment.

In pyvinecopulib the structure object is
:class:`pyvinecopulib.core.RVineStructure`, with the two well-known
specializations

* :class:`pyvinecopulib.core.DVineStructure` — every tree is a
  path; convenient for ordered-conditioning regression.
* :class:`pyvinecopulib.core.CVineStructure` — every tree is a
  star; convenient when one variable drives the others.

Both can be passed anywhere an R-vine is accepted, and
:meth:`~pyvinecopulib.core.RVineStructure.simulate` draws structures uniformly at
random (Joe, 2011) — the basis of the
:ref:`concepts-structure-selection` section below.


.. _concepts-pcc:

Pair-copula construction
------------------------

With the structure fixed, the vine-copula density factorizes as

.. math::

   c(\mathbf u)
   \;=\;
   \prod_{t = 1}^{d - 1} \prod_{e \in E_t}
   c_{j_e, k_e \mid D_e}\!\bigl(
   u_{j_e \mid D_e}, \; u_{k_e \mid D_e} \,\big|\, \mathbf u_{D_e}
   \bigr),

where the conditional pseudo-observations
:math:`u_{j_e \mid D_e} = C_{j_e \mid D_e}(u_{j_e} \mid \mathbf u_{D_e})`
are obtained tree-by-tree from h-functions of the lower-tree pair
copulas. Concretely, if an edge in tree :math:`t + 1` has
conditioned set :math:`\{j, k\}` and conditioning set
:math:`D \cup \{\ell\}`, then

.. math::

   u_{j \mid D \cup \{\ell\}}
   \;=\;
   h_1\!\bigl(u_{j \mid D}, \; u_{\ell \mid D};\;
              c_{j, \ell \mid D}\bigr),

i.e. the conditional pseudo-obs at level :math:`t + 1` is an
h-function of the pair-copula one level down. This recursive
structure is what
:meth:`pyvinecopulib.core.Vinecop.rosenblatt` evaluates forward
(data :math:`\to` independent uniforms) and what
:meth:`pyvinecopulib.core.Vinecop.inverse_rosenblatt` evaluates
backward (uniforms :math:`\to` data — used by
:meth:`pyvinecopulib.core.Vinecop.simulate`).

CDF evaluation is harder: there is no closed form for
:math:`C(\mathbf u)` in general, so
:meth:`pyvinecopulib.core.Vinecop.cdf` uses Monte-Carlo integration
with the quasi-random uniforms produced by
:func:`pyvinecopulib.utils.simulate_uniform`
(:func:`~pyvinecopulib.utils.sobol` or
:func:`~pyvinecopulib.utils.ghalton` sequences). Increase ``N`` to
trade compute for accuracy.

The same factorization backs the PyTorch port
:class:`pyvinecopulib.torch.TorchVinecop` (every pair copula is a
:class:`pyvinecopulib.torch.TorchBicop`); its cascade matches the C++
evaluator byte-for-byte and additionally offers a ``batched=True``
fast path (one stacked bicop call per tree level).


.. _concepts-simplifying:

Simplifying assumption
----------------------

The conditional pair-copulas in the previous section may, in
principle, depend on the conditioning value
:math:`\mathbf u_{D_e}`. The **simplifying assumption** states
that they do not:

.. math::

   c_{j_e, k_e \mid D_e}(u, v \mid \mathbf u_{D_e})
   \;\equiv\;
   c_{j_e, k_e \mid D_e}(u, v),
   \qquad \forall \mathbf u_{D_e}.

Under this assumption the model reduces to a collection of
two-dimensional copulas, which is what makes vines practical.
The choice of vine structure becomes load-bearing — different
structures yield different approximations of the same true
density. See Stoeber, Joe & Czado (2013), Spanhel & Kurz (2019)
and Nagler (2025) for the theoretical discussion. Every fit in
:class:`pyvinecopulib.core.Vinecop` (and hence the sklearn /
torch wrappers) uses the simplified model;
:class:`pyvinecopulib.core.FitControlsVinecop` controls *which*
families to consider, *which* structures to search, and *how* to
truncate the model in higher dimensions.

The backend-neutral :class:`~pyvinecopulib.core.VinecopBase` can go
further: a :class:`~pyvinecopulib.core.ConditioningContext` lets each
pair copula also depend on its conditioning value
:math:`\mathbf u_{D_e}` (and on external covariates), giving a
genuinely **non-simplified**, conditional vine. This is an advanced
extension point (see :ref:`concepts-extending`) — the built-in fit
stays simplified.


.. _concepts-families:

Available families
------------------

Every pair-copula in pyvinecopulib belongs to one of the families
below. The first column links the family constant, which lives on
:mod:`pyvinecopulib.families` and can be passed to
:attr:`pyvinecopulib.core.FitControlsBicop.family_set`
(or the same attribute on
:class:`pyvinecopulib.core.FitControlsVinecop`) to restrict the
fit-time search space.

Parameter ranges below are the conventional textbook ones; the
exact bounds the C++ library enforces are visible via
:attr:`pyvinecopulib.core.Bicop.parameters_lower_bounds` and
:attr:`pyvinecopulib.core.Bicop.parameters_upper_bounds`. The
"Kendall's :math:`\tau`" column says how
:meth:`pyvinecopulib.core.Bicop.parameters_to_tau` obtains the value:
in closed form, or by quadrature.

:meth:`~pyvinecopulib.core.Bicop.parameters_to_tau` is available for every
family. Its inverse :meth:`~pyvinecopulib.core.Bicop.tau_to_parameters` is
not: it needs :math:`\tau` to determine the parameters completely, which
holds only for the one-parameter families ``indep``, ``gaussian``,
``clayton``, ``gumbel``, ``frank`` and ``joe``. In particular it raises for
``student``, even though ``student`` belongs to
:data:`pyvinecopulib.families.itau` — that group means "fittable by
inverting :math:`\tau`", not "invertible here", because :math:`\tau` pins
the correlation and leaves the degrees of freedom free.

.. list-table::
   :header-rows: 1
   :widths: 14 22 15 6 16 10 13 16

   * - Family
     - Identifier
     - Type
     - Pars
     - Range
     - Rotations
     - Tail dep.
     - Kendall's :math:`\tau`
   * - Independence
     - :data:`pyvinecopulib.families.indep`
     - —
     - 0
     - —
     - none
     - none
     - :math:`0`
   * - Gaussian
     - :data:`pyvinecopulib.families.gaussian`
     - elliptical
     - 1
     - :math:`\rho \in (-1, 1)`
     - rotationless
     - none
     - :math:`(2/\pi) \arcsin \rho`
   * - Student
     - :data:`pyvinecopulib.families.student`
     - elliptical
     - 2
     - :math:`\rho \in (-1, 1)`, :math:`\nu > 2`
     - rotationless
     - symmetric
     - :math:`(2/\pi) \arcsin \rho`
   * - Clayton
     - :data:`pyvinecopulib.families.clayton`
     - Archimedean
     - 1
     - :math:`\theta > 0`
     - 0° / 90° / 180° / 270°
     - lower
     - :math:`\theta / (\theta + 2)`
   * - Gumbel
     - :data:`pyvinecopulib.families.gumbel`
     - Arch. / EV
     - 1
     - :math:`\theta \ge 1`
     - 0° / 90° / 180° / 270°
     - upper
     - :math:`1 - 1/\theta`
   * - Frank
     - :data:`pyvinecopulib.families.frank`
     - Archimedean
     - 1
     - :math:`\theta \in \mathbb{R} \setminus \{0\}`
     - rotationless
     - none
     - via Debye function
   * - Joe
     - :data:`pyvinecopulib.families.joe`
     - Archimedean
     - 1
     - :math:`\theta \ge 1`
     - 0° / 90° / 180° / 270°
     - upper
     - series expansion
   * - BB1
     - :data:`pyvinecopulib.families.bb1`
     - Arch. (2-par)
     - 2
     - :math:`\theta > 0`, :math:`\delta \ge 1`
     - 0° / 90° / 180° / 270°
     - lower + upper
     - closed form
   * - BB6
     - :data:`pyvinecopulib.families.bb6`
     - Arch. (2-par)
     - 2
     - :math:`\theta \ge 1`, :math:`\delta \ge 1`
     - 0° / 90° / 180° / 270°
     - upper
     - by quadrature
   * - BB7
     - :data:`pyvinecopulib.families.bb7`
     - Arch. (2-par)
     - 2
     - :math:`\theta \ge 1`, :math:`\delta > 0`
     - 0° / 90° / 180° / 270°
     - lower + upper
     - by quadrature
   * - BB8
     - :data:`pyvinecopulib.families.bb8`
     - Arch. (2-par)
     - 2
     - :math:`\theta \ge 1`, :math:`\delta \in (0, 1]`
     - 0° / 90° / 180° / 270°
     - upper
     - by quadrature
   * - Tawn
     - :data:`pyvinecopulib.families.tawn`
     - extreme-value
     - 3
     - bounded; see C++ bounds
     - 0° / 90° / 180° / 270°
     - upper (asymmetric)
     - by quadrature
   * - TLL
     - :data:`pyvinecopulib.families.tll`
     - nonparametric
     - —
     - —
     - data-driven
     - data-driven
     - rank-based

The non-elliptical, non-radially-symmetric parametric families are
labeled with a *rotation* in
:math:`\{0°, 90°, 180°, 270°\}`; rotation 0° is the base form
(positive dependence, lower-tail-heavy for Clayton, upper-tail-heavy
for Gumbel / Joe, …), 180° is the *survival* copula (covers the
opposite tail), and 90° / 270° provide negative dependence variants.
:meth:`pyvinecopulib.core.Bicop.flip` flips between rotations
0° :math:`\leftrightarrow` 180° and 90° :math:`\leftrightarrow` 270°.
Each fitted pair copula also reports its tail-dependence coefficients in the
four corners of the unit square as the read-only
:attr:`~pyvinecopulib.core.Bicop.taildep` matrix (NaN for ``tll``) and
Blomqvist's :math:`\beta = 4\,C(0.5, 0.5) - 1` as
:attr:`~pyvinecopulib.core.Bicop.beta`; the maps
:meth:`~pyvinecopulib.core.Bicop.parameters_to_taildep` and
:meth:`~pyvinecopulib.core.Bicop.parameters_to_beta` evaluate them at arbitrary
parameters.

Family-group constants (also in :mod:`pyvinecopulib.families`) are
pre-built lists you can pass directly to
:attr:`~pyvinecopulib.core.FitControlsBicop.family_set`:

* ``all`` — every family listed above.
* ``parametric`` — every family except ``tll``.
* ``nonparametric`` — ``indep`` and ``tll``.
* ``one_par`` / ``two_par`` / ``three_par`` — grouped by parameter
  count (one-parameter parametric / two-parameter parametric /
  three-parameter parametric).
* ``elliptical`` — ``gaussian``, ``student``.
* ``archimedean`` — ``clayton``, ``gumbel``, ``frank``, ``joe``,
  ``bb1``, ``bb6``, ``bb7``, ``bb8``.
* ``extreme_value`` — ``tawn``, ``gumbel``.
* ``bb`` — ``bb1``, ``bb6``, ``bb7``, ``bb8``.
* ``rotationless`` — families that already cover both positive
  and negative dependence (``indep``, ``gaussian``, ``student``,
  ``frank``, ``tll``).
* ``lt`` / ``ut`` — families with lower- / upper-tail dependence.
* ``itau`` — families that support estimation by Kendall's-:math:`\tau`
  inversion (``indep``, ``gaussian``, ``student``, ``clayton``,
  ``gumbel``, ``frank``, ``joe``).

The notebook ``examples/01_bivariate_copulas.ipynb`` walks through
a fit on synthetic data for several of these families.
:func:`pyvinecopulib.utils.benchmark` compares several families on
standard test problems.


.. _concepts-estimation:

Estimation
----------

Vine fitting is a two-step procedure inherited from
:ref:`concepts-sklar`:

1. **Marginals.** Each :math:`F_j` is estimated independently —
   :class:`pyvinecopulib.utils.Kde1d` (a boundary-corrected 1-d
   KDE) is the default both for the sklearn estimators and the
   notebook examples. ``Kde1d`` supports continuous,
   ordered-discrete, and unordered-categorical input via its
   ``type`` argument.
2. **Copula.** Given pseudo-observations
   :math:`\hat U_{i \cdot} = (\hat F_1(X_{i,1}), \ldots, \hat F_d(X_{i,d}))`,
   the joint copula is fit by
   :meth:`pyvinecopulib.core.Vinecop.from_data` (or
   :meth:`pyvinecopulib.core.Vinecop.select` for in-place
   re-fitting).

The pair-copula estimator on each edge depends on the chosen
family. Three regimes are available via
:class:`pyvinecopulib.core.FitControlsBicop`:

* **Maximum likelihood** (``parametric_method="mle"``) — the
  default for parametric families; numerically optimizes the
  per-edge log-likelihood under the family's parameter
  constraints.
* **Kendall's :math:`\tau` inversion**
  (``parametric_method="itau"``) — restricted to families in the
  ``itau`` group; uses
  :meth:`~pyvinecopulib.core.Bicop.tau_to_parameters` to back out the parameters from
  the empirical :math:`\hat\tau`. Cheaper than MLE and the
  conventional choice for very high-dimensional vines.
* **Nonparametric TLL** (family ``tll`` —
  :data:`pyvinecopulib.families.tll`) — *Transformed Local
  Likelihood* (Geenens, 2014; Nagler, 2018).
  The copula density is estimated on a grid in the
  inverse-normal-transformed space
  :math:`(z_1, z_2) = (\Phi^{-1}(u_1), \Phi^{-1}(u_2))`, where
  local-likelihood machinery is well-behaved at the boundary of
  the unit square. Bandwidth selection is automatic; the
  ``nonparametric_method`` and ``nonparametric_mult`` knobs on
  ``FitControlsBicop`` tune the kernel order and the
  bandwidth multiplier respectively.

TLL is the default family for both the C++ and PyTorch backends
because it captures arbitrary non-Gaussian-like dependence (heavy
tails, asymmetry) without committing to a parametric form, and
because its density-grid representation is exactly what
:class:`pyvinecopulib.torch.TorchBicop` consumes for GPU and
autograd evaluation.

Family selection across the parametric set runs by AIC / BIC /
mBIC inside :meth:`pyvinecopulib.core.Bicop.select`; choose the
selection criterion via
:attr:`pyvinecopulib.core.FitControlsBicop.selection_criterion`.

For asymptotics on a *parametric* fit, both ``Bicop`` and ``Vinecop`` expose the
log-likelihood score (``scores``), its observation-average (``gradient``), the
Hessian (``hessian``), and the score covariance (``scores_cov``) at the fitted
parameters; on ``Vinecop`` these accept an optional per-observation
``parameters`` matrix for evaluation off the fitted point.


.. _concepts-distributions:

Vine distributions
------------------

A vine copula is a model for dependence: everything on the preceding
sections lives on the copula scale, :math:`u \in [0, 1]^d`. A **vine
distribution** is Sklar's theorem as an object —
:class:`pyvinecopulib.core.Vinedist` pairs any vine with one univariate
margin per variable, and evaluates on the scale of the data:

.. code-block:: python

   dist = pv.Vinedist.from_data(x)      # Kde1d margins by default
   dist.logpdf(x)                       # log f(x), not log c(u)
   dist.simulate(1000, seeds=[1])       # draws on the x scale
   dist.cdf(x)

Without it, using a fitted vine as a distribution means recomputing
:math:`\hat F_j(x_j)` by hand, assembling the copula-scale matrix,
evaluating :math:`\log c`, and adding :math:`\sum_j \log \hat f_j(x_j)`
back — a chain that is easy to get subtly wrong, and one whose discrete
case needs a second column per variable with atoms.

The marginal contract
~~~~~~~~~~~~~~~~~~~~~

A margin is anything with ``pdf``, ``cdf`` and ``icdf`` (the inverse
cdf). That is the whole required surface — the
:class:`pyvinecopulib.core.MarginLike` protocol — and it is deliberately
small, because every member added is one a foreign distribution object
must happen to have. ``scipy.stats`` frozen distributions,
``torch.distributions`` objects and ``Kde1d`` are all accepted;
:func:`pyvinecopulib.margins.as_margin` adapts what needs adapting and
is idempotent, and
:func:`pyvinecopulib.margins.register_margin_adapter` teaches it about
an ecosystem it does not know.

One definition does most of the work. ``pdf`` means *the density with
respect to the margin's own reference measure*: a Lebesgue density for a
continuous margin, a probability mass at an atom, and whichever applies
pointwise for a mixed one. Write :math:`w_j` for that quantity,

.. math::

   w_j \;=\;
   \begin{cases}
     f_j(x_j) & \text{$x_j$ in the continuous part,} \\[2pt]
     F_j(x_j) - F_j(x_j^-) & \text{$x_j$ an atom,}
   \end{cases}

and the log-density identity of :ref:`concepts-sklar` holds verbatim in
every case:

.. math::

   \log f(\mathbf x)
   \;=\;
   \log c(\mathbf u) \;+\; \sum_{j=1}^{d} \log w_j .

There is no discreteness branch in the likelihood path — the discrete
pair-copula density already divides by the marginal mass
(:ref:`concepts-discrete`), so :math:`w_j` cancels exactly where it
should. This is why ``pdf`` and not ``pmf`` is the name on the protocol,
and it is worth knowing when adapting a foreign object: in SciPy's newer
distribution API ``pdf`` at an atom is :math:`+\infty`, and the mass is
``pmf``, so a discrete SciPy object must be routed through
``as_margin`` rather than passed straight through.

Everything past ``pdf`` / ``cdf`` / ``icdf`` is optional and discovered
at run time: ``var_type`` (``"c"``, ``"d"``, ``"zi"``), ``cdf_left``
for :math:`F(x^-)`, ``logpdf``, ``simulate``, ``support``. Each has a
correct continuous default, so a margin that declares none of them
behaves as a continuous margin.

Choosing margins
~~~~~~~~~~~~~~~~

``margins=`` accepts a string alias, one margin broadcast across
columns, a sequence of length :math:`d`, or a dict keyed by column:

.. code-block:: python

   pv.Vinedist.from_data(x)                       # "kde" (the default)
   pv.Vinedist.from_data(x, margins=EmpiricalMargin())
   pv.Vinedist.from_data(df, margins={"income": MarginSelector(),
                                      "score": st.norm(0, 1)})

Margins follow the same construct-then-``fit`` pattern as ``Bicop``,
``Vinecop`` and ``Kde1d``, with ``fit`` returning ``self``. One class is
therefore both the specification and the fitted object, which is what
lets a single ``margins=`` argument mix the two: ``from_data`` fits the
margins that are not yet fitted and leaves the already-fitted ones
alone. So ``st.norm(0, 1)`` above stays exactly :math:`N(0, 1)` while
``MarginSelector`` estimates its family from the ``income`` column.

:class:`pyvinecopulib.margins.MarginSelector` fits every admissible
candidate and keeps the best by AIC, BIC or AICc, reporting the rest.
Two choices in it are deliberate. The candidate set is **curated and
partitioned by support**, not "every family in SciPy": an unconstrained
sweep is actively misleading, because a family whose reported support is
wider than its density integrates over can win on likelihood without
being a candidate any statistician would accept. And a candidate that
fails is **reported with a reason** rather than silently skipped — a
column where everything fails falls back to a KDE margin with a
warning, never to a normal, since marginal misspecification distorts
the pseudo-observations and biases the copula.

The status quo as a special case
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:func:`pyvinecopulib.utils.to_pseudo_obs` is the scaled empirical cdf,
so fitting a ``Vinecop`` on ``to_pseudo_obs(x)`` *is* fitting a
``Vinedist`` with :class:`pyvinecopulib.margins.EmpiricalMargin`
margins. The familiar workflow did not become a different thing; it
became the configuration in which the margins are the rank transform.

What two-step estimation costs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``from_data`` estimates the margins first and then the copula on the
resulting pseudo-observations — the inference-functions-for-margins
estimator (Joe & Xu, 1996). It is consistent and it is what makes
high-dimensional vines tractable, but it is not full maximum
likelihood, and three consequences are worth stating:

* **Standard errors from the copula step alone are too small**, because
  they condition on :math:`\hat F_j` as if it were :math:`F_j`. Honest
  inference needs the sandwich form that accounts for the marginal
  estimation (Godambe, 1991).
* **Marginal error propagates into the copula, not the reverse.** A
  misspecified margin distorts every pseudo-observation and therefore
  the fitted dependence; a misspecified copula leaves the marginal fits
  untouched. Spend the modeling effort accordingly.
* **Family selection is itself an estimation step.** A likelihood or an
  interval computed at the selected margin ignores the selection, so
  ``report_`` is there to be read — a winner that beat the runner-up by
  a fraction of an AIC unit is not an established family.

``examples/11_vine_distributions.ipynb`` works through the whole
surface, including a mixed continuous / count example and a custom
margin written against :class:`pyvinecopulib.core.MarginBase`.


.. _concepts-discrete:

Discrete data
-------------

Sklar's theorem does not identify a unique copula when a margin has
atoms, so a discrete variable cannot be summarized by
:math:`\hat F_j(x)` alone: the estimator also needs the left limit
:math:`\hat F_j(x^-)`, the value the cdf takes just below the jump.
Mark which variables are discrete with ``var_types`` — ``"c"`` for
continuous, ``"d"`` for discrete, one entry per variable — and supply
the extra column each discrete variable needs.

There are two ways to lay those columns out, and every method taking
``u`` accepts either:

* **Expanded**, :math:`n \times 2d`: the :math:`d` cdf values, then
  :math:`d` left-limit values. The left-limit column of a continuous
  variable is simply equal to its cdf column.
* **Compact**, :math:`n \times (d + k)`: the :math:`d` cdf values, then
  one left-limit column per discrete variable, in variable order. With
  :math:`k` discrete variables out of :math:`d`, this is the expanded
  layout with the redundant continuous columns dropped.

The two agree exactly — the compact form only omits columns that carry
no information. For an all-continuous model :math:`k = 0` and both
collapse to the familiar :math:`n \times d` matrix.

You only assemble those columns when you are working with a copula
directly. :class:`pyvinecopulib.core.Vinedist` builds the compact layout
itself from each margin's ``cdf`` and ``cdf_left``, so a model with
atom-carrying margins takes the plain :math:`n \times d` data matrix and
the left limits never surface in user code (:ref:`concepts-distributions`).
The reason the same likelihood formula covers continuous, discrete and
mixed variables without a branch is that the pair-copula density below
divides by the marginal mass, which is exactly the factor the marginal
term contributes.

The same rule applies pairwise to :class:`pyvinecopulib.core.Bicop`
(:math:`n \times 4` expanded, :math:`n \times (2 + k)` compact) and to
the conditioning values passed to
:meth:`pyvinecopulib.core.Vinecop.simulate_conditional`, where a
discrete conditioning variable likewise contributes two columns.

The ``examples/04_discrete_variables.ipynb`` notebook works an example
end to end, from raw counts to a fitted vine;
``examples/11_vine_distributions.ipynb`` shows the same model as a
distribution, with the layout handled for you.

Two arguments elsewhere in the API are consequences of the same
non-uniqueness, and both change what you get rather than only how fast
you get it:

* ``randomize_discrete`` (on
  :meth:`pyvinecopulib.core.Vinecop.rosenblatt`, default ``True``)
  decides what the transform returns at an atom. The conditional
  distribution there is an interval, not a point, so the transform draws
  uniformly within :math:`[F(x^-), F(x)]`. That is what makes the output
  genuinely uniform — the randomized transform is the one whose
  distributional statement holds — but it also makes the call
  non-deterministic. Pass ``seeds`` to reproduce it, or
  ``randomize_discrete=False`` to take the upper endpoint instead and
  accept output that is not uniform.
* ``step_wise`` (on the ``Vinecop`` score family: ``scores``,
  ``gradient``, ``hessian``, ``scores_cov``) selects which likelihood is
  differentiated. ``True`` differentiates the step-wise (sequential,
  tree-by-tree) estimator that was actually fitted; ``False``
  differentiates the joint likelihood. They answer different questions:
  the step-wise gradient vanishes at the fitted model by construction,
  the joint one does not. Sandwich standard errors for a sequentially
  fitted vine want the step-wise form.


.. _concepts-serialization:

Saving and loading models
-------------------------

``Bicop``, ``Vinecop`` and ``RVineStructure`` all serialize the same
way. ``to_file`` / ``from_file`` write and read a file; ``to_json`` /
``from_json`` do the same through a string:

.. code-block:: python

   vine.to_file("model.json")             # JSON text
   vine.to_file("model.cbor")             # binary CBOR
   reloaded = pv.Vinecop.from_file("model.cbor")

The format follows the filename: a name ending in ``.cbor`` selects
binary `CBOR <https://cbor.io>`_, anything else JSON. CBOR is smaller
and faster to parse and is the better choice for large vines; JSON stays
the default because it is readable and diffable.

All three classes also pickle, which serializes through the JSON form,
so a fitted model round-trips through anything that speaks pickle
(``copy.deepcopy``, ``joblib``, a multiprocessing queue). The
:mod:`pyvinecopulib.sklearn` estimators pickle as ordinary estimators.


.. _concepts-structure-selection:

Structure selection
-------------------

The vine structure :math:`\mathcal V` is rarely known in advance,
and the number of regular vines on :math:`d` variables grows as
:math:`2^{(d-3)(d-2)/2 - 1} d!` (Morales-Napoles, 2011; Joe, 2011)
— super-exponential, so exhaustive search is infeasible beyond a
handful of variables. Two algorithms are exposed via
:attr:`pyvinecopulib.core.FitControlsVinecop.tree_algorithm`:

* ``"mst_prim"`` — Dissmann's greedy heuristic (Dissmann et al.,
  2013), which is the default. Builds the trees one at a time as
  maximum-spanning-trees weighted by absolute Kendall's :math:`\tau`
  (or the criterion of your choice — see
  :attr:`~pyvinecopulib.core.FitControlsVinecop.tree_criterion`).
  It is fast and remains the de-facto standard (Czado & Nagler,
  2022).
* ``"random_weighted"`` — Wilson-weighted random spanning trees;
  draws a random structure with edge probabilities proportional to
  the absolute dependence on each candidate edge. Use it when you
  want a randomized alternative to the greedy structure — for
  example to draw several candidates and compare their held-out
  likelihoods.

Selection is further tunable through
:class:`~pyvinecopulib.core.FitControlsVinecop`: supply a custom edge weight via
:attr:`~pyvinecopulib.core.FitControlsVinecop.tree_criterion_function` (any
callable ``(data, weights) -> float`` returning a scalar dependence measure),
make random searches reproducible with
:attr:`~pyvinecopulib.core.FitControlsVinecop.seeds`, and toggle whether rotated
families are considered with
:attr:`~pyvinecopulib.core.FitControlsVinecop.allow_rotations`.

Vine truncation
~~~~~~~~~~~~~~~

For sparse models in higher dimensions, set
:attr:`pyvinecopulib.core.FitControlsVinecop.trunc_lvl` to an
integer :math:`T < d - 1`: pair copulas in trees
:math:`T + 1, \ldots, d - 1` are forced to ``indep``, reducing
both fit time and statistical degrees of freedom. The same effect
can be triggered automatically via the mBIC criterion (Nagler,
2019).

Conditional sampling and conditioning-aware vines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A fitted vine can sample from the *conditional* distribution of a subset of
variables given fixed values of the rest, via
:meth:`pyvinecopulib.core.Vinecop.simulate_conditional`. The conditioning
variables are the last ``k`` of the vine order
(:attr:`~pyvinecopulib.core.Vinecop.order`): each row of ``u_cond`` is one
conditioning point and the corresponding output row is drawn from the remaining
variables' distribution conditional on that point (implemented as a Rosenblatt
transform of the conditioning variables followed by an inverse Rosenblatt
transform, so discrete conditioning variables are supported too).

The conditioning variables can also be named explicitly, with
``simulate_conditional(u_cond, conditioning_set=[...])`` — 1-based indices.
Note that the two forms map ``u_cond``'s columns differently: without the
argument, column ``i`` is the ``i``-th variable of the order tail; with it,
column ``i`` is ``conditioning_set[i]``. The same argument is available on
:meth:`pyvinecopulib.core.Vinecop.rosenblatt` and
:meth:`~pyvinecopulib.core.Vinecop.inverse_rosenblatt`, which then hold those
variables fixed rather than the order tail. In every case the vine is
evaluated through a reoriented view and the model is left unchanged.

Not every set is admissible as a sampling-order tail; when it is not, a
``RuntimeError`` says so. Conditioning is cheapest when the set already sits at
the tail, and two tools arrange that:

* :attr:`pyvinecopulib.core.FitControlsVinecop.conditioning_set` — select a vine
  whose order ends with a chosen set of variables (their own optimal sub-vine is
  fit first, then placed at the tail). Requires an MST ``tree_algorithm``.
* :meth:`pyvinecopulib.core.Vinecop.reorient` — relabel an already-fitted vine
  to an equivalent one whose order tail equals a given set, without refitting.
  This is value-preserving: ``pdf`` and ``loglik`` are invariant.

These operate on the *simplified* vine's exact conditional. For fully
non-simplified / conditional pair copulas on a custom (e.g. neural) backend, see
the extending guide (example notebook 10) instead. A fitted vine's tree-by-tree
decomposition — the conditioned pairs, conditioning sets, and pair-copulas — is
available as nested lists through
:meth:`pyvinecopulib.core.Vinecop.get_trees`, and a bare structure round-trips
through :meth:`pyvinecopulib.core.RVineStructure.get_trees` /
:meth:`pyvinecopulib.core.RVineStructure.from_trees`.


.. _concepts-regression:

Vine regression
---------------

For conditional inference, fix one variable as the response
:math:`Y` and stack it with the predictors to model the joint
distribution :math:`(Y, \mathbf X)`. Any conditional summary
of interest is then solved out of the estimating equation
(Nagler, 2018)

.. math::

   \int \psi_\beta(y) \,
   \hat f(y \mid \mathbf x)
   \, dy \;=\; 0,

with :math:`\psi_\beta(y) = y - \beta` recovering the
conditional mean :math:`\beta(\mathbf x) = \mathbb E[Y \mid \mathbf X = \mathbf x]`
and :math:`\psi_\beta(y) = \mathbb 1\{y < \beta(\mathbf x)\} - \tau`
recovering the conditional :math:`\tau`-quantile. In practice the
integral is taken on the probability scale, substituting
:math:`p = \hat F_Y(y)`, and replaced by a weighted sum over a
fixed grid of levels :math:`\{p_1, \ldots, p_G\} \subset (0, 1)`:

.. math::

   \sum_{g = 1}^G
   \psi_\beta\!\bigl(\hat F_Y^{-1}(p_g)\bigr) \,
   \hat c_{\mathcal V, \mathcal D}\!\bigl(
   p_g, \, \hat F_{\mathbf X}(\mathbf x)
   \bigr) \, \Delta_g
   \;\approx\; 0,

with :math:`\Delta_g` the spacing of the nodes. The substitution
cancels the marginal density :math:`\hat f_Y` that the
response-scale form of the same sum carries, so the rule asks the
response margin for nothing but its inverse CDF and places every
node inside its support by construction.

:class:`pyvinecopulib.sklearn.VineRegressor` solves this
numerically. Pass the quantile levels you want via the
``quantiles=`` constructor argument; the predicted conditional
mean is always returned when ``mean=True``. The node count is
``n_nodes=``, and ``use_grid=False`` swaps the quadrature for the
exact weighted sum over the training responses. Batching over test
rows is controlled by ``batch_size=``.


Where to next
-------------

* :class:`pyvinecopulib.core.Bicop` and
  :class:`pyvinecopulib.core.Vinecop` — the C++/nanobind classes
  that implement everything above. The notebooks
  ``examples/01_bivariate_copulas.ipynb``,
  ``examples/02_vine_copulas.ipynb``, and
  ``examples/03_vine_copulas_fit_sample.ipynb`` walk through
  end-to-end use.
* :mod:`pyvinecopulib.sklearn` — scikit-learn-compatible
  estimators :class:`~pyvinecopulib.sklearn.VineDensity` and
  :class:`~pyvinecopulib.sklearn.VineRegressor`. The notebook
  ``examples/08_sklearn_estimators.ipynb`` demonstrates them. Both
  estimators accept a backend (default C++, optional PyTorch) via
  :mod:`pyvinecopulib.sklearn.backends`.
* :mod:`pyvinecopulib.torch` — PyTorch evaluators
  :class:`~pyvinecopulib.torch.TorchBicop` and
  :class:`~pyvinecopulib.torch.TorchVinecop` for GPU placement and
  autograd. Notebook ``examples/09_torch_backend.ipynb``.
* :mod:`pyvinecopulib.utils` —
  :class:`~pyvinecopulib.utils.Kde1d` for the marginals (notebook
  ``examples/07_kde1d.ipynb``);
  :func:`~pyvinecopulib.utils.wdm` for weighted dependence
  measures (notebook ``examples/06_weighted_dependence_measures.ipynb``);
  :func:`~pyvinecopulib.utils.sobol`,
  :func:`~pyvinecopulib.utils.ghalton`,
  :func:`~pyvinecopulib.utils.simulate_uniform` for the
  low-discrepancy sequences that back Monte-Carlo CDF evaluation;
  :func:`~pyvinecopulib.utils.to_pseudo_obs` and
  :func:`~pyvinecopulib.utils.pairs_copula_data` for input
  preparation and pair-plot diagnostics.
* The :doc:`features` page is the autogenerated API reference; the
  :doc:`examples` toctree lists all worked notebooks.

The ``examples/04_discrete_variables.ipynb`` notebook covers the
discrete-margin extension (Panagiotelis, Czado & Joe, 2012; Funk,
Nagler & Czado, 2025), which replaces the marginal CDF derivatives in
:ref:`concepts-sklar` by finite differences (transparent to the
user — pass ``var_types=["d", ...]`` to
:meth:`pyvinecopulib.core.Vinecop.from_data` or set
``type="d"`` on :class:`pyvinecopulib.utils.Kde1d`).

A discrete variable needs its left limit :math:`F(x^-)` alongside
:math:`F(x)`, and there are two ways to supply them. The **expanded** layout is
``(n, 2d)``: the first :math:`d` columns are :math:`F(x)`, the next :math:`d`
are :math:`F(x^-)` for the same variables in the same order, with continuous
columns simply repeated. Its shape does not depend on ``var_types``, which is
what makes it the easier convention to write against. The **compact** layout is
``(n, d + k)``: the same first :math:`d` columns, then left limits for the
:math:`k` discrete variables only, in the order they appear. Both are accepted
wherever discrete data is; for an all-continuous model they coincide at
``(n, d)``. The same applies per pair to
:meth:`pyvinecopulib.core.Bicop.from_data`, with :math:`d = 2`.


.. _concepts-extending:

Extending: custom and conditional pair copulas
-----------------------------------------------

The evaluators :class:`pyvinecopulib.core.Bicop` /
:class:`~pyvinecopulib.core.Vinecop` and their PyTorch counterparts
:class:`pyvinecopulib.torch.TorchBicop` /
:class:`~pyvinecopulib.torch.TorchVinecop` are concrete implementations
of two backend-neutral contracts, evaluated on either NumPy or PyTorch
arrays:

* :class:`~pyvinecopulib.core.BicopLike` — a pair copula, exposing
  ``pdf`` / ``cdf`` / ``hfunc1`` / ``hfunc2`` / ``hinv1`` / ``hinv2`` /
  ``simulate``;
* :class:`~pyvinecopulib.core.VinecopLike` — a fitted vine, exposing
  ``pdf`` / ``cdf`` / ``rosenblatt`` / ``inverse_rosenblatt`` /
  ``simulate`` on an :class:`~pyvinecopulib.core.RVineStructure`.

You can plug your **own** pair copula into a vine by implementing the
contract — most easily by subclassing the canonical partial
implementations :class:`~pyvinecopulib.core.BicopBase` /
:class:`~pyvinecopulib.core.VinecopBase`, which fill in almost
everything from a few primitives. A ``BicopBase`` subclass need only
define ``pdf`` / ``hfunc1`` / ``hfunc2`` and inherits numerical
``hinv1`` / ``hinv2``, ``simulate``, ``loglik``, and ``plot``; a
``VinecopBase`` subclass need only return its pairs from
``_get_pair_copula`` and inherits the whole tree-by-tree cascade. The
bases are pure Python (no PyTorch), so custom pairs also work in a
torch-less environment.

Every method carries an optional trailing conditioning matrix ``x``.
For the common **simplified, unconditional** vine it is ``None``
everywhere (the default
:class:`~pyvinecopulib.core.SimplifiedContext`). To lift the
:ref:`simplifying assumption <concepts-simplifying>`, host the pairs
under a :class:`~pyvinecopulib.core.NonSimplifiedContext`: the cascade
then assembles each edge's conditioning-set values
:math:`\mathbf u_{D_e}` (and any external covariates) into ``x`` and
threads them to the pair copula, giving a **non-simplified /
conditional** vine. :meth:`pyvinecopulib.core.VinecopBase.fit`
is the seam for *fitting* such a vine edge by edge.

The ``examples/10_extending_pyvinecopulib.ipynb`` notebook is a
worked, end-to-end walk-through: a custom Gaussian pair copula hosted
first in a simplified vine (matching
:meth:`pyvinecopulib.core.Vinecop.from_structure`) and then made
non-simplified and conditional.


References
----------

* **Sklar (1959).** *Fonctions de répartition à n dimensions et
  leurs marges.* Publ. Inst. Statist. Univ. Paris 8, 229–231.
* **Bedford & Cooke (2001, 2002).** *Probability density
  decomposition for conditionally dependent random variables
  modeled by vines* / *Vines — a new graphical model for dependent
  random variables.* Annals of Mathematics and Artificial
  Intelligence 32, 245–268 / Annals of Statistics 30(4),
  1031–1068.
* **Joe (1996).** *Families of m-variate distributions with given
  margins and m(m-1)/2 bivariate dependence parameters.* In:
  *Distributions with Fixed Marginals and Related Topics* (IMS
  Lecture Notes 28), 120–141.
* **Joe & Xu (1996).** *The estimation method of inference
  functions for margins for multivariate models.* Technical
  Report 166, Department of Statistics, University of British
  Columbia.
* **Godambe (1991).** *Estimating Functions.* Oxford Statistical
  Science Series 7, Oxford University Press.
* **Aas, Czado, Frigessi & Bakken (2009).** *Pair-copula
  constructions of multiple dependence.* Insurance: Mathematics
  and Economics 44(2), 182–198.
* **Dissmann, Brechmann, Czado & Kurowicka (2013).** *Selecting
  and estimating regular vine copulae and application to financial
  returns.* Computational Statistics & Data Analysis 59, 52–69.
* **Joe (2011).** *Dependence comparisons of vine copulae with four
  or more variables.* In: *Dependence Modeling: Vine Copula
  Handbook*, 139–164.
* **Geenens (2014).** *Probit Transformation for Kernel Density
  Estimation on the Unit Interval.* JASA 109(505), 346–358.
* **Nagler (2018).** *A Generic Approach to Nonparametric Function
  Estimation with Mixed Data.* Statistics & Probability Letters
  137, 326–330.
* **Nagler, Schepsmeier, Stoeber, Brechmann, Graeler & Erhardt
  (2018).** *VineCopula: Statistical inference of vine copulas.*
  R package.
* **Nagler (2019).** *Model selection in sparse high-dimensional
  vine copula models with an application to portfolio risk.*
  Journal of Multivariate Analysis 172, 180–198.
* **Spanhel & Kurz (2019).** *Simplified vine copula models:
  approximations based on the simplifying assumption.* Electronic
  Journal of Statistics 13(1), 1254–1291.
* **Czado & Nagler (2022).** *Vine Copula Based Modeling.* Annual
  Review of Statistics and Its Application 9, 453–477.
* **Panagiotelis, Czado & Joe (2012).** *Pair Copula Constructions
  for Multivariate Discrete Data.* JASA 107(499), 1063–1072.
* **Funk, Nagler & Czado (2025).** *Discrete and mixed
  pair-copula constructions revisited.* (In press.)
