# Publication plan: pyvinecopulib in JMLR MLOSS

Planning document for a JMLR Machine Learning Open Source Software (MLOSS)
paper on pyvinecopulib, written on 2026-10-29 against `main` at 1.1.0.dev0
(#369 merged). Every count, date and quotation below is as of that day.

The manuscript skeleton is `pyvinecopulib.tex`, beside this file, in the
official JMLR format; `./build.sh` fetches the style file and builds the PDF.
Each of its sections holds the plan from §3 in place of prose.

One question is open: whether the GPU methods paper is the JSS "companion"
paper. §10 plans for both cases.

## 1. Executive recommendation

**Venue: MLOSS remains right.** What pyvinecopulib adds over the shared C++
engine is software architecture and ecosystem integration, which is what MLOSS
reviews. The statistical methods and the engine are claimed by the JSS
manuscript and its companion paper, so a second JSS paper would overlap, and
JOSS would undersell it. MLOSS admits prior publication of the *methods*
provided the *software* has not been published. pyvinecopulib has not been;
because the engine it shares will be, the repository should carry a page
listing what pyvinecopulib adds (MLOSS asks for exactly that; Appendix C).

**The repository meets most MLOSS criteria already.**

- Tests are extensive, with coverage close to 100% in the maintainers' runs.
- Documentation is complete, and there are ten executed notebooks.
- CI covers five platforms.
- Adoption is strong for a specialized library. pyvinecopulib's downloads lead
  every other Python vine package: about 3× torchvinecopulib, which itself
  depends on pyvinecopulib. Dependents span finance, atmospheric physics,
  bioinformatics and ML research.
- Two core maintainers lead development, with external contributors to both
  pyvinecopulib and the vinecopulib engine.

**Fixable gaps:**

- CI does not publish the coverage figure.
- No source file carries a license notice.
- CUDA and the declared dependency floors are untested in CI.
- conda-forge lags behind PyPI.

**Framing: composable vine distributions.** The thesis is Candidate B: four
array-agnostic contracts plus algorithms written once against them. Candidate
C, level-wise execution, is the mechanism that shows the contracts do real
work, because one cascade hosts compiled, NumPy, learnable-PyTorch, stacked-TLL,
discrete and conditional pairs. Candidates A (integration) and D
(differentiability) are consequences shown in the listing and in Table 1.

D should not lead. Differentiable vines and torchvinecopulib were already
presented in the authors' own arXiv:2506.13318 (Cheng, Vatter, Nagler, Chen,
2025). pyvinecopulib's PyTorch lane first implemented that computational-graph
formulation, then replaced it with level-wise execution because that is more
efficient. What is new here is the architecture that made the replacement a
change to one execution model rather than to every algorithm.

**Why it does not duplicate JSS.** The JSS paper owns vine theory, the engine
(computation index, families, TLL, selection, derivatives), validation against
VineCopula, R benchmarks and two applications. The MLOSS paper owns the Python
object model, how one cascade hosts any implementation of a part, scikit-learn
and PyTorch integration, and fidelity *between implementations inside the
package*. Background shrinks to three sentences plus citations. Four
cross-paper fixes are needed first (§9.3), including a JSS table that
currently overclaims for pyvinecopulib.

**Version: freeze on pyvinecopulib 1.1.0** (tagged, with its Zenodo version
DOI), not 1.0.1. #369 is what makes the architectural claim true for every
vine: in 1.0.1, level batching was a `TorchVinecop` special case behind
`NotBatchable`, and nothing in CI held a vine hosting `Bicop` to `Vinecop`. Tag
1.1.0 after the P0 items in §11, cite only what 1.1.0 ships, and allow 1.1.x
patch releases during review.

**Main risks.**

- **Prior work by the same authors (arXiv:2506.13318).** Cite it as the
  predecessor of the PyTorch lane. Say in one sentence what replaced it and
  why, and disclose it in the cover letter (§6, §10, §12).
- **"Is this ML?"** Frame vines as density-estimation and generative
  components, and cite ML uses (§12).
- **Consistency with the JSS manuscript.** Fix §9.3 first.

---

## 2. Paper thesis

### 2.1 Framing decision

| Candidate | Role | Reason |
| --- | --- | --- |
| A. Scientific-Python integration | consequence | "Integration" reads as glue; a reviewer asks what idea makes it more than adapters |
| **B. Composable architecture (four contracts)** | **thesis** | Python-only (in R, pair-copula families are a closed set, JSS §5.3; non-simplified vines are not provided); A, C and D follow from it |
| C. Level-wise execution | mechanism for B | Shows the contracts do real work: one cascade, any pair class. As the whole thesis it slides into methods-paper territory |
| D. Differentiable vine distributions | consequence (Listing 1) | The authors' own arXiv:2506.13318 already presented differentiable vines; the lane implemented that formulation first and replaced it with level-wise execution for efficiency. Leading with D would re-present the preprint's idea rather than what is new. The torch lane's built-in pair family is also TLL only |

**Central phrase: *composable vine distributions*.** It is stated once in the
introduction as: *"the vine algorithms are written once, against four
contracts; every implementation of a part — compiled, NumPy, PyTorch or the
user's own — inherits them, and decides how its tree levels execute."*

Two vocabulary rules:

- AGENTS.md forbids "backend" for the NumPy-vs-PyTorch distinction. The paper
  should say "array library", "lane" or "execution strategy", so that it
  matches the documentation it points to.
- The repository's banned-word list (`.codespell-prose.txt`) applies to the
  manuscript sources too.

### 2.2 Titles

1. **pyvinecopulib: Composable Vine Copula Distributions for Scientific Python
   and Machine Learning** *(preferred)*
2. pyvinecopulib: A Library for Vine Copula Distributions as Composable NumPy,
   scikit-learn and PyTorch Components
3. pyvinecopulib: Vine Copula Models from Compiled Estimation to Trainable
   PyTorch Modules

Why title 1:

- It keeps the user's wording.
- Adding "Copula" makes it findable: readers search for "vine copula", while
  "vine distribution" is the library's own term.
- It names the thesis.
- It avoids "Scalable", the working title in the JSS bib, which invites
  GPU-speed claims that belong to the methods paper and to 2506.13318.
- Running title (≤ 50 characters): "Composable Vine Copula Distributions".

### 2.3 Abstract-level thesis (≈ 135 words; MLOSS allows ≤ 200)

> Vine copulas build flexible multivariate distributions from bivariate
> building blocks, and the C++ library vinecopulib provides mature estimation
> for them. pyvinecopulib turns these models into components of Python
> machine-learning workflows. Its core is four array-agnostic contracts — for
> pair copulas, vines, univariate margins and joint distributions on the data
> scale — with base classes that implement structure selection, the
> evaluation and simulation cascades and two-step estimation once, against
> those contracts. Compiled vinecopulib objects, NumPy and PyTorch
> implementations and user-defined parts all satisfy the same contracts, so
> each inherits the same algorithms, while the pair copulas hosted in a vine
> decide how each tree level executes: as compiled calls, on a thread pool, or
> as one stacked GPU operation. On this foundation the library provides
> scikit-learn estimators and differentiable PyTorch modules, validated
> against the compiled engine.

### 2.4 Contributions

1. **Four stable contracts for vine models.** `BicopLike`, `VinecopLike`,
   `MarginLike` and `VinedistLike`, with canonical bases, make pair copulas,
   vines, margins and full distributions interchangeable. This includes
   user-defined, learnable, discrete-edge and non-simplified (conditional)
   pair copulas, which plug into the same selection, evaluation, simulation
   and estimation code as the compiled ones.
2. **One level-wise execution model.** Each vine algorithm is expressed once,
   as tree levels and, for the inverse transform, dependency waves. The hosted
   pair class decides whether a level runs as per-pair compiled calls, on a
   thread pool, or as one stacked tensor operation. This gives GPU execution
   and autograd without a second implementation of any algorithm.
3. **Vine distributions in the Python ML ecosystem.** The package provides
   scikit-learn estimators and PyTorch modules for full data-scale
   distributions with continuous, discrete and mixed margins (kernel, SciPy,
   `torch.distributions` or user-defined). These are validated against
   vinecopulib (identical selected structures, densities to 1e-12 relative)
   and tested on five platforms. The scikit-learn and PyTorch suites run on
   three of them: Linux, macOS arm64 and Windows.

---

## 3. Four-page manuscript architecture

Constraints taken from the MLOSS rules and the three precedents (Appendix A):

- 4 pages of main text from the title page, plus references.
- Lettered appendices are common practice (9 of 10 volume-27 MLOSS papers;
  median total length 10 pages). Use a modest one: ≤ 3 pages, as OptunaHub
  does.
- Full-width Figure 1 at the top of page 2.
- One `§2` paragraph per figure box.
- At most one listing of 15 lines or fewer.
- No equations.
- A one-paragraph conclusion ending on page 4.
- Acknowledgments at the top of page 5.

| § | Pages | Argument | Evidence | Display | Omitted |
| --- | --- | --- | --- | --- | --- |
| Title, abstract, keywords | 0.45 | — | — | — | — |
| 1 Introduction | 0.7 | Vines are useful ML density models (tabular density estimation, synthetic data, conditional simulation, learned dependence). Python has narrow vine packages and differentiable libraries without vines (Table 1). A binding alone cannot host user or learnable parts, run on accelerators, or sit in sklearn/PyTorch pipelines. Contribution paragraph and three bullets. Methods, engine and R interface are cited to [JSS]; the earlier computational-graph formulation to [2506.13318] | §6 | Table 1 (≈ 0.3 p; too wide to wrap) | vine theory beyond three sentences; families; selection algorithms |
| Figure 1 | 0.4 | — | — | Figure 1 (§4) | — |
| 2 Design | 1.15 | **2.1 Contracts.** Protocol + base: compiled classes satisfy the protocols structurally (a nanobind class cannot inherit a Python ABC), bases derive everything else, unserved members have defaults, fitting by declaration (`bicop_class`, `vinecop_class`, `margin_class`, `controls_class`). **2.2 Algorithms written once.** Cascades; selection that reproduces the engine's matrix; IFM; conditional sampling; discrete edges through `var_types`; conditional vines through a context; placement inferred for NumPy/PyTorch, device and dtype. **2.3 Execution.** The engine's computation index [JSS §3.1] wires whole levels; the inverse runs in dependency waves; a level runs as its pairs' stack or per pair, and is fitted stacked, on a thread pool, or edge by edge. One sentence: an earlier version of the lane followed the computational-graph formulation of [2506.13318], and level-wise execution replaced it for efficiency | parity tests; code structure | — | index internals; wave algorithm; kernels |
| 3 Vine distributions in ML workflows | 0.55 | scikit-learn (`VineDensity`, `VineRegressor`, `distribution=` a class, mixed DataFrames, `check_estimator`); PyTorch (`nn.Module` throughout, `.to`, `state_dict`, autograd, discrete and trainable margins); SciPy / `torch.distributions` / user margins via adapters | notebooks 08–10 | Listing 1 (§5, 17 lines) | API enumeration; controls |
| 4 Fidelity and execution | 0.4 | Implementations agree; each execution strategy wins somewhere; one sentence on tests, coverage, CI, typing, documentation | Table 2 (§7) | Table 2 (≈ 0.34 p) | R comparisons; GPU analysis |
| 5 Related software and availability | 0.25 | How Table 1 is read; complementary libraries; MIT; PyPI/conda-forge; docs; version 1.1.0 and its DOI; pointers to [JSS] and the methods paper | §6 | — | roadmap promises |

Total ≈ 3.9 pages: one figure, two tables, one listing. If the pages overflow,
move Table 1 to Appendix A and keep a two-sentence summary in §1, as OptunaHub
and MarkDiffusion do.

**Measured by the skeleton** (`pyvinecopulib.tex`, built with `jmlr2e`; page
fractions of the 9-inch text height):

| Display | Page fraction |
| --- | --- |
| Title block (title, two authors, abstract, keywords) | 0.55 |
| Figure 1 with its caption | 0.36 |
| Table 1 with its caption | 0.30 |
| Listing 1 with its caption (17 lines) | 0.43 |
| Table 2 with its caption | 0.34 |

The four displays take about 1.4 pages, which leaves about 2.0 pages for
prose. The table above implies 2.5 pages of prose, so half a page has to go:

- move Table 1 to Appendix A and keep two sentences in §1 (0.30), as OptunaHub
  and MarkDiffusion do;
- shorten §2 to 1.0 page (0.15);
- if needed, move the listing's imports into its caption (0.05).

The skeleton runs to page 5 now only because its plan boxes and section heads
stand in for that prose.

**Appendix (≤ 3 pages, before the references).**

- A. Extended comparison with Python software: the full verified matrix with
  sources.
- B. Fidelity and timing setup: hardware, versions, full tables including
  memory.
- C. The complete code for Listing 1, and the hosted-`Bicop` reproduction of
  `Vinecop`.
- D. What changed from the 0.x bindings to 1.x, mirroring the repository's
  scope page and following torchgfn's "What's New in v2".

---

## 4. Architecture figure (Figure 1)

**Goal.** In ten seconds the reader should see three things:

- a narrow, stable layer of four contracts, with many consumers above and many
  implementations below;
- which implementations come from C++ and which Python adds;
- that one vine algorithm runs over all of them, each tree level executed the
  way its pairs choose.

It is an hourglass, not a file diagram. Precedent: every Figure 1 in the
three papers is full width, about a third of a page, and vector art.

**Layout, top to bottom.** Full text width, about 0.4 page, greyscale-safe,
one accent color per implementation family.

1. **Consumers.** Three boxes in a row:
   - *scikit-learn*: `VineDensity`, `VineRegressor` (`distribution=` a class);
   - *PyTorch*: `nn.Module` composition, `.to(device)`, `state_dict`,
     autograd, optimizers;
   - *NumPy / SciPy user code*.

   One arrow down, labeled "depend only on".
2. **The waist: four contracts.** One rounded box labeled "stable since 1.0":
   - `MarginLike` and `VinecopLike` side by side, both with arrows into
     `VinedistLike` ("d margins + one copula");
   - `BicopLike` under `VinecopLike`, with an arrow labeled "hosted at each
     (tree, edge)";
   - a one-line note inside the box: "protocol (structural) + base (derives
     the rest)".
3. **Implementations.** Four columns, each aligned under the contracts it
   implements:
   - *compiled, from vinecopulib* (color A): `Bicop`, `Vinecop`, `Kde1d`;
   - *generic bases* (color B): `BicopBase`, `VinecopBase`, `MarginBase`, and
     `VinedistBase` → `Vinedist`;
   - *PyTorch lane* (color C): `TorchTllBicop`, `TorchVinecop`, `TorchKde1d` /
     `TorchDistributionMargin`, `TorchVinedist`;
   - *adapters and user parts* (color D, dashed border): `SciPyMargin`,
     `as_margin(...)`, custom / learnable / non-simplified pairs, custom
     margins.
4. **Inset beside `VinecopBase`: "one cascade, run level by level".**
   - Tree levels T1–T3, plus a small tag reading "inverse: dependency waves".
   - From one level, a decision diamond: "does the pairs' class supply a
     stack?"
     - *yes*: "one stacked call, e.g. all TLL pairs of a level on the GPU";
     - *no*: "each pair's own call, identical to edge by edge; fitted on a
       thread pool when the pair releases the GIL (`Bicop`)".
5. **Foundation.** A thin color-A band: "vinecopulib · kde1d · wdm (C++17,
   Eigen), shared with rvinecopulib [JSS]", with "nanobind" arrows from the
   color-A classes.

**Rules.**

- No file names and no private names: `LoopedPairs` and `_stack_pairs` stay
  out.
- The caption is three lines at most and states the reading order.
- Draw it in TikZ (as the JSS figures are), or as SVG exported to PDF.
- Reuse the same figure on a docs architecture page (§13).

---

## 5. Code example (Listing 1)

| Candidate | Shows | Misses | Lines |
| --- | --- | --- | --- |
| A. `TorchVinedist.from_data(y_cuda)`, then `log_prob(y).backward()` | data scale, GPU, autograd | composition; why it is more than a binding | 4 |
| B. `VineDensity(distribution=TorchVinedist)` on a mixed DataFrame, then `est.distribution_` on CUDA with autograd | sklearn and torch, mixed types | custom parts | 8 |
| **C. A user-defined learnable pair copula, hosted by the stock `TorchVinecop` and composed with fitted kernel margins into a `TorchVinedist`, trained with Adam on the GPU** | all five requested properties | sklearn (stated in the prose) | 17 |

**Selected: C**: 17 lines in the repository's ruff style (2-space indentation,
one statement per line), as the notebooks are written. The precedents stop at
15; if space runs short, the three imports move into the caption. The pair's three formulas are abridged in the main text;
the complete class, a Gaussian with `ndtr` / `ndtri`, already runs as notebook
10 §1.1 and goes in Appendix C.

```python
import torch
from pyvinecopulib.core import BicopBase
from pyvinecopulib.torch import TorchVinecop, TorchVinedist


class Gaussian(BicopBase[torch.Tensor], torch.nn.Module):  # learnable pair
  def _pdf_raw(self, u): ...  # density and the two conditional
  def _hfunc1_raw(self, u): ...  # distribution functions; parameter
  def _hfunc2_raw(self, u): ...  # and constructor in Appendix C


y = torch.as_tensor(data, device="cuda")  # raw data, on the GPU
fit = TorchVinedist.from_data(y)  # kernel margins + vine
pairs = [[Gaussian() for _ in range(d - 1 - t)] for t in range(d - 1)]
copula = TorchVinecop(pairs, fit.vinecop.structure)  # stock vine, user pairs
dist = TorchVinedist(copula, fit.margins)  # an nn.Module
opt = torch.optim.Adam(dist.parameters(), lr=0.05)
for _ in range(200):  # copula MLE, data scale
  opt.zero_grad()
  (-dist.log_prob(y).mean()).backward()
  opt.step()
```

The API names were checked against the source: `VinedistBase.from_data(y, controls, *, ...)`,
the `.vinecop` and `.margins` properties, `TorchVinedist(vinecop, margins)`,
`TorchVinecop(pair_copulas, structure)` registering the pairs as nested
`ModuleList`s, and `log_prob`. The full listing must still be executed
before submission (§13).

**What the listing shows:**

- **Composition through the contracts.** A pair writes three continuous
  leaves. It inherits `hinv*` (bisection), the log-likelihood, discrete-edge
  dispatch and plotting, and the stock vine and distribution classes host it
  unchanged. Sampling needs one more hook, `_sample_uniform`, because NumPy
  and PyTorch draw random numbers differently.
- **Data scale.** The loss is the joint log-density of raw `y` (kernel
  margins × copula), not of pseudo-observations.
- **Ordinary Python conventions.** `nn.Module`, `.parameters()`, a stock
  optimizer.
- **The ML ecosystem.** It is PyTorch end to end, on the GPU.
- **Absent from the JSS paper.** R's pair families are closed and nothing
  there is differentiable.

The surrounding prose adds two sentences:

- A network of covariates `x` in place of the constant correlation gives a
  conditional, non-simplified vine (notebook 10 §3).
- Naming `bicop_class` and `vinecop_class` on subclasses makes `from_data`
  select such a vine, and lets `VineDensity(distribution=...)` fit it through
  scikit-learn.

Fallback, if even 14 lines does not fit: candidate B.

---

## 6. Related software

Source-level check, 2026-10-29: 30 Python packages and 2 R packages. Each
package was checked against its code at a pinned commit, its PyPI metadata and
its documentation.

**Classes.**

1. **True alternatives.** These are the only maintained Python vine
   implementations besides pyvinecopulib:
   - **VineCopulas** (42★, GPL-3.0, JOSS 2024): pure Python, parametric.
   - **torchvinecopulib** (252★, MIT): PyTorch, nonparametric, continuous
     only, and it depends on pyvinecopulib. It accompanies arXiv:2506.13318,
     which the authors co-wrote; an earlier version of pyvinecopulib's PyTorch
     lane implemented the same computational-graph formulation. Present it as
     related work by the same group, and **state the relationship in the
     paper and the cover letter.**
   - **mixedvines** (50★, GPL-3.0): C-vines only, exact mixed likelihood.
   - **torchvine** (0★, created 2026-02) makes claims that contradict its own
     docstrings, so mention it at most in passing.
2. **General copula libraries with limited or no vines.**
   - Copulas (SDV, 658★) has a `VineCopula` that warns it is "not fully
     tested", fits 3 trees by default and supports only Clayton, Frank and
     Gumbel. Its license is **BUSL-1.1, which is not OSI-approved**.
   - The others have no vines: copulae (low activity), statsmodels copulas
     (experimental), OpenTURNS, copulax (JAX), chaospy, UQpy.
3. **Complementary libraries.**
   - Normalizing flows: zuko, nflows, normflows.
   - NumPyro's `GaussianCopula`.
   - scikit-learn's `KernelDensity` and `GaussianMixture`.
   - `torch.distributions` and TFP, which ship no copula class.
   - There is no JAX vine package.
4. **The R interface to the same engine.** rvinecopulib is not an
   independent competitor. VineCopula (R) is context only.

**Table 1 (compact, for the paper).** Every cell traces to a source in the
list below.

| | pyvinecopulib 1.1 | VineCopulas 2.0.3 | torchvinecopulib 1.1.2 | mixedvines 1.3.5 | Copulas 0.14.1 |
| --- | --- | --- | --- | --- | --- |
| General R-vines with selection | ✓ | ✓ | ✓ | C-vine only | 3 trees by default, "not fully tested" |
| Pair copulas: parametric / nonparametric | 11 / TLL | 6 / — | — / TLL, KDE | 3 / — | 3 / — |
| Discrete and mixed data | exact | jittered | — | exact (C-vine) | — |
| Joint distribution: density and sampling | ✓ | — | ✓ (KDE margins) | ✓ | sampling only |
| scikit-learn estimators | ✓ | — | — | — | — |
| Autograd and GPU | ✓ | — | ✓ | — | — |
| User pair copulas hosted in a vine | ✓ incl. learnable, non-simplified | — | — | ✓ subclass | — |
| License | MIT | GPL-3.0 | MIT | GPL-3.0 | BUSL-1.1 |

Statistical capability rows (truncation, Rosenblatt, conditional sampling,
scores, and so on) stay in JSS Table 6, which also covers VineCopulas. The
paper refers to it there instead of repeating them.

**Sources (pinned):**

- VineCopulas: `VU-IVM/VineCopulas@60f0207`, `vinecopula.py`, `bivariate.py#L24`,
  `marginals.py#L257-L274`; JOSS 10.21105/joss.06728.
- torchvinecopulib: `TY-Cheng/torchvinecopulib@7673ae4`, `bicop/__init__.py#L4-L5`,
  `vinecop/__init__.py` (L585-589 truncation test, L708-714 Rosenblatt raises);
  the `dev_refactor` README says it is not designed for discrete or mixed
  margins.
- mixedvines: `asnelt/mixedvines@b7b712a`, `mixedvine.py#L227-L276`,
  `copula.py#L39`.
- Copulas: `sdv-dev/Copulas@176705b`, `LICENSE`, `multivariate/vine.py`
  (L45-46, L84, L158, L328), `bivariate/__init__.py#L104-L145`.
- pyvinecopulib: the docs pages for `FitControlsVinecop`, `Vinedist`,
  `TorchVinecop`, `NonSimplifiedContext`, and notebooks 04, 05, 08 and 10.

**Runtime and memory comparison with other packages** (MLOSS criterion 4).
None of the three precedents reports one, and JSS refuses to time VineCopulas
because the language is a confound. Recommendation: in the main text, the
feature table plus within-package execution (Table 2). Optionally, in
Appendix B, two matched comparisons:

- against VineCopulas: CPU, parametric families both packages share, d ∈ {5, 10};
- against torchvinecopulib: GPU, TLL, continuous. This one measures the change
  of execution model (computational graph → levels), so its full analysis
  belongs in the methods paper; MLOSS quotes at most one number from it.

Run them only if a pilot shows the settings can be matched; otherwise say why
in one sentence, as JSS does.

---

## 7. Validation and performance plan (minimum evidence)

**Principle.** MLOSS needs evidence that the *software claims* hold:
interchangeable implementations agree, and each execution strategy is useful
in some regime. Algorithmic speedups, launch counts and crossover analysis
belong to the methods paper. MarkDiffusion's Table 3, which reports "ours vs
official" reproduction fidelity, is the model for Table 2.

**Setup.**

- Use the **v1.1.0 wheels** in float64.
- Hosts:
  - CPU: Intel Core Ultra 9 285H, 16 threads.
  - GPU: RTX 5070 Laptop, 8 GB.
  - Optionally the A100 VM used for #369.
- Time with `scripts/bench_torch_runner.py`: paired repetitions, medians with
  95% CIs, and peak memory.
- Scripts and raw JSON go in `paper/mloss/` and are archived with the
  submission.

**Table 2: one model, seven implementations.** One `tll` vine, d = 10, on the
1.1.0 wheels. Each row is an implementation. Its columns are its maximum
relative error against `Vinecop` on held-out points, then its time for
`logpdf` at n = 10^3 and n = 2·10^4, and for `sample` and `from_data` (with
selection) at n = 2·10^4. Keeping both sizes for `logpdf` shows both regimes:
C++ wins the small one.

| Row | Implementation | What it shows |
| --- | --- | --- |
| 1 | `Vinecop`, 1 thread | the reference |
| 2 | `Vinecop`, all cores | C++ threads |
| 3 | `VinecopBase` hosting `Bicop`, all cores | the cost of the generic layer; same structure, densities to `rtol=1e-12` |
| 4 | `TorchVinecop`, CPU | the torch lane without a GPU |
| 5 | `TorchVinecop`, CUDA, pair by pair | the loop over pairs |
| 6 | `TorchVinecop`, CUDA, stacked | agrees with row 5 to rounding; its inverse exactly |
| 7 | row 6, compiled | evaluation only |

Two checks go in the caption: `TorchVinedist` against `Vinedist` with mixed
margins (data-scale `logpdf`), and a threaded level fit equal to the serial
one. The skeleton's Table 2 has exactly this shape.

Earlier measurements already show both directions, for example "16-thread C++
still wins at d = 20, n = 2000 on the selection path" (#316); the paper says so
in one sentence and leaves the analysis to the methods paper.

**Software-quality numbers** (for §4 and the cover letter):

- the coverage figure for the release, as published by the CI coverage job
  (§11). The CUDA host adds its branches;
- collected test count;
- counts of CUDA- and `compile`-marked tests run on the release;
- notebook wall time.

**Not measured for MLOSS:**

- comparisons with R packages (JSS covers them);
- TLL accuracy against the truth (`bench_torch_precision.py`, methods paper);
- memory scaling curves;
- kernel-launch counts;
- CUDA-graph replay gains.

---

## 8. MLOSS readiness audit

Criteria are quoted from the MLOSS page (Appendix A). Status values: ready,
needs evidence, needs work.

| Criterion | Status | Evidence or gap | Action |
| --- | --- | --- | --- |
| Impact: novelty, breadth, significance | ready | The only Python package combining general R-vines, TLL, exact discrete data, sklearn, torch and hostable user pairs (Table 1) | Sharpen the claim in the paper |
| Community: an active user community, shown in the cover letter | ready | Downloads lead every Python vine package; ~55 dependent repositories; 38 external issue authors; external contributors to both repositories (§12) | Take a dated snapshot and verify the named dependents (§12) |
| Openness: repository, tracker, forum | ready, small gap | GitHub repo and issues; email to info@vinecopulib.org; no forum | Enable GitHub Discussions |
| Comparison with previous implementations (run time, memory, features) | needs work | Features verified (Table 1); timings exist only in PR bodies | Tables 1–2; optional Appendix B runtime comparison |
| Installation instructions | ready, with gaps | README covers pip, conda and source; wheels for 5 platforms. conda-forge stops at **1.0.0** and has **no osx-arm64** | conda-forge 1.1.0 with osx-arm64 |
| Tutorials | ready | Ten notebooks, executed in CI | Fix notebook 09's stale text; remove the network fetch in 08/09 |
| Non-trivial examples | ready | Notebooks 03, 08, 10 | — |
| API documentation as a website | ready | Read the Docs; autosummary; numpydoc-validated; nitpicky `-W` | — |
| Ease of contribution | mostly ready | CONTRIBUTING is clear; 8 external contributors have landed code across pyvinecopulib and vinecopulib. The design spec is a 150 KB `AGENTS.md` | A short human-oriented architecture page in front of it; optionally labels and Discussions |
| Clarity of software design | ready | Four layers pinned by `test_import_surface.py`; four contracts | Architecture page reusing Figure 1 |
| Developer documentation | mostly ready | CONTRIBUTING, AGENTS.md, notebook 10, the "Extending" part of concepts | Architecture page; fix the "each of each of" typo in CONTRIBUTING.md:8 |
| Unit and integration tests | ready | 1,070 test functions in 44 modules; parity suites; `check_estimator`; pickling; notebooks; docstring examples; import-surface and prose tests | — |
| Coverage "close to 100%" | ready, needs reporting | Close to 100% in the maintainers' runs, but CI publishes no figure: every leg passes `--no-cov`, and nothing is uploaded | Add a CI coverage job with a Codecov badge (rvinecopulib already uses Codecov); quote the release figure in the cover letter |
| CI on all platforms, several dependency versions | mostly ready, needs work | 5 platforms; Python 3.11–3.14; Eigen 3.4 and 5.0.1; sdist install. Dependency floors are never resolved; CUDA never runs in CI | Add a `--resolution lowest-direct` leg; attach CUDA logs to each release (or add a self-hosted GPU runner) |
| Platform breadth, including an open-source OS | ready | Linux glibc and musl, macOS x86-64 and arm64, Windows | Linux aarch64 wheels are optional |
| Freedom of the code | ready | MIT; open-source dependencies; CUDA optional through PyTorch | — |
| Checklist: license stated in the source package **and in each source file** | **needs work** | 0 of 65 files in `src/` carry a notice (the vinecopulib headers do) | SPDX/MIT header in every source file; a third-party notice for Eigen (MPL-2.0), Boost (BSL-1.0), nanobind (BSD-3) |
| Checklist: compiles and runs on every supported platform | ready | CI | — |
| Checklist: no extra files in the source package | ready | Allowlisted sdist plus a CI assertion | — |
| Exact version to review | needs work | 1.1.0 not tagged | §11 |
| Archive | ready | Zenodo concept DOI 10.5281/zenodo.4288292; version DOIs are automatic (v1.0.1 → 10.5281/zenodo.23184037). The post-acceptance `*-code.tar.gz` is a manual step | — |

---

## 9. JSS overlap audit (`vinecopulib/rvinecopulib#348`, head `13323078`)

**What the manuscript is.**

- *"rvinecopulib: A Fast and Comprehensive Implementation of Vine Copula
  Models in C++ and R"*, Nagler & Vatter.
- 40 A4 pages, with the body ending on page 35.
- Versions: rvinecopulib 1.0.0.1.0 and vinecopulib 1.0.0.
- Its division of labor (L126–135): *"Every algorithm for copulas and vines
  discussed below is implemented in vinecopulib. rvinecopulib contributes the
  R object model, the marginal-modeling layer, and the validation of both."*
- It never mentions PyTorch, GPU, scikit-learn, autograd or NumPy.
- Its only pointer to this paper (L2184–2185): *"A Python interface to the
  same engine, with machine-learning ecosystem integration, is described
  separately"*.
- It cites a **companion** paper, *"Algorithms for Vine Copula Models:
  Conditional Simulation, Nonparametric Estimation with Mixed Data, and
  Differentiable Likelihoods"* (in preparation).

### 9.1 Ownership

| Topic | JSS owns | MLOSS owns | Methods / companion |
| --- | --- | --- | --- |
| Vine theory, Sklar, R-vines, simplifying assumption | §2, Fig. 1 | three sentences and a citation | — |
| Engine: natural order, triangular array, **computation index** (min-array, two h-function masks), O(nd) memory, nested threading | §3, §3.1, Fig. 2 | cites it; `VinecopBase` *reads* this index to wire whole levels | — |
| Families (11 + independence + TLL), rotations, the TLL estimator | §3.2–3.3, Table 2 | one clause; `TorchTllBicop` is presented as a reimplementation, not explained | TLL-in-torch numerics → methods |
| Structure selection, truncation, thresholding, mBICV | §3.4, §4.3, §6.1 | only that `VinecopBase.select` re-expresses it over any pair class and reproduces the engine's matrix | — |
| Discrete / mixed / zero-inflated data | §2.2, §4.4, §5.1 | how a *custom* pair joins a discrete edge (`with_var_types`, a `_cdf_raw`); that torch batches discrete levels | discrete TLL formulation → methods / companion |
| Vine distributions and margins (R: `vine()`, univariateML, S3 margin protocols) | §2.3, §4.5, §5.1 | `Vinedist` as a composite of contracts: array-agnostic, a torch module, an sklearn target; SciPy and `torch.distributions` adapters | — |
| Conditional simulation, conditioning-aware selection, Rosenblatt | §4.6 (algorithm in the companion) | one clause: available on the data scale and in torch | companion |
| Scores, Hessians, `vcov`, coverage study, derivative cache | §3.5, §4.7 | nothing | companion |
| Custom tree criteria | §5.2 | nothing | — |
| Extensibility | margins only; **pair families are a closed set (§5.3)**; no non-simplified vines | **the contribution** | — |
| Applications (S&P 500, insurance claims) | §6 | none | — |
| Validation against VineCopula (Table 7) | §7.2 | fidelity *between implementations inside the package* | — |
| Runtime against VineCopula / kdecopula (Tables 8–9) | §7.3–7.4 | execution strategies inside the package | crossover analysis → methods |
| Capability matrix (Table 6) | §7.1 | Python table on ML dimensions only | — |
| R object model, S3, Rcpp | §4 | — | — |

### 9.2 Material not to repeat from #348

- R-vine theory beyond three sentences, and any vine drawing.
- The computation index, its masks and the skip fraction φ.
- Any family table or list of families.
- TLL estimator details.
- Dissmann/MST/Wilson selection, truncation and mBICV.
- Derivatives, scores, `vcov`, and the coverage study.
- The discrete-data algebra.
- IFM theory.
- Either application.
- Any comparison against VineCopula, kdecopula, CDVine or MATVines.
- Any statistical-capability row of JSS Table 6.
- The R/C++ division of labor.
- Custom tree criteria.

### 9.3 Cross-paper fixes (before either submission)

1. **JSS Table 6 overclaims for pyvinecopulib.** Its caption calls the
   pyvinecopulib column "identical by construction", but it ticks *Standard
   errors* and *Robust (sandwich) covariance*. Those come from R-level
   `vcov()` / `confint()`, which pyvinecopulib lacks: it has `scores` /
   `hessian` / `scores_cov`, but no covariance or interval. Fix the two cells
   and reword the caption: the engine features are identical, the interface
   layers are not.
2. **Title.** The JSS bib cites *"pyvinecopulib: Scalable Vine Copula Modeling
   for Python and Machine Learning"*. Change it to the title chosen in §2.2.
3. **Engine version.** JSS uses vinecopulib 1.0.0; pyvinecopulib ≥ 1.0.1 pins
   v1.0.1, whose TLL and discrete fixes change fits. The JSS claim of
   "identical numbers whenever they ship the same engine version" is true but
   conditional, and both papers should say so.
4. **Post the JSS manuscript to arXiv before the MLOSS submission.** The MLOSS
   paper can then cite it for everything it hands off, and the editors can
   check the boundary themselves.

---

## 10. Methods-paper boundary

Two facts constrain this section:

- The JSS companion already claims "differentiable likelihoods", meaning
  analytic derivatives.
- **arXiv:2506.13318**, by the same authors, presents the vine computational
  graph (a DAG of vine computations), sampling-order scheduling,
  conditional-sampling algorithms and GPU torchvinecopulib. pyvinecopulib's
  PyTorch lane first implemented that formulation, then replaced it with
  level-wise execution because the latter is more efficient.

So the methods paper presents level-wise execution as the successor to the
computational-graph implementation and quantifies the gain against it. The
MLOSS paper cites 2506.13318 once, as the earlier formulation.

The "Methods" column below reads "the GPU methods paper". **If that paper is
the companion** (undecided), the same rows move there, and the companion's
title should drop "Differentiable Likelihoods" or say "analytic derivatives",
so that it is not confused with 2506.13318.

| Contribution | MLOSS main text | MLOSS docs / appendix | Methods |
| --- | :---: | :---: | :---: |
| Four contracts and bases; fitting by declaration | ✓ core | API docs, notebook 10 | — |
| Algorithms written once (cascades, selection port, IFM, serialization) | ✓ one paragraph | concepts page | — |
| Level-wise execution (stack, loop or thread pool, chosen by the pair class) | ✓ concept, Fig. 1 inset | architecture page | formal cost model |
| Dependency waves of the inverse Rosenblatt | named | docs | ✓, with the gain over the 2506.13318 implementation |
| Stacked TLL fitting (`from_data_batched`, dropping converged lanes, fused bandwidth search, O(n) window smoother) | capability, one number | — | ✓ |
| TLL in torch (closed-form `hinv`, exact prefix-table integrals, compensated discrete masses, transposition-equivariant normalization) | not mentioned | — | ✓ |
| `torch.compile` / CUDA graphs | a control that exists | docs | ✓ |
| CPU vs GPU crossover | Table 2, both regimes, no analysis | raw data archived | ✓ |
| Discrete/mixed TLL numerical formulation | — | — | ✓ or companion |
| Conditional simulation on accelerators | capability | notebooks 05, 09 | ✓, compared with the 2506.13318 implementation |
| Parity methodology (normalization identity, bit-identity of loop and per-edge evaluation) | one sentence, Table 2 | tests | optional |
| Neural / covariate-dependent pair copulas | one sentence (extension point) | notebook 10 §3 | separate methodology, if pursued |

---

## 11. Review version: freeze on 1.1.0

**Why 1.1.0 and not 1.0.1.** #369 (merged 2026-10-29) is what makes
contribution 2 true for every vine:

- one cascade per transform, shared by simplified and conditional vines;
- level stacking decided by the pair class (`NotBatchable` removed);
- threaded level fits for compiled pairs, bit-identical to the serial fit;
- hosted-`Bicop` ↔ `Vinecop` parity in CI.

**Correction to the brief.** That parity is **not** bit-for-bit on the release
wheels. FMA contraction (x86-64-v3 wheels, clang on arm64) compiles the pair
functions differently inside `Vinecop`'s cascade. CI guarantees the same
selected structure matrix and densities to `rtol=1e-12`, over 16 cases:

- `tll` and `itau` families;
- continuous and mixed data;
- given and selected structures;
- serial and threaded fits.

The bit-exact claims that do hold:

- threaded = serial level fit;
- per-pair loop = per-edge walk;
- `TorchVinecop.inverse_rosenblatt`, batched = per edge;
- `VinecopBase.select` matrix = `Vinecop.select` matrix;
- the array-agnostic `sample_conditional` / `reorient` = their compiled
  versions.

**Before tagging 1.1.0 (P0):**

1. License headers in all 65 source files, plus a third-party notice.
2. Notebook 09: delete "the device / dtype targets" (removed in #339).
3. Notebooks 08 and 09: stop downloading the UCI Concrete data at run time
   (vendor it, CC BY 4.0, or cache it).
4. A coverage job in CI that publishes the figure (a Codecov badge), so
   reviewers can see it without running the suite.
5. `docs/citing.rst` (Appendix C) and an architecture page.
6. From issue #348:
   - trusted publishing;
   - Read the Docs badges pointing at `stable`;
   - the `find_libclang.py` pattern, if it is still missing.
7. The mandatory CUDA rehearsal from `docs/releasing.md` step 5, with its raw
   JSON attached to the release. It doubles as the paper's GPU evidence.
8. Small fixes:
   - the comment at `torch/vinedist.py:202-204` (it says `Vinedist`; the base
     is `VinedistBase`);
   - CONTRIBUTING.md:8, "each of each of".

**Before submission, in 1.1.0 or 1.1.1 (P1):**

- the dependency-floor CI leg;
- conda-forge 1.1.0 including osx-arm64.

Optional: GitHub Discussions and "good first issue" labels.

**Freeze rule.** The cover letter names 1.1.0 and its version DOI. 1.1.x patch
releases are fine during review. No 1.2 feature is cited, even if it lands
before the paper is accepted.

---

## 12. Cover letter and community evidence

**Required content.**

- **MLOSS:**
  - the track statement;
  - the license (MIT);
  - the URLs (GitHub, Read the Docs);
  - the version (1.1.0 and its DOI);
  - evidence of an active user community;
  - suggested reviewers.
- **JMLR in general:**
  - disclosure of overlapping publications, with the earlier papers attached:
    the JSS manuscript, arXiv:2506.13318, and the companion if it is public;
  - co-author consent;
  - conflicts of interest;
  - 3–5 action editors;
  - 3–5 reviewers (not co-authors of 2506.13318);
  - keywords.
- **Point to `docs/citing.rst`** as the repository document "listing the
  software package's improvements and extensions".

**Evidence collected 2026-10-29, robust first:**

| Signal | Value | Use |
| --- | --- | --- |
| Named dependents | ~55 repositories and ~23 packages (GitHub's estimate), including `hudson-and-thames/arbitragelab` (703★, finance), `TY-Cheng/torchvinecopulib` (252★), `open-atmos/PySDM` (93★, cloud microphysics), `PriorLabs/relarena` (72★, relational-learning benchmark), `dmey/synthia` (70★, synthetic data), `KempnerInstitute/vine-denoising-copula` (NeurIPS 2026 VDC paper, arXiv:2604.20568), `maiziezhoulab/FEAST` (spatial transcriptomics), `Techtonique/nnetsauce`, normalizing-flow research repos (`bossemel/Flow_PC`, `CMFlows`) | lead with 6–8 of these; read each one's usage first |
| ML research use | Kempner Institute VDC (NeurIPS 2026); arXiv:2605.03061 cites the Zenodo record | yes, after a Scholar search |
| Downloads | PyPI (non-mirror): 29,581 in the last 30 days; 121,030 over 2026-04-08…10-28 (13.8k–17.8k per month from May to August, ~25–26k in September and October). conda-forge: 161,632 cumulative. Same 30 days for comparison: torchvinecopulib 10,141 (which installs pyvinecopulib), copulae 6,069, vinecopulas 194, mixedvines 56 | yes; lead with the comparison, not the raw count |
| External issues | 52 issues from 38 distinct non-maintainer accounts, 2019–2026 (17 in 2024) | yes |
| External code contributions | pyvinecopulib: slayoo, dmey, lorepirri, jondo (10 pull requests, 8 merged). vinecopulib engine: slayoo (72 commits), jschueller (OpenTURNS; 9 merged pull requests), n5 (6 commits), bodashtart (1) | yes |
| Core maintenance | Vatter and Nagler (squash-merged commits: pyvinecopulib 102 / 8; vinecopulib 361 / 271) | state as is |
| Stars and forks | pyvinecopulib 130 / 23, starred steadily from 2019-08 to 2026-10; vinecopulib 36 / 14; rvinecopulib 40 / 11 | one line |

**How to present it.** pyvinecopulib is a specialized library with a broad
user base:

- it is downloaded more than all other dedicated Python vine packages
  combined (torchvinecopulib, VineCopulas, mixedvines);
- its dependents span finance, atmospheric physics, bioinformatics, synthetic
  data and ML research;
- both its Python layer and its engine have external contributors.

Two core maintainers lead development; say so alongside the contributor list.

**Evidence that can be collected legitimately before submission:**

- Google Scholar and Semantic Scholar counts for the Zenodo DOIs and the
  repository URL.
- Twelve months of PyPI downloads from the BigQuery `pypi.file_downloads`
  table, split by version to show 1.x adoption.
- conda downloads per version.
- A verified dependents list, with one line of use per project.
- External issues opened since 1.0.0.
- Repository traffic from the GitHub API, which needs push access.
- Whatever the measures in §8 produce — Discussions, labels, a 1.1
  announcement. Opening those channels is fine; they must not be engineered
  to generate numbers.

---

## 13. Submission worklist (prioritized)

**1. Code and release**

- P0: Everything in §11, then tag v1.1.0 and confirm PyPI, Read the Docs
  `stable`, the GitHub release and the Zenodo DOI.
- P0: Turn Listing 1's full version into a tested notebook cell or test.
- P1: Dependency-floor leg (`uv sync --resolution lowest-direct` plus the
  extras on one Linux leg), or raise the floors to what is tested.
- P1: CUDA evidence per release: logs as release assets, or a self-hosted GPU
  runner.
- P1: conda-forge feedstock at 1.1.0 with osx-arm64 (and linux-aarch64).
- P2: Linux aarch64 wheels.
- Cross-repository, P0: fix JSS Table 6 and the bib title in
  rvinecopulib#348 (§9.3).

**2. Documentation and repository**

- P0: `docs/citing.rst` plus a "Citing" section in the README (Appendix C).
- P1: Align the author order of `CITATION.cff` (Vatter, Nagler) and
  `.zenodo.json` (Nagler, Vatter), which the Zenodo record follows, with the
  paper's; and replace `.zenodo.json`'s link to
  `vinecopulib.github.io/pyvinecopulib`, which returns 404, with the Read the
  Docs URL.
- P1: An architecture page (Figure 1, the four contracts, level-wise
  execution) as the human entry point ahead of AGENTS.md.
- P2: GitHub Discussions; "good first issue" labels.
- P2: A performance guidance page ("which lane when"), linking to the
  archived benchmarks.

**3. Measurements and data**

- P0: Tables 2a and 2b on the 1.1.0 wheels (§7).
- P0: The release's coverage figure, from the CI job (§11), plus the collected
  test count and the CUDA and `compile` counts from the CUDA rehearsal.
- P1: Re-check every Table 1 cell against each package's latest release, with
  the date.
- P1: A dated community snapshot (§12).
- P2: A pilot for the optional cross-package runtime comparison (§6).

**4. Manuscript** (starting from the skeleton this branch adds)

- P0: Figure 1 in TikZ, replacing its placeholder; Table 2 filled in; Listing 1
  executed.
- P1: Replace each plan box with prose to §3's page budget. Have an ML reader
  with no copula background check that the four pages stand alone. Then
  write the ≤ 3-page appendix.
- P1: Post the JSS manuscript to arXiv first; align the titles.

**5. Cover letter and community**

- P0: Draft the letter to §12's structure.
- P1: Shortlist action editors and reviewers: people with copula, density
  estimation or probabilistic-software background, and none who co-authored
  with the authors recently.
- P1: Collect the evidence listed in §12.

---

## Appendix A: MLOSS rules and precedents

**Rules** (https://www.jmlr.org/mloss/mloss-info.html, read 2026-10-29).

- **Length:** "An up to four page description based on the JMLR format, with
  additional pages for references."
- **Review criteria**, grouped under four headings:
  - *The four page summary*: its quality.
  - *Impact*:
    - "novelty, breadth, and significance";
    - "Evidence of an active user community (documented in the cover letter,
      not the paper)";
    - openness ("public source code repository, bug tracker, mailing
      list/forum");
    - "comparison to previous (if any) related implementations, w.r.t.
      run-time, memory requirements, features".
  - *User documentation*: installation instructions; tutorials; non-trivial
    examples; full API documentation, preferably as a website.
  - *Openness to new developers*: ease of contributing; clarity of design;
    developer documentation.
  - *Best practices*:
    - "Extensive unit and integration testing, and report of code coverage
      of tests. It is expected that test coverage is close to 100%."
    - "Continuous integration, ideally on all supported platforms and with
      multiple versions of potential dependencies."
    - Platform breadth, including an open-source OS.
    - "freedom of the code".
- **Submission checklist:** the license must be "clearly stated in the source
  package and mentioned in each source file"; the source package must contain
  no extra files.
- **Prior publication:** "prior publication of the method is admissible, as
  long as the software has not been published elsewhere. If the software has
  been the main content of a paper appearing in a peer reviewed conference or
  journal, then there should be a document in the code repository (referred to
  in the cover letter of the submission), listing the software package's
  improvements and extensions."
- **JMLR in general:**
  - disclose overlapping publications and attach them;
  - no simultaneous submission;
  - preprints are fine;
  - `jmlr2e.sty`, 11 pt, single column, `latexmk` / pdflatex;
  - a running title of ≤ 50 characters; five keywords; an abstract of
    ≤ 200 words;
  - after acceptance, a `*-code.tar.gz` source archive.

**Precedents** (JMLR volume 27).

| | 26-0434 NeuralOperator | 25-2424 OptunaHub | 25-2553 MarkDiffusion |
| --- | --- | --- | --- |
| Pages (main / total) | 4 / 6 | 4 / 10 | 4 / 22 |
| Figures / tables / listings in the 4 pages | 1 / 0 / 0 | 2 / 0 / 1 (14 lines) | 2 / 3 / 0 |
| Figure 1 | Five-column pipeline (Datasets → Transforms → Layers → Models → Meta-algos), full width | Actor/component diagram (Registry, Web, Module, Developer, User) | Three panels (implementation, visualization, evaluation) |
| Competitors | none | prose, plus a 1.75-page appendix of positioning | appendix table, 2 tools |
| Evidence | none | adoption count (> 100 packages) | reproduction table, "ours / official" |
| Quality claims | "covered by unit tests", CI/CD, docs site, MIT | review process | almost none |
| Appendix | none | positioning; package development | design, background, comparison, usage, algorithms, visuals |

Common shape:

- Introduction (~1 page) with the gap and "we present X".
- Full-width Figure 1 at the top of page 2.
- A design section with one paragraph per figure box.
- An optional evaluation.
- A one-paragraph conclusion ending on page 4.

None of the three compares run time or memory with other libraries.
torchsom (also volume 27) puts "90% test coverage" in its abstract.

**What makes them read as research rather than README text:**

- the code is framed as the reference implementation of a research
  program;
- a taxonomy of the field;
- explicitly defended design trade-offs;
- fidelity against the official implementations.

## Appendix B: verified repository facts (2026-10-29)

**Level-wise execution.** Checked in `core/_vinecop_levels.py` and
`core/vinecop_base.py`.

- `level_wiring` derives each tree level's input columns, and the outputs the
  next tree reads, from the structure alone, using the engine's
  `min_array` / `needed_hfunc*`.
- `inverse_waves` levels the inverse's static `(var, tree)` dependency graph
  by longest path. On a conditional vine, a cell also waits for its
  conditioning columns.
- `stack_of` returns a pair class's own `_stack_pairs` stack when every pair
  in the level is of that one class. Otherwise `LoopedPairs` calls each pair's
  dispatchers on exactly the input of the per-edge walk.
- `batched=` chooses per level. A level with covariates or a conditioning
  context takes the loop.
- Fitting is resolved by `_resolve_fit_level` (`vinecop_base.py:1885`), in
  this order:
  1. the caller's `fit_level`;
  2. `bicop_class.from_data_batched`, when the fit batches;
  3. `bicop_class.from_data` on a `ThreadPoolExecutor`, when
     `num_threads > 1`;
  4. edge by edge.

  A caller's `fit_edge` disables both of the middle options.

**Size.**

- `src/pyvinecopulib`: ~27.7k lines of Python. The binding is ~3.8k lines of
  C++.
- Tests: 44 modules, ~26.5k lines, 1,070 test functions before
  parametrization.
- The `device` fixture (`tests/conftest.py:371`) runs the torch tests on `cpu`
  and `cuda`; CUDA is skipped without a GPU.
- 4 `compile` tests, run only by the `build_sdist` job.
- Tolerances asserted in the parity suites:
  - `TorchTllBicop` ↔ `Bicop`: mostly `1e-12`;
  - batched ↔ per edge: `1e-13`–`1e-15`, with the inverse exact;
  - CPU ↔ CUDA: `1e-9`–`1e-12`.

**CI** (`.github/workflows/pypi.yml`).

| Job | What it covers |
| --- | --- |
| `lint` | ruff, codespell, bandit, version check |
| `build` (10 legs) | cibuildwheel on manylinux / musllinux x86-64, macOS x86-64 / arm64, Windows × cp311 / cp312-abi3; pytest **without extras**, so torch and sklearn tests skip |
| `eigen5` (2 legs) | Ubuntu and Windows against Eigen 5.0.1 |
| `install_and_unit_test` (12 legs) | Ubuntu, macOS arm64, Windows × Python 3.11–3.14, with every extra; full pytest plus the ten notebooks |
| `build_sdist` | editable install; `make check` (ruff, ty); `compile` tests; sdist build, install and tests; content assertions |
| `verify_docs_build` | Sphinx with `-W`, nitpicky |
| upload | the tag must be on `main`; version metadata must match the tag |

Not covered:

- **coverage reporting**: every leg passes `--no-cov`, so CI publishes no
  figure;
- **CUDA**: no GPU runner;
- **dependency floors**: numpy ≥ 2.0, scikit-learn ≥ 1.4, torch ≥ 2.2,
  scipy ≥ 1.16 and pandas ≥ 2.0 are never resolved at their minimum;
- macOS x86-64 with extras (there is no torch wheel);
- Linux aarch64.

**Distribution.**

- PyPI: 10 wheels plus the sdist.
- conda-forge: newest build **1.0.0**, on linux-64, osx-64 and win-64 only
  (checked through api.anaconda.org).
- Zenodo: the GitHub integration works; v1.0.1 was archived on 2026-10-06.
- `CITATION.cff` already reads 1.1.0.

**Documentation.**

- Pages: the README; `docs/concepts.rst` (1,323 lines); the autosummary API
  reference; the migration guide; `docs/releasing.md`; CONTRIBUTING; AGENTS.md;
  Code of Conduct; SECURITY; the CHANGELOG.
- Ten executed notebooks. Notebook 10 covers a custom pair, a learnable torch
  pair trained with Adam, hosting `Bicop` to reproduce `Vinecop`, discrete
  edges, conditioning, non-simplified vines, and a custom margin and
  distribution.
- Defects found:
  - notebook 09's stale controls sentence;
  - notebooks 08 and 09 fetch data over the network;
  - no "how to cite" anywhere (issue #128 asked for one);
  - the `torch/vinedist.py:202-204` comment;
  - the CONTRIBUTING.md typo.

**Performance evidence already in PR bodies** (to be re-measured, not
reused):

- #307: batched CUDA `pdf` / `rosenblatt` at 1.3×–11.9× over 16-thread C++.
- #316: the batched level fit at 1.07×–2.7× over per edge, with C++ winning at
  d = 20, n = 2000 on selection (0.75×).
- #354: discrete batching at 2–3.5×.
- #369: hosted `Bicop` at 0.94–1.06× of the code before it.

## Appendix C: the repository scope page (`docs/citing.rst`)

**Why a "Citing" page rather than `docs/publication-scope.md`.** It renders on
Read the Docs, the `-W` docs build checks it, it answers users (issue #128),
and it is the document MLOSS asks the cover letter to point to. Add a short
"Citing" section to the README that links to it.

Outline:

1. **How to cite.** The MLOSS paper once it is accepted; until then the Zenodo
   concept DOI plus the version DOI. The JSS paper for methods, the engine and
   its validation. kde1d and wdm where those components matter.
2. **What each publication documents.** One row each for JSS, JMLR MLOSS, the
   companion, and the GPU methods paper (if it is a separate paper). Also
   arXiv:2506.13318 and torchvinecopulib: the authors' earlier
   computational-graph formulation, which a first version of the PyTorch lane
   implemented before level-wise execution replaced it.
3. **What is shared.** The same vinecopulib / kde1d / wdm headers, used
   unmodified; the same numbers *when the engine version is the same*. A table
   of the engine version each release pins.
4. **What pyvinecopulib adds, by release.**
   - 1.0: the four contracts and bases, `Vinedist` and `margins`, `sklearn`,
     `torch`, and the conditional / non-simplified extension points.
   - 1.0.1: discrete and mixed vines on the batched torch path; faster batched
     fits.
   - 1.1: one level-wise cascade for every vine; stacking chosen by the pair
     class; threaded level fits for compiled pairs.
5. **Why the papers do not overlap.** One paragraph restating §9.1.
