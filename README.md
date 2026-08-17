# pyvinecopulib

[![Documentation](https://readthedocs.org/projects/pyvinecopulib/badge/?version=latest)](https://pyvinecopulib.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/vinecopulib/pyvinecopulib/actions/workflows/pypi.yml/badge.svg)](https://github.com/vinecopulib/pyvinecopulib/actions/workflows/pypi.yml)
[![DOI](https://zenodo.org/badge/196999069.svg)](https://zenodo.org/badge/latestdoi/196999069)

## Introduction

### What are vine copulas?

Sklar's theorem factorizes every joint distribution into one-dimensional
**marginals** and a **copula** that carries the dependence between
variables. A vine copula decomposes that copula into bivariate building
blocks — **pair copulas** — arranged on a sequence of trees called an
**R-vine** (Bedford & Cooke, 2002;
[Aas et al., 2009](https://mediatum.ub.tum.de/doc/1083600/1083600.pdf)).
The decomposition makes pair-by-pair estimation scale gracefully into
high dimensions and gives a natural place to drop in non-parametric
pair-copula estimators like Transformed Local Likelihood (TLL).

A short primer is available on the
[concepts page](https://pyvinecopulib.readthedocs.io/en/latest/concepts.html);
a comprehensive list of publications lives on
[vine-copula.org](http://vine-copula.org).

### What is pyvinecopulib?

[pyvinecopulib](https://pyvinecopulib.readthedocs.io) is the
Python interface to vinecopulib, a header-only C++ library for vine
copula models based on
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). It
provides high-performance implementations of the core features of the
popular [VineCopula R library](https://github.com/tnagler/VineCopula),
in particular inference algorithms for both vine copula and bivariate
copula models. Advantages over VineCopula are

* a stand-alone C++ library with interfaces to both R and Python,
* a sleeker and more modern API,
* shorter runtimes and lower memory consumption, especially in high
  dimensions,
* nonparametric and multi-parameter families.

### First core fit

The core package is enough to fit, inspect, evaluate, and sample a model:

```python
import numpy as np
import pyvinecopulib as pv

rng = np.random.default_rng(0)
cov = [[1.0, 0.7, 0.3], [0.7, 1.0, 0.5], [0.3, 0.5, 1.0]]
x = rng.multivariate_normal([0, 0, 0], cov, size=500)

u = pv.to_pseudo_obs(x)              # ranks, on the copula scale
vine = pv.Vinecop.from_data(u)       # selects structure and families
print(vine)                          # the fitted trees, pair by pair
vine.loglik(u), vine.bic()           # fit diagnostics
draws = vine.sample(100, seeds=[1])  # new copula-scale observations
```

For a distribution on the original data scale, pair the fitted copula with
one margin per variable through `pv.Vinedist`. Notebooks 03, 07, and 11 build
out these core workflows.

### Optional subpackages

Three opt-in subpackages extend the core library:

* `pyvinecopulib.margins` — parametric margins and family selection
  (`ParametricMargin`, `MarginSelector`) to pair with `Vinedist` when a
  kernel-density margin is not what you want:

  ```python
  from pyvinecopulib.core import Vinedist
  from pyvinecopulib.margins import MarginSelector
  dist = Vinedist.from_data(x, margins=MarginSelector())
  print(dist.margins[0].selected_)
  ```

  Install with `pip install pyvinecopulib[scipy]` (or `[openturns]`).

* `pyvinecopulib.sklearn` — scikit-learn-compatible estimators
  (`VineDensity`, `VineRegressor`). Drop a vine
  into any sklearn pipeline:

  ```python
  from pyvinecopulib.sklearn import VineDensity
  density = VineDensity().fit(X)             # default backend (C++)
  density.score_samples(X[:3]); density.cdf(X[:3])
  ```

  Install with `pip install pyvinecopulib[sklearn]`.

* `pyvinecopulib.torch` — pure-PyTorch evaluators and data-scale modules
  (`TorchBicop`, `TorchVinecop`, `TorchKde1d`, `TorchMargin`, and
  `TorchVinedist`) for GPU placement and autograd:

  ```python
  import torch
  from pyvinecopulib.torch import TorchVinedist, FitControlsTorchVinecop
  dist = TorchVinedist.from_data(
    torch.as_tensor(x), controls=FitControlsTorchVinecop(device="cuda")
  )
  y = torch.as_tensor(x[:5], device="cuda").requires_grad_(True)
  dist.log_prob(y).sum().backward()   # autograd through the whole vine
  print(y.grad)                       # d log f / dy, on the GPU
  ```

  The same evaluator backs the sklearn estimators through
  `pyvinecopulib.sklearn.backends.TorchVinecopBackend`.

  Install with `pip install pyvinecopulib[torch]`.

### API stability

`pyvinecopulib.core`, `pyvinecopulib.families` and `pyvinecopulib.utils` are
stable: changes there follow semantic versioning, with a deprecation cycle
before anything is removed.

`pyvinecopulib.margins`, `pyvinecopulib.sklearn` and `pyvinecopulib.torch` are
**provisional in 1.x**. Their contracts are new in this release and may still
change in a minor version as they meet real data -- the margin contract
(`MarginLike` / `MarginBase`) and the torch/C++ evaluation parity are the parts
already treated as load-bearing. Pin an exact version if you depend on their
surface.

### Custom and conditional pair copulas

The core evaluators (`Bicop` / `Vinecop`, and their torch counterparts)
implement two backend-neutral contracts, `BicopLike` and `VinecopLike`.
Subclass the canonical, pure-Python `BicopBase` / `VinecopBase` (NumPy or
PyTorch) to plug your **own** pair copula into a vine. A pair may depend on its
vine conditioning-set values (a **non-simplified** vine), on row-aligned
external covariates, or on both. `Vinedist` can compose covariate-dependent
margins and such a copula into a full data-scale distribution `Y | X`.

This joint conditional model is an extension seam, not a built-in fitter:
`Vinedist.from_data(y, x=...)` can fit custom conditional margin
specifications, but fits an `x`-independent compiled `Vinecop` for the copula
half. Fit custom conditional pairs through `VinecopBase.fit` and compose the
parts explicitly when dependence must also vary with `X`. See the
[concepts page](https://pyvinecopulib.readthedocs.io/en/latest/concepts.html)
and notebooks `examples/10_extending_pyvinecopulib.ipynb` and
`examples/11_vine_distributions.ipynb`.

### Conditional sampling and likelihood diagnostics

A fitted `Vinecop` can draw from the conditional distribution of a subset of
variables given the rest (`sample_conditional`), select or `reorient` a
structure so a chosen conditioning set sits at the order tail, and expose its
tree-by-tree decomposition with the fitted pair copulas (`get_trees`).
Parametric fits additionally provide analytic log-likelihood scores, gradient,
Hessian, and score covariance (`scores` / `gradient` / `hessian` /
`scores_cov`, on both `Bicop` and `Vinecop`) for gradient-based inference. See
the `examples/05_conditional_sampling_and_vines.ipynb` notebook.

### License

pyvinecopulib is provided under an MIT license that can be found in the LICENSE
file. By using, distributing, or contributing to this project, you agree to the
terms and conditions of this license.

### Contact

If you have any questions regarding the library, feel free to
[open an issue](https://github.com/vinecopulib/pyvinecopulib/issues/new) or
send a mail to <info@vinecopulib.org>.

## Installation

On x86-64, the distributed wheels require the x86-64-v3 ISA baseline (AVX2 and
FMA). The package checks this before loading its native extension and explains
how to use a source build when a CPU or VM masks those features.

### With pip

The latest release can be installed using `pip`:

```bash
pip install pyvinecopulib
```

### With conda

Similarly, it can be installed with `conda`:

```bash
conda install conda-forge::pyvinecopulib
```

Or with `mamba`:

```bash
mamba install conda-forge::pyvinecopulib
```

### From source

Start by cloning this repository, noting the `--recursive` option which is needed for the `vinecopulib`, `wdm`, and `kde1d` submodules:

```bash
git clone --recursive https://github.com/vinecopulib/pyvinecopulib.git
cd pyvinecopulib
```

The main build time prerequisites are:

* scikit-build-core (>=0.5.0),
* nanobind (>=2.7.0),
* libclang (>=18) — used to regenerate `src/include/docstr.hpp` from the C++ headers as a step of the build,
* numpy / matplotlib / networkx — imported by the post-build stub-generation step,
* a compiler with C++17 support.

When installing via `pip install .` (the default), all of these are pulled into an isolated build environment automatically via `[build-system] requires` in `pyproject.toml`; you don't need to install them yourself.

To install from source, `Eigen` and `Boost` also need to be available, and CMake will try to find suitable versions automatically.

The recommended way to install `pyvinecopulib` from source is to use `conda`/`mamba` for the native build prerequisites and [`uv`](https://docs.astral.sh/uv/) for the Python side:

```bash
mamba create -n pyvinecopulib python=3.11 boost eigen 'python-clang=18.*' uv
mamba activate pyvinecopulib
make sync
```

See the [contributing guide](https://github.com/vinecopulib/pyvinecopulib/blob/main/CONTRIBUTING.md)
for the full developer workflow.

Alternatively, you can specify manually the location of `Eigen` and `Boost` using the environment variables `EIGEN3_INCLUDE_DIR` and `Boost_INCLUDE_DIR` respectively.
On Linux, you can install the required packages and set the environment variables as follows:

```bash
sudo apt-get install libeigen3-dev libboost-all-dev
export Boost_INCLUDE_DIR=/usr/include
export EIGEN3_INCLUDE_DIR=/usr/include/eigen3
```

Finally, you can build and install `pyvinecopulib` using `pip`:

```bash
pip install .
```

The build automatically regenerates `src/include/docstr.hpp` (from the C++ headers via libclang) and `src/pyvinecopulib/__init__.pyi` (from the freshly built extension). Both files are gitignored — they're pure build artifacts.

For an editable install (recommended for development), use `--no-build-isolation` so the conda env's libclang is reused and `editable.rebuild = true` regenerates everything on each `import`:

```bash
pip install -e . --no-build-isolation
```

## Documentation

Stable docs are published at <https://pyvinecopulib.readthedocs.io>. They are
rebuilt automatically by Read the Docs whenever a new release is tagged on
`main` and published to PyPI.

To build the documentation locally:

```bash
make docs           # one-shot HTML build → docs/_build/html/
```

## Contributing

Development setup, the build pipeline, the Makefile + pre-commit conventions,
the CI workflow, and the release flow are all documented in the
[contributing guide](https://github.com/vinecopulib/pyvinecopulib/blob/main/CONTRIBUTING.md).
