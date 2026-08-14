# pyvinecopulib

[![Documentation](https://readthedocs.org/projects/pyvinecopulib/badge/?version=latest)](https://pyvinecopulib.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/vinecopulib/pyvinecopulib/actions/workflows/pypi.yml/badge.svg)](https://github.com/vinecopulib/pyvinecopulib/actions/workflows/pypi.yml)
[![DOI](https://zenodo.org/badge/196999069.svg)](https://zenodo.org/badge/latestdoi/196999069)

## Introduction

### What are vine copulas?

Sklar's theorem factorises every joint distribution into one-dimensional
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

[pyvinecopulib](https://vinecopulib.github.io/pyvinecopulib/) is the
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

### Optional backends

Two opt-in subpackages extend the core library:

* `pyvinecopulib.sklearn` — scikit-learn-compatible estimators
  (`VineDensity`, `VineRegressor`). Drop a vine
  into any sklearn pipeline:

  ```python
  from pyvinecopulib.sklearn import VineDensity
  density = VineDensity().fit(X)             # default backend (C++)
  density.score_samples(X[:3]); density.cdf(X[:3])
  ```

  Install with `pip install pyvinecopulib[sklearn]`.

* `pyvinecopulib.torch` — pure-PyTorch evaluators (`TorchBicop`,
  `TorchVinecop`) for GPU placement and autograd:

  ```python
  from pyvinecopulib.sklearn import VineDensity
  from pyvinecopulib.sklearn.backends import TorchVinecopBackend
  from pyvinecopulib.torch import FitControlsTorchVinecop
  controls = FitControlsTorchVinecop(device="cuda")
  density_gpu = VineDensity(backend=TorchVinecopBackend(controls=controls)).fit(X)
  ```

  Install with `pip install pyvinecopulib[torch]`.

### Custom and conditional pair copulas

The core evaluators (`Bicop` / `Vinecop`, and their torch counterparts)
implement two backend-neutral contracts, `BicopLike` and `VinecopLike`.
Subclass the canonical, pure-Python `BicopBase` / `VinecopBase` (NumPy or
PyTorch) to plug your **own** pair copula into a vine — including
**non-simplified** vines where each pair copula conditions on its
conditioning set. See the
[concepts page](https://pyvinecopulib.readthedocs.io/en/latest/concepts.html)
and the `examples/10_extending_pyvinecopulib.ipynb` notebook.

### Conditional sampling and likelihood diagnostics

A fitted `Vinecop` can draw from the conditional distribution of a subset of
variables given the rest (`simulate_conditional`), select or `reorient` a
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
[open an issue](https://github.com/pyvinecopulib/pyvinecopulib/issues/new) or
send a mail to <info@vinecopulib.org>.

## Installation

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

See the [contributing guide](https://pyvinecopulib.readthedocs.io/en/latest/CONTRIBUTING.html)
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
[contributing guide](https://pyvinecopulib.readthedocs.io/en/latest/CONTRIBUTING.html).
