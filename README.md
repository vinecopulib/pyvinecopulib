# pyvinecopulib

[![Documentation](https://readthedocs.org/projects/pyvinecopulib/badge/?version=latest)](https://pyvinecopulib.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/vinecopulib/pyvinecopulib/actions/workflows/pypi.yml/badge.svg)](https://github.com/vinecopulib/pyvinecopulib/actions/workflows/pypi.yml)
[![DOI](https://zenodo.org/badge/196999069.svg)](https://zenodo.org/badge/latestdoi/196999069)

## Introduction

### What are vine copulas?

Vine copulas are a flexible class of dependence models consisting of bivariate
building blocks (see e.g.,
[Aas et al., 2009](https://mediatum.ub.tum.de/doc/1083600/1083600.pdf)).
You can find a comprehensive list of publications and other materials on
[vine-copula.org](http://vine-copula.org).

### What is pyvinecopulib?

[pyvinecopulib](https://vinecopulib.github.io/pyvinecopulib/) is the python interface to vinecopulib, a header-only C++ library for vine copula models based on
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). It provides
high-performance implementations of the core features of the popular
[VineCopula R library](https://github.com/tnagler/VineCopula), in particular
inference algorithms for both vine copula and bivariate copula models.
Advantages over VineCopula are

* a stand-alone C++ library with interfaces to both R and Python,
* a sleaker and more modern API,
* shorter runtimes and lower memory consumption, especially in high dimensions,
* nonparametric and multi-parameter families.

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

Start by cloning this repository, noting the `--recursive` option which is needed for the `vinecopulib` and `wdm` submodules:

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

See [CONTRIBUTING.md](CONTRIBUTING.md) for the full developer workflow.

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
make docs-serve     # live-reload dev server (sphinx-autobuild)
```

## Contributing

Development setup, the build pipeline, the Makefile + pre-commit conventions,
the CI workflow, and the release flow are all documented in
[CONTRIBUTING.md](CONTRIBUTING.md).
