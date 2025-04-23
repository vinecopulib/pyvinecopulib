# pyvinecopulib

[![Documentation](https://img.shields.io/website/http/vinecopulib.github.io/pyvinecopulib.svg)](https://vinecopulib.github.io/pyvinecopulib/)
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

* scikit-build-core (>=0.4.3),
* nanobind (>=2.7.0),
* a compiler with C++17 support.

To install from source, `Eigen` and `Boost` also need to be available, and CMake will try to find suitable versions automatically.

The recommended way to install `pyvinecopulib` from source is to use `conda` or `mamba`.
A reproducible environment, also including requirements for the `pyvinecopulib`'s development and documentation, can be created using:

```bash
python scripts/generate_requirements.py --format yml # from pyvinecopulib's root
mamba env create -f environment.yml
mamba activate pyvinecopulib
```

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

Note that the `generate_requirements.py` script can also be used to generate a `requirements.txt` file for use with `pip` via the `--format` option:

```bash
python scripts/generate_requirements.py --format txt
```

### Building the documentation

Documentation for the example project is generated using Sphinx and the "Read the Docs" theme.
The following command generates HTML-based reference documentation; for other
formats please refer to the Sphinx manual:

```bash
cd docs
python serve_sphinx.py
```
