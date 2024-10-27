# pyvinecopulib

[![Build Status](https://github.com/vinecopulib/pyvinecopulib/workflows/Build%20Status/badge.svg?branch=main)](https://github.com/vinecopulib/pyvinecopulib/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/3c0056d3ca5244a5ba6a2b32f87be4cf)](https://www.codacy.com/gh/vinecopulib/pyvinecopulib?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=vinecopulib/pyvinecopulib&amp;utm_campaign=Badge_Grade)
[![Documentation](https://img.shields.io/website/http/vinecopulib.github.io/pyvinecopulib.svg)](https://vinecopulib.github.io/pyvinecopulib/)
[![DOI](https://zenodo.org/badge/196999069.svg)](https://zenodo.org/badge/latestdoi/196999069)

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

### Prerequisites

* numpy (>=1.14)
* matplotlib (>=3.0),
* networkx (>=3.0),
* pydot (>=3.0),
* To install from source:
    * pybind11 (>=2.4)
    * setuptools (>=30.3.0)
    * setuptools_scm (>=2.0.0)
    * a compiler with C++11 support (Linux, OS X) or Visual Studio 2015 (required for all Python versions, see notes below)
    * Eigen (the environment variable `EIGEN3_INCLUDE_DIR` must be set to the directory containing the Eigen headers)
    * boost (the environment variable `Boost_INCLUDE_DIR` must be set to the directory containing the boost headers)

### Installation

#### With pip

The latest release can be installed using `pip`:

```
pip install pyvinecopulib
```

#### With conda

Installing `pyvinecopulib` from the `conda-forge` channel can be achieved by adding `conda-forge` to your channels with:

```
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Once the `conda-forge` channel has been enabled, `pyvinecopulib` can be installed with `conda`:

```
conda install pyvinecopulib
```

or with `mamba`:

```
mamba install pyvinecopulib
```

It is possible to list all of the versions of `pyvinecopulib` available on your platform with `conda`:

```
conda search pyvinecopulib --channel conda-forge
```

or with `mamba`:

```
mamba search pyvinecopulib --channel conda-forge
```

Alternatively, `mamba repoquery` may provide more information:

```
# Search all versions available on your platform:
mamba repoquery search pyvinecopulib --channel conda-forge

# List packages depending on `pyvinecopulib`:
mamba repoquery whoneeds pyvinecopulib --channel conda-forge

# List dependencies of `pyvinecopulib`:
mamba repoquery depends pyvinecopulib --channel conda-forge
```

#### From source

To install from source, Eigen and Boost need to be available on your system for the build to succeed, using the environment variables `EIGEN3_INCLUDE_DIR` and `Boost_INCLUDE_DIR` respectively.
On Linux, you can install the required packages and set the environment variables as follows:

```bash
sudo apt-get install libeigen3-dev libboost-all-dev
export Boost_INCLUDE_DIR=/usr/include
export EIGEN3_INCLUDE_DIR=/usr/include/eigen3
```

Then, just clone this repository and do `pip install`.
Note the `--recursive` option which is needed for the `vinecopulib` and `wdm` submodules:

```bash
git clone --recursive https://github.com/vinecopulib/pyvinecopulib.git
pip install ./pyvinecopulib
```

If the required dependencies are not installed, a reproducible environment, which also include stuff requirement for the library's development and documentation, can be created using:

```bash
mamba create -n pyvinecopulib numpy mypy ruff pytest sphinx-rtd-theme sphinx-autodoc-typehints pydot networkx matplotlib pybind11 setuptools-scm python=3.11
mamba activate pyvinecopulib
```

### Examples

Jupyter notebooks with examples can be found in the examples folder.

### Documentation

For documentation of the `pyvinecopulib`'s functionality and
instructions how to use it, check out our
[website](https://vinecopulib.github.io/pyvinecopulib/) or the `docs/` folder
in this repository.

#### Building the documentation

Documentation for the example project is generated using Sphinx and the "Read the Docs" theme.
The following command generates HTML-based reference documentation; for other
formats please refer to the Sphinx manual:

* `pip install sphinx-rtd-theme sphinx-autodoc-typehints`
* `cd pyvinecopulib/docs`
* `python serve_sphinx.py`

### License

pyvinecopulib is provided under an MIT license that can be found in the LICENSE
file. By using, distributing, or contributing to this project, you agree to the
terms and conditions of this license.

### Special notes for Windows

**Compiler requirements**

This package requires a C++11 compliant compiler, i.e Visual Studio 2015 on Windows.
Unlike regular C extension modules, it's perfectly fine to compile a pyvinecopulib module with a VS version newer than the target Python's VS version.

**Runtime requirements**

The Visual C++ 2015 redistributable packages are a runtime requirement for this
project.

### Contact

If you have any questions regarding the library, feel free to
[open an issue](https://github.com/pyvinecopulib/pyvinecopulib/issues/new) or
send a mail to <info@vinecopulib.org>.
