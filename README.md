# pyvinecopulib

[![Build Status](https://travis-ci.org/vinecopulib/pyvinecopulib.svg?branch=master)](https://travis-ci.org/vinecopulib/pyvinecopulib)
[![Build status](https://ci.appveyor.com/api/projects/status/2fn1v67sxdrmp2po/branch/master?svg=true)](https://ci.appveyor.com/project/vinecopulib/pyvinecopulib-x6s0i/branch/master)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

### What are vine copulas?

Vine copulas are a flexible class of dependence models consisting of bivariate
building blocks (see e.g.,
[Aas et al., 2009](https://mediatum.ub.tum.de/doc/1083600/1083600.pdf)).
You can find a comprehensive list of publications and other materials on
[vine-copula.org](http://www.statistics.ma.tum.de/en/research/vine-copula-models/).

### What is pyvinecopulib?

pyvinecopulib is the python interface to vinecopulib, a header-only C++ library for vine copula models based on
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

* NumPy
* To install from source: a compiler with C++11 support (Linux, OS X) or Visual Studio 2015 (required for all Python versions, see notes below)

### Installation

The easiest way to install the latest release is to use `pip`:

```
pip install pyvinecopulib
```

Alternatively, `pyvinecopulib` is also available via `conda` (linux and macos only):

```
conda install -c vinecopulib pyvinecopulib
```

To install from source, just clone this repository and do `pip install`.
Note the `--recursive` option which is needed for the `pybind11`, `vinecopulib` and `wdm` submodules:

```bash
git clone --recursive https://github.com/vinecopulib/pyvinecopulib.git
pip install ./pyvinecopulib
```

### Examples 

Jupyter notebooks with examples can be found in the examples folder.

### Building the documentation

Documentation for the example project is generated using Sphinx.
The following command generates HTML-based reference documentation; for other
formats please refer to the Sphinx manual:

 - `cd pyvinecopulib/docs`
 - `make html`

### License

pyvinecopulib is provided under an MIT license that can be found in the LICENSE
file. By using, distributing, or contributing to this project, you agree to the
terms and conditions of this license.

### Special notes for Windows

**Compiler requirements**

Pyvinecopulib requires a C++11 compliant compiler, i.e Visual Studio 2015 on Windows.
This applies to all Python versions, including 2.7. Unlike regular C extension
modules, it's perfectly fine to compile a pyvinecopulib module with a VS version newer
than the target Python's VS version.

**Runtime requirements**

The Visual C++ 2015 redistributable packages are a runtime requirement for this
project. If you use the Anaconda Python
distribution, you can add `vs2015_runtime` as a platform-dependent runtime
requirement for you package: see the `conda.recipe/meta.yaml` file in this example.


### Contact

If you have any questions regarding the library, feel free to
[open an issue](https://github.com/pyvinecopulib/pyvinecopulib/issues/new) or
send a mail to <info@vinecopulib.org>.

