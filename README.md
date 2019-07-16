# pyvinecopulib

[![Build Status](https://travis-ci.org/vinecopulib/pyvinecopulib.svg?branch=master)](https://travis-ci.org/vinecopulib/pyvinecopulib)
[![Build status](https://ci.appveyor.com/api/projects/status/2fn1v67sxdrmp2po/branch/master?svg=true)](https://ci.appveyor.com/project/vinecopulib/pyvinecopulib-x6s0i/branch/master)

## Prerequisites

**On Unix (Linux, OS X)**

* A compiler with C++11 support
* CMake >= 2.8.12

**On Windows**

* Visual Studio 2015 (required for all Python versions, see notes below)
* CMake >= 3.1

## Installation

Just clone this repository and pip install. Note the `--recursive` option which is
needed for the pybind11 submodule:

```bash
git clone --recursive https://github.com/pybind/pyvinecopulib.git
pip install ./pyvinecopulib
```

With the `setup.py` file included in this example, the `pip install` command will
invoke CMake and build the pybind11 module as specified in `CMakeLists.txt`.


## Special notes for Windows

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


## Building the documentation

Documentation for the example project is generated using Sphinx. Sphinx has the
ability to automatically inspect the signatures and documentation strings in
the extension module to generate beautiful documentation in a variety formats.
The following command generates HTML-based reference documentation; for other
formats please refer to the Sphinx manual:

 - `cd pyvinecopulib/docs`
 - `make html`


## License

pyvinecopulib is provided under an MIT license that can be found in the LICENSE
file. By using, distributing, or contributing to this project, you agree to the
terms and conditions of this license.


## Test call

```python
import pyvinecopulib
pyvinecopulib.add(1, 2)
```

