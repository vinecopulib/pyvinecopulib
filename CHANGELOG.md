# Changelog

## 0.7.1

## Bug fixes in `pyvinecopulib`

* Upgrade nanobind to allow for single row matrices (fix #169 and #170)

## New features in `pyvinecopulib`

* Add pickle support for all classes (#168)
* Add `allow_rotation` option to `FitControlsBicop` and `FitControlsVinecop` (#168)

## Changes in `vinecopulib`

These changes originate from the underlying C++ library, [`vinecopulib`](https://github.com/vinecopulib/vinecopulib), which powers `pyvinecopulib`.

### NEW FEATURES

* add `allow_rotation` option to `FitControlsBicop` and `FitControlsVinecop`
  to allow for the rotation of the pair copulas ([#628](https://github.com/vinecopulib/vinecopulib/pull/628)).
* add a `FitControlsConfig` struct to create flexible and yet safe constructors
  for `FitControlsBicop` and `FitControlsVinecop` ([#629](https://github.com/vinecopulib/vinecopulib/pull/629)).

### BUG FIXES

* restrict parameter range for fitting Tawn copulas; fix handling of their
  shape/argument order ([#620](https://github.com/vinecopulib/vinecopulib/pull/620)).
* compute and save loglik/nobs in `Vinecop::fit()` ([#623](https://github.com/vinecopulib/vinecopulib/pull/623))
* disable unwanted compiler output related to BOOST_CONCEPT checks ([#624](https://github.com/vinecopulib/vinecopulib/pull/624))

## 0.7.0

This version introduces a switch to nanobind as a backend (#160): i.e., the C++ bindings, now use [nanobind](https://nanobind.readthedocs.io/) instead of [pybind11](https://pybind11.readthedocs.io/).
It allows for considerable performance improvements (~8x speedup in our latest benchmarks) and smaller binaries.

### Breaking API changes in `pyvinecopulib`

* Removal of the overloaded constructors:
    * For all classes, only one constructor is now available.
    The reason is that the overloaded constructors were un-Pythonic, error-prone, and could not be properly documented with Sphinx.
    They have been replaced by a single constructor for each class, along with factory `from_xzy` methods.
    * For the ``Bicop`` class:
        * ``Bicop.from_family()``: Instantiate from a family, rotation, parameters, and variable types.
        * ``Bicop.from_data()``: Instantiate from data, as well as optional controls and variable types.
        * ``Bicop.from_file()``: Instantiate from a file.
        * ``Bicop.from_json()``: Instantiate from a JSON-like string.
    * For the ``Vinecop`` class:
        * ``Vinecop.from_dimension()``: Instantiate an empty vine copula of a given dimension.
        * ``Vinecop.from_data()``: Instantiate from data, as well as an optional ``FitControlsVinecop``, an ``RVineStructure`` or matrix, and variable types.
        * ``Vinecop.from_structure()``: Instantiate from an ``RVineStructure`` or matrix, as well as optional pair-copulas and variable types.
        * ``Vinecop.from_file()``: Instantiate from a file.
        * ``Vinecop.from_json()``: Instantiate from a JSON-like string.
    * For the ``RVineStructure`` class:
        * ``RVineStructure.from_dimension()``: Instantiate a default structure of a given dimension and truncation level.
        * ``RVineStructure.from_order()``: Instantiate from an order vector.
        * ``RVineStructure.from_matrix()``: Instantiate from a matrix.
        * ``RVineStructure.from_file()``: Instantiate from a file.
        * ``RVineStructure.from_json()``: Instantiate from a JSON-like string.

### New features in `pyvinecopulib`

* Expose more structure methods to python (#157)
* Switch to nanobind as a backend (#160)
* New IO methods for `Bicop` and `Vinecop` classes to use JSON-like strings (#160)
* Extensive documentation revamp (#160)
* Adding a benchmark example (#160)
* Convertion of all examples to Jupyter notebooks (#160)

### Bug fixes in `pyvinecopulib`

* Install and test source distribution (#164)

### Changes in `vinecopulib`

These changes originate from the underlying C++ library, [`vinecopulib`](https://github.com/vinecopulib/vinecopulib), which powers `pyvinecopulib`.

#### New features

* Use analytical derivatives in discrete pdf/hfuncs ([#572](https://github.com/vinecopulib/vinecopulib/pull/572))
* Allow for alternative for `"prim"` vs `"kruskal"` in MST-based model selection ([#577](https://github.com/vinecopulib/vinecopulib/pull/577))
* Improve the dependencies install script to use it in other projects ([#576](https://github.com/vinecopulib/vinecopulib/pull/576))
* Add tawn copula ([#579](https://github.com/vinecopulib/vinecopulib/pull/579))
* Improve doc ([#580](https://github.com/vinecopulib/vinecopulib/pull/580), [#585](https://github.com/vinecopulib/vinecopulib/pull/585), [#607](https://github.com/vinecopulib/vinecopulib/pull/607))
* Allow for the discrete Rosenblatt transform ([#581](https://github.com/vinecopulib/vinecopulib/pull/581))
* Add `Vinecop::fit()` ([#584](https://github.com/vinecopulib/vinecopulib/pull/584))
* Improve `Bicop::str()` ([#588](https://github.com/vinecopulib/vinecopulib/pull/588)) and `Vinecop::str()` ([#589](https://github.com/vinecopulib/vinecopulib/pull/589))
* Properly handle discrete variables for the TLL family ([#597](https://github.com/vinecopulib/vinecopulib/pull/597))
* Weighted pseudo-observations ([#602](https://github.com/vinecopulib/vinecopulib/pull/602))
* Cross-platform random numbers and add seeds options to `to_pseudo_obs` ([#603](https://github.com/vinecopulib/vinecopulib/pull/603))
* Improve performance by
    * aligning with the `R` defaults (e.g., `BOOST_NO_AUTO_PTR`, `BOOST_ALLOW_DEPRECATED_HEADERS`, `BOOST_MATH_PROMOTE_DOUBLE_POLICY=false`, `std::string nonparametric_method = "constant"` for the TLL instead of `"quadratic"`, `-O3 -march=native` compiler flags) and add benchmarking example ([#592](https://github.com/vinecopulib/vinecopulib/pull/592), [#611](https://github.com/vinecopulib/vinecopulib/pull/611), [#613](https://github.com/vinecopulib/vinecopulib/pull/613)),
    * using `Eigen` element-wise operations instead of `boost` whenever possible ([#598](https://github.com/vinecopulib/vinecopulib/pull/598), [#612](https://github.com/vinecopulib/vinecopulib/pull/612)),
    * using binary search in the TLL for `get_indices` ([#613](https://github.com/vinecopulib/vinecopulib/pull/613)).

#### Bug fixes

* Improve stability in BB7 PDF ([#573](https://github.com/vinecopulib/vinecopulib/pull/573))
* Revamped CI/CD pipeline, tests discoverable by CTest, boost version on windows (([66cf8b0](https://github.com/vinecopulib/vinecopulib/commit/66cf8b0)))
* Fix ASAN issues ([#583](https://github.com/vinecopulib/vinecopulib/pull/583))
* Fix interface includes and other CMake issue ([#586](https://github.com/vinecopulib/vinecopulib/pull/586), [#599](https://github.com/vinecopulib/vinecopulib/pull/599), [#601](https://github.com/vinecopulib/vinecopulib/pull/601), [#608](https://github.com/vinecopulib/vinecopulib/pull/608)), by @jschueller
