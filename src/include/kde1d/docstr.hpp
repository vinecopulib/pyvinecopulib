#pragma once

//! Custom docstrings for kde1d Python bindings
//!
//! This file contains manually written docstrings for the kde1d bindings,
//! based on the R interface documentation (kde1d.R and kde1d-methods.R).
//! These docstrings follow the same style as those in the main docstr.hpp
//! file but are adapted for the Python/C++ interface differences.

// Docstrings for Kde1d class bindings
// Based on the R interface documentation but adapted for Python/C++
// functionality

namespace kde1d_docstrings {

constexpr const char* kde1d_class_doc = R"""(
A class for univariate kernel density estimation.

The ``Kde1d`` class provides methods for univariate kernel density estimation
using local polynomial fitting. It can handle data with bounded, unbounded,
and discrete support.

The estimator uses a Gaussian kernel in all cases. A log-transform is used if
there is only one boundary; a probit transform is used if there are two
boundaries. Discrete variables are handled via jittering.

Zero-inflated densities are estimated by a hurdle-model with discrete mass at 0
and the remainder estimated as for continuous data.

Examples
--------
>>> import numpy as np
>>> import pyvinecopulib as pv
>>>
>>> # Unbounded data
>>> x = np.random.normal(0, 1, 500)
>>> fit = pv.Kde1d()
>>> fit.fit(x)
>>> pdf_vals = fit.pdf(np.array([0.0]))
>>> fit.plot(x)
>>>
>>> # Bounded data
>>> x = np.random.gamma(1, size=500)
>>> fit = pv.Kde1d(xmin=0.0, degree=1)
>>> fit.fit(x)
>>> fit.plot(x)
>>>
>>> # Discrete data
>>> x = np.random.binomial(5, 0.5, 500)
>>> fit = pv.Kde1d(xmin=0, xmax=5, type="discrete")
>>> fit.fit(x)
>>> fit.plot(x)

References
----------
Geenens, G. (2014). *Probit transformation for kernel density estimation on
the unit interval.* Journal of the American Statistical Association,
109(505), 346–358.
[arXiv:1303.4121](https://arxiv.org/abs/1303.4121)

Geenens, G., & Wang, C. (2018). *Local-likelihood transformation kernel
density estimation for positive random variables.* Journal of Computational
and Graphical Statistics, 27(4), 822–835.
[arXiv:1602.04862](https://arxiv.org/abs/1602.04862)

Loader, C. (2006). *Local Regression and Likelihood.* Springer Science &
Business Media.

Nagler, T. (2018a). *A generic approach to nonparametric function estimation
with mixed data.* Statistics & Probability Letters, 137, 326–330.
[arXiv:1704.07457](https://arxiv.org/abs/1704.07457)

Nagler, T. (2018b). *Asymptotic analysis of the jittering kernel density
estimator.* Mathematical Methods of Statistics, 27, 32–46.
[arXiv:1705.05431](https://arxiv.org/abs/1705.05431)
)""";

constexpr const char* kde1d_constructor_doc = R"""(
Constructor for the ``Kde1d`` class.

Parameters
----------
xmin : float, optional
    Lower bound for the support of the density. ``NaN`` means no boundary.
    Default is ``NaN``.
xmax : float, optional
    Upper bound for the support of the density. ``NaN`` means no boundary.
    Default is ``NaN``.
type : str, optional
    Variable type. Must be one of ``"continuous"``, ``"discrete"``, or
    ``"zero_inflated"``. Default is ``"continuous"``.
multiplier : float, optional
    Bandwidth multiplier. The actual bandwidth used is ``bandwidth * multiplier``.
    Default is 1.0.
bandwidth : float, optional
    Bandwidth parameter. Must be a positive number or ``NaN`` for automatic
    selection using the plug-in methodology. Default is ``NaN``.
degree : int, optional
    Degree of the local polynomial. Either 0, 1, or 2 for log-constant,
    log-linear, and log-quadratic fitting, respectively. Default is 2.
grid_size : int, optional
    Number of grid points for the interpolation grid. Must be at least 4.
    Default is 400.
)""";

constexpr const char* kde1d_from_params_doc = R"""(
Create a Kde1d object from parameters.

This is a factory method that creates a Kde1d object with specified parameters.
The object needs to be fitted to data using the ``fit()`` method.

Parameters
----------
xmin : float, optional
    Lower bound for the support of the density. Default is ``NaN``.
xmax : float, optional
    Upper bound for the support of the density. Default is ``NaN``.
type : str, optional
    Variable type (``"continuous"``, ``"discrete"``, or ``"zero_inflated"``).
    Default is ``"continuous"``.
multiplier : float, optional
    Bandwidth multiplier. Default is 1.0.
bandwidth : float, optional
    Bandwidth parameter (``NaN`` for automatic selection). Default is ``NaN``.
degree : int, optional
    Degree of the local polynomial (0, 1, or 2). Default is 2.
grid_size : int, optional
    Number of grid points for the interpolation grid. Must be at least 4.
    Default is 400.

Returns
-------
Kde1d
    A Kde1d object ready for fitting.
)""";

constexpr const char* kde1d_from_grid_doc = R"""(
Create a Kde1d object from grid points and density values.

This factory method creates a Kde1d object from pre-computed grid points and
corresponding density values. This is useful for loading previously fitted
models or creating models from externally computed densities.

Parameters
----------
grid_points : array_like
    Vector of grid points where the density was evaluated.
values : array_like
    Vector of density values corresponding to the grid points.
    Must have the same length as ``grid_points``.
xmin : float, optional
    Lower bound for the support of the density. Default is ``NaN``.
xmax : float, optional
    Upper bound for the support of the density. Default is ``NaN``.
type : str, optional
    Variable type. Default is ``"continuous"``.
prob0 : float, optional
    Point mass at 0 (for zero-inflated models). Default is 0.0.

Returns
-------
Kde1d
    A fitted Kde1d object ready for evaluation.
)""";

constexpr const char* fit_doc = R"""(
Fit the kernel density estimate to data.

Parameters
----------
x : array_like
    Vector of observations to fit the density to.
weights : array_like, optional
    Vector of weights for individual observations. If not provided,
    all observations are weighted equally.

Notes
-----
After calling this method, the object will be fitted and can be used
for density evaluation, sampling, etc.
)""";

constexpr const char* pdf_doc = R"""(
Evaluate the probability density function.

Computes the pdf of the kernel density estimate by interpolation.

Parameters
----------
x : array_like
    Vector of evaluation points.
check_fitted : bool, optional
    Whether to check if the model is fitted before evaluation.
    Default is ``True``.

Returns
-------
array_like
    Vector of pdf values at the evaluation points.
)""";

constexpr const char* cdf_doc = R"""(
Evaluate the cumulative distribution function.

Computes the cdf of the kernel density estimate by numerical integration.

Parameters
----------
x : array_like
    Vector of evaluation points.
check_fitted : bool, optional
    Whether to check if the model is fitted before evaluation.
    Default is ``True``.

Returns
-------
array_like
    Vector of cdf values at the evaluation points.
)""";

constexpr const char* quantile_doc = R"""(
Evaluate the quantile function.

Computes quantiles of the kernel density estimate by numerical inversion
of the cumulative distribution function.

Parameters
----------
x : array_like
    Vector of probabilities (between 0 and 1).
check_fitted : bool, optional
    Whether to check if the model is fitted before evaluation.
    Default is ``True``.

Returns
-------
array_like
    Vector of quantiles corresponding to the input probabilities.
)""";

constexpr const char* simulate_doc = R"""(
Simulate data from the fitted density.

Generates random samples from the kernel density estimate.

Parameters
----------
n : int
    Number of observations to simulate.
seeds : list of int, optional
    Optional vector of random seeds for reproducibility.
check_fitted : bool, optional
    Whether to check if the model is fitted before simulation.
    Default is ``True``.

Returns
-------
array_like
    Vector of simulated observations from the kernel density.
)""";

constexpr const char* set_xmin_xmax_doc = R"""(
Set the boundary parameters.

Parameters
----------
xmin : float, optional
    Lower bound for the support. ``NaN`` means no boundary.
    Default is ``NaN``.
xmax : float, optional
    Upper bound for the support. ``NaN`` means no boundary.
    Default is ``NaN``.
)""";

// Property docstrings
constexpr const char* xmin_doc = R"""(Lower bound of the density support.)""";
constexpr const char* xmax_doc = R"""(Upper bound of the density support.)""";
constexpr const char* type_doc = R"""(Variable type as VarType enum.)""";
constexpr const char* type_str_doc = R"""(Variable type as string.)""";
constexpr const char* prob0_doc =
    R"""(Point mass at 0 (for zero-inflated models).)""";
constexpr const char* multiplier_doc = R"""(Bandwidth multiplier.)""";
constexpr const char* bandwidth_doc = R"""(Bandwidth parameter.)""";
constexpr const char* degree_doc = R"""(Degree of the local polynomial.)""";
constexpr const char* grid_size_doc =
    R"""(Number of grid points for interpolation.)""";
constexpr const char* loglik_doc = R"""(Log-likelihood of the fitted model.)""";
constexpr const char* edf_doc = R"""(Effective degrees of freedom.)""";
constexpr const char* grid_points_doc =
    R"""(Grid points used for interpolation.)""";
constexpr const char* values_doc = R"""(Density values at grid points.)""";

}  // namespace kde1d_docstrings
