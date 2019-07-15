// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/parametric.hpp>

namespace vinecopulib {

namespace constant {
    constexpr double pi = 3.141592653589793238462643383279502884;
}

//! @brief An abstract class for elliptical copula families
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! @literature
//! Joe, Harry. Dependence modeling with copulas. CRC Press, 2014.
class EllipticalBicop : public ParBicop
{
private:
    // hfunction and its inverse
    Eigen::VectorXd hfunc2(
        const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
    );

    Eigen::VectorXd hinv2(
        const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
    );

    // link between Kendall's tau and the par_bicop parameter
    double parameters_to_tau(const Eigen::MatrixXd &parameters);
};
}

#include <vinecopulib/bicop/implementation/elliptical.ipp>
