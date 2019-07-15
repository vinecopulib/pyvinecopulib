// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/archimedean.hpp>

namespace vinecopulib {
//! @brief The Joe copula
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! @literature
//! Joe, Harry. Dependence modeling with copulas. CRC Press, 2014.
class JoeBicop : public ArchimedeanBicop
{
public:
    // constructor
    JoeBicop();

private:
    // generator, its inverse and derivatives for the archimedean copula
    double generator(const double &u);

    double generator_inv(const double &u);

    double generator_derivative(const double &u);

    double generator_derivative2(const double &u);

    // pdf
    Eigen::VectorXd pdf_raw(const Eigen::Matrix<double, Eigen::Dynamic, 2> &u);

    // inverse hfunction
    Eigen::VectorXd hinv1(
        const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
    );

    // link between Kendall's tau and the par_bicop parameter
    Eigen::MatrixXd tau_to_parameters(const double &tau);

    double parameters_to_tau(const Eigen::MatrixXd &par);

    Eigen::VectorXd get_start_parameters(const double tau);
};
}

double qcondjoe(double *q, double *u, double *de);

#include <vinecopulib/bicop/implementation/joe.ipp>
