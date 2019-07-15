// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <vinecopulib/bicop/kernel.hpp>

namespace vinecopulib {
//! @brief The transformation local-constant likelihood estimator
//!
//! This class is used in the implementation underlying the Bicop class.
//! Users should not use AbstractBicop or derived classes directly, but
//! always work with the Bicop interface.
//!
//! @literature
//! Nagler, Thomas. *kdecopula: An R Package for the Kernel Estimation of
//! Copula Densities*. arXiv:1603.04229 [stat.CO], 2016
class TllBicop : public KernelBicop
{
public:
    TllBicop();

private:
    Eigen::VectorXd gaussian_kernel_2d(
        const Eigen::Matrix<double, Eigen::Dynamic, 2> &x);

    Eigen::Matrix2d select_bandwidth(
        const Eigen::Matrix<double, Eigen::Dynamic, 2> &x,
        std::string method,
        const Eigen::VectorXd& weights);

    Eigen::MatrixXd fit_local_likelihood(
        const Eigen::Matrix<double, Eigen::Dynamic, 2> &x,
        const Eigen::Matrix<double, Eigen::Dynamic, 2> &x_data,
        const Eigen::Matrix2d &B,
        std::string method,
        const Eigen::VectorXd& weights);

    double calculate_infl(const size_t &n,
                          const double &f0,
                          const Eigen::Vector2d &b,
                          const Eigen::Matrix2d &B,
                          const double &det_irB,
                          const Eigen::Matrix2d &S,
                          const std::string &method,
                          const double& weight);

    void fit(const Eigen::Matrix<double, Eigen::Dynamic, 2> &data,
             std::string method, 
             double mult, 
             const Eigen::VectorXd& weights);
};
}

#include <vinecopulib/bicop/implementation/tll.ipp>
