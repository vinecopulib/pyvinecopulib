// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <boost/math/distributions.hpp>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

namespace tools_stats {


//! @brief Density function of the Standard normal distribution.
//!
//! @param x evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated densities.
inline Eigen::MatrixXd dnorm(const Eigen::MatrixXd &x)
{
    boost::math::normal dist;
    auto f = [&dist](double y) { return boost::math::pdf(dist, y); };
    return tools_eigen::unaryExpr_or_nan(x, f);
}

//! @brief Distribution function of the Standard normal distribution.
//!
//! @param x evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated probabilities.
inline Eigen::MatrixXd pnorm(const Eigen::MatrixXd &x)
{
    boost::math::normal dist;
    auto f = [&dist](double y) { return boost::math::cdf(dist, y); };
    return tools_eigen::unaryExpr_or_nan(x, f);
}

//! @brief Quantile function of the Standard normal distribution.
//!
//! @param x evaluation points.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated quantiles.
inline Eigen::MatrixXd qnorm(const Eigen::MatrixXd &x)
{
    boost::math::normal dist;
    auto f = [&dist](double y) { return boost::math::quantile(dist, y); };
    return tools_eigen::unaryExpr_or_nan(x, f);
}

//! @brief Density function of the Student t distribution.
//!
//! @param x evaluation points.
//! @param nu degrees of freedom parameter.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated densities.
inline Eigen::MatrixXd dt(const Eigen::MatrixXd &x, double nu)
{
    boost::math::students_t dist(nu);
    auto f = [&dist](double y) { return boost::math::pdf(dist, y); };
    return tools_eigen::unaryExpr_or_nan(x, f);
}

//! @brief Distribution function of the Student t distribution.
//!
//! @param x evaluation points.
//! @param nu degrees of freedom parameter.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated probabilities.
inline Eigen::MatrixXd pt(const Eigen::MatrixXd &x, double nu)
{
    boost::math::students_t dist(nu);
    auto f = [&dist](double y) { return boost::math::cdf(dist, y); };
    return tools_eigen::unaryExpr_or_nan(x, f);
}

//! @brief Quantile function of the Student t distribution.
//!
//! @param x evaluation points.
//! @param nu degrees of freedom parameter.
//!
//! @return An \f$ n \times d \f$ matrix of evaluated quantiles.
inline Eigen::MatrixXd qt(const Eigen::MatrixXd &x, double nu)
{
    boost::math::students_t dist(nu);
    auto f = [&dist](double y) { return boost::math::quantile(dist, y); };
    return tools_eigen::unaryExpr_or_nan(x, f);
}

Eigen::MatrixXd simulate_uniform(const size_t& n, const size_t& d,
                                 std::vector<int> seeds = std::vector<int>());

Eigen::MatrixXd simulate_uniform(const size_t& n, const size_t& d, bool qrng,
                                  std::vector<int> seeds = std::vector<int>());

Eigen::VectorXd to_pseudo_obs_1d(Eigen::VectorXd x,
                                 std::string ties_method = "average");

Eigen::MatrixXd to_pseudo_obs(Eigen::MatrixXd x,
                              std::string ties_method = "average");

double pairwise_mcor(const Eigen::Matrix<double, Eigen::Dynamic, 2>& x,
                     const Eigen::VectorXd &weights = Eigen::VectorXd());

Eigen::MatrixXd dependence_matrix(const Eigen::MatrixXd &x,
                                  const std::string &measure);

Eigen::MatrixXd ghalton(const size_t& n, const size_t& d,
                        std::vector<int> seeds = std::vector<int>());

Eigen::MatrixXd sobol(const size_t& n, const size_t& d,
                      std::vector<int> seeds = std::vector<int>());

Eigen::VectorXd pbvt(const Eigen::Matrix<double, Eigen::Dynamic, 2> &z,
                     int nu, double rho);

Eigen::VectorXd pbvnorm(const Eigen::Matrix<double, Eigen::Dynamic, 2> &z,
                        double rho);
}

}

#include <vinecopulib/misc/implementation/tools_stats.ipp>
