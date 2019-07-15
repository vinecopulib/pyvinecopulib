// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <boost/math/special_functions/log1p.hpp>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {
inline ClaytonBicop::ClaytonBicop()
{
    family_ = BicopFamily::clayton;
    parameters_ = Eigen::VectorXd(1);
    parameters_lower_bounds_ = Eigen::VectorXd(1);
    parameters_upper_bounds_ = Eigen::VectorXd(1);
    parameters_ << 0;
    parameters_lower_bounds_ << 0;
    parameters_upper_bounds_ << 28;
}

inline double ClaytonBicop::generator(const double &u)
{
    double theta = double(this->parameters_(0));
    return (std::pow(u, -theta) - 1) / theta;
}

inline double ClaytonBicop::generator_inv(const double &u)
{
    double theta = double(this->parameters_(0));
    return std::pow(1 + theta * u, -1 / theta);
}

inline double ClaytonBicop::generator_derivative(const double &u)
{
    return (-1) * std::pow(u, -1 - this->parameters_(0));
}

//inline double ClaytonBicop::generator_derivative2(const double &u)
//{
//    double theta = double(this->parameters_(0));
//    return (1 + theta) * std::pow(u, -2 - theta);
//}

inline Eigen::VectorXd ClaytonBicop::pdf_raw(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
)
{
    double theta = static_cast<double>(parameters_(0));
    // avoid numerical issues when copula is too close to independence
    if (theta < 1e-10) {
        auto f = [](const double &, const double &) { return 1.0; };
        return tools_eigen::binaryExpr_or_nan(u, f);
    }

    auto f = [theta](const double &u1, const double &u2) {
        double temp = boost::math::log1p(theta)-(1.0+theta)*std::log(u1*u2);
        temp = temp - (2.0+1.0/(theta))*std::log(std::pow(u1,-theta) +
                                                     std::pow(u2,-theta)-1.0);
        return std::exp(temp);
    };
    return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd ClaytonBicop::hinv1(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
)
{
    double theta = double(this->parameters_(0));
    Eigen::VectorXd hinv = u.col(0).array().pow(theta + 1.0);
    if (theta < 75) {
        hinv = u.col(1).cwiseProduct(hinv);
        hinv = hinv.array().pow(-theta / (theta + 1.0));
        Eigen::VectorXd x = u.col(0);
        x = x.array().pow(-theta);
        hinv = hinv - x + Eigen::VectorXd::Ones(x.size());
        hinv = hinv.array().pow(-1 / theta);
    } else {
        hinv = hinv1_num(u);
    }
    return hinv;
}

inline Eigen::MatrixXd ClaytonBicop::tau_to_parameters(const double &tau)
{
    Eigen::VectorXd parameters(1);
    parameters(0) = 2 * std::fabs(tau) / (1 - std::fabs(tau));
    return parameters;
}

inline double
ClaytonBicop::parameters_to_tau(const Eigen::MatrixXd &parameters)
{
    return parameters(0) / (2 + std::fabs(parameters(0)));
}

inline Eigen::VectorXd ClaytonBicop::get_start_parameters(const double tau)
{
    Eigen::VectorXd par = tau_to_parameters(tau);
    par = par.cwiseMax(parameters_lower_bounds_);
    par = par.cwiseMin(parameters_upper_bounds_);
    return par;
}
}
