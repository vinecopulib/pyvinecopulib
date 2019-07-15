// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {
inline Bb1Bicop::Bb1Bicop()
{
    family_ = BicopFamily::bb1;
    parameters_ = Eigen::VectorXd(2);
    parameters_lower_bounds_ = Eigen::VectorXd(2);
    parameters_upper_bounds_ = Eigen::VectorXd(2);
    parameters_ << 0, 1;
    parameters_lower_bounds_ << 0, 1;
    parameters_upper_bounds_ << 7, 7;
}

inline double Bb1Bicop::generator(const double &u)
{
    return std::pow(std::pow(u, -parameters_(0)) - 1, parameters_(1));
}

inline double Bb1Bicop::generator_inv(const double &u)
{
    return std::pow(std::pow(u, 1 / parameters_(1)) + 1, -1 / parameters_(0));
}

inline double Bb1Bicop::generator_derivative(const double &u)
{
    double theta = double(parameters_(0));
    double delta = double(parameters_(1));
    double res = -delta * theta * std::pow(u, -(1 + theta));
    return res * std::pow(std::pow(u, -theta) - 1, delta - 1);
}

//inline double Bb1Bicop::generator_derivative2(const double &u)
//{
//    double theta = double(parameters_(0));
//    double delta = double(parameters_(1));
//    double res = delta * theta * std::pow(std::pow(u, -theta) - 1, delta);
//    res /= (std::pow(std::pow(u, theta) - 1, 2) * std::pow(u, 2));
//    return res * (1 + delta * theta - (1 + theta) * std::pow(u, theta));
//}

inline Eigen::VectorXd Bb1Bicop::pdf_raw(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
)
{
    double theta = static_cast<double>(parameters_(0));
    double delta = static_cast<double>(parameters_(1));

    auto f = [theta, delta](const double &u1, const double &u2) {

        double t1 = std::pow(u1,-theta);
        double t2 = t1 - 1.0;
        double t3 = std::pow(t2, delta);
        double t16 = 1.0 / u1;
        double t17 = 1.0 / t2;
        double t38 = t1 * t16;
        double t39 = t38 * t17;
        double t4 = std::pow(u2,-theta);
        double t5 = t4 - 1.0;
        double t6 = std::pow(t5,delta);
        double t7 = t3 + t6;
        double t9 = std::pow(t7,1.0/delta);
        double t10 = 1.0 + t9;
        double t12 = std::pow(t10,-1.0/theta);
        double t13 = t12 * t9;
        double t20 = 1.0 / t10;
        double t24 = t9 * t9;
        double t25 = t12 * t24;
        double t27 = 1.0 / u2;
        double t29 = 1.0 / t5;
        double t32 = t7 * t7;
        double t33 = 1.0 / t32;
        double t34 = t10 * t10;
        double t36 = t33 / t34;
        double t43 = t4 * theta;
        double t59 = t43 * t27 * t29;

        return t25 * t6 * t27 * t4 * t29 * t36 * t3 * t39
               - t13 * t6 * t43 * t27 * t29 * t33 * t3 * t38 * t17 * t20
               + t13 * t3 * t38 * t17 * t33 * t20 * t6 * delta * t59
               + t25 * t3 * t39 * t36 * t6 * t59;
    };
    return tools_eigen::binaryExpr_or_nan(u, f);
}

inline double Bb1Bicop::parameters_to_tau(const Eigen::MatrixXd &parameters)
{
    return 1 - 2 / (parameters(1) * (parameters(0) + 2));
}

inline Eigen::MatrixXd Bb1Bicop::tau_to_parameters(const double &tau)
{
    return no_tau_to_parameters(tau);
}
}
