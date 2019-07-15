// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_integration.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {
inline Bb6Bicop::Bb6Bicop()
{
    family_ = BicopFamily::bb6;
    parameters_ = Eigen::VectorXd(2);
    parameters_lower_bounds_ = Eigen::VectorXd(2);
    parameters_upper_bounds_ = Eigen::VectorXd(2);
    parameters_ << 1, 1;
    parameters_lower_bounds_ << 1, 1;
    parameters_upper_bounds_ << 6, 8;
}

inline double Bb6Bicop::generator(const double &u)
{
    double res = tools_stl::log1p(-std::pow(1 - u, parameters_(0)));
    return std::pow((-1) * res, parameters_(1));
}

inline double Bb6Bicop::generator_inv(const double &u)
{
    double res = std::expm1(-std::pow(u, 1 / parameters_(1)));
    return 1 - std::pow(-res, 1 / parameters_(0));
}

inline double Bb6Bicop::generator_derivative(const double &u)
{
    double theta = double(parameters_(0));
    double delta = double(parameters_(1));
    double res = tools_stl::log1p(-std::pow(1 - u, theta));
    res = delta * theta * std::pow((-1) * res, delta - 1);
    return res * std::pow(1 - u, theta - 1) / (std::pow(1 - u, theta) - 1);
}

//inline double Bb6Bicop::generator_derivative2(const double &u)
//{
//    double theta = double(parameters_(0));
//    double delta = double(parameters_(1));
//    double tmp = std::pow(1 - u, theta);
//    double tmp2 = tools_stl::log1p(-tmp);
//    double res = std::pow((-1) * tmp2, delta - 2);
//    res *= ((delta - 1) * theta * tmp - (tmp + theta - 1) * tmp2);
//    return res * delta * theta * std::pow(1 - u, theta - 2) /
//           std::pow(tmp - 1, 2);
//}

inline Eigen::VectorXd Bb6Bicop::pdf_raw(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
)
{
    double theta = static_cast<double>(parameters_(0));
    double delta = static_cast<double>(parameters_(1));

    double t12 = 1 / delta;
    double t16 = 1 / theta;
    double t32 = delta - 1.0;
    double t38 = 2.0 * delta;
    double t39 = -1.0 + t38;
    double t47 = 3.0 * delta - 1.0;

    auto f = [theta, delta, t12, t16, t32, t38, t39, t47](const double &u1,
                                                          const double &u2) {

        double t1 = 1.0 - u1;
        double t2 = std::pow(t1,theta);
        double t3 = 1.0 - t2;
        double t4 = std::log(t3);
        double t5 = std::pow(-t4,delta);
        double t40 = std::pow(-t4,t39);
        double t50 = std::pow(-t4,t32);
        double t61 = std::pow(-t4,t47);
        double t90 = std::pow(-t4,t38);

        double t6 = 1.0 - u2;
        double t7 = std::pow(t6,theta);
        double t8 = 1.0 - t7;
        double t9 = std::log(t8);
        double t10 = std::pow(-t9,delta);
        double t11 = t5 + t10;
        double t13 = std::pow(t11,t12);
        double t14 = std::exp(-t13);
        double t35 = std::pow(t11,-2.0 * t32 * t12);
        double t36 = t35 * theta;
        double t37 = std::exp(t13);
        double t42 = std::pow(-t9,t39);
        double t48 = std::pow(-t9,t47);
        double t53 = t13 * delta;
        double t56 = std::pow(-t9,t32);
        double t57 = t37 * t50 * t56;
        double t59 = t13 * theta;
        double t78 = t37 - 1.0;
        double t80 = std::pow(t78 * t14, t16);
        double t87 = t78 * t78;
        double t93 = std::pow(-t9,t38);

        return (2.0 * t36 * t37 * t40 * t42
                + t36 * t37 * t48 * t50
                + t53 * theta * t57 - t59 * t57
                + t36 * t37 * t61 * t56 - 2.0 * t35 * t40 * t42
                - t35 * t61 * t56 - t53 * theta * t50 * t56+t59 * t50 * t56
                - t35 * t48 * t50)
               * t80 * t7 * t2 / t3 / t8 / t87 /
            (t90 + 2.0 * t5 * t10 + t93) / t1 / t6;
    };
    return tools_eigen::binaryExpr_or_nan(u, f);
}

inline double Bb6Bicop::parameters_to_tau(const Eigen::MatrixXd &parameters)
{
    double theta = parameters(0);
    double delta = parameters(1);
    auto f = [&theta, &delta](const double &v) {
        double res = -4 * (1 - v - std::pow(1 - v, -theta) +
                           std::pow(1 - v, -theta) * v);
        return 1 / (delta * theta) *
               tools_stl::log1p(-std::pow(1 - v, theta)) * res;
    };
    return 1 + tools_integration::integrate_zero_to_one(f);
}

inline Eigen::MatrixXd Bb6Bicop::tau_to_parameters(const double &tau)
{
    return no_tau_to_parameters(tau);
}
}
