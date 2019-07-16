// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_eigen.hpp>
#include <vinecopulib/misc/tools_integration.hpp>

namespace vinecopulib {
inline Bb7Bicop::Bb7Bicop()
{
  family_ = BicopFamily::bb7;
  parameters_ = Eigen::VectorXd(2);
  parameters_lower_bounds_ = Eigen::VectorXd(2);
  parameters_upper_bounds_ = Eigen::VectorXd(2);
  parameters_ << 1, 1;
  parameters_lower_bounds_ << 1, 0;
  parameters_upper_bounds_ << 6, 25;
}

inline double
Bb7Bicop::generator(const double& u)
{
  double theta = double(parameters_(0));
  double delta = double(parameters_(1));
  return std::pow(1 - std::pow(1 - u, theta), -delta) - 1;
}

inline double
Bb7Bicop::generator_inv(const double& u)
{
  double theta = double(parameters_(0));
  double delta = double(parameters_(1));
  return 1 - std::pow(1 - std::pow(1 + u, -1 / delta), 1 / theta);
}

inline double
Bb7Bicop::generator_derivative(const double& u)
{
  double theta = double(parameters_(0));
  double delta = double(parameters_(1));
  double res = delta * theta * std::pow(1 - std::pow(1 - u, theta), -1 - delta);
  return -res * std::pow(1 - u, theta - 1);
}

// inline double Bb7Bicop::generator_derivative2(const double &u)
//{
//    double theta = double(parameters_(0));
//    double delta = double(parameters_(1));
//    double tmp = std::pow(1 - u, theta);
//    double res = delta * theta * std::pow(1 - tmp, -2 - delta) *
//                 std::pow(1 - u, theta - 2);
//    return res * (theta - 1 + (1 + delta * theta) * tmp);
//}

inline Eigen::VectorXd
Bb7Bicop::pdf_raw(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  double theta = static_cast<double>(parameters_(0));
  double delta = static_cast<double>(parameters_(1));

  auto f = [theta, delta](const double& u1, const double& u2) {
    double t1 = 1.0 - u1;
    double t2 = std::pow(t1, theta);
    double t3 = 1.0 - t2;
    double t4 = std::pow(t3, -delta);
    double t5 = 1.0 - u2;
    double t6 = std::pow(t5, theta);
    double t7 = 1.0 - t6;
    double t8 = std::pow(t7, -delta);
    double t9 = t4 + t8 - 1.0;
    double t11 = std::pow(t9, -1.0 / delta);
    double t12 = 1.0 - t11;
    double t14 = std::pow(t12, 1.0 / theta);
    double t15 = t11 * t11;
    double t16 = t14 * t15;
    double t18 = 1.0 / t5;
    double t20 = 1.0 / t7;
    double t23 = t9 * t9;
    double t24 = 1.0 / t23;
    double t25 = t12 * t12;
    double t27 = t24 / t25;
    double t30 = t2 / t1;
    double t31 = 1.0 / t3;
    double t32 = t30 * t31;
    double t35 = t14 * t11;
    double t37 = t6 * theta;
    double t42 = 1.0 / t12;
    double t54 = t37 * t18 * t20;

    return -t16 * t8 * t6 * t18 * t20 * t27 * t4 * t32 +
           t35 * t8 * t37 * t18 * t20 * t24 * t4 * t30 * t31 * t42 +
           t35 * t4 * t30 * t31 * t24 * t42 * t8 * delta * t54 +
           t16 * t4 * t32 * t27 * t8 * t54;
  };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline double
Bb7Bicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  double theta = parameters(0);
  double delta = parameters(1);
  auto f = [&theta, &delta](const double& v) {
    double tmp = std::pow(1 - v, theta);
    double res = -4 * (std::pow(1 - tmp, -delta) - 1) / (theta * delta);
    return res / (std::pow(1 - v, theta - 1) * std::pow(1 - tmp, -delta - 1));
  };
  return 1 + tools_integration::integrate_zero_to_one(f);
}

inline Eigen::MatrixXd
Bb7Bicop::tau_to_parameters(const double& tau)
{
  return no_tau_to_parameters(tau);
}
}
