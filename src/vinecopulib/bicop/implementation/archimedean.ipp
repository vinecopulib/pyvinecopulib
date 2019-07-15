// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {
// inline Eigen::VectorXd ArchimedeanBicop::pdf(
//    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
//)
//{
//    auto f = [this](const double &u1, const double &u2) {
//        double temp = generator_inv(generator(u1) + generator(u2));
//        temp = log(std::abs(generator_derivative2(temp))) -
//            3.0 * log(std::abs(generator_derivative(temp)));
//        temp += std::log(std::abs(generator_derivative(u1)));
//        temp += std::log(std::abs(generator_derivative(u2)));
//        return std::exp(temp);
//    };
//    return tools_eigen::binaryExpr_or_nan(u, f);
//}

inline Eigen::VectorXd
ArchimedeanBicop::cdf(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  auto f = [this](const double& u1, const double& u2) {
    return generator_inv(generator(u1) + generator(u2));
  };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd
ArchimedeanBicop::hfunc1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  auto f = [this](const double& u1, const double& u2) {
    double temp = generator_inv(generator(u1) + generator(u2));
    temp = generator_derivative(u1) / generator_derivative(temp);
    return std::min(temp, 1.0);
  };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd
ArchimedeanBicop::hfunc2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  return hfunc1(tools_eigen::swap_cols(u));
}

inline Eigen::VectorXd
ArchimedeanBicop::hinv1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  Eigen::VectorXd hinv = hinv1_num(u);
  return hinv;
}

inline Eigen::VectorXd
ArchimedeanBicop::hinv2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  return hinv1(tools_eigen::swap_cols(u));
}

inline Eigen::VectorXd
ArchimedeanBicop::get_start_parameters(const double)
{
  Eigen::MatrixXd lb = this->get_parameters_lower_bounds();
  Eigen::VectorXd parameters = lb + Eigen::VectorXd::Constant(2, 0.1);
  return parameters;
}
}
