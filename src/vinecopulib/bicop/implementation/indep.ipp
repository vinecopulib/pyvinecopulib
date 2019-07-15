// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

inline IndepBicop::IndepBicop()
{
  family_ = BicopFamily::indep;
  parameters_ = Eigen::MatrixXd();
}

inline Eigen::VectorXd
IndepBicop::pdf_raw(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  auto f = [](double, double) { return 1.0; };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd
IndepBicop::cdf(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  return u.rowwise().prod();
}

inline Eigen::VectorXd
IndepBicop::hfunc1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  auto f = [](double, double u2) { return u2; };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd
IndepBicop::hfunc2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  auto f = [](double u1, double) { return u1; };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd
IndepBicop::hinv1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  auto f = [](double, double u2) { return u2; };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd
IndepBicop::hinv2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  auto f = [](double u1, double) { return u1; };
  return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::MatrixXd
IndepBicop::tau_to_parameters(const double&)
{
  return Eigen::VectorXd();
}

inline double
IndepBicop::parameters_to_tau(const Eigen::MatrixXd&)
{
  return 0.0;
}

inline Eigen::VectorXd
IndepBicop::get_start_parameters(const double tau)
{
  return tau_to_parameters(tau);
}

inline void
IndepBicop::flip()
{
  // nothing to do because independence copula is radially syemmetric
}
}
