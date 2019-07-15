// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_stats.hpp>

namespace vinecopulib {
inline GaussianBicop::GaussianBicop()
{
  family_ = BicopFamily::gaussian;
  parameters_ = Eigen::VectorXd(1);
  parameters_lower_bounds_ = Eigen::VectorXd(1);
  parameters_upper_bounds_ = Eigen::VectorXd(1);
  parameters_ << 0;
  parameters_lower_bounds_ << -1;
  parameters_upper_bounds_ << 1;
}

inline Eigen::VectorXd
GaussianBicop::pdf_raw(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  // Inverse Cholesky of the correlation matrix
  double rho = double(this->parameters_(0));
  Eigen::Matrix2d L;
  L(0, 0) = 1;
  L(1, 1) = 1 / sqrt(1.0 - pow(rho, 2.0));
  L(0, 1) = -rho * L(1, 1);
  L(1, 0) = 0;

  // Compute copula density
  Eigen::VectorXd f = Eigen::VectorXd::Ones(u.rows());
  Eigen::Matrix<double, Eigen::Dynamic, 2> tmp = tools_stats::qnorm(u);
  f = f.cwiseQuotient(tools_stats::dnorm(tmp).rowwise().prod());
  tmp = tmp * L;
  f = f.cwiseProduct(tools_stats::dnorm(tmp).rowwise().prod());
  return f / sqrt(1.0 - pow(rho, 2.0));
}

inline Eigen::VectorXd
GaussianBicop::cdf(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  return tools_stats::pbvnorm(tools_stats::qnorm(u),
                              double(this->parameters_(0)));
}

inline Eigen::VectorXd
GaussianBicop::hfunc1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  double rho = double(this->parameters_(0));
  Eigen::VectorXd h = Eigen::VectorXd::Zero(u.rows());
  Eigen::Matrix<double, Eigen::Dynamic, 2> tmp = tools_stats::qnorm(u);
  h = (tmp.col(1) - rho * tmp.col(0)) / sqrt(1.0 - pow(rho, 2.0));
  return tools_stats::pnorm(h);
}

inline Eigen::VectorXd
GaussianBicop::hinv1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  double rho = double(this->parameters_(0));
  Eigen::VectorXd hinv = Eigen::VectorXd::Zero(u.rows());
  Eigen::Matrix<double, Eigen::Dynamic, 2> tmp = tools_stats::qnorm(u);
  hinv = tmp.col(1) * sqrt(1.0 - pow(rho, 2.0)) + rho * tmp.col(0);
  return tools_stats::pnorm(hinv);
}

inline Eigen::VectorXd
GaussianBicop::get_start_parameters(const double tau)
{
  return tau_to_parameters(tau);
}

inline Eigen::MatrixXd
GaussianBicop::tau_to_parameters(const double& tau)
{
  Eigen::VectorXd parameters = this->parameters_;
  parameters(0) = std::sin(tau * constant::pi / 2);
  return parameters;
}
}
