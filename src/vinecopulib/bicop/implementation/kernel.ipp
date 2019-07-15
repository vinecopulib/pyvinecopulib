// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_interpolation.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <wdm/eigen.hpp>

namespace vinecopulib {
inline KernelBicop::KernelBicop()
{
  // construct default grid (equally spaced on Gaussian scale)
  size_t m = 30;
  Eigen::VectorXd grid_points(m);
  for (size_t i = 0; i < m; ++i)
    grid_points(i) = -3.25 + i * (6.5 / static_cast<double>(m - 1));
  interp_grid_ = std::make_shared<tools_interpolation::InterpolationGrid>(
    tools_stats::pnorm(grid_points),
    Eigen::MatrixXd::Constant(m, m, 1.0) // independence
  );
}

inline Eigen::VectorXd
KernelBicop::pdf_raw(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  return interp_grid_->interpolate(u);
}

inline Eigen::VectorXd
KernelBicop::cdf(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  return interp_grid_->integrate_2d(u);
}

inline Eigen::VectorXd
KernelBicop::hfunc1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  return interp_grid_->integrate_1d(u, 1);
}

inline Eigen::VectorXd
KernelBicop::hfunc2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  return interp_grid_->integrate_1d(u, 2);
}

inline Eigen::VectorXd
KernelBicop::hinv1(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  return hinv1_num(u);
}

inline Eigen::VectorXd
KernelBicop::hinv2(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u)
{
  return hinv2_num(u);
}

inline double
KernelBicop::parameters_to_tau(const Eigen::MatrixXd& parameters)
{
  set_parameters(parameters);
  std::vector<int> seeds = {
    204967043, 733593603, 184618802, 399707801, 290266245
  };
  auto u = tools_stats::ghalton(1000, 2, seeds);
  u.col(1) = hinv1(u);
  return wdm::wdm(u, "tau")(0, 1);
}

inline double
KernelBicop::calculate_npars()
{
  return npars_;
}

inline Eigen::MatrixXd
KernelBicop::get_parameters() const
{
  return interp_grid_->get_values();
}

inline Eigen::MatrixXd
KernelBicop::get_parameters_lower_bounds() const
{
  return Eigen::MatrixXd::Constant(30, 30, 0.0);
}

inline Eigen::MatrixXd
KernelBicop::get_parameters_upper_bounds() const
{
  return Eigen::MatrixXd::Constant(30, 30, 1e4);
}

inline void
KernelBicop::set_parameters(const Eigen::MatrixXd& parameters)
{
  if (parameters.minCoeff() < 0) {
    std::stringstream message;
    message << "density should be larger than 0. ";
    throw std::runtime_error(message.str().c_str());
  }
  interp_grid_->set_values(parameters);
}

inline void
KernelBicop::flip()
{
  interp_grid_->flip();
}

inline Eigen::MatrixXd
KernelBicop::tau_to_parameters(const double& tau)
{
  return no_tau_to_parameters(tau);
}
}
