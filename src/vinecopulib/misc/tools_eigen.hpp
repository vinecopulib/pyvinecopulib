// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <boost/math/special_functions/fpclassify.hpp> // isnan
#include <functional>

namespace vinecopulib {

//! Tools for working with Eigen types
namespace tools_eigen {
//! An `Eigen::Matrix` containing `bool`s (similar to `Eigen::MatrixXd`).
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;

template<typename T>
Eigen::MatrixXd
unaryExpr_or_nan(const Eigen::MatrixXd& x, const T& func)
{
  return x.unaryExpr([&func](const double& y) {
    if ((boost::math::isnan)(y)) {
      return std::numeric_limits<double>::quiet_NaN();
    } else {
      return func(y);
    }
  });
}

template<typename T>
Eigen::VectorXd
binaryExpr_or_nan(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u,
                  const T& func)
{
  auto func_or_nan = [&func](const double& u1, const double& u2) {
    if ((boost::math::isnan)(u1) | (boost::math::isnan)(u2)) {
      return std::numeric_limits<double>::quiet_NaN();
    } else {
      return func(u1, u2);
    }
  };
  return u.col(0).binaryExpr(u.col(1), func_or_nan);
}

void
remove_nans(Eigen::MatrixXd& x);

void
remove_nans(Eigen::MatrixXd& x, Eigen::VectorXd& weights);

bool
check_if_in_unit_cube(const Eigen::MatrixXd& u);

Eigen::Matrix<double, Eigen::Dynamic, 2>
swap_cols(Eigen::Matrix<double, Eigen::Dynamic, 2> u);

Eigen::VectorXd
invert_f(const Eigen::VectorXd& x,
         std::function<Eigen::VectorXd(const Eigen::VectorXd&)> f,
         const double lb = 1e-20,
         const double ub = 1 - 1e-20,
         int n_iter = 35);

Eigen::Matrix<double, Eigen::Dynamic, 2>
expand_grid(const Eigen::VectorXd& grid_points);

Eigen::MatrixXd
read_matxd(const char* filename, int max_buffer_size = static_cast<int>(1e6));

Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>
read_matxs(const char* filename, int max_buffer_size = static_cast<int>(1e6));
}
}

#include <vinecopulib/misc/implementation/tools_eigen.ipp>
