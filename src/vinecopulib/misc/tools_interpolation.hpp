// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>

namespace vinecopulib {

namespace tools_interpolation {
//! A class for cubic spline interpolation of bivariate copulas
//!
//! The class is used for implementing kernel estimators. It makes storing the
//! observations obsolete and allows for fast numerical integration.
class InterpolationGrid
{
public:
    InterpolationGrid()
    {
    }

    InterpolationGrid(const Eigen::VectorXd &grid_points,
                      const Eigen::MatrixXd &values,
                      int norm_times = 3);

    Eigen::MatrixXd get_values() const;

    void set_values(const Eigen::MatrixXd &values, int norm_times = 3);

    void flip();

    void normalize_margins(int times);

    Eigen::VectorXd interpolate(const Eigen::MatrixXd &x);

    Eigen::VectorXd integrate_1d(const Eigen::MatrixXd &u, size_t cond_var);

    Eigen::VectorXd integrate_2d(const Eigen::MatrixXd &u);

private:

    Eigen::Matrix<ptrdiff_t, 1, 2> get_indices(double x0, double x1);
    double bilinear_interpolation(double z11, double z12, double z21, double z22,
                                  double x1, double x2, double y1, double y2,
                                  double x, double y);
    double int_on_grid(const double &upr, const Eigen::VectorXd &vals,
                       const Eigen::VectorXd &grid);

    Eigen::VectorXd grid_points_;
    Eigen::MatrixXd values_;
};
}

}

#include <vinecopulib/misc/implementation/tools_interpolation.ipp>
