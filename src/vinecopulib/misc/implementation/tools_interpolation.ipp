// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <stdexcept>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {

namespace tools_interpolation {
//! Constructor
//!
//! @param grid_points an ascending sequence of grid_points; used in both
//! dimensions.
//! @param values a dxd matrix of copula density values evaluated at
//! (grid_points_i, grid_points_j).
//! @param norm_times how many times the normalization routine should run.
inline InterpolationGrid::InterpolationGrid(const Eigen::VectorXd &grid_points,
                                            const Eigen::MatrixXd &values,
                                            int norm_times)
{
    if (values.cols() != values.rows()) {
        throw std::runtime_error("values must be a quadratic matrix");
    }
    if (grid_points.size() != values.rows()) {
        throw std::runtime_error(
            "number of grid_points must equal dimension of values");
    }

    grid_points_ = grid_points;
    values_ = values;
    normalize_margins(norm_times);
}

inline Eigen::MatrixXd InterpolationGrid::get_values() const
{
    return values_;
}

inline void InterpolationGrid::set_values(const Eigen::MatrixXd &values,
                                          int norm_times)
{
    if (values.size() != values_.size()) {
        if (values.rows() != values_.rows()) {
            std::stringstream message;
            message <<
                    "values have has wrong number of rows; " <<
                    "expected: " << values_.rows() << ", " <<
                    "actual: " << values.rows() << std::endl;
            throw std::runtime_error(message.str().c_str());
        }
        if (values.cols() != values_.cols()) {
            std::stringstream message;
            message <<
                    "values have wrong number of columns; " <<
                    "expected: " << values_.cols() << ", " <<
                    "actual: " << values.cols() << std::endl;
            throw std::runtime_error(message.str().c_str());
        }
    }

    values_ = values;
    normalize_margins(norm_times);
}

inline void InterpolationGrid::flip()
{
    values_.transposeInPlace();
}

//! renormalizes the estimate to uniform margins
//!
//! @param times how many times the normalization routine should run.
inline void InterpolationGrid::normalize_margins(int times)
{
    size_t m = grid_points_.size();
    for (int k = 0; k < times; ++k) {
        for (size_t i = 0; i < m; ++i) {
            values_.row(i) /= int_on_grid(1.0, values_.row(i), grid_points_);
        }
        for (size_t j = 0; j < m; ++j) {
            values_.col(j) /= int_on_grid(1.0, values_.col(j), grid_points_);
        }
    }
}

inline Eigen::Matrix<ptrdiff_t, 1, 2> InterpolationGrid::get_indices(
    double x0, double x1)
{
    Eigen::Matrix<ptrdiff_t, 1, 2> out;
    out << 0, 0;
    bool found_i = false;
    bool found_j = false;
    for (ptrdiff_t k = 1; k < (grid_points_.size() - 1); ++k) {
        if ((x0 >= grid_points_(k))) {
            out(0) = k;
        } else {
            found_i = true;
        }
        if ((x1 >= grid_points_(k))) {
            out(1) = k;
        } else {
            found_j = true;
        }
        if (found_i & found_j) {
            break;
        }
    }
    return out;
}

//! Interpolate linearly in two dimensions
//!
//! @param z11 value corresponding to (x1, y1)
//! @param z12 value corresponding to (x1, y2)
//! @param z21 value corresponding to (x2, y1)
//! @param z22 value corresponding to (x2, y2)
//! @param x1 first cell value for the first dimension
//! @param x2 second cell value for the first dimension
//! @param y1 first cell value for the second dimension
//! @param y2 second cell value for the second dimension
//! @param x evaluation point for the first dimension
//! @param y evaluation point for the second dimension
inline double InterpolationGrid::bilinear_interpolation(
    double z11, double z12, double z21, double z22,
    double x1, double x2, double y1, double y2,
    double x, double y)
{
    double x2x1, y2y1, x2x, y2y, yy1, xx1;
    x2x1 = x2 - x1;
    y2y1 = y2 - y1;
    x2x = x2 - x;
    y2y = y2 - y;
    yy1 = y - y1;
    xx1 = x - x1;
    return (z11 * x2x * y2y +
            z21 * xx1 * y2y +
            z12 * x2x * yy1 +
            z22 * xx1 * yy1) / (x2x1 * y2y1);
}

//! Interpolation in two dimensions
//!
//! @param x mx2 matrix of evaluation points.
//! @return a vector of resulting interpolated values
inline Eigen::VectorXd
InterpolationGrid::interpolate(const Eigen::MatrixXd &x)
{

    auto f = [this](double x0, double x1) {

        auto indices = this->get_indices(x0, x1);
        return fmax(bilinear_interpolation(
            this->values_(indices(0), indices(1)),
            this->values_(indices(0), indices(1) + 1),
            this->values_(indices(0) + 1, indices(1)),
            this->values_(indices(0) + 1, indices(1) + 1),
            this->grid_points_(indices(0)),
            this->grid_points_(indices(0) + 1),
            this->grid_points_(indices(1)),
            this->grid_points_(indices(1) + 1),
            x0, x1), 1e-15);
    };

    return tools_eigen::binaryExpr_or_nan(x, f);
}

//! Integrate the grid along one axis
//!
//! @param u mx2 matrix of evaluation points
//! @param cond_var either 1 or 2; the axis considered fixed.
//! @return a vector of resulting integral values
inline Eigen::VectorXd
InterpolationGrid::integrate_1d(const Eigen::MatrixXd &u,
                                 size_t cond_var)
{
    ptrdiff_t m = grid_points_.size();
    Eigen::VectorXd tmpvals(m);
    Eigen::MatrixXd tmpgrid(m, 2);

    auto f = [this, m, cond_var, &tmpvals, &tmpgrid](double u1, double u2) {
        double upr = 0.0;
        double tmpint = 0.0, int1;
        if (cond_var == 1) {
            upr = u2;
            tmpgrid.col(0) = Eigen::VectorXd::Constant(m, u1);
            tmpgrid.col(1) = grid_points_;
        } else if (cond_var == 2) {
            upr = u1;
            tmpgrid.col(0) = grid_points_;
            tmpgrid.col(1) = Eigen::VectorXd::Constant(m, u2);
        }
        tmpvals = interpolate(tmpgrid).array().max(1e-4);
        tmpint = int_on_grid(upr, tmpvals, grid_points_);
        int1 = int_on_grid(1.0, tmpvals, grid_points_);

        return fmin(fmax(tmpint / int1, 1e-10), 1 - 1e-10);
    };

    return tools_eigen::binaryExpr_or_nan(u, f);
}

//! Integrate the grid along the two axis
//!
//! @param u mx2 matrix of evaluation points
//! @return a vector of resulting integral values
inline Eigen::VectorXd
InterpolationGrid::integrate_2d(const Eigen::MatrixXd &u)
{
    ptrdiff_t m = grid_points_.size();
    Eigen::VectorXd tmpvals(m), tmpvals2(m);
    Eigen::MatrixXd tmpgrid(m, 2);
    tmpgrid.col(1) = grid_points_;

    auto f = [this, m, &tmpvals, &tmpvals2, &tmpgrid](double u1, double u2) {
        double upr, tmpint, tmpint1;
        upr = u2;
        for (ptrdiff_t k = 0; k < m; ++k) {
            tmpgrid.col(0) = Eigen::VectorXd::Constant(m, grid_points_(k));
            tmpvals = interpolate(tmpgrid);
            tmpint = int_on_grid(upr, tmpvals, grid_points_);
            tmpvals2(k) = tmpint;
        }
        upr = u1;
        tmpint = int_on_grid(upr, tmpvals2, grid_points_);
        tmpint1 = int_on_grid(1.0, tmpvals2, grid_points_);
        return fmin(fmax(tmpint / tmpint1, 1e-10), 1 - 1e-10);
    };

    return tools_eigen::binaryExpr_or_nan(u, f);
}


// ---------------- Utility functions for integration ----------------


//! Integrate using a trapezoid rule
//!
//! @param upr upper limit of integration (lower is 0).
//! @param vals vector of values to be interpolated and integrated.
//! @param grid vector of grid points on which vals has been computed.
//!
//! @return integral of a piecewise linear function defined by (grid_i, vals_i).
inline double InterpolationGrid::int_on_grid(const double &upr,
                                             const Eigen::VectorXd &vals,
                                             const Eigen::VectorXd &grid)
{
    double tmpint = 0.0;

    if (upr > grid(0)) {
        // go up the grid and integrate
        for (ptrdiff_t k = 0; k <  (grid.size() - 1); ++k) {
            // stop loop if fully integrated
            if (upr < grid(k))
                break;

            // don't integrate over full cell if upr is in interior
            if (upr < grid(k + 1)) {
                tmpint += (2 * vals(k) +
                    (vals(k + 1) - vals(k)) *
                        (upr - grid(k)) / (grid(k + 1) - grid(k))) *
                    (upr - grid(k)) / 2.0;
            } else {
                tmpint += (vals(k + 1) + vals(k)) *
                    (grid(k + 1) - grid(k)) / 2.0;
            }
        }
    }

    return tmpint;
}
}

}
