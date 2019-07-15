// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/bicop/family.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_interpolation.hpp>
#include <boost/math/special_functions/fpclassify.hpp> // isnan
#include <wdm/eigen.hpp>

namespace vinecopulib {
inline TllBicop::TllBicop()
{
    family_ = BicopFamily::tll;
}

inline Eigen::VectorXd TllBicop::gaussian_kernel_2d(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &x
)
{
    return tools_stats::dnorm(x).rowwise().prod();
}

//! selects the bandwidth matrix for local líkelihood estimator (covariance
//! times appropriate factor).
inline Eigen::Matrix2d TllBicop::select_bandwidth(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &x,
    std::string method,
    const Eigen::VectorXd& weights
)
{
    size_t n = x.rows();
    double cor = wdm::wdm(x, "cor", weights)(0, 1);
    Eigen::Matrix2d cov = Eigen::MatrixXd::Identity(2, 2);
    cov(0, 1) = cor;
    cov(1, 0) = cor;

    double mult;
    if (method == "constant") {
        mult = std::pow(n, -1.0 / 3.0);
    } else {
        double degree;
        if (method == "linear") {
            degree = 1.0;
        } else {
            degree = 2.0;
        }
        mult = 1.5 * std::pow(n, -1.0 / (2.0 * degree + 1.0));
    }
    double mcor = tools_stats::pairwise_mcor(x, weights);
    double scale = std::pow(std::fabs(cor / mcor), 0.5 * mcor);

    return mult * cov * scale;
}

//! calculates the cholesky root of a 2x2 matrix.
inline Eigen::Matrix2d chol22(const Eigen::Matrix2d &B)
{

    Eigen::Matrix2d rB;

    rB(0, 0) = std::sqrt(B(0, 0));
    rB(0, 1) = 0.0;
    rB(1, 0) = B(1, 0) / rB(0, 0);
    rB(1, 1) = std::sqrt(B(1, 1) - rB(1, 0) * rB(1, 0));

    return rB;
}

//! evaluates local likleihood density estimate.
//!
//! @param x evaluation points.
//! @param x_data observations.
//! @param B bandwidth matrix.
//! @param method order of local polynomial approximation; either `"constant"`,
//!   `"linear"`, or `"quadratic"`.
//! @param weights vector of weights for the observations
//! @return a two-column matrix; first column is estimated density, second
//!    column is influence of evaluation point.
inline Eigen::MatrixXd TllBicop::fit_local_likelihood(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &x,
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &x_data,
    const Eigen::Matrix2d &B,
    std::string method,
    const Eigen::VectorXd& weights)
{
    size_t m = x.rows();       // number of evaluation points
    size_t n = x_data.rows();  // number of observations

    // pre-calculate inverse root of bandwidth matrix and determinant
    Eigen::Matrix2d irB = chol22(B).inverse();
    double det_irB = irB.determinant();

    // de-correlate data by applying B^{-1/2}
    Eigen::MatrixXd z = (irB * x.transpose()).transpose();
    Eigen::MatrixXd z_data = (irB * x_data.transpose()).transpose();


    Eigen::MatrixXd res(m, 2);
    res.col(0) = Eigen::VectorXd::Ones(m);  // result will be a product
    Eigen::VectorXd kernels(n);
    double f0;
    Eigen::Vector2d f1;
    Eigen::Vector2d b;
    Eigen::Matrix2d S(B);
    Eigen::MatrixXd zz(n, 2), zz2(n, 2);
    for (size_t k = 0; k < m; ++k) {
        zz = z_data - z.row(k).replicate(n, 1);
        kernels = gaussian_kernel_2d(zz) * det_irB;
        if (weights.size() > 0)
            kernels = kernels.cwiseProduct(weights);
        f0 = kernels.mean();
        if (method != "constant") {
            zz = (irB * zz.transpose()).transpose();
            f1 = zz.cwiseProduct(kernels.replicate(1, 2)).colwise().mean();
            b = f1 / f0;
            if (method == "quadratic") {
                zz2 = zz.cwiseProduct(kernels.replicate(1, 2)) /
                      (f0 * static_cast<double>(n));
                b = B * b;
                S = (B * (zz.transpose() * zz2) * B -
                     b * b.transpose()).inverse();
                res(k) *= std::sqrt(S.determinant()) / det_irB;
            }
            res(k) *= std::exp(-0.5 * double(b.transpose() * S * b));
            if ((boost::math::isnan)(res(k)) |
                (boost::math::isinf)(res(k))) {
                // inverse operation might go wrong due to rounding when
                // true value is equal or close to zero
                res(k) = 0.0;
            }
        }
        res(k, 0) *= f0;
        if (weights.size() > 0) {
            res(k, 1) = calculate_infl(n, f0, b, B, det_irB, S, method, weights(k));
        } else {
            res(k, 1) = calculate_infl(n, f0, b, B, det_irB, S, method, 1.0);
        }
    }

    return res;
}

//! calculate influence for data point for density estimate based on
//! quantities pre-computed in `fit_local_likelihood()`.
inline double TllBicop::calculate_infl(const size_t &n,
                                       const double &f0,
                                       const Eigen::Vector2d &b,
                                       const Eigen::Matrix2d &B,
                                       const double &det_irB,
                                       const Eigen::Matrix2d &S,
                                       const std::string &method,
                                       const double& weight)
{
    Eigen::MatrixXd M;
    if (method == "constant") {
        M = Eigen::MatrixXd::Constant(1, 1, f0);
    } else if (method == "linear") {
        M = Eigen::MatrixXd(3, 3);
        M(0, 0) = f0;
        M.col(0).tail(2) = B * b * f0;
        M.row(0).tail(2) = M.col(0).tail(2);
        M.block(1, 1, 2, 2) = f0 * B + f0 * B * b * b.transpose() * B;
    } else if (method == "quadratic") {
        M = Eigen::MatrixXd::Zero(6, 6);
        M(0, 0) = f0;
        M.col(0).segment(1, 2) = f0 * b;
        M.row(0).segment(1, 2) = M.col(0).segment(1, 2);
        M.block(1, 1, 2, 2) = f0 * B + f0 * b * b.transpose();
        M(3, 0) = 0.5 * M(1, 1);
        M(4, 0) = 0.5 * M(2, 2);
        M(5, 0) = M(1, 2);
        M.row(0).tail(3) = M.col(0).tail(3);
        Eigen::MatrixXd Si = S.inverse();
        M(3, 1) = 0.5 * f0 * (3.0 * Si(0, 0) * b(0) + std::pow(b(0), 3));
        M(4, 2) = 0.5 * f0 * (3.0 * Si(1, 1) * b(1) + std::pow(b(1), 3));
        M(4, 1) = 0.5 * f0;
        M(4, 1) *= 2.0 * Si(0, 1) * b(1) + Si(1, 1) * b(0) + b(0) * b(1) * b(1);
        M(3, 2) = 0.5 * f0;
        M(3, 2) *= 2.0 * Si(0, 1) * b(0) + Si(0, 0) * b(1) + b(1) * b(0) * b(0);
        M(5, 1) = 2.0 * M(3, 2);
        M(5, 2) = 2.0 * M(4, 1);
        M.block(1, 3, 2, 3) = M.block(3, 1, 3, 2).transpose();
        M(3, 3) = 0.25 * f0;
        M(3, 3) *= 3.0 * Si(0, 0) * Si(0, 0) + 6.0 * Si(0, 0) * b(0) * b(0) +
                   std::pow(b(0), 4);
        M(4, 4) = 0.25 * f0;
        M(4, 4) *= 3.0 * Si(1, 1) * Si(1, 1) + 6.0 * Si(1, 1) * b(1) * b(1) +
                   std::pow(b(1), 4);
        M(5, 5) =
            Si(0, 0) * Si(1, 1) + 2.0 * S(0, 1) + b(0) * b(0) * b(1) * b(1);
        M(5, 5) += 4.0 * Si(0, 1) * b(0) * b(1);
        M(5, 5) += Si(0, 0) * b(1) * b(1) + Si(1, 1) * b(0) * b(0);
        M(5, 5) *= f0;
        M(4, 3) = M(5, 5) * 0.25;
        M(3, 4) = M(4, 3);
        M(5, 3) = 3.0 * Si(0, 0) * Si(0, 1) + 3.0 * Si(0, 1) * b(0) * b(0);
        M(5, 3) += 3.0 * Si(0, 0) * b(0) * b(1) + b(1) * std::pow(b(0), 3);
        M(5, 3) *= 0.5 * f0;
        M(3, 5) = M(5, 3);
        M(5, 4) = 3.0 * Si(1, 1) * Si(0, 1) + 3.0 * Si(0, 1) * b(1) * b(1);
        M(5, 4) += 3.0 * Si(1, 1) * b(0) * b(1) + b(0) * std::pow(b(1), 3);
        M(5, 4) *= 0.5 * f0;
        M(4, 5) = M(5, 4);
    }

    double infl = gaussian_kernel_2d(Eigen::MatrixXd::Zero(1, 2))(0) * det_irB;
    infl *= M.inverse()(0, 0) * weight / static_cast<double>(n);
    return infl;
}


inline void TllBicop::fit(const Eigen::Matrix<double, Eigen::Dynamic, 2> &data,
                          std::string method,
                          double mult,
                          const Eigen::VectorXd& weights)
{
    using namespace tools_interpolation;

    // construct default grid (equally spaced on Gaussian scale)
    size_t m = 30;
    Eigen::VectorXd grid_points(m);
    for (size_t i = 0; i < m; ++i)
        grid_points(i) = -3.25 + i * (6.5 / static_cast<double>(m - 1));
    grid_points = tools_stats::pnorm(grid_points);

    // expand the interpolation grid; a matrix with two columns where each row
    // contains one combination of the grid points
    auto grid_2d = tools_eigen::expand_grid(grid_points);

    // transform evaluation grid and data by inverse Gaussian cdf
    Eigen::Matrix<double, Eigen::Dynamic, 2> z = tools_stats::qnorm(grid_2d);
    Eigen::Matrix<double, Eigen::Dynamic, 2> z_data = tools_stats::qnorm(data);

    // find bandwidth matrix
    Eigen::Matrix2d B = select_bandwidth(z_data, method, weights);
    B *= mult;

    // compute the density estimator (first column estimate, second influence)
    Eigen::MatrixXd ll_fit = fit_local_likelihood(z, z_data, B, method, weights);

    // transform density estimate to copula scale
    Eigen::VectorXd c =
        ll_fit.col(0).cwiseQuotient(tools_stats::dnorm(z).rowwise().prod());
    // store values in mxm grid
    Eigen::MatrixXd values(m, m);
    values = Eigen::Map<Eigen::MatrixXd>(c.data(), m, m).transpose();

    // for interpolation, we shift the limiting gridpoints to 0 and 1
    grid_points(0) = 0.0;
    grid_points(m - 1) = 1.0;
    interp_grid_ = std::make_shared<InterpolationGrid>(grid_points, values);

    // compute effective degrees of freedom via interpolation ---------
    Eigen::VectorXd infl_vec = ll_fit.col(1);
    // stabilize interpolation by restricting to plausible range
    infl_vec = infl_vec.array().min(1.0).max(-1.0);
    Eigen::MatrixXd infl(m, m);
    infl = Eigen::Map<Eigen::MatrixXd>(infl_vec.data(), m, m).transpose();
    // don't normalize margins of the EDF! (norm_times = 0)
    auto infl_grid = InterpolationGrid(grid_points, infl, 0);
    npars_ = infl_grid.interpolate(data).sum();
    set_loglik(loglik(data, weights));
}
}
