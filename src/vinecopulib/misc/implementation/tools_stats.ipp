// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_stats_ghalton.hpp>
#include <vinecopulib/misc/tools_stats_sobol.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <random>
#include <unsupported/Eigen/FFT>
#include <wdm/eigen.hpp>

namespace vinecopulib {

//! Utilities for statistical analysis
namespace tools_stats {

//! @brief simulates from the multivariate uniform distribution.
//!
//! @param n number of observations.
//! @param d dimension.
//! @param seeds seeds of the random number generator; if empty (default),
//!   the random number generator is seeded randomly.
//!
//! @return An \f$ n \times d \f$ matrix of independent
//! \f$ \mathrm{U}[0, 1] \f$ random variables.
inline Eigen::MatrixXd simulate_uniform(const size_t& n, const size_t& d,
                                        std::vector<int> seeds)
{
    if ((n < 1) | (d < 1)) {
        throw std::runtime_error("both n and d must be at least 1.");
    }

    if (seeds.size() == 0) {
        // no seeds provided, seed randomly
        std::random_device rd{};
        seeds = std::vector<int>(5);
        for (auto& s : seeds)
            s = static_cast<int>(rd());
    }

    // initialize random engine and uniform distribution
    std::seed_seq seq(seeds.begin(), seeds.end());
    std::mt19937 generator(seq);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    Eigen::MatrixXd U(n, d);
    return U.unaryExpr([&](double) { return distribution(generator); });
}

//! simulates from the multivariate uniform distribution
//!
//! @param n number of observations.
//! @param d dimension.
//! @param qrng if true, quasi-numbers are generated.
//! @param seeds seeds of the random number generator; if empty (default),
//!   the random number generator is seeded randomly.
//!
//! If `qrng = TRUE`, generalized Halton sequences (see `ghalton()`) are used
//! for `d <= 300` and Sobol sequences otherwise (see `sobol()`).
//!
//! @return An \f$ n \times d \f$ matrix of independent
//! \f$ \mathrm{U}[0, 1] \f$ random variables.
inline Eigen::MatrixXd simulate_uniform(const size_t& n, const size_t& d,
                                        bool qrng, std::vector<int> seeds)
{
    Eigen::MatrixXd u(n, d);
    if (qrng) {
        if (d > 300) {
            u = tools_stats::sobol(n, d, seeds);
        } else {
            u = tools_stats::ghalton(n, d, seeds);
        }
    } else {
        u = tools_stats::simulate_uniform(n, d, seeds);
    }
    return u;
}

//! @brief applies the empirical probability integral transform to a data matrix.
//!
//! Gives pseudo-observations from the copula by applying the empirical
//! distribution function (scaled by n + 1) to each margin/column.
//!
//! @param x a matrix of real numbers.
//! @param ties_method indicates how to treat ties; same as in R, see
//! https://stat.ethz.ch/R-manual/R-devel/library/base/html/rank.html.
//! @return Psuedo-observations of the copula, i.e. F_X(X) (column-wise)
inline Eigen::MatrixXd
to_pseudo_obs(Eigen::MatrixXd x, std::string ties_method)
{
    for (int j = 0; j < x.cols(); ++j)
        x.col(j) = to_pseudo_obs_1d(static_cast<Eigen::VectorXd>(x.col(j)),
                                    ties_method);

    return x;
}

//! @brief applies the empirical probability integral transform to a data vector.
//!
//! Gives pseudo-observations from the copula by applying the empirical
//! distribution function (scaled by n + 1) to each margin/column.
//!
//! @param x a vector of real numbers.
//! @param ties_method indicates how to treat ties; same as in R, see
//! https://stat.ethz.ch/R-manual/R-devel/library/base/html/rank.html.
//! @return Psuedo-observations of the copula, i.e. F_X(X) (column-wise)
inline Eigen::VectorXd
to_pseudo_obs_1d(Eigen::VectorXd x, std::string ties_method)
{
    size_t n = x.size();
    std::vector<double> xvec(x.data(), x.data() + n);
    auto order = tools_stl::get_order(xvec);
    if (ties_method == "first") {
        for (auto i : order)
            x[order[i]] = static_cast<double>(i + 1);
    } else if (ties_method == "average") {
        for (size_t i = 0, reps; i < n; i += reps) {
            // find replications
            reps = 1;
            while ((i + reps < n) && (x[order[i]] == x[order[i + reps]]))
                ++reps;
            // assign average rank of the tied values
            for (size_t k = 0; k < reps; ++k)
                x[order[i + k]] = i + 1 + (reps - 1) / 2.0;
        }
    } else if (ties_method == "random") {
        // set up random number generator
        std::random_device rd;
        std::default_random_engine gen(rd());
        auto sim = [&](int m) {
            std::uniform_int_distribution<> distr(0, m - 1);
            return distr(gen);
        };
        for (size_t i = 0, reps; i < n; i += reps) {
            // find replications
            reps = 1;
            while ((i + reps < n) && (x[order[i]] == x[order[i + reps]]))
                ++reps;
            // assign random rank between ties
            std::vector <size_t> rvals(reps);
            std::iota(rvals.begin(), rvals.end(), 0);  // 0, 1, 2, ...
            std::random_shuffle(rvals.begin(), rvals.end(), sim);
            for (size_t k = 0; k < reps; ++k)
                x[order[i + k]] = static_cast<double>(i + 1 + rvals[k]);
        }
    } else {
        std::stringstream msg;
        msg << "unknown ties method (" << ties_method << ")";
        throw std::runtime_error(msg.str().c_str());
    }

    return x / (x.size() + 1.0);
}

//! window smoother
inline Eigen::VectorXd win(const Eigen::VectorXd &x, size_t wl = 5)
{
    size_t n = x.size();
    Eigen::VectorXd xx = Eigen::VectorXd::Zero(n + 2 * wl);
    Eigen::VectorXd yy = Eigen::VectorXd::Zero(n + 2 * wl);
    xx.block(2 * wl, 0, n, 1) = x;
    yy.block(0, 0, 2 * wl + 1, 1) = Eigen::VectorXd::Ones(2 * wl + 1);

    Eigen::FFT<double> fft;
    Eigen::VectorXcd tmp1 = fft.fwd(xx);
    Eigen::VectorXcd tmp2 = fft.fwd(yy);
    tmp2 = tmp2.conjugate();
    tmp1 = tmp1.cwiseProduct(tmp2);
    tmp2 = fft.inv(tmp1);
    Eigen::VectorXd result = tmp2.real().block(wl, 0, n, 1);
    result /= 2.0 * static_cast<double>(wl) + 1.0;
    result.block(0, 0, wl, 1) = Eigen::VectorXd::Constant(wl, result(wl));
    result.block(n - wl, 0, wl, 1) = Eigen::VectorXd::Constant(wl, result(n - wl- 1));
    return result;
}

//! helper routine for ace (In R, this would be win(x[ind], wl)[ranks])
inline Eigen::VectorXd cef(const Eigen::VectorXd &x,
                    const Eigen::Matrix<size_t, Eigen::Dynamic, 1> &ind,
                    const Eigen::Matrix<size_t, Eigen::Dynamic, 1> &ranks,
                    size_t wl = 5)
{
    size_t n_ind = ind.size();
    Eigen::VectorXd cey(n_ind);
    for (size_t i = 0; i < n_ind; i++) {
        cey(i) = x(ind(i));
    }
    cey = win(cey, wl);

    size_t n_ranks = ranks.size();
    Eigen::VectorXd result(n_ranks);
    for (size_t i = 0; i < n_ranks; i++) {
        result(i) = cey(ranks(i));
    }
    return result;
}

//! alternating conditional expectation algorithm
inline Eigen::Matrix<double, Eigen::Dynamic, 2> ace(
        const Eigen::Matrix<double, Eigen::Dynamic, 2>& data, // data
        const Eigen::VectorXd& weights = Eigen::VectorXd(), // weights
        size_t wl = 0, // window length for the smoother
        size_t outer_iter_max = 100, // max number of outer iterations
        size_t inner_iter_max = 10, // max number of inner iterations
        double outer_abs_tol = 2e-15, // outer stopping criterion
        double inner_abs_tol = 1e-4) // inner stopping criterion
{
    // sample size and memory allocation for the outer/inner loops
    size_t n = data.rows();
    Eigen::VectorXd tmp(n);

    size_t nw = weights.size();
    Eigen::VectorXd w(n);
    if (nw == 0) {
        w = Eigen::VectorXd::Ones(n);
    } else {
        if (nw != n) {
            throw std::runtime_error("weights should have a length equal to "
                                         "the number of rows in data");
        }
        w = weights;
    }

    // default window size
    double n_dbl = static_cast<double>(n);
    if (wl == 0) {
        wl = static_cast<size_t>(std::ceil(n_dbl/5));
    }

    // assign order/ranks to ind/ranks
    Eigen::Matrix<size_t, Eigen::Dynamic, 2> ind(n, 2);
    Eigen::Matrix<size_t, Eigen::Dynamic, 2> ranks(n, 2);
    for (size_t i = 0; i < 2; i++) {
        std::vector<double> xvec(data.data() + n * i, data.data() + n * (i + 1));
        auto order = tools_stl::get_order(xvec);
        for (auto j : order) {
            ind(j, i) = order[j];
            ranks(order[j], i) = j;
        }
    }

    // initialize output
    Eigen::Matrix<double, Eigen::Dynamic, 2> phi = ranks.cast<double>();
    phi.array() -= (n_dbl - 1.0)/2.0 - 1.0;
    phi /= std::sqrt(n_dbl * (n_dbl - 1.0) / 12.0);
    if (nw > 0) {
        phi.col(0) = phi.col(0).cwiseProduct(w);
        phi.col(1) = phi.col(1).cwiseProduct(w);
    }

    // initialize variables for the outer loop
    size_t outer_iter = 1;
    double outer_eps = 1.0;
    double outer_abs_err = 1.0;

    // outer loop (expectation of the first variable given the second)
    while (outer_iter <= outer_iter_max && outer_abs_err > outer_abs_tol) {
        // initialize variables for the inner loop
        size_t inner_iter = 1;
        double inner_eps = 1.0;
        double inner_abs_err = 1.0;

        // inner loop (expectation of the second variable given the first)
        while(inner_iter <= inner_iter_max && inner_abs_err > inner_abs_tol) {
            // conditional expectation
            phi.col(1) = cef(phi.col(0).cwiseProduct(w) , ind.col(1),
                             ranks.col(1), wl);


            // center and standardize
            double m1 = phi.col(1).sum() / n_dbl;
            phi.col(1).array() -= m1;
            double s1 = std::sqrt(phi.col(1).cwiseAbs2().sum() / (n_dbl - 1));
            phi.col(1) /= s1;

            // compute error and increase step
            inner_abs_err = inner_eps;
            tmp = (phi.col(1) - phi.col(0));
            inner_eps = tmp.cwiseAbs2().sum() / n_dbl;
            inner_abs_err = std::fabs(inner_abs_err - inner_eps);
            inner_iter = inner_iter + 1;
        }

        // conditional expectation
        phi.col(0) = cef(phi.col(1).cwiseProduct(w) , ind.col(0),
                         ranks.col(0), wl);

        // center and standardize
        double m0 = phi.col(0).sum() / n_dbl;
        phi.col(0).array() -= m0;
        double s0 = std::sqrt(phi.col(0).cwiseAbs2().sum() / (n_dbl - 1));
        phi.col(0) /= s0;

        // compute error and increase step
        outer_abs_err = outer_eps;
        tmp = (phi.col(1) - phi.col(0));
        outer_eps = tmp.cwiseAbs2().sum() / n_dbl;
        outer_abs_err = std::fabs(outer_abs_err - outer_eps);
        outer_iter = outer_iter + 1;
    }

    // return result
    return phi;
}

//! calculates the pairwise maximum correlation coefficient.
inline double pairwise_mcor(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &x,
    const Eigen::VectorXd &weights)
{
    Eigen::Matrix<double, Eigen::Dynamic, 2> phi = ace(x, weights);
    return wdm::wdm(phi, "cor", weights)(0, 1);
}
//! @}

//! @brief simulates from the multivariate Generalized Halton Sequence.
//!
//! For more information on Generalized Halton Sequence, see
//! Faure, H., Lemieux, C. (2009). Generalized Halton Sequences in 2008:
//! A Comparative Study. ACM-TOMACS 19(4), Article 15.
//!
//! @param n number of observations.
//! @param d dimension.
//! @param seeds seeds to scramble the quasi-random numbers; if empty (default),
//!   the quasi-random number generator is seeded randomly.
//!
//! @return An \f$ n \times d \f$ matrix of quasi-random
//! \f$ \mathrm{U}[0, 1] \f$ variables.
inline Eigen::MatrixXd ghalton(const size_t& n, const size_t& d,
                               std::vector<int> seeds)
{

    Eigen::MatrixXd res(d, n);

    // Coefficients of the shift
    Eigen::MatrixXi shcoeff(d, 32);
    Eigen::VectorXi base = tools_ghalton::primes.block(0, 0, d, 1);
    Eigen::MatrixXd u = Eigen::VectorXd::Zero(d, 1);
    auto U = simulate_uniform(d, 32, seeds);
    for (int k = 31; k >= 0; k--) {
        shcoeff.col(k) = (base.cast<double>()).cwiseProduct(
            U.block(0, k, d, 1)).cast<int>();
        u = (u + shcoeff.col(k).cast<double>()).cwiseQuotient(
            base.cast<double>());
    }
    res.block(0, 0, d, 1) = u;

    Eigen::VectorXi perm = tools_ghalton::permTN2.block(0, 0, d, 1);
    Eigen::MatrixXi coeff(d, 32);
    Eigen::VectorXi tmp(d);
    int k;
    auto mod = [](const int &u1, const int &u2) { return u1 % u2; };
    for (size_t i = 1; i < n; i++) {

        // Find i in the prime base
        tmp = Eigen::VectorXi::Constant(d, static_cast<int>(i));
        coeff = Eigen::MatrixXi::Zero(d, 32);
        k = 0;
        while ((tmp.maxCoeff() > 0) && (k < 32)) {
            coeff.col(k) = tmp.binaryExpr(base, mod);
            tmp = tmp.cwiseQuotient(base);
            k++;
        }

        u = Eigen::VectorXd::Zero(d);
        k = 31;
        while (k >= 0) {
            tmp = perm.cwiseProduct(coeff.col(k)) + shcoeff.col(k);
            u = u + tmp.binaryExpr(base, mod).cast<double>();
            u = u.cwiseQuotient(base.cast<double>());
            k--;
        }
        res.block(0, i, d, 1) = u;
    }

    return res.transpose();
}

//! @brief simulates from the multivariate Sobol sequence.
//!
//! For more information on the Sobol sequence, see S. Joe and F. Y. Kuo
//! (2008), Constructing Sobol  sequences with better two-dimensional
//! projections, SIAM J. Sci. Comput. 30, 2635–2654.
//!
//! @param n number of observations.
//! @param d dimension.
//! @param seeds seeds to scramble the quasi-random numbers; if empty (default),
//!   the quasi-random number generator is seeded randomly.
//!
//! @return An \f$ n \times d \f$ matrix of quasi-random
//! \f$ \mathrm{U}[0, 1] \f$ variables.
inline Eigen::MatrixXd sobol(const size_t& n, const size_t& d,
                             std::vector<int> seeds)
{

    // output matrix
    Eigen::MatrixXd output = Eigen::MatrixXd::Zero(n, d);

    // L = max number of bits needed
    size_t L = static_cast<size_t>(std::ceil(log(static_cast<double>(n))/log(2.0)));

    // Vector of scrambling factors
    Eigen::MatrixXd scrambling = simulate_uniform(d, 1, seeds);

    // C(i) = index from the right of the first zero bit of i + 1
    Eigen::Matrix<size_t, Eigen::Dynamic, 1> C(n);
    C(0) = 1;
    for (size_t i = 1; i < n ; i++) {
        C(i) = 1;
        size_t value = i;
        while (value & 1) {
            value >>= 1;
            C(i)++;
        }
    }

    // Compute the first dimension

    // Compute direction numbers scaled by pow(2,32)
    Eigen::Matrix<size_t, Eigen::Dynamic, 1> V(L);
    for (size_t i = 0; i < L; i++) {
        V(i) = static_cast<size_t>(std::pow(2, 32 - (i + 1))); // all m's = 1
    }

    // Evalulate X scaled by pow(2,32)
    Eigen::Matrix<size_t, Eigen::Dynamic, 1> X(n);
    X(0) = scrambling(0) * std::pow(2.0, 32);
    for (size_t i = 1; i < n; i++) {
        X(i) = X(i - 1) ^ V(C(i - 1) - 1);
    }
    output.block(0, 0, n, 1) = X.cast<double> ();

    // Compute the remaining dimensions
    for (size_t j = 0; j < d - 1; j++) {

        // Get parameters from static vectors
        size_t a = tools_sobol::a_sobol[j];
        size_t s = tools_sobol::s_sobol[j];

        Eigen::Map<Eigen::Matrix<size_t,Eigen::Dynamic,1> >
            m(tools_sobol::minit_sobol[j], s);

        // Compute direction numbers scaled by pow(2,32)
        for (size_t i = 0; i < std::min(L, s); i++)
            V(i) = m(i) << (32 - (i + 1));

        if (L > s) {
            for (size_t i = s; i < L; i++) {
                V(i) = V(i - s) ^ (V(i - s) >> s);
                for (size_t k = 0; k < s - 1; k++)
                    V(i) ^= (((a >> (s - 2 - k)) & 1) * V(i - k - 1));
            }
        }

        // Evalulate X
        X(0) = scrambling(j + 1) * std::pow(2.0, 32);
        for (size_t i = 1; i < n; i++)
            X(i) = X(i - 1) ^ V(C(i - 1) - 1);
        output.block(0, j + 1, n, 1) = X.cast <double> ();
    }

    // Scale output by pow(2,32)
    output /= std::pow(2.0, 32);

    return output;
}

//! @brief computes bivariate t probabilities.
//!
//! Based on the method described by
//! Dunnett, C.W. and M. Sobel, (1954),
//! A bivariate generalization of Student's t-distribution
//! with tables for certain special cases,
//! Biometrika 41, pp. 153-169. Translated from the Fortran routines of
//! Alan Genz (www.math.wsu.edu/faculty/genz/software/fort77/mvtdstpack.f).
//!
//! @param z an \f$ n \times 2 \f$ matrix of evaluation points.
//! @param nu number of degrees of freedom.
//! @param rho correlation.
//!
//! @return An \f$ n \times 1 \f$ vector of probabilities.
inline Eigen::VectorXd pbvt(const Eigen::Matrix<double, Eigen::Dynamic, 2> &z,
                            int nu, double rho)
{
    double snu = sqrt(static_cast<double>(nu));
    double ors = 1 - pow(rho, 2.0);

    auto f = [snu, nu, ors, rho](double h, double k) {

        double d1, d2, d3, hkn, hpk, bvt, gmph, gmpk, hkrn,
            qhrk, xnkh, xnhk, btnckh, btnchk, btpdkh, btpdhk;
        int hs, ks;

        double hrk = h - rho * k;
        double krh = k - rho * h;
        if (std::fabs(hrk) + ors > 0.) {
            /* Computing 2nd power */
            d1 = hrk;
            /* Computing 2nd power */
            d2 = hrk;
            /* Computing 2nd power */
            d3 = k;
            xnhk = d1 * d1 / (d2 * d2 + ors * (nu + d3 * d3));
            /* Computing 2nd power */
            d1 = krh;
            /* Computing 2nd power */
            d2 = krh;
            /* Computing 2nd power */
            d3 = h;
            xnkh = d1 * d1 / (d2 * d2 + ors * (nu + d3 * d3));
        } else {
            xnhk = 0.;
            xnkh = 0.;
        }
        d1 = h - rho * k;
        hs = static_cast<int>(d1 >= 0 ? 1 : -1);//d_sign(&c_b91, &d1);
        d1 = k - rho * h;
        ks = static_cast<int>(d1 >= 0 ? 1 : -1);
        if (nu % 2 == 0) {
            bvt = atan2(sqrt(ors), -rho) / 6.2831853071795862;
            /* Computing 2nd power */
            d1 = h;
            gmph = h / sqrt((nu + d1 * d1) * 16);
            /* Computing 2nd power */
            d1 = k;
            gmpk = k / sqrt((nu + d1 * d1) * 16);
            btnckh = atan2(sqrt(xnkh), sqrt(1 - xnkh)) * 2 /
                     3.14159265358979323844;
            btpdkh = sqrt(xnkh * (1 - xnkh)) * 2 / 3.14159265358979323844;
            btnchk = atan2(sqrt(xnhk), sqrt(1 - xnhk)) * 2 /
                     3.14159265358979323844;
            btpdhk = sqrt(xnhk * (1 - xnhk)) * 2 / 3.14159265358979323844;
            size_t i1 = static_cast<size_t>(nu / 2);
            for (size_t j = 1; j <= i1; ++j) {
                bvt += gmph * (ks * btnckh + 1);
                bvt += gmpk * (hs * btnchk + 1);
                btnckh += btpdkh;
                btpdkh = (j << 1) * btpdkh * (1 - xnkh) / ((j << 1) + 1);
                btnchk += btpdhk;
                btpdhk = (j << 1) * btpdhk * (1 - xnhk) / ((j << 1) + 1);
                /* Computing 2nd power */
                d1 = h;
                gmph = gmph * ((j << 1) - 1) / ((j << 1) * (d1 * d1 / nu + 1)
                );
                /* Computing 2nd power */
                d1 = k;
                gmpk = gmpk * ((j << 1) - 1) / ((j << 1) * (d1 * d1 / nu + 1)
                );
            }
        } else {
            /* Computing 2nd power */
            d1 = h;
            /* Computing 2nd power */
            d2 = k;
            qhrk = sqrt(
                d1 * d1 + d2 * d2 - rho * 2 * h * k + nu * ors);
            hkrn = h * k + rho * nu;
            hkn = h * k - nu;
            hpk = h + k;
            bvt =
                atan2(-snu * (hkn * qhrk + hpk * hkrn), hkn * hkrn - nu * hpk *
                                                                     qhrk) /
                6.2831853071795862;
            if (bvt < -1e-15) {
                bvt += 1;
            }
            /* Computing 2nd power */
            d1 = h;
            gmph = h / (snu * 6.2831853071795862 * (d1 * d1 / nu + 1));
            /* Computing 2nd power */
            d1 = k;
            gmpk = k / (snu * 6.2831853071795862 * (d1 * d1 / nu + 1));
            btnckh = sqrt(xnkh);
            btpdkh = btnckh;
            btnchk = sqrt(xnhk);
            btpdhk = btnchk;
            size_t i1 = static_cast<size_t>((nu - 1) / 2);
            for (size_t j = 1; j <= i1; ++j) {
                bvt += gmph * (ks * btnckh + 1);
                bvt += gmpk * (hs * btnchk + 1);
                btpdkh = ((j << 1) - 1) * btpdkh * (1 - xnkh) / (j << 1);
                btnckh += btpdkh;
                btpdhk = ((j << 1) - 1) * btpdhk * (1 - xnhk) / (j << 1);
                btnchk += btpdhk;
                /* Computing 2nd power */
                d1 = h;
                gmph = (j << 1) * gmph / (((j << 1) + 1) * (d1 * d1 / nu + 1)
                );
                /* Computing 2nd power */
                d1 = k;
                gmpk = (j << 1) * gmpk / (((j << 1) + 1) * (d1 * d1 / nu + 1)
                );
            }
        }
        return bvt;
    };

    return tools_eigen::binaryExpr_or_nan(z, f);
}

//! @brief compute bivariate normal probabilities.
//!
//! A function for computing bivariate normal probabilities;
//! developed using Drezner, Z. and Wesolowsky, G. O. (1989),
//! On the Computation of the Bivariate Normal Integral,
//! J. Stat. Comput. Simul.. 35 pp. 101-107.
//! with extensive modications for double precisions by
//! Alan Genz and Yihong Ge. Translated from the Fortran routines of
//! Alan Genz (www.math.wsu.edu/faculty/genz/software/fort77/mvtdstpack.f).
//!
//! @param z an \f$ n \times 2 \f$ matrix of evaluation points.
//! @param rho correlation.
//!
//! @return An \f$ n \times 1 \f$ vector of probabilities.
inline Eigen::VectorXd
pbvnorm(const Eigen::Matrix<double, Eigen::Dynamic, 2> &z,
        double rho)
{

    static struct
    {
        double e_1[3];
        double fill_2[7];
        double e_3[6];
        double fill_4[4];
        double e_5[10];
    } equiv_112 = {{.1713244923791705,  .3607615730481384, .4679139345726904},
                   {0},
                   {.04717533638651177, .1069393259953183,
                                                           .1600783285433464,  .2031674267230659,
                       .2334925365383547, .2491470458134029},
                   {0},
                   {.01761400713915212, .04060142980038694,
                                                           .06267204833410906, .08327674157670475,
                       .1019301198172404, .1181945319615184,
                       .1316886384491766, .1420961093183821,
                       .1491729864726037, .1527533871307259}};
    auto w = reinterpret_cast<double *>(&equiv_112);

    static struct
    {
        double e_1[3];
        double fill_2[7];
        double e_3[6];
        double fill_4[4];
        double e_5[10];
    } equiv_113 = {{-.9324695142031522, -.6612093864662647, -.238619186083197},
                   {0},
                   {-.9815606342467191, -.904117256370475,  -.769902674194305,
                                                                                -.5873179542866171, -.3678314989981802, -.1252334085114692},
                   {0},
                   {-.9931285991850949, -.9639719272779138,
                                                            -.9122344282513259, -.8391169718222188,
                                                                                                    -.7463319064601508, -.636053680726515,
                       -.5108670019508271, -.3737060887154196,
                       -.2277858511416451, -.07652652113349733}};
    auto x = reinterpret_cast<double *>(&equiv_113);


    boost::math::normal dist;
    auto phi = [&dist](double y) { return boost::math::cdf(dist, y); };

    size_t lg, ng;
    if (std::fabs(rho) < .3f) {
        ng = 1;
        lg = 3;
    } else if (std::fabs(rho) < .75f) {
        ng = 2;
        lg = 6;
    } else {
        ng = 3;
        lg = 10;
    }

    auto f = [ng, lg, rho, x, w, phi](double h, double k) {
        size_t i1;
        double a, b, c, d, d1, d2, as, bs, hk, hs, sn, rs, xs, bvn, asr;
        h = -h;
        k = -k;
        hk = h * k;
        bvn = 0.;
        if (std::fabs(rho) < .925f) {
            hs = (h * h + k * k) / 2;
            asr = asin(rho);
            i1 = lg;
            for (size_t i = 1; i <= i1; ++i) {
                sn = sin(asr * (x[i + ng * 10 - 11] + 1) / 2);
                bvn +=
                    w[i + ng * 10 - 11] * exp((sn * hk - hs) / (1 - sn * sn));
                sn = sin(asr * (-x[i + ng * 10 - 11] + 1) / 2);
                bvn +=
                    w[i + ng * 10 - 11] * exp((sn * hk - hs) / (1 - sn * sn));
            }
            d1 = -h;
            d2 = -k;
            bvn = bvn * asr / 12.566370614359172 + phi(d1) * phi(d2);
        } else {
            if (rho < 0.) {
                k = -k;
                hk = -hk;
            }
            if (std::fabs(rho) < 1.) {
                as = (1 - rho) * (rho + 1);
                a = sqrt(as);
                /* Computing 2nd power */
                d1 = h - k;
                bs = d1 * d1;
                c = (4 - hk) / 8;
                d = (12 - hk) / 16;
                bvn = a * exp(-(bs / as + hk) / 2) * (1 - c * (bs - as) * (1 -
                                                                           d *
                                                                           bs /
                                                                           5) /
                                                          3 +
                                                      c * d * as * as / 5);
                if (hk > -160.) {
                    b = sqrt(bs);
                    d1 = -b / a;
                    bvn -= exp(-hk / 2) * sqrt(6.283185307179586) * phi(d1)
                           * b * (1 - c * bs * (1 - d * bs / 5) / 3);
                }
                a /= 2;
                i1 = lg;
                for (size_t i = 1; i <= i1; ++i) {
                    /* Computing 2nd power */
                    d1 = a * (x[i + ng * 10 - 11] + 1);
                    xs = d1 * d1;
                    rs = sqrt(1 - xs);
                    bvn += a * w[i + ng * 10 - 11] * (exp(-bs / (xs * 2) - hk /
                                                                           (rs +
                                                                            1)) /
                                                      rs -
                                                      exp(-(bs / xs + hk) / 2) *
                                                      (c * xs
                                                       * (d * xs + 1) + 1));
                    /* Computing 2nd power */
                    d1 = -x[i + ng * 10 - 11] + 1;
                    xs = as * (d1 * d1) / 4;
                    rs = sqrt(1 - xs);
                    /* Computing 2nd power */
                    d1 = rs + 1;
                    bvn += a * w[i + ng * 10 - 11] * exp(-(bs / xs + hk) / 2) *
                           (exp(-hk * xs / (d1 * d1 * 2)) / rs - (c * xs *
                                                                  (d * xs +
                                                                   1) + 1));
                }
                bvn = -bvn / 6.283185307179586;
            }
            if (rho > 0.) {
                d1 = -std::max(h, k);
                bvn += phi(d1);
            } else {
                bvn = -bvn;
                if (k > h) {
                    if (h < 0.) {
                        bvn = bvn + phi(k) - phi(h);
                    } else {
                        d1 = -h;
                        d2 = -k;
                        bvn = bvn + phi(d1) - phi(d2);
                    }
                }
            }
        }
        return bvn;
    };

    return tools_eigen::binaryExpr_or_nan(z, f);
}
}

}
