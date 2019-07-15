// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp> // isnan
#include <boost/math/special_functions/log1p.hpp>
#include <vinecopulib/misc/tools_eigen.hpp>

namespace vinecopulib {
inline GumbelBicop::GumbelBicop()
{
    family_ = BicopFamily::gumbel;
    parameters_ = Eigen::VectorXd(1);
    parameters_lower_bounds_ = Eigen::VectorXd(1);
    parameters_upper_bounds_ = Eigen::VectorXd(1);
    parameters_ << 1;
    parameters_lower_bounds_ << 1;
    parameters_upper_bounds_ << 50;
}

inline double GumbelBicop::generator(const double &u)
{
    return std::pow(std::log(1 / u), this->parameters_(0));
}

inline double GumbelBicop::generator_inv(const double &u)
{
    return std::exp(-std::pow(u, 1 / this->parameters_(0)));
}

inline double GumbelBicop::generator_derivative(const double &u)
{
    double theta = double(this->parameters_(0));
    return std::pow(std::log(1 / u), theta - 1) * (-theta / u);
}

//inline double GumbelBicop::generator_derivative2(const double &u)
//{
//    double theta = double(this->parameters_(0));
//    return (theta - 1 - std::log(u)) * std::pow(std::log(1 / u), theta - 2) *
//           (theta / std::pow(u, 2));
//}

inline Eigen::VectorXd GumbelBicop::pdf_raw(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
)
{
    double theta = static_cast<double>(parameters_(0));
    double thetha1 = 1.0/theta;
    auto f = [theta,thetha1](const double &u1, const double &u2) {
        double t1 = std::pow(-std::log(u1),theta)+std::pow(-std::log(u2),theta);
        double temp = -std::pow(t1,thetha1)+(2*thetha1-2.0)*std::log(t1)
                      + (theta-1.0)*std::log(std::log(u1)*std::log(u2))
                      - std::log(u1*u2)
                      + boost::math::log1p((theta-1.0)*std::pow(t1,-thetha1));
        return std::exp(temp);
    };
    return tools_eigen::binaryExpr_or_nan(u, f);
}

inline Eigen::VectorXd GumbelBicop::hinv1(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u
)
{
    double theta = double(this->parameters_(0));
    double u1, u2;
    Eigen::VectorXd hinv = Eigen::VectorXd::Zero(u.rows());
    for (int j = 0; j < u.rows(); ++j) {
        u1 = u(j, 1);
        u2 = u(j, 0);
        if ((boost::math::isnan)(u1) | (boost::math::isnan)(u2)) {
            hinv(j) = std::numeric_limits<double>::quiet_NaN();
        } else {
            hinv(j) = qcondgum(&u1, &u2, &theta);
        }
    }

    return hinv;
}

inline Eigen::MatrixXd GumbelBicop::tau_to_parameters(const double &tau)
{
    return Eigen::VectorXd::Constant(1, 1.0 / (1 - std::fabs(tau)));
}

inline double
GumbelBicop::parameters_to_tau(const Eigen::MatrixXd &parameters)
{
    return (parameters(0) - 1) / parameters(0);
}

inline Eigen::VectorXd GumbelBicop::get_start_parameters(const double tau)
{
    Eigen::VectorXd par = tau_to_parameters(tau);
    par = par.cwiseMax(parameters_lower_bounds_);
    par = par.cwiseMin(parameters_upper_bounds_);
    return par;
}
}

// This is copy&paste from the VineCopula package
inline double qcondgum(double *q, double *u, double *de)
{
    double a, p, g, gp, z1, z2, con, de1, dif;
    double mxdif;
    int iter;

    p = 1 - *q;
    z1 = -log(*u);
    con = log(1. - p) - z1 + (1. - *de) * log(z1);
    de1 = *de - 1.;
    a = pow(2. * pow(z1, *de), 1. / (*de));
    mxdif = 1;
    iter = 0;
    dif = .1;  // needed in case first step leads to NaN
    while ((mxdif > 1.e-6) && (iter < 20)) {
        g = a + de1 * log(a) + con;
        gp = 1. + de1 / a;
        if ((boost::math::isnan)(g) ||
            (boost::math::isnan)(gp) ||
            (boost::math::isnan)(g / gp)) {
            // added for de>50
            dif /= -2.;
        } else {
            dif = g / gp;
        }
        a -= dif;
        iter++;
        int it = 0;
        while ((a <= z1) && (it < 20)) {
            dif /= 2.;
            a += dif;
            ++it;
        }
        mxdif = fabs(dif);
    }
    z2 = pow(pow(a, *de) - pow(z1, *de), 1. / (*de));
    return (exp(-z2));
}
