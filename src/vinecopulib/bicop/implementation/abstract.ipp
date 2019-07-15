// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <stdexcept>

#include <vinecopulib/misc/tools_eigen.hpp>
#include <vinecopulib/bicop/bb1.hpp>
#include <vinecopulib/bicop/bb6.hpp>
#include <vinecopulib/bicop/bb7.hpp>
#include <vinecopulib/bicop/bb8.hpp>
#include <vinecopulib/bicop/clayton.hpp>
#include <vinecopulib/bicop/frank.hpp>
#include <vinecopulib/bicop/gaussian.hpp>
#include <vinecopulib/bicop/gumbel.hpp>
#include <vinecopulib/bicop/indep.hpp>
#include <vinecopulib/bicop/joe.hpp>
#include <vinecopulib/bicop/student.hpp>
#include <vinecopulib/bicop/tll.hpp>

namespace vinecopulib {

//! virtual destructor
inline AbstractBicop::~AbstractBicop()
{
}

//! Create a bivariate copula using the default contructor
//!
//! @param family the copula family.
//! @param parameters the copula parameters (optional, must be compatible
//!     with family).
//! @return A pointer to an object that inherits from AbstractBicop.
//! @{
inline BicopPtr AbstractBicop::create(BicopFamily family,
                                      const Eigen::MatrixXd &parameters)
{
    BicopPtr new_bicop;
    switch (family) {
        case BicopFamily::indep:
            new_bicop = BicopPtr(new IndepBicop());
            break;
        case BicopFamily::gaussian:
            new_bicop = BicopPtr(new GaussianBicop());
            break;
        case BicopFamily::student:
            new_bicop = BicopPtr(new StudentBicop());
            break;
        case BicopFamily::clayton:
            new_bicop = BicopPtr(new ClaytonBicop());
            break;
        case BicopFamily::gumbel:
            new_bicop = BicopPtr(new GumbelBicop());
            break;
        case BicopFamily::frank:
            new_bicop = BicopPtr(new FrankBicop());
            break;
        case BicopFamily::joe:
            new_bicop = BicopPtr(new JoeBicop());
            break;
        case BicopFamily::bb1:
            new_bicop = BicopPtr(new Bb1Bicop());
            break;
        case BicopFamily::bb6:
            new_bicop = BicopPtr(new Bb6Bicop());
            break;
        case BicopFamily::bb7:
            new_bicop = BicopPtr(new Bb7Bicop());
            break;
        case BicopFamily::bb8:
            new_bicop = BicopPtr(new Bb8Bicop());
            break;
        case BicopFamily::tll:
            new_bicop = BicopPtr(new TllBicop());
            break;

        default:
            throw std::runtime_error(std::string("Family not implemented"));
    }

    if (parameters.size() > 0) {
        new_bicop->set_parameters(parameters);
    }

    return new_bicop;
}

//!@}

inline Eigen::MatrixXd AbstractBicop::no_tau_to_parameters(const double &)
{
    throw std::runtime_error("Method not implemented for this family");
}

//! Getters and setters.
//! @{
inline BicopFamily AbstractBicop::get_family() const
{
    return family_;
}

inline std::string AbstractBicop::get_family_name() const
{
    return vinecopulib::get_family_name(family_);
}

inline double AbstractBicop::get_loglik() const
{
    return loglik_;
}

inline void AbstractBicop::set_loglik(const double loglik)
{
    loglik_ = loglik;
}
//! @}


//! evaluates the pdf, but truncates it's value by DBL_MIN and DBL_MAX.
//! @param u matrix of evaluation points.
inline Eigen::VectorXd AbstractBicop::pdf(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u)
{
    auto trim = [] (const double &x) {
        return std::min(DBL_MAX, std::max(x, DBL_MIN));
    };
    return tools_eigen::unaryExpr_or_nan(pdf_raw(u), trim);
}

//! evaluates the log-likelihood.
//! @param u data matrix.
//! @param weights optional weights for each observation.
inline double AbstractBicop::loglik(
    const Eigen::Matrix<double, Eigen::Dynamic, 2> &u,
    const Eigen::VectorXd weights)
{
    Eigen::MatrixXd log_pdf = this->pdf(u).array().log();
    if (weights.size() > 0) {
        log_pdf = log_pdf.cwiseProduct(weights);
    }
    tools_eigen::remove_nans(log_pdf);
    return log_pdf.sum();
}

//! Numerical inversion of h-functions
//!
//! These are generic functions to invert the hfunctions numerically.
//! They can be used in derived classes to define \c hinv1 and \c hinv2.
//!
//! @param u \f$m \times 2\f$ matrix of evaluation points.
//! @return The numerical inverse of h-functions.
//! @{
inline Eigen::VectorXd
AbstractBicop::hinv1_num(const Eigen::Matrix<double, Eigen::Dynamic, 2> &u)
{
    Eigen::Matrix<double, Eigen::Dynamic, 2> u_new = u;
    auto h1 = [&](const Eigen::VectorXd &v) {
        u_new.col(1) = v;
        return hfunc1(u_new);
    };

    return tools_eigen::invert_f(u.col(1), h1);
}

inline Eigen::VectorXd
AbstractBicop::hinv2_num(const Eigen::Matrix<double, Eigen::Dynamic, 2> &u)
{
    Eigen::Matrix<double, Eigen::Dynamic, 2> u_new = u;
    auto h1 = [&](const Eigen::VectorXd &x) {
        u_new.col(0) = x;
        return hfunc2(u_new);
    };

    return tools_eigen::invert_f(u.col(0), h1);
}
//! @}
}
