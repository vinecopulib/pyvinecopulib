// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_bobyqa.hpp>
#include <vinecopulib/misc/tools_optimization.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <wdm/eigen.hpp>

namespace vinecopulib {
inline Eigen::MatrixXd
ParBicop::get_parameters() const
{
  return parameters_;
}

inline Eigen::MatrixXd
ParBicop::get_parameters_lower_bounds() const
{
  return parameters_lower_bounds_;
}

inline Eigen::MatrixXd
ParBicop::get_parameters_upper_bounds() const
{
  return parameters_upper_bounds_;
}

inline void
ParBicop::set_parameters(const Eigen::MatrixXd& parameters)
{
  check_parameters(parameters);
  parameters_ = parameters;
}

inline void
ParBicop::flip()
{
  // Most parametric families can be flipped by changing the rotation.
  // This is done in Bicop::flip() directly. All other families need to
  // override this method.
}

// calculate number of parameters
inline double
ParBicop::calculate_npars()
{
  // indepence copula has no parameters
  if (family_ == BicopFamily::indep) {
    return 0.0;
  }
  // otherwise, return length of parameter vector
  return static_cast<double>(parameters_.size());
}

// fit
inline void
ParBicop::fit(const Eigen::Matrix<double, Eigen::Dynamic, 2>& data,
              std::string method,
              double,
              const Eigen::VectorXd& weights)
{
  // for independence copula we don't have to do anything
  if (family_ == BicopFamily::indep) {
    set_loglik(0.0);
    return;
  }

  check_fit_method(method);
  double tau = wdm::wdm(data, "tau", weights)(0, 1);

  // for method itau and one-parameter families we don't need to optimize
  int npars = static_cast<int>(calculate_npars()) - (method == "itau");
  if (npars == 0) {
    set_parameters(tau_to_parameters(tau));
    set_loglik(loglik(data, weights));
    return;
  }

  // Set bounds and starting values
  auto lb = get_parameters_lower_bounds();
  auto ub = get_parameters_upper_bounds();
  adjust_parameters_bounds(lb, ub, tau, method);
  auto initial_parameters = get_start_parameters(winsorize_tau(tau));

  // find (pseudo-) mle
  std::function<double(const Eigen::VectorXd&)> objective;
  if (method == "mle") {
    objective = [&data, &weights, this](const Eigen::VectorXd& pars) {
      this->set_parameters(pars);
      return this->loglik(data, weights);
    };
  } else {
    // profile likelihood
    set_parameters(initial_parameters);
    initial_parameters(0) = initial_parameters(1);
    initial_parameters.conservativeResize(1);
    objective = [&data, &weights, this](const Eigen::VectorXd& pars) {
      Eigen::VectorXd newpars(2);
      newpars(0) = this->get_parameters()(0);
      newpars(1) = pars(0);
      this->set_parameters(newpars);
      return this->loglik(data, weights);
    };
  }

  tools_optimization::Optimizer optimizer;
  auto newpars = optimizer.optimize(initial_parameters, lb, ub, objective);

  // check if fit is reasonable, otherwise increase search interval
  // and refit
  if (tools_stl::is_member(family_, bicop_families::one_par) &&
      (optimizer.get_objective_max() < -0.1)) {
    newpars = optimizer.optimize(initial_parameters,
                                 get_parameters_lower_bounds(),
                                 get_parameters_upper_bounds(),
                                 objective);
  }

  // finalize fitted model
  if (method == "itau") {
    // only the second parameter has been optimized
    newpars.conservativeResize(2);
    std::swap(newpars(0), newpars(1));
    newpars(0) = get_parameters()(0);
  }

  set_parameters(newpars);
  set_loglik(optimizer.get_objective_max());
}

//! ensures that starting values are sufficiently separated from bounds
//! @param tau kendall's tau
inline double
ParBicop::winsorize_tau(double tau) const
{
  double sign = 1.0;
  if (tau < 0) {
    sign = -1.0;
  }
  if (std::abs(tau) < 0.01) {
    tau = 0.01 * sign;
  } else if (std::abs(tau) > 0.9) {
    tau = 0.9 * sign;
  }
  return tau;
}

//! adjusts parameter bounds for better search intervals.
inline void
ParBicop::adjust_parameters_bounds(Eigen::MatrixXd& lb,
                                   Eigen::MatrixXd& ub,
                                   const double& tau,
                                   const std::string method)
{
  if (method == "itau") {
    // for pseudo mle, we can ignore the first parameter
    lb(0) = lb(1);
    ub(0) = ub(1);
    lb.conservativeResize(1, 1);
    ub.conservativeResize(1, 1);
    if (family_ == BicopFamily::student) {
      // the df parameter doesn't need to be estimated as accurately
      ub(0) = 15;
    }
  }

  // refine search interval for Brent algorithm
  if (tools_stl::is_member(family_, bicop_families::one_par)) {
    auto lb2 = lb;
    auto ub2 = ub;
    if (tools_stl::is_member(family_, bicop_families::rotationless)) {
      lb = tau_to_parameters(std::max(tau - 0.1, -0.99));
      ub = tau_to_parameters(std::min(tau + 0.1, 0.99));
    } else {
      lb = tau_to_parameters(std::max(std::fabs(tau) - 0.1, 1e-10));
      ub = tau_to_parameters(std::min(std::fabs(tau) + 0.1, 0.95));
    }
    // make sure that parameter bounds are respected
    lb = lb2.cwiseMax(lb);
    ub = ub2.cwiseMin(ub);
  }
}

//! Sanity checks
//! @{
inline void
ParBicop::check_parameters(const Eigen::MatrixXd& parameters)
{
  check_parameters_size(parameters);
  check_parameters_lower(parameters);
  check_parameters_upper(parameters);
}

inline void
ParBicop::check_parameters_size(const Eigen::MatrixXd& parameters)
{
  if (parameters.size() != parameters_.size()) {
    if (parameters.rows() != parameters_.rows()) {
      std::stringstream message;
      message << "parameters have has wrong number of rows "
              << "for " << get_family_name() << " copula; "
              << "expected: " << parameters_.rows() << ", "
              << "actual: " << parameters.rows() << std::endl;
      throw std::runtime_error(message.str().c_str());
    }
    if (parameters.cols() != parameters_.cols()) {
      std::stringstream message;
      message << "parameters have wrong number of columns "
              << "for " << get_family_name() << " copula; "
              << "expected: " << parameters_.cols() << ", "
              << "actual: " << parameters.cols() << std::endl;
      throw std::runtime_error(message.str().c_str());
    }
  }
}

inline void
ParBicop::check_parameters_lower(const Eigen::MatrixXd& parameters)
{
  if (parameters_lower_bounds_.size() > 0) {
    std::stringstream message;
    if ((parameters.array() < parameters_lower_bounds_.array()).any()) {
      message << "parameters exceed lower bound "
              << "for " << get_family_name() << " copula; " << std::endl
              << "bound:" << std::endl
              << parameters_lower_bounds_ << std::endl
              << "actual:" << std::endl
              << parameters << std::endl;
      throw std::runtime_error(message.str().c_str());
    }
  }
}

inline void
ParBicop::check_parameters_upper(const Eigen::MatrixXd& parameters)
{
  if (parameters_upper_bounds_.size() > 0) {
    std::stringstream message;
    if ((parameters.array() > parameters_upper_bounds_.array()).any()) {
      message << "parameters exceed upper bound "
              << "for " << get_family_name() << " copula; " << std::endl
              << "bound:" << std::endl
              << parameters_upper_bounds_ << std::endl
              << "actual:" << std::endl
              << parameters << std::endl;
      throw std::runtime_error(message.str().c_str());
    }
  }
}

inline void
ParBicop::check_fit_method(const std::string& method)
{
  if (!tools_stl::is_member(method, { "itau", "mle" })) {
    throw std::runtime_error("Method not implemented.");
  }

  if (method == "itau") {
    if (!tools_stl::is_member(family_, bicop_families::itau)) {
      throw std::runtime_error("itau method is not available for this family.");
    }
  }
}

//! @}
}
