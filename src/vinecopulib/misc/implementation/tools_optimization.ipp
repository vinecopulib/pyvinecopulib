// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <boost/math/tools/minima.hpp>
#include <vinecopulib/bicop/parametric.hpp>
#include <vinecopulib/misc/tools_bobyqa.hpp>

namespace vinecopulib {

//! Utilities for numerical optimization (based on Bobyqa)
namespace tools_optimization {

//! creates an Optimizer using the default controls, see BobyqaControls.
//!
//! @param n_parameters Number of parameters to optimize
//! @param lower_bounds
//! @param upper_bounds
//! @param objective The optimizer's objective function
inline Optimizer::Optimizer()
  : controls_(BobyqaControls())
{}

//! set the optimizer's controls.
//!
//! @param initial_trust_region initial trust region.
//! @param final_trust_region final trust region.
//! @param maxeval maximal number of evaluations of the objective.
inline void
Optimizer::set_controls(double initial_trust_region,
                        double final_trust_region,
                        int maxeval)
{
  controls_ = BobyqaControls(initial_trust_region, final_trust_region, maxeval);
}

//! @brief solve the maximization problem.
//!
//! @param initial_parameters of starting values for the optimization
//!     algorithm.
//! @param lower_bounds lower bounds for the parameters.
//! @param upper_bounds upper bounds for the parameters.
//! @param the objective function to maximize.
//! @return the optimal parameters.
inline Eigen::VectorXd
Optimizer::optimize(const Eigen::VectorXd& initial_parameters,
                    const Eigen::VectorXd& lower_bounds,
                    const Eigen::VectorXd& upper_bounds,
                    std::function<double(const Eigen::VectorXd&)> objective)
{
  check_parameters_size(initial_parameters, lower_bounds, upper_bounds);
  size_t n_parameters = initial_parameters.size();
  auto optimal_parameters = initial_parameters;
  if (n_parameters > 1) {
    // const int number_interpolation_conditions = (n_parameters + 1) *
    //        (n_parameters + 2)/2;
    int number_interpolation_conditions = n_parameters + 3;
    std::function<double(long, const double*)> f =
      [objective, this](long n, const double* x) {
        Eigen::Map<const Eigen::VectorXd> par(x, n);
        this->objective_calls_++;
        return -objective(par);
      };
    auto result = tools_bobyqa::bobyqa(f,
                                       n_parameters,
                                       number_interpolation_conditions,
                                       initial_parameters,
                                       lower_bounds,
                                       upper_bounds,
                                       controls_.get_initial_trust_region(),
                                       controls_.get_final_trust_region(),
                                       controls_.get_maxeval());
    optimal_parameters = result.first;
    objective_max_ = -result.second;
  } else {
    double eps = 1e-6;
    std::function<double(double)> f = [objective, this](double x) {
      Eigen::Map<const Eigen::VectorXd> par(&x, 1);
      this->objective_calls_++;
      return -objective(par);
    };
    auto result = boost::math::tools::brent_find_minima(
      f,
      static_cast<double>(lower_bounds(0)) + eps,
      static_cast<double>(upper_bounds(0)) - eps,
      20);
    optimal_parameters(0) = result.first;
    objective_max_ = -result.second;
  }

  return optimal_parameters;
}

//! @brief returns how often the objective function was called.
inline size_t
Optimizer::get_objective_calls() const
{
  return objective_calls_;
}

//! @brief returns the objective value at the maximum.
inline double
Optimizer::get_objective_max() const
{
  return objective_max_;
}

//! checks whether sizes of parameters and bounds match.
inline void
Optimizer::check_parameters_size(const Eigen::VectorXd& initial_parameters,
                                 const Eigen::VectorXd& lower_bounds,
                                 const Eigen::VectorXd& upper_bounds) const
{
  if (initial_parameters.size() != upper_bounds.size()) {
    throw std::runtime_error(
      "initial parameters and and bounds must have same size.");
  }
  if (lower_bounds.size() != upper_bounds.size()) {
    throw std::runtime_error("lower and upper bounds must have same size.");
  }
  if (lower_bounds.size() < 1) {
    throw std::runtime_error("n_parameters should be larger than 0.");
  }
}

//! Create controls using the default contructor
//!
//! The defaults are
//! ```
//! initial_trust_region_ = 1e-4;
//! final_trust_region_ = 1e3;
//! maxeval_ = 1000;
//! ```
inline BobyqaControls::BobyqaControls()
{
  initial_trust_region_ = 1e-4;
  final_trust_region_ = 1e3;
  maxeval_ = 1000;
}

//! Create controls by passing the arguments
//!
//! @param initial_trust_region initial trust region.
//! @param final_trust_region final trust region.
//! @param maxeval maximal number of evaluations of the objective.
inline BobyqaControls::BobyqaControls(double initial_trust_region,
                                      double final_trust_region,
                                      int maxeval)
{
  check_parameters(initial_trust_region, final_trust_region, maxeval);
  initial_trust_region_ = initial_trust_region;
  final_trust_region_ = final_trust_region;
  maxeval_ = maxeval;
}

inline void
BobyqaControls::check_parameters(double initial_trust_region,
                                 double final_trust_region,
                                 int maxeval)
{
  if (initial_trust_region <= 0) {
    throw std::runtime_error("initial_trust_region should be larger than 0");
  }
  if (final_trust_region <= 0) {
    throw std::runtime_error("final_trust_region should be larger than 0");
  }
  if (maxeval <= 0) {
    throw std::runtime_error("maxeval should be larger than 0");
  }
}

//! @name Getters and setters
//! @{

//! @return the initial trust region.
inline double
BobyqaControls::get_initial_trust_region()
{
  return initial_trust_region_;
}

//! @return the final trust region.
inline double
BobyqaControls::get_final_trust_region()
{
  return final_trust_region_;
}

//! @return the maximal number of evaluations of the objective.
inline int
BobyqaControls::get_maxeval()
{
  return maxeval_;
}

//! @}
}
}
