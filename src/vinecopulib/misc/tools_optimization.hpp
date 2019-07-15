// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>

namespace vinecopulib {

// forward declaration of parametric subclass
class ParBicop;

namespace tools_optimization {

//! @brief A class for the controls to Bobyqa
class BobyqaControls
{
public:

    BobyqaControls();

    BobyqaControls(double initial_trust_region,
                   double final_trust_region,
                   int maxeval);

    double get_initial_trust_region();

    double get_final_trust_region();

    int get_maxeval();

private:
    double initial_trust_region_; //! Initial trust region
    double final_trust_region_; //! Final trust region
    int maxeval_; //! Maximal number of evaluations of the objective

    //! Sanity checks
    //! @{
    void check_parameters(double initial_trust_region,
                          double final_trust_region,
                          int maxeval);
    //! @}
};

//! @brief A class for optimization (wrapping Bobyqa).
class Optimizer
{
public:
    Optimizer();

    void set_controls(double initial_trust_region,
                      double final_trust_region,
                      int maxeval);

    Eigen::VectorXd optimize(
        const Eigen::VectorXd &initial_parameters,
        const Eigen::VectorXd &lower_bounds,
        const Eigen::VectorXd &upper_bounds,
        std::function<double(const Eigen::VectorXd&)> objective);

    size_t get_objective_calls() const;
    double get_objective_max() const;

private:
    void check_parameters_size(const Eigen::VectorXd &initial_parameters,
                           const Eigen::VectorXd &lower_bounds,
                           const Eigen::VectorXd &upper_bounds) const;

    BobyqaControls controls_;
    size_t objective_calls_{0};
    double objective_max_{0};
};
}

}

#include <vinecopulib/misc/implementation/tools_optimization.ipp>
