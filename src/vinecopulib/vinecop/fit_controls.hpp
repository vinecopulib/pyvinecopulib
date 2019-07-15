// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <limits>
#include <vinecopulib/bicop/fit_controls.hpp>

#if defined(__GNUC__) || defined(__clang__)
#define DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED
#endif

namespace vinecopulib {
//! @brief A class for controlling fit of bivariate copula models.
//!
class FitControlsVinecop : public FitControlsBicop
{
public:
  // Constructor
  FitControlsVinecop();

  FitControlsVinecop(std::vector<BicopFamily> family_set,
                     std::string parametric_method = "mle",
                     std::string nonparametric_method = "quadratic",
                     double nonparametric_mult = 1.0,
                     size_t trunc_lvl = std::numeric_limits<size_t>::max(),
                     std::string tree_criterion = "tau",
                     double threshold = 0.0,
                     std::string selection_criterion = "bic",
                     const Eigen::VectorXd& weights = Eigen::VectorXd(),
                     double psi0 = 0.9,
                     bool preselect_families = true,
                     bool select_trunc_lvl = false,
                     bool select_threshold = false,
                     bool show_trace = false,
                     size_t num_threads = 1);

  FitControlsVinecop(const FitControlsBicop& controls,
                     size_t trunc_lvl = std::numeric_limits<size_t>::max(),
                     std::string tree_criterion = "tau",
                     double threshold = 0.0,
                     bool select_trunc_lvl = false,
                     bool select_threshold = false,
                     bool show_trace = false,
                     size_t num_threads = 1);

  // Getters
  DEPRECATED size_t get_truncation_level() const;
  size_t get_trunc_lvl() const;

  std::string get_tree_criterion() const;

  double get_threshold() const;

  bool get_show_trace() const;

  DEPRECATED bool get_select_truncation_level() const;
  bool get_select_trunc_lvl() const;

  bool get_select_threshold() const;

  bool needs_sparse_select() const;

  FitControlsBicop get_fit_controls_bicop() const;

  // Setters
  DEPRECATED void set_truncation_level(size_t trunc_lvl);
  void set_trunc_lvl(size_t trunc_lvl);

  void set_tree_criterion(std::string tree_criterion);

  void set_threshold(double threshold);

  void set_show_trace(bool show_trace);

  DEPRECATED void set_select_truncation_level(bool select_trunc_lvl);
  void set_select_trunc_lvl(bool select_trunc_lvl);

  void set_select_threshold(bool select_threshold);

  void set_fit_controls_bicop(FitControlsBicop controls);

private:
  size_t trunc_lvl_;
  std::string tree_criterion_;
  double threshold_;
  bool show_trace_;
  bool select_trunc_lvl_;
  bool select_threshold_;

  void check_tree_criterion(std::string tree_criterion);

  void check_threshold(double threshold);
};
}

#include <vinecopulib/vinecop/implementation/fit_controls.ipp>
