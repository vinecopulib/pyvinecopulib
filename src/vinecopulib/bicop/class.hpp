// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <boost/property_tree/ptree.hpp>
#include <vinecopulib/bicop/fit_controls.hpp>

namespace vinecopulib {

// forward declaration of Abstract class
class AbstractBicop;
using BicopPtr = std::shared_ptr<AbstractBicop>;

//! @brief A class for bivariate copula models.
//!
//! The copula model is fully characterized by the family, rotation,
//! and parameters.
class Bicop
{

public:
  // Constructors
  Bicop(const BicopFamily family = BicopFamily::indep,
        const int rotation = 0,
        const Eigen::MatrixXd& parameters = Eigen::MatrixXd());

  Bicop(const Eigen::Matrix<double, Eigen::Dynamic, 2>& data,
        const FitControlsBicop& controls = FitControlsBicop());

  Bicop(const char* filename);

  Bicop(const boost::property_tree::ptree input);

  // Serialize
  boost::property_tree::ptree to_ptree() const;

  void to_json(const char* filename) const;

  // Getters and setters
  BicopFamily get_family() const;

  std::string get_family_name() const;

  int get_rotation() const;

  Eigen::MatrixXd get_parameters() const;

  double get_tau() const;

  double get_loglik() const;
  size_t get_nobs() const;
  double get_aic() const;
  double get_bic() const;
  double get_mbic(const double psi0) const;

  void set_rotation(const int rotation);

  void set_parameters(const Eigen::MatrixXd& parameters);

  // Stats methods
  Eigen::VectorXd pdf(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u) const;

  Eigen::VectorXd cdf(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u) const;

  Eigen::VectorXd hfunc1(
    const Eigen::Matrix<double, Eigen::Dynamic, 2>& u) const;

  Eigen::VectorXd hfunc2(
    const Eigen::Matrix<double, Eigen::Dynamic, 2>& u) const;

  Eigen::VectorXd hinv1(
    const Eigen::Matrix<double, Eigen::Dynamic, 2>& u) const;

  Eigen::VectorXd hinv2(
    const Eigen::Matrix<double, Eigen::Dynamic, 2>& u) const;

  Eigen::Matrix<double, Eigen::Dynamic, 2> simulate(
    const size_t& n,
    const bool qrng = false,
    const std::vector<int>& seeds = std::vector<int>()) const;

  // Methods modifying the family/rotation/parameters
  void fit(const Eigen::Matrix<double, Eigen::Dynamic, 2>& data,
           const FitControlsBicop& controls = FitControlsBicop());

  void select(const Eigen::Matrix<double, Eigen::Dynamic, 2>& data,
              FitControlsBicop controls = FitControlsBicop());

  // Fit statistics
  double loglik(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u =
                  Eigen::Matrix<double, Eigen::Dynamic, 2>()) const;

  double aic(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u =
               Eigen::Matrix<double, Eigen::Dynamic, 2>()) const;

  double bic(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u =
               Eigen::Matrix<double, Eigen::Dynamic, 2>()) const;

  double mbic(const Eigen::Matrix<double, Eigen::Dynamic, 2>& u =
                Eigen::Matrix<double, Eigen::Dynamic, 2>(),
              const double psi0 = 0.9) const;

  // Misc
  std::string str() const;

  double calculate_npars() const;

  double parameters_to_tau(const Eigen::MatrixXd& parameters) const;

  Eigen::MatrixXd tau_to_parameters(const double& tau) const;

  void flip();

  Eigen::MatrixXd get_parameters_lower_bounds() const;

  Eigen::MatrixXd get_parameters_upper_bounds() const;

private:
  Eigen::Matrix<double, Eigen::Dynamic, 2> cut_and_rotate(
    const Eigen::Matrix<double, Eigen::Dynamic, 2>& u) const;

  void check_rotation(int rotation) const;

  void check_weights_size(const Eigen::VectorXd& weights,
                          const Eigen::MatrixXd& data) const;

  void check_fitted() const;

  double compute_mbic_penalty(const size_t nobs, const double psi0) const;

  BicopPtr get_bicop() const;

  BicopPtr bicop_;
  int rotation_;
  size_t nobs_;
};
}

#include <vinecopulib/bicop/implementation/class.ipp>
