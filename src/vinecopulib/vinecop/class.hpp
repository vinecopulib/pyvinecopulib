// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <boost/property_tree/ptree.hpp>
#include <vinecopulib/vinecop/fit_controls.hpp>
#include <vinecopulib/vinecop/rvine_structure.hpp>

namespace vinecopulib {

// forward declarations
class Bicop;
namespace tools_select {
class VinecopSelector;
}

//! @brief A class for vine copula models
//!
//! A vine copula model is characterized by the structure matrix (see
//! TriangularArray) and the pair-copulas.
class Vinecop
{
public:
  // default constructors
  Vinecop() {}

  Vinecop(size_t d);

  // Constructors with structure only
  Vinecop(const RVineStructure& vine_struct);

  Vinecop(const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix,
          const bool check_matrix = true);

  Vinecop(const std::vector<size_t>& order,
          const TriangularArray<size_t>& struct_array,
          const bool check_array = true);

  // Constructors with pair_copulas + structure
  Vinecop(const std::vector<std::vector<Bicop>>& pair_copulas,
          const RVineStructure& vine_struct);

  Vinecop(const std::vector<std::vector<Bicop>>& pair_copulas,
          const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix,
          const bool check_matrix = true);

  Vinecop(const std::vector<std::vector<Bicop>>& pair_copulas,
          const std::vector<size_t>& order,
          const TriangularArray<size_t>& struct_array,
          const bool check_array = true);

  // Constructors from data
  Vinecop(const Eigen::MatrixXd& data,
          const FitControlsVinecop& controls = FitControlsVinecop());

  Vinecop(const Eigen::MatrixXd& data,
          const RVineStructure& vine_struct,
          FitControlsVinecop controls = FitControlsVinecop());

  Vinecop(const Eigen::MatrixXd& data,
          const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& matrix,
          FitControlsVinecop controls = FitControlsVinecop(),
          const bool check_matrix = true);

  Vinecop(const Eigen::MatrixXd& data,
          const std::vector<size_t>& order,
          const TriangularArray<size_t>& struct_array,
          FitControlsVinecop controls = FitControlsVinecop(),
          const bool check_array = true);

  // Constructors from files/serialized objects
  Vinecop(const char* filename, const bool check_matrix = true);

  Vinecop(const boost::property_tree::ptree input,
          const bool check_matrix = true);

  // Serialize
  boost::property_tree::ptree to_ptree() const;

  void to_json(const char* filename) const;

  // Methods modifying structure and/or families and parameters
  void select_all(const Eigen::MatrixXd& data,
                  const FitControlsVinecop& controls = FitControlsVinecop());

  void select_families(
    const Eigen::MatrixXd& data,
    const FitControlsVinecop& controls = FitControlsVinecop());

  // Getters for a single pair copula
  Bicop get_pair_copula(size_t tree, size_t edge) const;

  BicopFamily get_family(size_t tree, size_t edge) const;

  int get_rotation(size_t tree, size_t edge) const;

  Eigen::MatrixXd get_parameters(size_t tree, size_t edge) const;

  double get_tau(size_t tree, size_t edge) const;

  size_t get_trunc_lvl() const;

  // Getters for all pair copulas
  std::vector<std::vector<Bicop>> get_all_pair_copulas() const;

  std::vector<std::vector<BicopFamily>> get_all_families() const;

  std::vector<std::vector<int>> get_all_rotations() const;

  std::vector<std::vector<Eigen::MatrixXd>> get_all_parameters() const;

  std::vector<std::vector<double>> get_all_taus() const;

  // Getters for the structure
  size_t get_dim() const;

  std::vector<size_t> get_order() const;

  RVineStructure get_rvine_structure() const;

  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> get_matrix() const;

  TriangularArray<size_t> get_struct_array() const;

  // getters for fit statistics
  double get_threshold() const;
  double get_loglik() const;
  size_t get_nobs() const;
  double get_aic() const;
  double get_bic() const;
  double get_mbicv(const double psi0) const;

  // Stats methods
  Eigen::VectorXd pdf(const Eigen::MatrixXd& u,
                      const size_t num_threads = 1) const;

  Eigen::VectorXd cdf(const Eigen::MatrixXd& u,
                      const size_t N = 1e4,
                      const size_t num_threads = 1,
                      std::vector<int> seeds = std::vector<int>()) const;

  Eigen::MatrixXd simulate(
    const size_t n,
    const bool qrng = false,
    const size_t num_threads = 1,
    const std::vector<int>& seeds = std::vector<int>()) const;

  Eigen::MatrixXd rosenblatt(const Eigen::MatrixXd& u,
                             const size_t num_threads = 1) const;
  Eigen::MatrixXd inverse_rosenblatt(const Eigen::MatrixXd& u,
                                     const size_t num_threads = 1) const;

  // Fit statistics
  double calculate_npars() const;

  double loglik(const Eigen::MatrixXd& u, const size_t num_threads = 1) const;

  double aic(const Eigen::MatrixXd& u, const size_t num_threads = 1) const;

  double bic(const Eigen::MatrixXd& u, const size_t num_threads = 1) const;

  double mbicv(const Eigen::MatrixXd& u,
               const double psi0,
               const size_t num_threads = 1) const;

  // Misc methods
  static std::vector<std::vector<Bicop>> make_pair_copula_store(
    const size_t d,
    const size_t trunc_lvl = std::numeric_limits<size_t>::max());
  void truncate(size_t trunc_lvl);

  std::string str() const;

protected:
  size_t d_;
  RVineStructure vine_struct_;
  std::vector<std::vector<Bicop>> pair_copulas_;
  double threshold_;
  double loglik_;
  size_t nobs_;

  void check_data_dim(const Eigen::MatrixXd& data) const;
  void check_pair_copulas_rvine_structure(
    const std::vector<std::vector<Bicop>>& pair_copulas) const;
  double calculate_mbicv_penalty(const size_t nobs, const double psi0) const;
  void finalize_fit(const tools_select::VinecopSelector& selector);
  void check_weights_size(const Eigen::VectorXd& weights,
                          const Eigen::MatrixXd& data) const;
  void check_enough_data(const Eigen::MatrixXd& data) const;
  void check_fitted() const;
};
}

#include <vinecopulib/vinecop/implementation/class.ipp>
