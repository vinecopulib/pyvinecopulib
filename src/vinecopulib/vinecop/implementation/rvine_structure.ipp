// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>

namespace vinecopulib {

//! @brief instantiates an RVineStructure object from a matrix representing an
//! R-vine array.
//!
//! The matrix must contain zeros in the lower right triangle and
//! the upper left triangle must be a valid R-vine array. Truncated vines can
//! be encoded by putting zeros above the digonal in all rows below the
//! truncation level. Example of a 1-truncated matrix:
//! ```
//! 4 4 4 4
//! 0 0 3 0
//! 0 2 0 0
//! 1 0 0 0
//! ```
//! @param mat a matrix representing a valid R-vine array.
//! @param check whether `mat` shall be checked for validity.
inline RVineStructure::RVineStructure(
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat,
  bool check)
{
  d_ = mat.cols();
  if (check) {
    check_if_quadratic(mat);
    check_lower_tri(mat);
  }

  order_ = get_order(mat);
  if (check)
    check_antidiagonal();

  trunc_lvl_ = find_trunc_lvl(mat);
  struct_array_ = to_rvine_array(mat);

  if (check)
    check_upper_tri();

  struct_array_ = to_natural_order();
  if (check)
    check_columns();

  min_array_ = compute_min_array();
  if (check)
    check_proximity_condition();

  needed_hfunc1_ = compute_needed_hfunc1();
  needed_hfunc2_ = compute_needed_hfunc2();
}

//! @brief instantiates an RVineStructure object to a D-vine with given ordering
//! of variables.
//! @param order the order of variables in the D-vine (diagonal entries in the
//!    R-vine array); must be a permutation of 1, ..., d.
//! @param check whether `order shall be checked for validity.
inline RVineStructure::RVineStructure(const std::vector<size_t>& order,
                                      bool check)
  : RVineStructure(order, order.size() - 1, check)
{}

inline RVineStructure::RVineStructure(const std::vector<size_t>& order,
                                      const size_t& trunc_lvl,
                                      bool check)
  : order_(order)
  , d_(order.size())
  , trunc_lvl_(std::min(trunc_lvl, d_ - 1))
{
  if (check)
    check_antidiagonal();

  if (trunc_lvl > 0) {
    struct_array_ = compute_dvine_struct_array();
    min_array_ = compute_min_array();
    needed_hfunc1_ = compute_needed_hfunc1();
    needed_hfunc2_ = compute_needed_hfunc2();
  } else {
    struct_array_ = TriangularArray<size_t>(d_, trunc_lvl);
    min_array_ = TriangularArray<size_t>(d_, trunc_lvl);
    needed_hfunc1_ = TriangularArray<size_t>(d_, trunc_lvl);
    needed_hfunc2_ = TriangularArray<size_t>(d_, trunc_lvl);
  }
}

//! @brief instantiates an RVineStructure object from the variable order
//! (diagonal elements of the R-vine array) and a triangular structure array
//! (all elements above the diagonal).
//!
//! @param order the order of variables (diagonal entries in the
//!    R-vine array); must be a permutation of 1, ..., d.
//! @param struct_array the structure array  (all elements
//!    above the diagonal in the R-vine array). For truncated vines, all rows
//!    below the truncation level are omitted.
//! @param is_natural_order whether `struct_array` is already in natural order.
//! @param check whether `order` and `struct_array` shall be checked for
//! validity.
inline RVineStructure::RVineStructure(
  const std::vector<size_t>& order,
  const TriangularArray<size_t>& struct_array,
  bool is_natural_order,
  bool check)
{
  d_ = order.size();
  if (check & (struct_array.get_dim() != d_)) {
    throw std::runtime_error("order and struct_array have "
                             "incompatible dimensions");
  }

  order_ = order;

  if (check)
    check_antidiagonal();

  trunc_lvl_ = struct_array.get_trunc_lvl();
  if (trunc_lvl_ > 0) {
    struct_array_ = struct_array;
    if (check)
      check_upper_tri();

    if (!is_natural_order)
      struct_array_ = to_natural_order();
    if (check)
      check_columns();

    min_array_ = compute_min_array();
    if (check)
      check_proximity_condition();

    needed_hfunc1_ = compute_needed_hfunc1();
    needed_hfunc2_ = compute_needed_hfunc2();
  } else {
    struct_array_ = TriangularArray<size_t>(d_, trunc_lvl_);
    min_array_ = TriangularArray<size_t>(d_, trunc_lvl_);
    needed_hfunc1_ = TriangularArray<size_t>(d_, trunc_lvl_);
    needed_hfunc2_ = TriangularArray<size_t>(d_, trunc_lvl_);
  }
}

//! extract the dimension of the vine.
inline size_t
RVineStructure::get_dim() const
{
  return d_;
}

//! extract the truncation level of the vine.
inline size_t
RVineStructure::get_trunc_lvl() const
{
  return trunc_lvl_;
}

//! @brief extract the order of variables in the vine (diagonal entries in the
//! R-vine array).
inline std::vector<size_t>
RVineStructure::get_order() const
{
  return order_;
}

//! @brief extract structure array (all elements above the diagonal in the
//! R-vine array).
inline TriangularArray<size_t>
RVineStructure::get_struct_array() const
{
  return struct_array_;
}

//! @brief extracts the minimum array.
//!
//! The minimum array is derived from an R-vine array by
//! iteratively computing the (elementwise) minimum of two subsequent rows
//! (starting from the top). It is used in estimation and evaluation algorithms
//! to find the two edges in the previous tree that are joined by the current
//! edge.
inline TriangularArray<size_t>
RVineStructure::get_min_array() const
{
  return min_array_;
}

//! @brief extracts an array indicating which of the first h-functions are
//! needed.
//!
//! (it is usually not necessary to compute both h-functions for each
//! pair-copula).
inline TriangularArray<size_t>
RVineStructure::get_needed_hfunc1() const
{
  return needed_hfunc1_;
}

//! @brief extracts an array indicating which of the second h-functions are
//! needed.
//!
//! (it is usually not necessary to compute both h-functions for each
//! pair-copula).
inline TriangularArray<size_t>
RVineStructure::get_needed_hfunc2() const
{
  return needed_hfunc2_;
}

//! @brief access elements of the structure array.
//! @param tree tree index.
//! @param edge edge index.
inline size_t
RVineStructure::struct_array(size_t tree, size_t edge) const
{
  return struct_array_(tree, edge);
}

//! @brief access elements of the minimum array.
//! @param tree tree index.
//! @param edge edge index.
inline size_t
RVineStructure::min_array(size_t tree, size_t edge) const
{
  return min_array_(tree, edge);
}

//! @brief truncates the R-vine structure.
//! @param trunc_lvl the truncation level.
//!
//! If the structure is already truncated at a level
//! less than `trunc_lvl`, the function does nothing.
inline void
RVineStructure::truncate(size_t trunc_lvl)
{
  if (trunc_lvl < trunc_lvl_) {
    struct_array_.truncate(trunc_lvl);
    min_array_.truncate(trunc_lvl);
    needed_hfunc1_.truncate(trunc_lvl);
    needed_hfunc2_.truncate(trunc_lvl);
    trunc_lvl_ = struct_array_.get_trunc_lvl();
  }
}

//! converts the structure to a string representation (most useful for
//! printing).
inline std::string
RVineStructure::str() const
{
  std::stringstream str;
  for (size_t i = 0; i < d_ - 1; i++) {
    for (size_t j = 0; j < d_ - i - 1; j++) {
      if (i < trunc_lvl_) {
        str << order_[struct_array_(i, j) - 1] << " ";
      } else {
        str << "  ";
      }
    }
    str << order_[d_ - 1 - i] << " " << std::endl;
  }
  str << order_[0] << " " << std::endl;

  return str.str();
}

//! @brief randomly sample a regular vine structure.
//! @param d the dimension.
//! @param natural_order should the sampled structure be in natural order?
//! @param seeds seeds of the random number generator; if empty (default),
//!   the random number generator is seeded randomly.
//! @note Implementation of Algorithm 13 in Harry Joe's 2014 book (p. 288),
//! but there's a typo: the end of line 6 in the book should be
//! 'column j' instead of 'column k'.
inline RVineStructure
RVineStructure::simulate(size_t d, bool natural_order, std::vector<int> seeds)
{
  auto U = tools_stats::simulate_uniform(d, d, false, seeds);

  // A is the R-vine matrix we want to create (upper right-triag format).
  // B is a random binary representation that we need to convert.
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> A(d, d), B(d, d);
  A.setZero();
  B = (U.leftCols(d).array() > 0.5).cast<size_t>();

  for (size_t i = 0; i < d; i++) {
    A(i, i) = i + 1;
    B(i, i) = 1;
    if (i > 0) {
      A(i - 1, i) = i;
      B(0, i) = 1;
      B(i - 1, i) = 1;
    }
  }
  if (d > 2) {
    A(0, 2) = 1;
  }

  for (size_t j = 3; j < d; j++) {
    int ac = j - 2;
    auto to_assign = tools_stl::seq_int(1, j - 1);
    for (ptrdiff_t k = j - 2; k >= 0; k--) {
      if (B(k, j) == 1) {
        A(k, j) = ac + 1;
        to_assign = tools_stl::set_diff(to_assign, { A(k, j) });
        if (k > 0) {
          // to_assign is always ordered ascendingly -> we pick largest
          ac = to_assign[to_assign.size() - 1] - 1;
        }
      } else {
        A(k, j) = A(k - 1, ac);
        to_assign = tools_stl::set_diff(to_assign, { A(k, j) });
      }
    }
  }

  // need to convert to upper left triangular form (our notation)
  auto rvm = RVineStructure(A.rowwise().reverse());

  // sampling the variable order randomly
  // the first column of U has not been used to construct B,
  // hence it is stochastically independent of B. Calling
  // pseudo_obs and rescaling gives us a permutation of (1, ..., d)
  // that is independent of B.
  if (!natural_order) {
    std::vector<size_t> order(d);
    U.col(0) = tools_stats::to_pseudo_obs_1d(U.col(0)) * (d + 1);
    for (size_t k = 0; k < d; k++) {
      order[k] = static_cast<size_t>(U(k, 0));
    }
    rvm = RVineStructure(order, rvm.get_struct_array(), true, false);
  }

  return rvm;
}

//! extract the R-vine matrix representation.
inline Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>
RVineStructure::get_matrix() const
{
  Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> array(d_, d_);
  array.fill(0);
  for (size_t i = 0; i < trunc_lvl_; ++i) {
    for (size_t j = 0; j < d_ - i - 1; ++j) {
      array(i, j) = order_[struct_array_(i, j) - 1];
    }
  }
  for (size_t i = 0; i < d_; ++i) {
    array(d_ - i - 1, i) = order_[i];
  }
  return array;
}

//! @brief find the truncation level in an R-vine array.
//!
//! The truncation level is
//! determined by the first row (starting from the bottom) that contains only
//! zeros above the diagonal.
//!
//! @param mat an array representing the R-vine array.
inline size_t
RVineStructure::find_trunc_lvl(
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
  size_t trunc_lvl;
  size_t d = mat.cols();

  std::stringstream problem;
  problem << "not a valid R-vine array: "
          << "a row with a 0 above the diagonal contains one or more "
          << "non-zero values.";

  for (trunc_lvl = d - 1; trunc_lvl > 0; trunc_lvl--) {
    std::vector<size_t> row_vec(d - trunc_lvl);
    Eigen::Matrix<size_t, Eigen::Dynamic, 1>::Map(&row_vec[0], d - trunc_lvl) =
      mat.row(trunc_lvl - 1).head(d - trunc_lvl);

    if (*(std::min_element(row_vec.begin(), row_vec.end())) != 0)
      break;
  }

  return trunc_lvl;
}

//! @brief find the order of an R-vine array.
//!
//! The order is contained in the counter-diagonal of the R-vine array.
//! @param mat a matrix representing the R-vine array.
inline std::vector<size_t>
RVineStructure::get_order(
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
  std::vector<size_t> order(d_);
  for (size_t i = 0; i < d_; i++)
    order[i] = mat(d_ - i - 1, i);

  return order;
}

//! @brief extracts the structure array (entries above the diagonal in R-vine
//! array).
//! @param mat a array representing the R-vine array.
inline TriangularArray<size_t>
RVineStructure::to_rvine_array(
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
  // copy upper triangle
  TriangularArray<size_t> struct_array(d_, trunc_lvl_);
  for (size_t j = 0; j < d_ - 1; j++) {
    for (size_t i = 0; i < std::min(d_ - 1 - j, trunc_lvl_); i++) {
      struct_array(i, j) = mat(i, j);
    }
  }

  return struct_array;
}

//! converts `struct_array_` to natural order.
inline TriangularArray<size_t>
RVineStructure::to_natural_order() const
{
  // create vector of new variable labels
  auto order = tools_stl::get_order(get_order());

  // relabel to natural order
  TriangularArray<size_t> struct_array(d_, trunc_lvl_);
  for (size_t j = 0; j < d_ - 1; j++) {
    for (size_t i = 0; i < std::min(d_ - 1 - j, trunc_lvl_); i++) {
      struct_array(i, j) = order[struct_array_(i, j) - 1] + 1;
    }
  }

  return struct_array;
}

//! creates a structure array corresponding to a D-vine (in natural order).
inline TriangularArray<size_t>
RVineStructure::compute_dvine_struct_array() const
{
  TriangularArray<size_t> struct_array(d_, trunc_lvl_);
  for (size_t j = 0; j < d_ - 1; j++) {
    for (size_t i = 0; i < std::min(d_ - 1 - j, trunc_lvl_); i++) {
      struct_array(i, j) = i + j + 2;
    }
  }

  return struct_array;
}

inline TriangularArray<size_t>
RVineStructure::compute_min_array() const
{
  TriangularArray<size_t> min_array = struct_array_;
  for (size_t j = 0; j < d_ - 1; j++) {
    for (size_t i = 1; i < std::min(d_ - 1 - j, trunc_lvl_); i++) {
      min_array(i, j) = std::min(struct_array_(i, j), min_array(i - 1, j));
    }
  }

  return min_array;
}

inline TriangularArray<size_t>
RVineStructure::compute_needed_hfunc1() const
{
  TriangularArray<size_t> needed_hfunc1(d_, trunc_lvl_);

  for (size_t i = 0; i < std::min(d_ - 2, trunc_lvl_ - 1); i++) {
    for (size_t j = 0; j < d_ - 2 - i; j++) {
      if (struct_array_(i + 1, j) != min_array_(i + 1, j))
        needed_hfunc1(i, min_array_(i + 1, j) - 1) = 1;
    }
  }

  return needed_hfunc1;
}

inline TriangularArray<size_t>
RVineStructure::compute_needed_hfunc2() const
{
  TriangularArray<size_t> needed_hfunc2(d_, trunc_lvl_);

  for (size_t i = 0; i < std::min(d_ - 2, trunc_lvl_ - 1); i++) {
    for (size_t j = 0; j < d_ - 2 - i; j++) {
      needed_hfunc2(i, j) = 1;
      if (struct_array_(i + 1, j) == min_array_(i + 1, j))
        needed_hfunc2(i, min_array_(i + 1, j) - 1) = 1;
    }
  }

  return needed_hfunc2;
}

inline void
RVineStructure::check_if_quadratic(
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
  std::string problem = "must be quadratic.";
  if (mat.rows() != mat.cols()) {
    throw std::runtime_error("not a valid R-vine array: " + problem);
  }
}

inline void
RVineStructure::check_lower_tri(
  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const
{
  std::string problem = "the lower right triangle must only contain zeros";
  size_t sum_lwr = 0;
  for (size_t j = 1; j < d_; ++j) {
    sum_lwr += mat.block(d_ - j, j, j, 1).array().sum();
    if (sum_lwr != 0) {
      throw std::runtime_error("not a valid R-vine array: " + problem);
    }
  }
}

inline void
RVineStructure::check_upper_tri() const
{
  std::string problem;
  problem += "the upper left triangle can only contain numbers ";
  problem += "between 1 and d (number of variables).";

  for (size_t j = 0; j < d_ - 1; ++j) {
    auto col_vec = struct_array_[j];
    auto minmax_in_col = std::minmax_element(col_vec.begin(), col_vec.end());
    if ((*(minmax_in_col.first) < 1) | (*(minmax_in_col.second) > d_)) {
      throw std::runtime_error("not a valid R-vine array: " + problem);
    }
  }
}

inline void
RVineStructure::check_columns() const
{
  std::string problem;
  problem += "the antidiagonal entry of a column must not be ";
  problem += "contained in any column further to the right; ";
  problem += "the entries of a column must be contained ";
  problem += "in all columns to the left.";

  // check that column j only contains unique indices in 1:(d - j).
  for (size_t j = 0; j < d_ - 1; ++j) {
    auto col_vec = struct_array_[j];
    std::sort(col_vec.begin(), col_vec.end());
    size_t unique_in_col =
      std::unique(col_vec.begin(), col_vec.end()) - col_vec.begin();
    if ((!tools_stl::is_member(col_vec, tools_stl::seq_int(1 + j, d_))) |
        (unique_in_col != col_vec.size())) {
      throw std::runtime_error("not a valid R-vine array: " + problem);
    }
  }
}

inline void
RVineStructure::check_antidiagonal() const
{
  std::string problem;
  problem += "the order/antidiagonal must contain the numbers ";
  problem += "1, ..., d (the number of variables)";
  if (!tools_stl::is_same_set(order_, tools_stl::seq_int(1, d_))) {
    throw std::runtime_error("not a valid R-vine array: " + problem);
  }
}

inline void
RVineStructure::check_proximity_condition() const
{
  for (size_t t = 1; t < trunc_lvl_; ++t) {
    for (size_t e = 0; e < d_ - t - 1; ++e) {
      std::vector<size_t> target_set(t + 1), test_set(t + 1);
      // conditioning set
      for (size_t i = 0; i < t; i++) {
        target_set[i] = struct_array_(i, e);
        test_set[i] = struct_array_(i, min_array_(t, e) - 1);
      }

      // non-diagonal conditioned variable
      target_set[t] = struct_array_(t, e);
      // diagonal conditioned variable in other column
      test_set[t] = min_array_(t, e);

      if (!tools_stl::is_same_set(target_set, test_set)) {
        std::stringstream problem;
        problem << "not a valid R-vine array: "
                << "proximity condition violated; "
                << "cannot extract conditional distribution (" << target_set[t]
                << " | ";
        for (size_t i = 0; i < t - 1; ++i) {
          problem << order_[target_set[i] - 1] << ", ";
        }
        problem << order_[target_set[t - 1] - 1] << ") from pair-copulas.";
        throw std::runtime_error(problem.str().c_str());
      }
    }
  }
}

//! @brief ostream method for RVineStructure, to be used with `std::cout`.
//! @param os output stream.
//! @param rvs r-vine structure array.
inline std::ostream&
operator<<(std::ostream& os, const RVineStructure& rvs)
{
  os << rvs.str();
  return os;
}
}
