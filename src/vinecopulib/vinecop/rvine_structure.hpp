// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <Eigen/Dense>
#include <vinecopulib/misc/triangular_array.hpp>

namespace vinecopulib {

//! @brief R-vine structures
//!
//! RVineStructure objects encode the tree structure of the vine, i.e. the
//! conditioned/conditioning variables of each edge. It is represented by a
//! triangular array. An exemplary array is
//! ```
//! 4 4 4 4
//! 3 3 3
//! 2 2
//! 1
//! ```
//! which encodes the following pair-copulas:
//! ```
//! | tree | edge | pair-copulas   |
//! |------|------|----------------|
//! | 0    | 0    | `(1, 4)`       |
//! |      | 1    | `(2, 4)`       |
//! |      | 2    | `(3, 4)`       |
//! | 1    | 0    | `(1, 3; 4)`    |
//! |      | 1    | `(2, 3; 4)`    |
//! | 2    | 0    | `(1, 2; 3, 4)` |
//! ```
//! Denoting by `M[i, j]` the array entry in row `i` and column `j`,
//! the pair-copula index for edge `e` in tree `t` of a `d` dimensional vine
//! is `(M[d - 1 - t, e], M[t, e]; M[t - 1, e], ..., M[0, e])`. Less
//! formally,
//! 1. Start with the counter-diagonal element of column `e` (first conditioned
//!    variable).
//! 2. Jump up to the element in row `t` (second conditioned variable).
//! 3. Gather all entries further up in column `e` (conditioning set).
//!
//! A valid R-vine array must satisfy several conditions which are checked
//! when `RVineStructure()` is called:
//! 1. It only contains numbers between 1 and d.
//! 2. The diagonal must contain the numbers 1, ..., d.
//! 3. The diagonal entry of a column must not be contained in any
//!    column further to the right.
//! 4. The entries of a column must be contained in all columns to the left.
//! 5. The proximity condition must hold: For all t = 1, ..., d - 2 and
//!    e = 0, ..., d - t - 1 there must exist an index j > d, such that
//!    `(M[t, e], {M[0, e], ..., M[t-1, e]})` equals either
//!    `(M[d-j-1, j], {M[0, j], ..., M[t-1, j]})` or
//!    `(M[t-1, j], {M[d-j-1, j], M[0, j], ..., M[t-2, j]})`.
//!
//! An R-vine array is said to be in natural order when the anti-diagonal
//! entries are \f$ 1, \dots, d \f$ (from left to right). The exemplary arrray
//! above is in natural order. Any R-vine array can be characterized by the
//! diagonal entries (called order) and the entries below the diagonal of the
//! corresponding R-vine array in natural order. Since most algorithms work
//! with the structure in natural order, this is how RVineStructure stores the
//! structure internally.
class RVineStructure {
public:
    RVineStructure() {}

    RVineStructure(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat,
        bool check = true);
    RVineStructure(const std::vector<size_t>& order,
                   bool check = true);
    RVineStructure(const std::vector<size_t>& order,
                   const size_t& trunc_lvl,
                   bool check = true);
    RVineStructure(const std::vector<size_t>& order,
                   const TriangularArray<size_t>& struct_array,
                   bool is_natural_order = false,
                   bool check = true);

    size_t get_dim() const;
    size_t get_trunc_lvl() const;
    std::vector<size_t> get_order() const;
    TriangularArray<size_t> get_struct_array() const;
    TriangularArray<size_t> get_min_array() const;
    TriangularArray<size_t> get_needed_hfunc1() const;
    TriangularArray<size_t> get_needed_hfunc2() const;
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> get_matrix() const;

    size_t struct_array(size_t tree, size_t edge) const;
    size_t min_array(size_t tree, size_t edge) const;

    void truncate(size_t trunc_lvl);
    std::string str() const;

protected:
    size_t find_trunc_lvl(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;
    std::vector<size_t> get_order(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;
    TriangularArray<size_t> to_rvine_array(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;

    TriangularArray<size_t> to_natural_order() const;
    TriangularArray<size_t> compute_dvine_struct_array() const;
    TriangularArray<size_t> compute_min_array() const;
    TriangularArray<size_t> compute_needed_hfunc1() const;
    TriangularArray<size_t> compute_needed_hfunc2() const;

    void check_if_quadratic(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;
    void check_lower_tri(
        const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat) const;
    void check_upper_tri() const;
    void check_columns() const;
    void check_antidiagonal() const;
    void check_proximity_condition() const;

    std::vector<size_t> order_;
    size_t d_;
    size_t trunc_lvl_;
    TriangularArray<size_t> struct_array_;
    TriangularArray<size_t> min_array_;
    TriangularArray<size_t> needed_hfunc1_;
    TriangularArray<size_t> needed_hfunc2_;
};

std::ostream& operator<<(std::ostream& os, const RVineStructure& rvs);

}

#include <vinecopulib/vinecop/implementation/rvine_structure.ipp>
