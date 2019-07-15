// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>

namespace vinecopulib {

//! @brief Triangular arrays.
//!
//! A triangular array behaves like a matrix with the structure
//! ```
//! x x x x x
//! x x x x
//! x x x
//! x x
//! x
//! ```
//! and all other elements omitted. This structure appears naturally in the
//! representation of a vine copula model and related algorithms. Each row
//! corresponds to one tree in the vine, starting from the top. In each row
//! (= tree), each column represents an edge in this tree.
//!
//! For truncated vine models the last few rows are omitted. For example, a
//! 3-truncated version of the above array contains the elements
//! ```
//! x x x x x
//! x x x x
//! x x x
//! ```
//! Only the elements indicated by `x`s are stored and can be accessed.
//!
//! The data structure is templated and any type or class can be used to fill
//! the entries (`x`s) of the triangular array.
template<typename T>
class TriangularArray
{
public:
  TriangularArray() = default;
  TriangularArray(size_t d);
  TriangularArray(size_t d, size_t trunc_lvl);

  T& operator()(size_t tree, size_t edge);
  T operator()(size_t tree, size_t edge) const;
  std::vector<T>& operator[](size_t column);
  std::vector<T> operator[](size_t column) const;
  bool operator==(const TriangularArray<T>& rhs) const;

  void set_column(size_t column, const std::vector<size_t>& new_column);
  void truncate(size_t trunc_lvl);

  size_t get_trunc_lvl() const;
  size_t get_dim() const;

  std::string str() const;

private:
  size_t d_;
  size_t trunc_lvl_;
  std::vector<std::vector<T>> mat_;
};

//! @brief construct a triangular array of dimension `d`.
//!
//! The array has `d-1` columns and `d-1` rows.
//! @param d the dimension of the underlying vine.
template<typename T>
TriangularArray<T>::TriangularArray(size_t d)
  : TriangularArray(d, d - 1)
{}

//! @brief construct a truncated triangular array
//!
//! The array has `d-1` columns and `min(trunv_lvl, d-1)` rows.
//! @param d the dimension of the vine.
//! @param trunc_lvl the truncation level.
template<typename T>
TriangularArray<T>::TriangularArray(size_t d, size_t trunc_lvl)
  : d_(d)
  , trunc_lvl_(std::min(d - 1, trunc_lvl))
{
  if (d < 2)
    throw std::runtime_error("d should be greater than 1");

  mat_ = std::vector<std::vector<T>>(d - 1);
  for (size_t i = 0; i < d - 1; i++)
    mat_[i] = std::vector<T>(std::min(d - i - 1, trunc_lvl));
}

//! @brief access one element of the trapezoid (writable).
//! @param tree the tree level.
//! @param edge the edge in this tree.
template<typename T>
T&
TriangularArray<T>::operator()(size_t tree, size_t edge)
{
  assert(tree < trunc_lvl_);
  assert(edge < d_ - 1 - tree);
  return mat_[edge][tree];
}

//! @brief access one element of the trapezoid (non-writable).
//! @param tree the tree level.
//! @param edge the edge in this tree.
template<typename T>
T
TriangularArray<T>::operator()(size_t tree, size_t edge) const
{
  assert(tree < trunc_lvl_);
  assert(edge < d_ - 1 - tree);
  return mat_[edge][tree];
}

//! @brief access one column of the trapezoid (writable).
//! @param column which column to extract.
template<typename T>
std::vector<T>& TriangularArray<T>::operator[](size_t column)
{
  assert(column < d_ - 1);
  return mat_[column];
}

//! @brief access one column of the trapezoid (non-writable).
//! @param column which column to extract.
template<typename T>
std::vector<T> TriangularArray<T>::operator[](size_t column) const
{
  assert(column < d_ - 1);
  return mat_[column];
}

//! @brief set one column of the trapezoid.
//! @param column which column to set.
//! @param new_column the column column to set.
template<typename T>
void
TriangularArray<T>::set_column(size_t column,
                               const std::vector<size_t>& new_column)
{
  if (column >= d_ - 1) {
    std::stringstream problem;
    problem << "column should be smaller than " << d_ - 1 << ".";
    throw std::runtime_error(problem.str());
  }
  if (new_column.size() != mat_[column].size()) {
    std::stringstream problem;
    problem << "column " << column << " should have size "
            << mat_[column].size() << ".";
    throw std::runtime_error(problem.str());
  }

  mat_[column] = new_column;
}

//! @brief truncates the trapezoid.
//! If the trapezoid is already truncated at a level
//! less than `trunc_lvl`, the function does nothing.
//! @param trunc_lvl the truncation level.
template<typename T>
void
TriangularArray<T>::truncate(size_t trunc_lvl)
{
  if (trunc_lvl < this->get_trunc_lvl()) {
    trunc_lvl_ = trunc_lvl;
    for (size_t column = 0; column < d_ - 1 - trunc_lvl; column++) {
      mat_[column].resize(trunc_lvl);
    }
  }
}

//! @brief equality operator to compare two TriangularArray objects.
//! @param rhs right-hand-side of the equality operator.
template<typename T>
bool
TriangularArray<T>::operator==(const TriangularArray<T>& rhs) const
{
  if ((d_ != rhs.get_dim()) | (trunc_lvl_ != rhs.get_trunc_lvl()))
    return false;

  for (size_t i = 0; i < d_ - 1; i++) {
    if (!((*this)[i] == rhs[i]))
      return false;
  }
  return true;
}

//! get the truncation level of the underlying vine.
template<typename T>
size_t
TriangularArray<T>::get_trunc_lvl() const
{
  return trunc_lvl_;
}

//! get the dimension of the underlying vine (the matrix has `d-1` columns and
//! `min(trunv_lvl, d-1)` rows).
template<typename T>
size_t
TriangularArray<T>::get_dim() const
{
  return d_;
}

//! represent RightTrapezoid as a string.
template<typename T>
std::string
TriangularArray<T>::str() const
{
  std::stringstream str;
  for (size_t i = 0; i < std::min(d_ - 1, trunc_lvl_); i++) {
    for (size_t j = 0; j < d_ - i - 1; j++) {
      str << (*this)(i, j) << " ";
    }
    str << std::endl;
  }
  return str.str();
}

//! @brief ostream method for RightTrapezoid, to be used with `std::cout`
//! @param os an output stream.
//! @param tri_array n triangular array.
template<typename T>
std::ostream&
operator<<(std::ostream& os, const TriangularArray<T>& tri_array)
{
  os << tri_array.str();
  return os;
}

} // end of namespace vinecopulib!
