// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <vector>

namespace vinecopulib {

namespace tools_stl {

template<typename T>
std::vector<size_t>
get_order(const std::vector<T>& x)
{
  std::vector<size_t> order(x.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](size_t i, size_t j) -> bool {
    return (x[i] < x[j]);
  });
  return order;
}

template<typename T>
bool
is_member(T element, std::vector<T> set)
{
  return std::find(set.begin(), set.end(), element) != set.end();
}

template<class T>
std::vector<T>
intersect(std::vector<T> x, std::vector<T> y)
{
  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end());
  std::vector<T> common;
  std::set_intersection(
    x.begin(), x.end(), y.begin(), y.end(), std::back_inserter(common));

  return common;
}

template<class T>
size_t
find_position(T x, std::vector<T> vec)
{
  return std::distance(vec.begin(), find(vec.begin(), vec.end(), x));
}

template<class T>
std::vector<T>
set_diff(std::vector<T> x, std::vector<T> y)
{
  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end());
  std::vector<T> different;
  set_difference(
    x.begin(), x.end(), y.begin(), y.end(), std::back_inserter(different));

  return different;
}

template<class T>
std::vector<T>
cat(std::vector<T> x, const std::vector<T>& y)
{
  x.reserve(x.size() + y.size());
  x.insert(x.end(), y.begin(), y.end());
  return x;
}

template<class T>
std::vector<T>
cat(T x, const std::vector<T>& y)
{
  std::vector<T> out(1);
  out[0] = x;
  out.reserve(1 + y.size());
  out.insert(out.end(), y.begin(), y.end());
  return out;
}

template<class T>
std::vector<T>
set_sym_diff(std::vector<T> x, std::vector<T> y)
{
  std::vector<T> dxy = set_diff(x, y);
  auto dyx = set_diff(y, x);
  return cat(dxy, dyx);
}

template<class T>
void
reverse(std::vector<T>& x)
{
  std::reverse(x.begin(), x.end());
}

template<class T>
bool
is_same_set(std::vector<T> x, std::vector<T> y)
{
  auto z = intersect(x, y);
  return ((z.size() == x.size()) & (z.size() == y.size()));
}

template<class T>
bool
is_member(std::vector<T> x, std::vector<T> y)
{
  auto z = intersect(x, y);
  return z.size() == x.size();
}

//! Integer sequence starting at 1
inline std::vector<size_t>
seq_int(size_t from, size_t length)
{
  std::vector<size_t> seq(length);
  std::iota(seq.begin(), seq.end(), from);
  return seq;
}

inline std::vector<size_t>
invert_permutation(const std::vector<size_t>& perm)
{
  auto inv_perm = seq_int(0, perm.size());
  std::sort(inv_perm.begin(), inv_perm.end(), [&](size_t i, size_t j) {
    return perm[i] < perm[j];
  });
  return inv_perm;
}

// safer version of log1p that avoids rounding errors if x ~ -1
inline double
log1p(const double& x)
{
  if (x < -0.5) {
    return std::log(1.0 + x);
  } else {
    return std::log1p(x);
  }
}
}
}
