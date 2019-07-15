// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <string>
#include <vector>

namespace vinecopulib {

//! @brief A bivariate copula family identifier.
enum class BicopFamily
{
  indep,    ///< Independence copula
  gaussian, ///< Gaussian copula
  student,  ///< Student t copula
  clayton,  ///< Clayton copula
  gumbel,   ///< Gumbel copula
  frank,    ///< Frank copula
  joe,      ///< Joe copula
  bb1,      ///< BB1 copula
  bb6,      ///< BB6 copula
  bb7,      ///< BB7 copula
  bb8,      ///< BB8 copula
  tll       ///< Transformation local likelihood kernel estimator
};

std::string
get_family_name(BicopFamily family);

BicopFamily
get_family_enum(std::string family);

//! Convenience definitions of sets of bivariate copula families
namespace bicop_families {

//! All implemented families
const std::vector<BicopFamily> all = {
  BicopFamily::indep,   BicopFamily::gaussian, BicopFamily::student,
  BicopFamily::clayton, BicopFamily::gumbel,   BicopFamily::frank,
  BicopFamily::joe,     BicopFamily::bb1,      BicopFamily::bb6,
  BicopFamily::bb7,     BicopFamily::bb8,      BicopFamily::tll
};

//! All parametric families
const std::vector<BicopFamily> parametric = {
  BicopFamily::indep,   BicopFamily::gaussian, BicopFamily::student,
  BicopFamily::clayton, BicopFamily::gumbel,   BicopFamily::frank,
  BicopFamily::joe,     BicopFamily::bb1,      BicopFamily::bb6,
  BicopFamily::bb7,     BicopFamily::bb8
};

//! All nonparametric families
const std::vector<BicopFamily> nonparametric = { BicopFamily::indep,
                                                 BicopFamily::tll };

//! All one-parameter families
const std::vector<BicopFamily> one_par = {
  BicopFamily::gaussian, BicopFamily::clayton, BicopFamily::gumbel,
  BicopFamily::frank,    BicopFamily::joe,
};

//! All two-parameter families
const std::vector<BicopFamily> two_par = { BicopFamily::student,
                                           BicopFamily::bb1,
                                           BicopFamily::bb6,
                                           BicopFamily::bb7,
                                           BicopFamily::bb8 };

//! All elliptical copulas
const std::vector<BicopFamily> elliptical = { BicopFamily::gaussian,
                                              BicopFamily::student };

//! All Archimedean copulas
const std::vector<BicopFamily> archimedean = {
  BicopFamily::clayton, BicopFamily::gumbel, BicopFamily::frank,
  BicopFamily::joe,     BicopFamily::bb1,    BicopFamily::bb6,
  BicopFamily::bb7,     BicopFamily::bb8
};

//! All BB copulas
const std::vector<BicopFamily> bb = { BicopFamily::bb1,
                                      BicopFamily::bb6,
                                      BicopFamily::bb7,
                                      BicopFamily::bb8 };

//! @brief All copulas that don't have a rotation
//!
//! (because they already cover positive and negative dependence)
const std::vector<BicopFamily> rotationless = { BicopFamily::indep,
                                                BicopFamily::gaussian,
                                                BicopFamily::student,
                                                BicopFamily::frank,
                                                BicopFamily::tll };

//! Families with stronger dependence in the lower tail
const std::vector<BicopFamily> lt = { BicopFamily::clayton,
                                      BicopFamily::bb1,
                                      BicopFamily::bb7 };

//! Families with stronger dependence in the upper tail
const std::vector<BicopFamily> ut = { BicopFamily::gumbel, BicopFamily::joe,
                                      BicopFamily::bb1,    BicopFamily::bb6,
                                      BicopFamily::bb7,    BicopFamily::bb8 };

//! Families for which `method = "itau"` is available in Bicop::fit()
const std::vector<BicopFamily> itau = {
  BicopFamily::indep,   BicopFamily::gaussian, BicopFamily::student,
  BicopFamily::clayton, BicopFamily::gumbel,   BicopFamily::frank,
  BicopFamily::joe
};

//! Families that can be flipped by adjusting the rotation.
const std::vector<BicopFamily> flip_by_rotation = {
  BicopFamily::clayton, BicopFamily::gumbel, BicopFamily::frank,
  BicopFamily::joe,     BicopFamily::bb1,    BicopFamily::bb6,
  BicopFamily::bb7,     BicopFamily::bb8
};

} // end of namespace BicopFamilies
} // end of namespace vinecopulib

#include <vinecopulib/bicop/implementation/family.ipp>
