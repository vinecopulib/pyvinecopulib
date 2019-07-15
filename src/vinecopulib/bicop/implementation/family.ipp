// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <boost/bimap.hpp>
#include <boost/assign.hpp>

namespace vinecopulib {

typedef boost::bimap <BicopFamily, std::string> family_bimap;
const family_bimap family_names =
    boost::assign::list_of<family_bimap::relation>
        (BicopFamily::indep, "Independence")
        (BicopFamily::gaussian, "Gaussian")
        (BicopFamily::student, "Student")
        (BicopFamily::clayton, "Clayton")
        (BicopFamily::gumbel, "Gumbel")
        (BicopFamily::frank, "Frank")
        (BicopFamily::joe, "Joe")
        (BicopFamily::bb1, "BB1")
        (BicopFamily::bb6, "BB6")
        (BicopFamily::bb7, "BB7")
        (BicopFamily::bb8, "BB8")
        (BicopFamily::tll, "TLL");

//! converts a BicopFamily into a string with its name.
//! @param family the family.
inline std::string get_family_name(BicopFamily family)
{
    return family_names.left.at(family);
}

//! converts a string name into a BicopFamily.
//! @param family the family name.
inline BicopFamily get_family_enum(std::string family)
{
    return family_names.right.at(family);
}
}
