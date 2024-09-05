#pragma once

#include "docstr.hpp"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vinecopulib.hpp>

namespace py = pybind11;
using namespace vinecopulib;

inline void
init_bicop_family(py::module_& module)
{

  constexpr auto& doc = pyvinecopulib_doc;
  constexpr auto& bicopfamily_doc = doc.vinecopulib.BicopFamily;

  py::enum_<BicopFamily>(module, "BicopFamily", py::arithmetic(), R"pbdoc(
   A bivariate copula family identifier.

   The following convenient sets of families are also provided:

   - ``all`` contains all the families,
   - ``parametric`` contains the parametric families (all except ``tll``),
   - ``nonparametric`` contains the nonparametric families
     (``indep`` and ``tll``)
   - ``one_par`` contains the parametric families with a single parameter,
     (``gaussian``, ``clayton``, ``gumbel``, ``frank``, and ``joe``),
   - ``two_par`` contains the parametric families with two parameters
     (``student``, ``bb1``, ``bb6``, ``bb7``, and ``bb8``),
   - ``elliptical`` contains the elliptical families,
   - ``archimedean`` contains the archimedean families,
   - ``bb`` contains the BB families,
   - ``itau`` families for which estimation by Kendall's tau inversion is
     available (``indep``, ``gaussian``, ``student``, ``clayton``,
     ``gumbel``, ``frank``, ``joe``),
   - ``lt`` contains the families that are lower-tail dependent,
   - ``ut`` contains the families that are upper-tail dependent.
   )pbdoc")
    .value("indep", BicopFamily::indep, bicopfamily_doc.indep.doc)
    .value("gaussian", BicopFamily::gaussian, bicopfamily_doc.gaussian.doc)
    .value("student", BicopFamily::student, bicopfamily_doc.student.doc)
    .value("clayton", BicopFamily::clayton, bicopfamily_doc.clayton.doc)
    .value("gumbel", BicopFamily::gumbel, bicopfamily_doc.gumbel.doc)
    .value("frank", BicopFamily::frank, bicopfamily_doc.frank.doc)
    .value("joe", BicopFamily::joe, bicopfamily_doc.joe.doc)
    .value("bb1", BicopFamily::bb1, bicopfamily_doc.bb1.doc)
    .value("bb6", BicopFamily::bb6, bicopfamily_doc.bb6.doc)
    .value("bb7", BicopFamily::bb7, bicopfamily_doc.bb7.doc)
    .value("bb8", BicopFamily::bb8, bicopfamily_doc.bb8.doc)
    .value("tll", BicopFamily::tll, bicopfamily_doc.tll.doc)
    .export_values();

  module.attr("all") = bicop_families::all;
  module.attr("parametric") = bicop_families::parametric;
  module.attr("nonparametric") = bicop_families::nonparametric;
  module.attr("one_par") = bicop_families::one_par;
  module.attr("two_par") = bicop_families::two_par;
  module.attr("elliptical") = bicop_families::elliptical;
  module.attr("archimedean") = bicop_families::archimedean;
  module.attr("bb") = bicop_families::bb;
  module.attr("lt") = bicop_families::lt;
  module.attr("ut") = bicop_families::ut;
  module.attr("itau") = bicop_families::itau;
}