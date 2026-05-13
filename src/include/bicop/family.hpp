#pragma once

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <vinecopulib.hpp>

#include "docstr.hpp"

namespace nb = nanobind;
using namespace vinecopulib;

inline void init_bicop_family(nb::module_& module) {
  constexpr auto& doc = pyvinecopulib_doc;
  constexpr auto& bicopfamily_doc = doc.vinecopulib.BicopFamily;

  nb::enum_<BicopFamily>(module, "BicopFamily", R"pbdoc(
    A bivariate copula family identifier.

    Members of this enum are the family tags used throughout the API
    (e.g. as values inside ``FitControlsBicop.family_set`` or returned
    by :meth:`pyvinecopulib.core.Bicop.family`):

      - ``indep``: Independent copula,
      - ``gaussian``: Gaussian copula,
      - ``student``: Student t copula,
      - ``clayton``: Clayton copula,
      - ``gumbel``: Gumbel copula,
      - ``frank``: Frank copula,
      - ``joe``: Joe copula,
      - ``bb1``: BB1 copula,
      - ``bb6``: BB6 copula,
      - ``bb7``: BB7 copula,
      - ``bb8``: BB8 copula,
      - ``tawn``: Tawn copula,
      - ``tll``: Transformation local-likelihood (nonparametric) copula.

    See the :mod:`pyvinecopulib.families` module documentation for the
    named family-group constants (``parametric``, ``elliptical``,
    ``archimedean``, ``bb``, ``itau``, ...).
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
      .value("tawn", BicopFamily::tawn, bicopfamily_doc.tawn.doc)
      .value("tll", BicopFamily::tll, bicopfamily_doc.tll.doc)
      .export_values();

  module.attr("all") = bicop_families::all;
  module.attr("parametric") = bicop_families::parametric;
  module.attr("nonparametric") = bicop_families::nonparametric;
  module.attr("one_par") = bicop_families::one_par;
  module.attr("two_par") = bicop_families::two_par;
  module.attr("three_par") = bicop_families::three_par;
  module.attr("elliptical") = bicop_families::elliptical;
  module.attr("archimedean") = bicop_families::archimedean;
  module.attr("extreme_value") = bicop_families::extreme_value;
  module.attr("bb") = bicop_families::bb;
  module.attr("rotationless") = bicop_families::rotationless;
  module.attr("lt") = bicop_families::lt;
  module.attr("ut") = bicop_families::ut;
  module.attr("itau") = bicop_families::itau;
}
