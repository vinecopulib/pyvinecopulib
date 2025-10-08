#pragma once

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <vinecopulib.hpp>
#include <wdm/eigen.hpp>

#include "docstr.hpp"

#include "docstr.hpp"

namespace nb = nanobind;
using namespace nb::literals;
using namespace vinecopulib;

const char* doc_wdm =
    R"""(Calculates (weighted) dependence measures.

This function computes various measures of dependence between two variables, optionally
using observation weights.

Parameters
----------
x, y :
    Input data vectors.
method :
    The dependence measure to compute. Possible values are:

    - ``"pearson"``, ``"prho"``, ``"cor"`` : Pearson correlation
    - ``"spearman"``, ``"srho"``, ``"rho"`` : Spearman’s :math:`\rho`
    - ``"kendall"``, ``"ktau"``, ``"tau"`` : Kendall’s :math:`\tau`
    - ``"blomqvist"``, ``"bbeta"``, ``"beta"`` : Blomqvist’s :math:`\beta`
    - ``"hoeffding"``, ``"hoeffd"``, ``"d"`` : Hoeffding’s :math:`D`
weights :
    Optional vector of observation weights.
remove_missing :
    If ``True``, all observations containing a ``NaN`` are removed. Otherwise, an error is raised
    if missing values are present.

Returns
-------
float
    The computed dependence measure.)""";

inline void init_stats(nb::module_& m) {
  constexpr auto& doc = pyvinecopulib_doc;
  constexpr auto& tools_stat_doc = doc.vinecopulib.tools_stats;

  m.def("simulate_uniform", &tools_stats::simulate_uniform,
        tools_stat_doc.simulate_uniform.doc, "n"_a, "d"_a, "qrng"_a = false,
        "seeds"_a = std::vector<int>(),
        nb::call_guard<nb::gil_scoped_release>());

  m.def("sobol", &tools_stats::sobol, tools_stat_doc.sobol.doc, "n"_a, "d"_a,
        "seeds"_a = std::vector<int>(),
        nb::call_guard<nb::gil_scoped_release>());

  m.def("ghalton", &tools_stats::ghalton, tools_stat_doc.ghalton.doc, "n"_a,
        "d"_a, "seeds"_a = std::vector<int>(),
        nb::call_guard<nb::gil_scoped_release>());

  m.def("to_pseudo_obs", &tools_stats::to_pseudo_obs,
        tools_stat_doc.to_pseudo_obs.doc, "x"_a, "ties_method"_a = "average",
        "weights"_a = Eigen::VectorXd(), "seeds"_a = std::vector<int>(),
        nb::call_guard<nb::gil_scoped_release>());
  m.def("wdm",
        static_cast<double (*)(const Eigen::VectorXd&, const Eigen::VectorXd&,
                               std::string, Eigen::VectorXd, bool)>(&wdm::wdm),
        doc_wdm, "x"_a, "y"_a, "method"_a, "weights"_a = Eigen::VectorXd(),
        "remove_missing"_a = true);
}
