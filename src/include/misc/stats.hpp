#pragma once

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <vinecopulib.hpp>
#include <vinecopulib/misc/tools_serialization.hpp>
#include <wdm/eigen.hpp>

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
    - ``"chatterjee"``, ``"cxi"``, ``"xi"`` : Chatterjee’s :math:`\xi`

    Every measure above is symmetric in ``x`` and ``y`` except Chatterjee’s
    :math:`\xi`, which measures how far ``y`` is a function of ``x``.
weights :
    Optional vector of observation weights.
remove_missing :
    If ``True``, all observations containing a ``NaN`` are removed. Otherwise, an error is raised
    if missing values are present.
seeds :
    Seeds for the random tie-breaking Chatterjee’s :math:`\xi` applies to ``x``; ignored by every
    other measure, and unused when ``x`` has no ties. Ordering tied predictors by ``y`` would
    manufacture dependence, so the ties are broken at random -- but from a fixed default, so
    :math:`\xi` is a function of its arguments. Pass seeds to vary that ordering or to average
    over it.

Returns
-------
float
    The computed dependence measure.)""";

inline void init_stats(nb::module_& m) {
  // The JSON-or-CBOR-by-extension file IO the compiled model classes use, so
  // the Python-side `Vinedist.to_file` / `from_file` write the same formats
  // through the same code rather than a second implementation.
  m.def(
      "_json_to_file",
      [](const std::string& filename, const std::string& json) {
        vinecopulib::tools_serialization::json_to_file(
            filename, nlohmann::json::parse(json));
      },
      "filename"_a, "json"_a,
      "Write a JSON string to a file, as CBOR when the filename ends in "
      "``.cbor``.\n\nParameters\n----------\nfilename : str\n"
      "    Path to write.\njson : str\n    The JSON string to write.");
  m.def(
      "_file_to_json",
      [](const std::string& filename) -> std::string {
        return vinecopulib::tools_serialization::file_to_json(filename).dump();
      },
      "filename"_a,
      "Read a file written by ``_json_to_file`` back into a JSON string."
      "\n\nParameters\n----------\nfilename : str\n    Path to read."
      "\n\nReturns\n-------\nstr\n    The payload as JSON.");

  constexpr auto& doc = pyvinecopulib_doc;
  constexpr auto& tools_stat_doc = doc.vinecopulib.tools_stats;

  m.def("sample_uniform", &tools_stats::simulate_uniform,
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
  m.def("find_latent_sample", &tools_stats::find_latent_sample,
        tools_stat_doc.find_latent_sample.doc, "u"_a, "b"_a, "niter"_a = 3,
        nb::call_guard<nb::gil_scoped_release>());

  m.def("wdm",
        static_cast<double (*)(const Eigen::VectorXd&, const Eigen::VectorXd&,
                               std::string, Eigen::VectorXd, bool,
                               std::vector<int>)>(&wdm::wdm),
        doc_wdm, "x"_a, "y"_a, "method"_a, "weights"_a = Eigen::VectorXd(),
        "remove_missing"_a = true, "seeds"_a = std::vector<int>(),
        nb::call_guard<nb::gil_scoped_release>());
}
