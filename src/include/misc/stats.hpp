#pragma once

#include "docstr.hpp"
#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <vinecopulib.hpp>

namespace nb = nanobind;
using namespace vinecopulib;

inline void
init_stats(nb::module_& m)
{

  constexpr auto& doc = pyvinecopulib_doc;
  constexpr auto& tools_stat_doc = doc.vinecopulib.tools_stats;

  m.def("simulate_uniform",
        &tools_stats::simulate_uniform,
        tools_stat_doc.simulate_uniform.doc,
        nb::arg("n"),
        nb::arg("d"),
        nb::arg("qrng") = false,
        nb::arg("seeds") = std::vector<int>());

  m.def("sobol",
        &tools_stats::sobol,
        tools_stat_doc.sobol.doc,
        nb::arg("n"),
        nb::arg("d"),
        nb::arg("seeds") = std::vector<int>());

  m.def("ghalton",
        &tools_stats::ghalton,
        tools_stat_doc.ghalton.doc,
        nb::arg("n"),
        nb::arg("d"),
        nb::arg("seeds") = std::vector<int>());

  m.def("to_pseudo_obs",
        &tools_stats::to_pseudo_obs,
        tools_stat_doc.to_pseudo_obs.doc,
        nb::arg("x"),
        nb::arg("ties_method") = "average");
}
