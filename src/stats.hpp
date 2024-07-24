#ifndef STATS_H
#define STATS_H

#include "docstr.hpp"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vinecopulib.hpp>

namespace py = pybind11;
using namespace vinecopulib;

inline void
init_stats(py::module_& m)
{

  constexpr auto& doc = pyvinecopulib_doc;
  constexpr auto& tools_stat_doc = doc.vinecopulib.tools_stats;

  m.def("simulate_uniform",
        &tools_stats::simulate_uniform,
        tools_stat_doc.simulate_uniform.doc,
        py::arg("n"),
        py::arg("d"),
        py::arg("qrng") = false,
        py::arg("seeds") = std::vector<int>());

  m.def("sobol",
        &tools_stats::sobol,
        tools_stat_doc.sobol.doc,
        py::arg("n"),
        py::arg("d"),
        py::arg("seeds") = std::vector<int>());

  m.def("ghalton",
        &tools_stats::ghalton,
        tools_stat_doc.ghalton.doc,
        py::arg("n"),
        py::arg("d"),
        py::arg("seeds") = std::vector<int>());

  m.def("to_pseudo_obs",
        &tools_stats::to_pseudo_obs,
        tools_stat_doc.to_pseudo_obs.doc,
        py::arg("x"),
        py::arg("ties_method") = "average");
}

#endif // STATS_H