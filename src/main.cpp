#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vinecopulib.hpp>

namespace py = pybind11;

PYBIND11_MODULE(pyvinecopulib, pv) {
  pv.doc() = R"pbdoc(
        Pyvinecopulib library
        -----------------------

        .. currentmodule:: pyvinecopulib

        .. autosummary::
           :toctree: _generate

           simulate_uniform
    )pbdoc";

  pv.def("simulate_uniform", &vinecopulib::tools_stats::simulate_uniform,
         R"pbdoc(
        Simulate uniform random numbers.

        Simulate a matrix of random numbers, with an option to get quasi random
        numbers or a seed.
    )pbdoc");
#ifdef VERSION_INFO
  pv.attr("__version__") = VERSION_INFO;
#else
  pv.attr("__version__") = "dev";
#endif
}
