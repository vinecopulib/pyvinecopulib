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

  py::enum_<vinecopulib::BicopFamily>(pv, "BicopFamily", py::arithmetic())
      .value("indep", vinecopulib::BicopFamily::indep)
      .value("gaussian", vinecopulib::BicopFamily::gaussian)
      .value("student", vinecopulib::BicopFamily::student)
      .value("clayton", vinecopulib::BicopFamily::clayton)
      .value("gumbel", vinecopulib::BicopFamily::gumbel)
      .value("frank", vinecopulib::BicopFamily::frank)
      .value("joe", vinecopulib::BicopFamily::joe)
      .value("bb1", vinecopulib::BicopFamily::bb1)
      .value("bb6", vinecopulib::BicopFamily::bb6)
      .value("bb7", vinecopulib::BicopFamily::bb7)
      .value("bb8", vinecopulib::BicopFamily::bb8)
      .value("tll", vinecopulib::BicopFamily::tll);

  py::class_<vinecopulib::Bicop> bicop(pv, "Bicop");
  bicop.def(py::init<const vinecopulib::BicopFamily, const int,
                     const Eigen::MatrixXd &>());
  bicop.def_property_readonly("rotation", &vinecopulib::Bicop::get_rotation);
  bicop.def_property_readonly("parameters",
                              &vinecopulib::Bicop::get_parameters);
  bicop.def_property_readonly("family", &vinecopulib::Bicop::get_family);
#ifdef VERSION_INFO
  pv.attr("__version__") = VERSION_INFO;
#else
  pv.attr("__version__") = "dev";
#endif
}
