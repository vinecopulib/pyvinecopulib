#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vinecopulib.hpp>

namespace py = pybind11;
using namespace vinecopulib;

PYBIND11_MODULE(pyvinecopulib, pv) {
  pv.doc() = R"pbdoc(
        Pyvinecopulib library
        -----------------------

        .. currentmodule:: pyvinecopulib

        .. autosummary::
           :toctree: _generate

           simulate_uniform
    )pbdoc";

  py::class_<RVineStructure>(pv, "RVineStructure")
      .def(py::init<
               const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> &,
               bool>(),
           "creates an array representation of regular vine structures.",
           py::arg("mat"), py::arg("check") = true)
      .def(py::init<const std::vector<size_t> &, bool>(),
           "creates a D-vine with given ordering of variables.",
           py::arg("order"), py::arg("check") = true)
      .def(py::init<const std::string, bool>(),
           "creates a structure from a JSON file.", py::arg("filename"),
           py::arg("check") = true)
      .def("to_json", &RVineStructure::to_json,
           "writes the file in a JSON file.", py::arg("filename"))
      .def_property_readonly("dim", &RVineStructure::get_dim, "The dimension.")
      .def_property_readonly("trunc_lvl", &RVineStructure::get_trunc_lvl,
                             "The truncation level.")
      .def_property_readonly("order",
                             (std::vector<size_t>(RVineStructure::*)() const) &
                                 RVineStructure::get_order,
                             "The variable order.")
      .def("struct_array", &RVineStructure::struct_array, py::arg("tree"),
           py::arg("edge"), py::arg("natural_order") = false)
      .def("truncate", &RVineStructure::truncate, "The truncation level.",
           py::arg("trunc_lvl"))
      .def_static("simulate", &RVineStructure::simulate,
                  "simulates a random R-vine array.", py::arg("d"),
                  py::arg("natural order") = false,
                  py::arg("seeds") = std::vector<size_t>())
      .def("__repr__", [](const RVineStructure &rvs) {
        return "<pyvinecopulib.RVineStructure>\n" + rvs.str();
      });

#ifdef VERSION_INFO
  pv.attr("__version__") = VERSION_INFO;
#else
  pv.attr("__version__") = "dev";
#endif
}
