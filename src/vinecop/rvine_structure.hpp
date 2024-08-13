#pragma once

#include "docstr.hpp"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vinecopulib.hpp>

namespace py = pybind11;
using namespace vinecopulib;

inline void
init_vinecop_rvine_structure(py::module_& module)
{

  constexpr auto& doc = pyvinecopulib_doc;
  constexpr auto& rvinestructure_doc = doc.vinecopulib.RVineStructure;
  constexpr auto& dvinestructure_doc = doc.vinecopulib.DVineStructure;
  constexpr auto& cvinestructure_doc = doc.vinecopulib.CVineStructure;

  py::class_<RVineStructure>(module, "RVineStructure", rvinestructure_doc.doc)
    .def(py::init<const size_t&, const size_t&>(),
         py::arg("d") = static_cast<size_t>(1),
         py::arg("trunc_lvl") = std::numeric_limits<size_t>::max(),
         rvinestructure_doc.ctor.doc_2args_d_trunc_lvl)
    .def(py::init<const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>&,
                  bool>(),
         py::arg("mat"),
         py::arg("check") = true,
         rvinestructure_doc.ctor.doc_2args_mat_check)
    .def(py::init<const std::vector<size_t>&, const size_t&, bool>(),
         py::arg("order"),
         py::arg("trunc_lvl") = std::numeric_limits<size_t>::max(),
         py::arg("check") = true,
         rvinestructure_doc.ctor.doc_3args_order_trunc_lvl_check)
    .def(py::init<const std::string, bool>(),
         py::arg("filename"),
         py::arg("check") = true,
         rvinestructure_doc.ctor.doc_2args_filename_check)
    .def("to_json",
         &RVineStructure::to_file,
         py::arg("filename"),
         rvinestructure_doc.to_file.doc)
    .def_property_readonly("dim", &RVineStructure::get_dim, "The dimension.")
    .def_property_readonly(
      "trunc_lvl", &RVineStructure::get_trunc_lvl, "The truncation level.")
    .def_property_readonly("order",
                           (std::vector<size_t>(RVineStructure::*)() const) &
                             RVineStructure::get_order,
                           "The variable order.")
    .def("struct_array",
         &RVineStructure::struct_array,
         py::arg("tree"),
         py::arg("edge"),
         py::arg("natural_order") = false,
         rvinestructure_doc.struct_array.doc)
    .def("truncate",
         &RVineStructure::truncate,
         py::arg("trunc_lvl"),
         rvinestructure_doc.truncate.doc)
    .def_static("simulate",
                &RVineStructure::simulate,
                py::arg("d"),
                py::arg("natural order") = false,
                py::arg("seeds") = std::vector<size_t>(),
                rvinestructure_doc.simulate.doc)
    .def("__repr__",
         [](const RVineStructure& rvs) {
           return "<pyvinecopulib.RVineStructure>\n" + rvs.str();
         })
    .def("str", &RVineStructure::str, rvinestructure_doc.str.doc);

  py::class_<DVineStructure, RVineStructure>(
    module, "DVineStructure", dvinestructure_doc.doc)
    .def(py::init<const std::vector<size_t>&>(),
         py::arg("order"),
         dvinestructure_doc.ctor.doc_1args)
    .def(py::init<const std::vector<size_t>&, size_t>(),
         py::arg("order"),
         py::arg("trunc_lvl"),
         dvinestructure_doc.ctor.doc_2args)
    .def("__repr__", [](const DVineStructure& rvs) {
      return "<pyvinecopulib.DVineStructure>\n" + rvs.str();
    });

  py::class_<CVineStructure, RVineStructure>(
    module, "CVineStructure", cvinestructure_doc.doc)
    .def(py::init<const std::vector<size_t>&>(),
         cvinestructure_doc.ctor.doc_1args,
         py::arg("order"))
    .def(py::init<const std::vector<size_t>&, size_t>(),
         py::arg("order"),
         py::arg("trunc_lvl"),
         cvinestructure_doc.ctor.doc_2args)
    .def("__repr__", [](const CVineStructure& rvs) {
      return "<pyvinecopulib.CVineStructure>\n" + rvs.str();
    });
}