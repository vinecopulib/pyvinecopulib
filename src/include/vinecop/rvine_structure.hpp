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
init_vinecop_rvine_structure(nb::module_& module)
{

  constexpr auto& doc = pyvinecopulib_doc;
  constexpr auto& rvinestructure_doc = doc.vinecopulib.RVineStructure;
  constexpr auto& dvinestructure_doc = doc.vinecopulib.DVineStructure;
  constexpr auto& cvinestructure_doc = doc.vinecopulib.CVineStructure;

  nb::class_<RVineStructure>(module, "RVineStructure", rvinestructure_doc.doc)
    .def(nb::init<const size_t&, const size_t&>(),
         nb::arg("d") = static_cast<size_t>(1),
         nb::arg("trunc_lvl") = std::numeric_limits<size_t>::max(),
         rvinestructure_doc.ctor.doc_2args_d_trunc_lvl)
    .def(nb::init<const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>&,
                  bool>(),
         nb::arg("mat"),
         nb::arg("check") = true,
         rvinestructure_doc.ctor.doc_2args_mat_check)
    .def(nb::init<const std::vector<size_t>&, const size_t&, bool>(),
         nb::arg("order"),
         nb::arg("trunc_lvl") = std::numeric_limits<size_t>::max(),
         nb::arg("check") = true,
         rvinestructure_doc.ctor.doc_3args_order_trunc_lvl_check)
    .def(nb::init<const std::string, bool>(),
         nb::arg("filename"),
         nb::arg("check") = true,
         rvinestructure_doc.ctor.doc_2args_filename_check)
    .def("to_json",
         &RVineStructure::to_file,
         nb::arg("filename"),
         rvinestructure_doc.to_file.doc)
    .def_prop_ro("dim", &RVineStructure::get_dim, "The dimension.")
    .def_prop_ro(
      "trunc_lvl", &RVineStructure::get_trunc_lvl, "The truncation level.")
    .def_prop_ro("order",
                 (std::vector<size_t>(RVineStructure::*)() const) &
                   RVineStructure::get_order,
                 "The variable order.")
    .def_prop_ro(
      "matrix", &RVineStructure::get_matrix, rvinestructure_doc.get_matrix.doc)
    .def("struct_array",
         &RVineStructure::struct_array,
         nb::arg("tree"),
         nb::arg("edge"),
         nb::arg("natural_order") = false,
         rvinestructure_doc.struct_array.doc)
    .def("min_array",
         &RVineStructure::min_array,
         nb::arg("tree"),
         nb::arg("edge"),
         rvinestructure_doc.min_array.doc)
    .def("needed_hfunc1",
         &RVineStructure::needed_hfunc1,
         nb::arg("tree"),
         nb::arg("edge"),
         rvinestructure_doc.needed_hfunc1.doc)
    .def("needed_hfunc2",
         &RVineStructure::needed_hfunc2,
         nb::arg("tree"),
         nb::arg("edge"),
         rvinestructure_doc.needed_hfunc2.doc)
    .def("truncate",
         &RVineStructure::truncate,
         nb::arg("trunc_lvl"),
         rvinestructure_doc.truncate.doc)
    .def_static("simulate",
                &RVineStructure::simulate,
                nb::arg("d"),
                nb::arg("natural order") = false,
                nb::arg("seeds") = std::vector<size_t>(),
                rvinestructure_doc.simulate.doc)
    .def("__repr__",
         [](const RVineStructure& rvs) {
           return "<pyvinecopulib.RVineStructure>\n" + rvs.str();
         })
    .def("str", &RVineStructure::str, rvinestructure_doc.str.doc);

  nb::class_<DVineStructure, RVineStructure>(
    module, "DVineStructure", dvinestructure_doc.doc)
    .def(nb::init<const std::vector<size_t>&>(),
         nb::arg("order"),
         dvinestructure_doc.ctor.doc_1args)
    .def(nb::init<const std::vector<size_t>&, size_t>(),
         nb::arg("order"),
         nb::arg("trunc_lvl"),
         dvinestructure_doc.ctor.doc_2args)
    .def("__repr__", [](const DVineStructure& rvs) {
      return "<pyvinecopulib.DVineStructure>\n" + rvs.str();
    });

  nb::class_<CVineStructure, RVineStructure>(
    module, "CVineStructure", cvinestructure_doc.doc)
    .def(nb::init<const std::vector<size_t>&>(),
         cvinestructure_doc.ctor.doc_1args,
         nb::arg("order"))
    .def(nb::init<const std::vector<size_t>&, size_t>(),
         nb::arg("order"),
         nb::arg("trunc_lvl"),
         cvinestructure_doc.ctor.doc_2args)
    .def("__repr__", [](const CVineStructure& rvs) {
      return "<pyvinecopulib.CVineStructure>\n" + rvs.str();
    });
}