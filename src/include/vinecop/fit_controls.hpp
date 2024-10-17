#pragma once

#include "docstr.hpp"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vinecopulib.hpp>

namespace py = pybind11;
using namespace vinecopulib;

inline void
init_vinecop_fit_controls(py::module_& module)
{

  constexpr auto& doc = pyvinecopulib_doc;
  constexpr auto& fitcontrolsvinecop_doc = doc.vinecopulib.FitControlsVinecop;

  py::class_<FitControlsVinecop>(
    module, "FitControlsVinecop", fitcontrolsvinecop_doc.doc)
    .def(py::init<std::vector<BicopFamily>,
                  std::string,
                  std::string,
                  double,
                  size_t,
                  std::string,
                  double,
                  std::string,
                  const Eigen::VectorXd&,
                  double,
                  bool,
                  bool,
                  bool,
                  bool,
                  size_t>(),
         py::arg("family_set") = bicop_families::all,
         py::arg("parametric_method") = "mle",
         py::arg("nonparametric_method") = "quadratic",
         py::arg("nonparametric_mult") = 1.0,
         py::arg("trunc_lvl") = std::numeric_limits<size_t>::max(),
         py::arg("tree_criterion") = "tau",
         py::arg("threshold") = 0.0,
         py::arg("selection_criterion") = "bic",
         py::arg("weights") = Eigen::VectorXd(),
         py::arg("psi0") = 0.9,
         py::arg("preselect_families") = true,
         py::arg("select_trunc_lvl") = false,
         py::arg("select_threshold") = false,
         py::arg("show_trace") = false,
         py::arg("num_threads") = 1,
         fitcontrolsvinecop_doc.ctor.doc_17args)
    .def_property("family_set",
                  &FitControlsVinecop::get_family_set,
                  &FitControlsVinecop::set_family_set,
                  "The family set.")
    .def_property("parametric_method",
                  &FitControlsVinecop::get_parametric_method,
                  &FitControlsVinecop::set_parametric_method,
                  "The fit method for parametric families.")
    .def_property("nonparametric_method",
                  &FitControlsVinecop::get_nonparametric_method,
                  &FitControlsVinecop::set_nonparametric_method,
                  "The fit method for nonparametric families.")
    .def_property("nonparametric_mult",
                  &FitControlsVinecop::get_nonparametric_mult,
                  &FitControlsVinecop::set_nonparametric_method,
                  "The multiplier for the smoothing parameters.")
    .def_property("trunc_lvl",
                  &FitControlsVinecop::get_trunc_lvl,
                  &FitControlsVinecop::set_trunc_lvl,
                  "The truncation level.")
    .def_property("tree_criterion",
                  &FitControlsVinecop::get_tree_criterion,
                  &FitControlsVinecop::set_tree_criterion,
                  "The tree criterion.")
    .def_property("threshold",
                  &FitControlsVinecop::get_threshold,
                  &FitControlsVinecop::set_threshold,
                  "The threshold.")
    .def_property("selection_criterion",
                  &FitControlsVinecop::get_selection_criterion,
                  &FitControlsVinecop::set_selection_criterion,
                  "The selection criterion.")
    .def_property("weights",
                  &FitControlsVinecop::get_weights,
                  &FitControlsVinecop::set_weights,
                  "The weights for the observations.")
    .def_property("psi0",
                  &FitControlsVinecop::get_psi0,
                  &FitControlsVinecop::set_psi0,
                  "The prior probability of non-independence.")
    .def_property("preselect_families",
                  &FitControlsVinecop::get_preselect_families,
                  &FitControlsVinecop::set_preselect_families,
                  "Whether to exclude families based on symmetry properties "
                  "of the data.")
    .def_property("select_trunc_lvl",
                  &FitControlsVinecop::get_select_trunc_lvl,
                  &FitControlsVinecop::set_select_trunc_lvl,
                  "Whether to select the truncation level.")
    .def_property("select_threshold",
                  &FitControlsVinecop::get_select_threshold,
                  &FitControlsVinecop::set_select_threshold,
                  "Whether to select the threshold.")
    .def_property("show_trace",
                  &FitControlsVinecop::get_show_trace,
                  &FitControlsVinecop::set_show_trace,
                  "Whether to show the trace.")
    .def_property("num_threads",
                  &FitControlsVinecop::get_num_threads,
                  &FitControlsVinecop::set_num_threads,
                  "The number of threads.")
    .def("__repr__",
         [](const FitControlsVinecop& controls) {
           return "<pyvinecopulib.FitControlsVinecop>\n" + controls.str();
         })
    .def("str", &FitControlsVinecop::str, fitcontrolsvinecop_doc.str.doc);
}