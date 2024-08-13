#pragma once

#include "docstr.hpp"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vinecopulib.hpp>

namespace py = pybind11;
using namespace vinecopulib;

inline void
init_bicop_fit_controls(py::module_& module)
{

  constexpr auto& doc = pyvinecopulib_doc;
  constexpr auto& fitcontrolsbicop_doc = doc.vinecopulib.FitControlsBicop;

  py::class_<FitControlsBicop>(module, "FitControlsBicop", fitcontrolsbicop_doc.doc)
    .def(py::init<std::vector<BicopFamily>,
                  std::string,
                  std::string,
                  double,
                  std::string,
                  const Eigen::VectorXd&,
                  double,
                  bool,
                  size_t>(),
         fitcontrolsbicop_doc.ctor.doc_9args,
         py::arg("family_set") = bicop_families::all,
         py::arg("parametric_method") = "mle",
         py::arg("nonparametric_method") = "quadratic",
         py::arg("nonparametric_mult") = 1.0,
         py::arg("selection_criterion") = "bic",
         py::arg("weights") = Eigen::VectorXd(),
         py::arg("psi0") = 0.9,
         py::arg("preselect_families") = true,
         py::arg("num_threads") = 1)
    /*      .def(py::init<std::string>(), */
    //      "creates default controls except for the parameteric method.",
    //      py::arg("parametric_method"))
    // .def(py::init<std::string, double>(),
    //      "creates default controls except for the nonparametric method.",
    /* py::arg("nonparametric_method"), py::arg("mult") = 1.0) */
    .def_property("family_set",
                  &FitControlsBicop::get_family_set,
                  &FitControlsBicop::set_family_set,
                  "The family set.")
    .def_property("parametric_method",
                  &FitControlsBicop::get_parametric_method,
                  &FitControlsBicop::set_parametric_method,
                  "The fit method for parametric families.")
    .def_property("nonparametric_method",
                  &FitControlsBicop::get_nonparametric_method,
                  &FitControlsBicop::set_nonparametric_method,
                  "The fit method for nonparametric families.")
    .def_property("nonparametric_mult",
                  &FitControlsBicop::get_nonparametric_mult,
                  &FitControlsBicop::set_nonparametric_method,
                  "The multiplier for the smoothing parameters.")
    .def_property("selection_criterion",
                  &FitControlsBicop::get_selection_criterion,
                  &FitControlsBicop::set_selection_criterion,
                  "The selection criterion.")
    .def_property("weights",
                  &FitControlsBicop::get_weights,
                  &FitControlsBicop::set_weights,
                  "The weights for the observations.")
    .def_property("psi0",
                  &FitControlsBicop::get_psi0,
                  &FitControlsBicop::set_psi0,
                  "The prior probability of non-independence.")
    .def_property("preselect_families",
                  &FitControlsBicop::get_preselect_families,
                  &FitControlsBicop::set_preselect_families,
                  "Whether to exclude families based on symmetry properties "
                  "(see ``FitControlsBicop``)")
    .def_property("num_threads",
                  &FitControlsBicop::get_num_threads,
                  &FitControlsBicop::set_num_threads,
                  "The number of threads.")
    .def("__repr__",
         [](const FitControlsBicop& controls) {
           return "<pyvinecopulib.FitControlsBicop>\n" + controls.str();
         })
    .def("str", &FitControlsBicop::str, fitcontrolsbicop_doc.str.doc);
}