#pragma once

#include "docstr.hpp"
#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <vinecopulib.hpp>

namespace nb = nanobind;
using namespace nb::literals;
using namespace vinecopulib;

inline void
init_bicop_fit_controls(nb::module_& module)
{

  constexpr auto& fitcontrolsbicop_doc =
    pyvinecopulib_doc.vinecopulib.FitControlsBicop;

  nb::class_<FitControlsBicop>(
    module, "FitControlsBicop", fitcontrolsbicop_doc.doc)
    .def(nb::init<std::vector<BicopFamily>,
                  std::string,
                  std::string,
                  double,
                  std::string,
                  const Eigen::VectorXd&,
                  double,
                  bool,
                  size_t>(),
         fitcontrolsbicop_doc.ctor.doc_9args,
         "family_set"_a = bicop_families::all,
         "parametric_method"_a = "mle",
         "nonparametric_method"_a = "constant",
         "nonparametric_mult"_a = 1.0,
         "selection_criterion"_a = "bic",
         "weights"_a = Eigen::VectorXd(),
         "psi0"_a = 0.9,
         "preselect_families"_a = true,
         "num_threads"_a = 1)
    /*      .def(nb::init<std::string>(), */
    //      "creates default controls except for the parameteric method.",
    //      "parametric_method"_a)
    // .def(nb::init<std::string, double>(),
    //      "creates default controls except for the nonparametric method.",
    /* "nonparametric_method"_a, "mult"_a = 1.0) */
    .def_prop_rw("family_set",
                 &FitControlsBicop::get_family_set,
                 &FitControlsBicop::set_family_set,
                 "The family set.")
    .def_prop_rw("parametric_method",
                 &FitControlsBicop::get_parametric_method,
                 &FitControlsBicop::set_parametric_method,
                 "The fit method for parametric families.")
    .def_prop_rw("nonparametric_method",
                 &FitControlsBicop::get_nonparametric_method,
                 &FitControlsBicop::set_nonparametric_method,
                 "The fit method for nonparametric families.")
    .def_prop_rw("nonparametric_mult",
                 &FitControlsBicop::get_nonparametric_mult,
                 &FitControlsBicop::set_nonparametric_method,
                 "The multiplier for the smoothing parameters.")
    .def_prop_rw("selection_criterion",
                 &FitControlsBicop::get_selection_criterion,
                 &FitControlsBicop::set_selection_criterion,
                 "The selection criterion.")
    .def_prop_rw("weights",
                 &FitControlsBicop::get_weights,
                 &FitControlsBicop::set_weights,
                 "The weights for the observations.")
    .def_prop_rw("psi0",
                 &FitControlsBicop::get_psi0,
                 &FitControlsBicop::set_psi0,
                 "The prior probability of non-independence.")
    .def_prop_rw("preselect_families",
                 &FitControlsBicop::get_preselect_families,
                 &FitControlsBicop::set_preselect_families,
                 "Whether to exclude families based on symmetry properties "
                 "of the data.")
    .def_prop_rw("num_threads",
                 &FitControlsBicop::get_num_threads,
                 &FitControlsBicop::set_num_threads,
                 "The number of threads.")
    .def("__repr__",
         [](const FitControlsBicop& controls) {
           return "<pyvinecopulib.FitControlsBicop>\n" + controls.str();
         })
    .def("str", &FitControlsBicop::str, fitcontrolsbicop_doc.str.doc);
}