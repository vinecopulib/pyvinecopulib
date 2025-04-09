#pragma once

#include "docstr.hpp"
#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
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
    .def(
      nb::init<std::vector<BicopFamily>,
               std::string,
               std::string,
               double,
               std::string,
               const Eigen::VectorXd&,
               double,
               bool,
               bool,
               size_t>(),
      fitcontrolsbicop_doc.ctor
        .doc_10args_family_set_parametric_method_nonparametric_method_nonparametric_mult_selection_criterion_weights_psi0_preselect_families_allow_rotations_num_threads,
      "family_set"_a = bicop_families::all,
      "parametric_method"_a = "mle",
      "nonparametric_method"_a = "constant",
      "nonparametric_mult"_a = 1.0,
      "selection_criterion"_a = "bic",
      "weights"_a = Eigen::VectorXd(),
      "psi0"_a = 0.9,
      "preselect_families"_a = true,
      "allow_rotations"_a = true,
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
    .def_prop_rw("allow_rotations",
                 &FitControlsBicop::get_allow_rotations,
                 &FitControlsBicop::set_allow_rotations,
                 "Whether to allow rotations for the families.")
    .def_prop_rw("num_threads",
                 &FitControlsBicop::get_num_threads,
                 &FitControlsBicop::set_num_threads,
                 "The number of threads.")
    .def("__repr__",
         [](const FitControlsBicop& controls) {
           return "<pyvinecopulib.FitControlsBicop>\n" + controls.str();
         })
    .def("str", &FitControlsBicop::str, fitcontrolsbicop_doc.str.doc)
    .def("__getstate__",
         [](const FitControlsBicop& controls) {
           return std::make_tuple(controls.get_family_set(),
                                  controls.get_parametric_method(),
                                  controls.get_nonparametric_method(),
                                  controls.get_nonparametric_mult(),
                                  controls.get_selection_criterion(),
                                  controls.get_weights(),
                                  controls.get_psi0(),
                                  controls.get_preselect_families(),
                                  controls.get_allow_rotations(),
                                  controls.get_num_threads());
         })
    .def("__setstate__",
         [](FitControlsBicop& controls,
            std::tuple<std::vector<BicopFamily>,
                       std::string,
                       std::string,
                       double,
                       std::string,
                       const Eigen::VectorXd&,
                       double,
                       bool,
                       bool,
                       size_t> state) {
           FitControlsConfig config;
           config.family_set = std::get<0>(state);
           config.parametric_method = std::get<1>(state);
           config.nonparametric_method = std::get<2>(state);
           config.nonparametric_mult = std::get<3>(state);
           config.selection_criterion = std::get<4>(state);
           config.weights = std::get<5>(state);
           config.psi0 = std::get<6>(state);
           config.preselect_families = std::get<7>(state);
           config.allow_rotations = std::get<8>(state);
           config.num_threads = std::get<9>(state);

           new (&controls) FitControlsBicop(config);
         });
}