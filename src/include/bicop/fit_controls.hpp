#pragma once

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <vinecopulib.hpp>

#include "docstr.hpp"

namespace nb = nanobind;
using namespace nb::literals;
using namespace vinecopulib;

inline void init_bicop_fit_controls(nb::module_& module) {
  constexpr auto& fitcontrolsbicop_doc =
      pyvinecopulib_doc.vinecopulib.FitControlsBicop;

  nb::class_<FitControlsBicop>(module, "FitControlsBicop",
                               fitcontrolsbicop_doc.doc)
      .def(
          nb::init<std::vector<BicopFamily>, std::string, std::string, double,
                   size_t, std::string, const Eigen::VectorXd&, double, bool,
                   bool, size_t>(),
          fitcontrolsbicop_doc.ctor
              .doc_11args_family_set_parametric_method_nonparametric_method_nonparametric_mult_nonparametric_grid_size_selection_criterion_weights_psi0_preselect_families_allow_rotations_num_threads,
          "family_set"_a = bicop_families::all, "parametric_method"_a = "mle",
          "nonparametric_method"_a = "constant", "nonparametric_mult"_a = 1.0,
          "nonparametric_grid_size"_a = 30, "selection_criterion"_a = "bic",
          "weights"_a = Eigen::VectorXd(), "psi0"_a = 0.9,
          "preselect_families"_a = true, "allow_rotations"_a = true,
          "num_threads"_a = 1, nb::call_guard<nb::gil_scoped_release>())
      /*      .def(nb::init<std::string>(), */
      //      "creates default controls except for the parameteric method.",
      //      "parametric_method"_a)
      // .def(nb::init<std::string, double>(),
      //      "creates default controls except for the nonparametric method.",
      /* "nonparametric_method"_a, "mult"_a = 1.0) */
      .def_prop_rw("family_set", &FitControlsBicop::get_family_set,
                   &FitControlsBicop::set_family_set, "The family set.",
                   nb::call_guard<nb::gil_scoped_release>())
      .def_prop_rw("parametric_method",
                   &FitControlsBicop::get_parametric_method,
                   &FitControlsBicop::set_parametric_method,
                   "The fit method for parametric families.")
      .def_prop_rw("nonparametric_method",
                   &FitControlsBicop::get_nonparametric_method,
                   &FitControlsBicop::set_nonparametric_method,
                   "The fit method for nonparametric families.")
      .def_prop_rw(
          "nonparametric_mult", &FitControlsBicop::get_nonparametric_mult,
          &FitControlsBicop::set_nonparametric_mult,
          "A factor with which the smoothing parameters are multiplied.")
      .def_prop_rw("nonparametric_grid_size",
                   &FitControlsBicop::get_nonparametric_grid_size,
                   &FitControlsBicop::set_nonparametric_grid_size,
                   "The grid size for the post-estimation interpolation in "
                   "nonparametric models.")
      .def_prop_rw("selection_criterion",
                   &FitControlsBicop::get_selection_criterion,
                   &FitControlsBicop::set_selection_criterion,
                   "The selection criterion.")
      .def_prop_rw("weights", &FitControlsBicop::get_weights,
                   &FitControlsBicop::set_weights,
                   "The weights for the observations.")
      .def_prop_rw("psi0", &FitControlsBicop::get_psi0,
                   &FitControlsBicop::set_psi0,
                   "The prior probability of non-independence.")
      .def_prop_rw("preselect_families",
                   &FitControlsBicop::get_preselect_families,
                   &FitControlsBicop::set_preselect_families,
                   "Whether to exclude families based on symmetry properties "
                   "of the data.")
      .def_prop_rw("allow_rotations", &FitControlsBicop::get_allow_rotations,
                   &FitControlsBicop::set_allow_rotations,
                   "Whether to allow rotations for the families.")
      .def_prop_rw("num_threads", &FitControlsBicop::get_num_threads,
                   &FitControlsBicop::set_num_threads, "The number of threads.")
      .def(
          "__repr__",
          [](const FitControlsBicop& controls) {
            return "<pyvinecopulib.FitControlsBicop>\n" + controls.str();
          },
          fitcontrolsbicop_doc.str.doc)
      .def(
          "__str__",
          [](const FitControlsBicop& controls) {
            return "<pyvinecopulib.FitControlsBicop>\n" + controls.str();
          },
          fitcontrolsbicop_doc.str.doc)
      .def("__getstate__",
           [](const FitControlsBicop& controls) {
             nb::dict state;
             state["family_set"] = controls.get_family_set();
             state["parametric_method"] = controls.get_parametric_method();
             state["nonparametric_method"] =
                 controls.get_nonparametric_method();
             state["nonparametric_mult"] = controls.get_nonparametric_mult();
             state["nonparametric_grid_size"] =
                 controls.get_nonparametric_grid_size();
             state["selection_criterion"] = controls.get_selection_criterion();
             state["weights"] = controls.get_weights();
             state["psi0"] = controls.get_psi0();
             state["preselect_families"] = controls.get_preselect_families();
             state["allow_rotations"] = controls.get_allow_rotations();
             state["num_threads"] = controls.get_num_threads();
             return state;
           })

      .def("__setstate__", [](FitControlsBicop& controls, nb::dict state) {
        FitControlsConfig config;
        config.family_set =
            nb::cast<std::vector<BicopFamily>>(state["family_set"]);
        config.parametric_method =
            nb::cast<std::string>(state["parametric_method"]);
        config.nonparametric_method =
            nb::cast<std::string>(state["nonparametric_method"]);
        config.nonparametric_mult =
            nb::cast<double>(state["nonparametric_mult"]);
        config.nonparametric_grid_size =
            nb::cast<std::size_t>(state["nonparametric_grid_size"]);
        config.selection_criterion =
            nb::cast<std::string>(state["selection_criterion"]);
        config.weights = nb::cast<Eigen::VectorXd>(state["weights"]);
        config.psi0 = nb::cast<double>(state["psi0"]);
        config.preselect_families = nb::cast<bool>(state["preselect_families"]);
        config.allow_rotations = nb::cast<bool>(state["allow_rotations"]);
        config.num_threads = nb::cast<std::size_t>(state["num_threads"]);

        new (&controls) FitControlsBicop(std::move(config));
      });
}
