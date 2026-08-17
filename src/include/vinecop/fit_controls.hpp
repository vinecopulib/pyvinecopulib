#pragma once

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <vinecopulib.hpp>

#include "docstr.hpp"

namespace nb = nanobind;
using namespace nb::literals;
using namespace vinecopulib;

// Factory function to create FitControlsVinecop from pair-copula controls
inline FitControlsVinecop fcv_from_bicop_controls(
    const FitControlsBicop& controls, size_t trunc_lvl,
    const std::string& tree_criterion, double threshold, bool select_trunc_lvl,
    bool select_threshold, bool select_families, bool show_trace,
    const std::string& tree_algorithm, const std::vector<int>& seeds) {
  return FitControlsVinecop(controls, trunc_lvl, tree_criterion, threshold,
                            select_trunc_lvl, select_threshold, select_families,
                            show_trace, tree_algorithm, seeds);
}

inline void init_vinecop_fit_controls(nb::module_& module) {
  constexpr auto& doc = pyvinecopulib_doc;
  constexpr auto& fitcontrolsvinecop_doc = doc.vinecopulib.FitControlsVinecop;
  // Inherited properties (`family_set`, `parametric_method`, …)
  // surface in the FitControlsBicop entry of `pyvinecopulib_doc`,
  // not under FitControlsVinecop. Reference the parent class's docs
  // for those.
  constexpr auto& fitcontrolsbicop_doc = doc.vinecopulib.FitControlsBicop;

  nb::class_<FitControlsVinecop>(module, "FitControlsVinecop",
                                 fitcontrolsvinecop_doc.doc)
      .def(nb::init<std::vector<BicopFamily>, std::string, std::string, double,
                    size_t, size_t, std::string, double, std::string,
                    const Eigen::VectorXd&, double, bool, bool, bool, bool,
                    bool, size_t, std::string, bool, std::vector<int>>(),
           "family_set"_a = bicop_families::all, "parametric_method"_a = "mle",
           "nonparametric_method"_a = "constant", "nonparametric_mult"_a = 1.0,
           "nonparametric_grid_size"_a = 30,
           "trunc_lvl"_a = std::numeric_limits<size_t>::max(),
           "tree_criterion"_a = "tau", "threshold"_a = 0.0,
           "selection_criterion"_a = "aic", "weights"_a = Eigen::VectorXd(),
           "psi0"_a = 0.9, "preselect_families"_a = true,
           "select_trunc_lvl"_a = false, "select_threshold"_a = false,
           "select_families"_a = true, "show_trace"_a = false,
           "num_threads"_a = 1, "tree_algorithm"_a = "mst_prim",
           "allow_rotations"_a = true, "seeds"_a = std::vector<int>(),
           fitcontrolsvinecop_doc.ctor.doc_20args,
           nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_bicop_controls", &fcv_from_bicop_controls, "controls"_a,
                  "trunc_lvl"_a = std::numeric_limits<size_t>::max(),
                  "tree_criterion"_a = "tau", "threshold"_a = 0.0,
                  "select_trunc_lvl"_a = false, "select_threshold"_a = false,
                  "select_families"_a = true, "show_trace"_a = false,
                  "tree_algorithm"_a = "mst_prim",
                  "seeds"_a = std::vector<int>(),
                  fitcontrolsvinecop_doc.ctor.doc_10args,
                  nb::call_guard<nb::gil_scoped_release>())
      .def_prop_rw("bicop_controls",
                   &FitControlsVinecop::get_fit_controls_bicop,
                   &FitControlsVinecop::set_fit_controls_bicop,
                   fitcontrolsvinecop_doc.get_fit_controls_bicop.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      // Inherited from FitControlsBicop — reference its docs.
      .def_prop_rw("family_set", &FitControlsVinecop::get_family_set,
                   &FitControlsVinecop::set_family_set,
                   fitcontrolsbicop_doc.get_family_set.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      .def_prop_rw("parametric_method",
                   &FitControlsVinecop::get_parametric_method,
                   &FitControlsVinecop::set_parametric_method,
                   fitcontrolsbicop_doc.get_parametric_method.doc)
      .def_prop_rw("nonparametric_method",
                   &FitControlsVinecop::get_nonparametric_method,
                   &FitControlsVinecop::set_nonparametric_method,
                   fitcontrolsbicop_doc.get_nonparametric_method.doc)
      .def_prop_rw("nonparametric_mult",
                   &FitControlsVinecop::get_nonparametric_mult,
                   &FitControlsVinecop::set_nonparametric_mult,
                   fitcontrolsbicop_doc.get_nonparametric_mult.doc)
      .def_prop_rw("nonparametric_grid_size",
                   &FitControlsVinecop::get_nonparametric_grid_size,
                   &FitControlsVinecop::set_nonparametric_grid_size,
                   fitcontrolsbicop_doc.get_nonparametric_grid_size.doc)
      // Vinecop-specific knobs.
      .def_prop_rw("trunc_lvl", &FitControlsVinecop::get_trunc_lvl,
                   &FitControlsVinecop::set_trunc_lvl,
                   fitcontrolsvinecop_doc.get_trunc_lvl.doc)
      .def_prop_rw("tree_criterion", &FitControlsVinecop::get_tree_criterion,
                   &FitControlsVinecop::set_tree_criterion,
                   fitcontrolsvinecop_doc.get_tree_criterion.doc)
      .def_prop_rw(
          "tree_criterion_function",
          [](const FitControlsVinecop& controls)
              -> std::optional<TreeCriterionFunction> {
            auto f = controls.get_tree_criterion_function();
            if (!f) {
              return std::nullopt;
            }
            return f;
          },
          [](FitControlsVinecop& controls,
             std::optional<TreeCriterionFunction> tree_criterion_function) {
            controls.set_tree_criterion_function(tree_criterion_function
                                                     ? *tree_criterion_function
                                                     : TreeCriterionFunction{});
          },
          fitcontrolsvinecop_doc.get_tree_criterion_function.doc)
      .def_prop_rw("threshold", &FitControlsVinecop::get_threshold,
                   &FitControlsVinecop::set_threshold,
                   fitcontrolsvinecop_doc.get_threshold.doc)
      // Inherited.
      .def_prop_rw("selection_criterion",
                   &FitControlsVinecop::get_selection_criterion,
                   &FitControlsVinecop::set_selection_criterion,
                   fitcontrolsbicop_doc.get_selection_criterion.doc)
      .def_prop_rw("weights", &FitControlsVinecop::get_weights,
                   &FitControlsVinecop::set_weights,
                   fitcontrolsbicop_doc.get_weights.doc)
      .def_prop_rw("psi0", &FitControlsVinecop::get_psi0,
                   &FitControlsVinecop::set_psi0,
                   fitcontrolsbicop_doc.get_psi0.doc)
      .def_prop_rw("preselect_families",
                   &FitControlsVinecop::get_preselect_families,
                   &FitControlsVinecop::set_preselect_families,
                   fitcontrolsbicop_doc.get_preselect_families.doc)
      // Vinecop-specific knobs.
      .def_prop_rw("select_trunc_lvl",
                   &FitControlsVinecop::get_select_trunc_lvl,
                   &FitControlsVinecop::set_select_trunc_lvl,
                   fitcontrolsvinecop_doc.get_select_trunc_lvl.doc)
      .def_prop_rw("select_threshold",
                   &FitControlsVinecop::get_select_threshold,
                   &FitControlsVinecop::set_select_threshold,
                   fitcontrolsvinecop_doc.get_select_threshold.doc)
      .def_prop_rw("select_families", &FitControlsVinecop::get_select_families,
                   &FitControlsVinecop::set_select_families,
                   fitcontrolsvinecop_doc.get_select_families.doc)
      .def_prop_rw("show_trace", &FitControlsVinecop::get_show_trace,
                   &FitControlsVinecop::set_show_trace,
                   fitcontrolsvinecop_doc.get_show_trace.doc)
      // Inherited.
      .def_prop_rw("num_threads", &FitControlsVinecop::get_num_threads,
                   &FitControlsVinecop::set_num_threads,
                   fitcontrolsbicop_doc.get_num_threads.doc)
      // Vinecop-specific.
      .def_prop_rw("tree_algorithm", &FitControlsVinecop::get_tree_algorithm,
                   &FitControlsVinecop::set_tree_algorithm,
                   fitcontrolsvinecop_doc.get_tree_algorithm.doc)
      // Inherited.
      .def_prop_rw("allow_rotations", &FitControlsVinecop::get_allow_rotations,
                   &FitControlsVinecop::set_allow_rotations,
                   fitcontrolsbicop_doc.get_allow_rotations.doc)
      // Vinecop-specific.
      .def_prop_rw("seeds", &FitControlsVinecop::get_seeds,
                   &FitControlsVinecop::set_seeds,
                   fitcontrolsvinecop_doc.get_seeds.doc)
      // Conditioning-aware structure selection: 1-based variable indices placed
      // at the tail of the vine order during ``Vinecop.select``. Not part of
      // the positional constructor; set it as a property.
      .def_prop_rw("conditioning_set",
                   &FitControlsVinecop::get_conditioning_set,
                   &FitControlsVinecop::set_conditioning_set,
                   fitcontrolsvinecop_doc.get_conditioning_set.doc)
      .def(
          "__repr__",
          [](const FitControlsVinecop& controls) {
            return "<pyvinecopulib.core.FitControlsVinecop>\n" + controls.str();
          },
          fitcontrolsvinecop_doc.str.doc)
      .def(
          "__str__",
          [](const FitControlsVinecop& controls) {
            return "<pyvinecopulib.core.FitControlsVinecop>\n" + controls.str();
          },
          fitcontrolsvinecop_doc.str.doc)
      .def("__getstate__",
           [](const FitControlsVinecop& controls) {
             nb::dict state;
             state["family_set"] = controls.get_family_set();
             state["parametric_method"] = controls.get_parametric_method();
             state["nonparametric_method"] =
                 controls.get_nonparametric_method();
             state["nonparametric_mult"] = controls.get_nonparametric_mult();
             state["nonparametric_grid_size"] =
                 controls.get_nonparametric_grid_size();
             state["trunc_lvl"] = controls.get_trunc_lvl();
             state["tree_criterion"] = controls.get_tree_criterion();
             // Empty std::function casts to None; a wrapped Python callable
             // casts back to the original object (picklable if module-level).
             state["tree_criterion_function"] =
                 controls.get_tree_criterion_function();
             state["threshold"] = controls.get_threshold();
             state["selection_criterion"] = controls.get_selection_criterion();
             state["weights"] = controls.get_weights();
             state["psi0"] = controls.get_psi0();
             state["preselect_families"] = controls.get_preselect_families();
             state["select_trunc_lvl"] = controls.get_select_trunc_lvl();
             state["select_threshold"] = controls.get_select_threshold();
             state["select_families"] = controls.get_select_families();
             state["show_trace"] = controls.get_show_trace();
             state["num_threads"] = controls.get_num_threads();
             state["tree_algorithm"] = controls.get_tree_algorithm();
             state["allow_rotations"] = controls.get_allow_rotations();
             state["seeds"] = controls.get_seeds();
             state["conditioning_set"] = controls.get_conditioning_set();
             return state;
           })

      .def("__setstate__", [](FitControlsVinecop& controls, nb::dict state) {
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
        config.trunc_lvl = nb::cast<std::size_t>(state["trunc_lvl"]);
        config.tree_criterion = nb::cast<std::string>(state["tree_criterion"]);
        // Guard the key so pre-feature pickles (without it) still load.
        if (state.contains("tree_criterion_function") &&
            !state["tree_criterion_function"].is_none()) {
          config.tree_criterion_function =
              nb::cast<TreeCriterionFunction>(state["tree_criterion_function"]);
        }
        config.threshold = nb::cast<double>(state["threshold"]);
        config.selection_criterion =
            nb::cast<std::string>(state["selection_criterion"]);
        config.weights = nb::cast<Eigen::VectorXd>(state["weights"]);
        config.psi0 = nb::cast<double>(state["psi0"]);
        config.preselect_families = nb::cast<bool>(state["preselect_families"]);
        config.select_trunc_lvl = nb::cast<bool>(state["select_trunc_lvl"]);
        config.select_threshold = nb::cast<bool>(state["select_threshold"]);
        config.select_families = nb::cast<bool>(state["select_families"]);
        config.show_trace = nb::cast<bool>(state["show_trace"]);
        config.num_threads = nb::cast<std::size_t>(state["num_threads"]);
        config.tree_algorithm = nb::cast<std::string>(state["tree_algorithm"]);
        config.allow_rotations = nb::cast<bool>(state["allow_rotations"]);
        config.seeds = nb::cast<std::vector<int>>(state["seeds"]);
        // Guard the key so pre-feature pickles (without it) still load.
        if (state.contains("conditioning_set")) {
          config.conditioning_set =
              nb::cast<std::vector<size_t>>(state["conditioning_set"]);
        }

        new (&controls) FitControlsVinecop(std::move(config));
      });
}
