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
init_vinecop_fit_controls(nb::module_& module)
{

  constexpr auto& doc = pyvinecopulib_doc;
  constexpr auto& fitcontrolsvinecop_doc = doc.vinecopulib.FitControlsVinecop;

  nb::class_<FitControlsVinecop>(
    module, "FitControlsVinecop", fitcontrolsvinecop_doc.doc)
    .def(nb::init<std::vector<BicopFamily>,
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
                  bool,
                  size_t,
                  std::string,
                  bool,
                  std::vector<int>>(),
         "family_set"_a = bicop_families::all,
         "parametric_method"_a = "mle",
         "nonparametric_method"_a = "constant",
         "nonparametric_mult"_a = 1.0,
         "trunc_lvl"_a = std::numeric_limits<size_t>::max(),
         "tree_criterion"_a = "tau",
         "threshold"_a = 0.0,
         "selection_criterion"_a = "bic",
         "weights"_a = Eigen::VectorXd(),
         "psi0"_a = 0.9,
         "preselect_families"_a = true,
         "select_trunc_lvl"_a = false,
         "select_threshold"_a = false,
         "select_families"_a = true,
         "show_trace"_a = false,
         "num_threads"_a = 1,
         "tree_algorithm"_a = "mst_prim",
         "allow_rotations"_a = true,
         "seeds"_a = std::vector<int>(),
         fitcontrolsvinecop_doc.ctor.doc_19args)
    .def_prop_rw("family_set",
                 &FitControlsVinecop::get_family_set,
                 &FitControlsVinecop::set_family_set,
                 "The family set.")
    .def_prop_rw("parametric_method",
                 &FitControlsVinecop::get_parametric_method,
                 &FitControlsVinecop::set_parametric_method,
                 "The fit method for parametric families.")
    .def_prop_rw("nonparametric_method",
                 &FitControlsVinecop::get_nonparametric_method,
                 &FitControlsVinecop::set_nonparametric_method,
                 "The fit method for nonparametric families.")
    .def_prop_rw("nonparametric_mult",
                 &FitControlsVinecop::get_nonparametric_mult,
                 &FitControlsVinecop::set_nonparametric_method,
                 "The multiplier for the smoothing parameters.")
    .def_prop_rw("trunc_lvl",
                 &FitControlsVinecop::get_trunc_lvl,
                 &FitControlsVinecop::set_trunc_lvl,
                 "The truncation level.")
    .def_prop_rw("tree_criterion",
                 &FitControlsVinecop::get_tree_criterion,
                 &FitControlsVinecop::set_tree_criterion,
                 "The tree criterion.")
    .def_prop_rw("threshold",
                 &FitControlsVinecop::get_threshold,
                 &FitControlsVinecop::set_threshold,
                 "The threshold.")
    .def_prop_rw("selection_criterion",
                 &FitControlsVinecop::get_selection_criterion,
                 &FitControlsVinecop::set_selection_criterion,
                 "The selection criterion.")
    .def_prop_rw("weights",
                 &FitControlsVinecop::get_weights,
                 &FitControlsVinecop::set_weights,
                 "The weights for the observations.")
    .def_prop_rw("psi0",
                 &FitControlsVinecop::get_psi0,
                 &FitControlsVinecop::set_psi0,
                 "The prior probability of non-independence.")
    .def_prop_rw("preselect_families",
                 &FitControlsVinecop::get_preselect_families,
                 &FitControlsVinecop::set_preselect_families,
                 "Whether to exclude families based on symmetry properties "
                 "of the data.")
    .def_prop_rw("select_trunc_lvl",
                 &FitControlsVinecop::get_select_trunc_lvl,
                 &FitControlsVinecop::set_select_trunc_lvl,
                 "Whether to select the truncation level.")
    .def_prop_rw("select_threshold",
                 &FitControlsVinecop::get_select_threshold,
                 &FitControlsVinecop::set_select_threshold,
                 "Whether to select the threshold.")
    .def_prop_rw("select_families",
                 &FitControlsVinecop::get_select_families,
                 &FitControlsVinecop::set_select_families,
                 "Whether to select the families.")
    .def_prop_rw("show_trace",
                 &FitControlsVinecop::get_show_trace,
                 &FitControlsVinecop::set_show_trace,
                 "Whether to show the trace.")
    .def_prop_rw("num_threads",
                 &FitControlsVinecop::get_num_threads,
                 &FitControlsVinecop::set_num_threads,
                 "The number of threads.")
    .def_prop_rw("tree_algorithm",
                 &FitControlsVinecop::get_tree_algorithm,
                 &FitControlsVinecop::set_tree_algorithm,
                 "The spanning tree algorithm.")
    .def_prop_rw("allow_rotations",
                 &FitControlsVinecop::get_allow_rotations,
                 &FitControlsVinecop::set_allow_rotations,
                 "Whether to allow rotations for the families.")
    .def_prop_rw("seeds",
                 &FitControlsVinecop::get_seeds,
                 &FitControlsVinecop::set_seeds,
                 "The seeds for the random number generator.")
    .def("__repr__",
         [](const FitControlsVinecop& controls) {
           return "<pyvinecopulib.FitControlsVinecop>\n" + controls.str();
         })
    .def("str", &FitControlsVinecop::str, fitcontrolsvinecop_doc.str.doc)
    .def("__getstate__",
         [](const FitControlsVinecop& controls) {
           return std::make_tuple(controls.get_family_set(),
                                  controls.get_parametric_method(),
                                  controls.get_nonparametric_method(),
                                  controls.get_nonparametric_mult(),
                                  controls.get_trunc_lvl(),
                                  controls.get_tree_criterion(),
                                  controls.get_threshold(),
                                  controls.get_selection_criterion(),
                                  controls.get_weights(),
                                  controls.get_psi0(),
                                  controls.get_preselect_families(),
                                  controls.get_select_trunc_lvl(),
                                  controls.get_select_threshold(),
                                  controls.get_select_families(),
                                  controls.get_show_trace(),
                                  controls.get_num_threads(),
                                  controls.get_tree_algorithm(),
                                  controls.get_allow_rotations(),
                                  controls.get_seeds());
         })
    .def("__setstate__",
         [](FitControlsVinecop& controls,
            std::tuple<std::vector<BicopFamily>,
                       std::string,
                       std::string,
                       double,
                       size_t,
                       std::string,
                       double,
                       std::string,
                       Eigen::VectorXd,
                       double,
                       bool,
                       bool,
                       bool,
                       bool,
                       bool,
                       size_t,
                       std::string,
                       bool,
                       std::vector<int>> state) {
           FitControlsConfig config;
           config.family_set = std::get<0>(state);
           config.parametric_method = std::get<1>(state);
           config.nonparametric_method = std::get<2>(state);
           config.nonparametric_mult = std::get<3>(state);
           config.trunc_lvl = std::get<4>(state);
           config.tree_criterion = std::get<5>(state);
           config.threshold = std::get<6>(state);
           config.selection_criterion = std::get<7>(state);
           config.weights = std::get<8>(state);
           config.psi0 = std::get<9>(state);
           config.preselect_families = std::get<10>(state);
           config.select_trunc_lvl = std::get<11>(state);
           config.select_threshold = std::get<12>(state);
           config.select_families = std::get<13>(state);
           config.show_trace = std::get<14>(state);
           config.num_threads = std::get<15>(state);
           config.tree_algorithm = std::get<16>(state);
           config.allow_rotations = std::get<17>(state);
           config.seeds = std::get<18>(state);

           new (&controls) FitControlsVinecop(config);
         });
}