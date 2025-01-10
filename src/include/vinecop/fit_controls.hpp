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
                  std::string>(),
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
         "mst_algorithm"_a = "prim",
         fitcontrolsvinecop_doc.ctor.doc_17args)
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
    .def_prop_rw("mst_algorithm",
                 &FitControlsVinecop::get_mst_algorithm,
                 &FitControlsVinecop::set_mst_algorithm,
                 "The minimum spanning tree algorithm.")
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
                                  controls.get_mst_algorithm());
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
                       std::string> state) {
           new (&controls) FitControlsVinecop(std::get<0>(state),
                                              std::get<1>(state),
                                              std::get<2>(state),
                                              std::get<3>(state),
                                              std::get<4>(state),
                                              std::get<5>(state),
                                              std::get<6>(state),
                                              std::get<7>(state),
                                              std::get<8>(state),
                                              std::get<9>(state),
                                              std::get<10>(state),
                                              std::get<11>(state),
                                              std::get<12>(state),
                                              std::get<13>(state),
                                              std::get<14>(state),
                                              std::get<15>(state),
                                              std::get<16>(state));
         });
}