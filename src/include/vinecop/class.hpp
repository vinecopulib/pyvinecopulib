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
init_vinecop_class(nb::module_& module)
{

  constexpr auto& vinecop_doc = pyvinecopulib_doc.vinecopulib.Vinecop;

  nb::class_<Vinecop>(module, "Vinecop", vinecop_doc.doc)
    .def(nb::init<const size_t>(), vinecop_doc.ctor.doc_1args_d, "d"_a)
    .def(nb::init<const RVineStructure&,
                  const std::vector<std::vector<Bicop>>&,
                  const std::vector<std::string>&>(),
         "structure"_a,
         "pair_copulas"_a = std::vector<size_t>(),
         "var_types"_a = std::vector<std::string>(),
         vinecop_doc.ctor.doc_2args_structure_constint)
    .def(nb::init<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>&,
                  const std::vector<std::vector<Bicop>>&,
                  const std::vector<std::string>&>(),
         "matrix"_a,
         "pair_copulas"_a = std::vector<size_t>(),
         "var_types"_a = std::vector<std::string>(),
         vinecop_doc.ctor.doc_2args_matrix_constint)
    .def(nb::init<const Eigen::MatrixXd&,
                  const RVineStructure&,
                  const std::vector<std::string>&,
                  const FitControlsVinecop&>(),
         "data"_a,
         "structure"_a = RVineStructure(),
         "var_types"_a = std::vector<std::string>(),
         "controls"_a.sig("FitControlsVinecop()") = FitControlsVinecop(),
         vinecop_doc.ctor.doc_4args_data_structure_var_types_controls)
    .def(nb::init<const Eigen::MatrixXd&,
                  const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>&,
                  const std::vector<std::string>&,
                  const FitControlsVinecop&>(),
         "data"_a,
         "matrix"_a = Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>(),
         "var_types"_a = std::vector<std::string>(),
         "controls"_a.sig("FitControlsVinecop()") = FitControlsVinecop(),
         vinecop_doc.ctor.doc_4args_data_matrix_var_types_controls)
    .def(nb::init<const std::string, bool>(),
         "filename"_a,
         "check"_a = true,
         vinecop_doc.ctor.doc_2args_filename_check)
    .def("to_json", &Vinecop::to_file, "filename"_a, vinecop_doc.to_file.doc)
    .def_prop_rw("var_types",
                 &Vinecop::get_var_types,
                 &Vinecop::set_var_types,
                 "The types of each variables.")
    .def_prop_ro("trunc_lvl", &Vinecop::get_trunc_lvl, "The truncation level.")
    .def_prop_ro("dim", &Vinecop::get_dim, "The dimension.")
    .def("get_pair_copula",
         &Vinecop::get_pair_copula,
         "Gets a pair-copula.",
         "tree"_a,
         "edge"_a)
    .def("get_family",
         &Vinecop::get_family,
         "Gets the family of a pair-copula.",
         "tree"_a,
         "edge"_a)
    .def("get_rotation",
         &Vinecop::get_rotation,
         "Gets the rotation of a pair-copula.",
         "tree"_a,
         "edge"_a)
    .def("get_parameters",
         &Vinecop::get_parameters,
         "Gets the parameters of a pair-copula.",
         "tree"_a,
         "edge"_a)
    .def("get_tau",
         &Vinecop::get_tau,
         "Gets the kendall's tau of a pair-copula.",
         "tree"_a,
         "edge"_a)
    .def_prop_ro(
      "pair_copulas", &Vinecop::get_all_pair_copulas, "All pair-copulas.")
    .def_prop_ro(
      "families", &Vinecop::get_all_families, "Families of all pair-copulas.")
    .def_prop_ro("rotations",
                 &Vinecop::get_all_rotations,
                 "The rotations of all pair-copulas.")
    .def_prop_ro("parameters",
                 &Vinecop::get_all_parameters,
                 "The parameters of all pair-copulas.")
    .def_prop_ro(
      "taus", &Vinecop::get_all_taus, "The Kendall's taus of all pair-copulas.")
    .def_prop_ro("order", &Vinecop::get_order, "The R-vine structure's order.")
    .def_prop_ro(
      "matrix", &Vinecop::get_matrix, "The R-vine structure's matrix.")
    .def_prop_ro(
      "structure", &Vinecop::get_rvine_structure, "The R-vine structure.")
    .def_prop_ro(
      "npars", &Vinecop::get_npars, "The total number of parameters.")
    .def_prop_ro("nobs",
                 &Vinecop::get_nobs,
                 "The number of observations (for fitted objects only).")
    .def_prop_ro("threshold",
                 &Vinecop::get_threshold,
                 "The threshold (for thresholded copulas only).")
    .def("select",
         &Vinecop::select,
         "data"_a,
         "controls"_a.sig("FitControlsVinecop()") = FitControlsVinecop(),
         vinecop_doc.select.doc)
    .def("fit",
         &Vinecop::fit,
         "data"_a,
         "controls"_a.sig("FitControlsBicop()") = FitControlsBicop(),
         "num_threads"_a = 1,
         vinecop_doc.fit.doc)
    .def("pdf", &Vinecop::pdf, "u"_a, "num_threads"_a = 1, vinecop_doc.pdf.doc)
    .def("cdf",
         &Vinecop::cdf,
         "u"_a,
         "N"_a = 10000,
         "num_threads"_a = 1,
         "seeds"_a = std::vector<int>(),
         vinecop_doc.cdf.doc)
    .def("simulate",
         &Vinecop::simulate,
         "n"_a,
         "qrng"_a = false,
         "num_threads"_a = 1,
         "seeds"_a = std::vector<int>(),
         vinecop_doc.simulate.doc)
    .def("rosenblatt",
         &Vinecop::rosenblatt,
         "u"_a,
         "num_threads"_a = 1,
         "randomize_discrete"_a = true,
         "seeds"_a = std::vector<int>(),
         vinecop_doc.rosenblatt.doc)
    .def("inverse_rosenblatt",
         &Vinecop::inverse_rosenblatt,
         "u"_a,
         "num_threads"_a = 1,
         vinecop_doc.inverse_rosenblatt.doc)
    .def("loglik",
         &Vinecop::loglik,
         "u"_a = Eigen::MatrixXd(),
         "num_threads"_a = 1,
         vinecop_doc.loglik.doc)
    .def("aic",
         &Vinecop::aic,
         "u"_a = Eigen::MatrixXd(),
         "num_threads"_a = 1,
         vinecop_doc.aic.doc)
    .def("bic",
         &Vinecop::bic,
         "u"_a = Eigen::MatrixXd(),
         "num_threads"_a = 1,
         vinecop_doc.bic.doc)
    .def("mbicv",
         &Vinecop::mbicv,
         "u"_a = Eigen::MatrixXd(),
         "psi0"_a = 0.9,
         "num_threads"_a = 1,
         vinecop_doc.mbicv.doc)
    .def(
      "__repr__",
      [](const Vinecop& cop) { return "<pyvinecopulib.Vinecop> " + cop.str(); })
    .def(
      "str",
      [](const Vinecop& cop, const std::vector<size_t>& trees = {}) {
        return "<pyvinecopulib.Vinecop> " + cop.str(trees);
      },
      "trees"_a = std::vector<size_t>{},
      vinecop_doc.str.doc)
    .def(
      "truncate", &Vinecop::truncate, "trunc_lvl"_a, vinecop_doc.truncate.doc)
    .def(
      "plot",
      [](const Vinecop& cop,
         nb::object tree = nb::none(),
         bool add_edge_labels = true,
         const std::string& layout = "graphviz",
         nb::object vars_names = nb::none()) {
        auto python_helpers_plotting =
          nb::module_::import_("pyvinecopulib._python_helpers.vinecop");

        // Import the Python plotting function
        nb::object vinecop_plot = python_helpers_plotting.attr("vinecop_plot");

        // Call the Python function with the C++ object and additional arguments
        vinecop_plot(nb::cast(cop), tree, add_edge_labels, layout, vars_names);
      },
      "tree"_a = nb::none(),
      "add_edge_labels"_a = true,
      "layout"_a = "graphviz",
      "vars_names"_a = nb::none(),
      nb::cast<std::string>(
        nb::module_::import_("pyvinecopulib._python_helpers.vinecop")
          .attr("VINECOP_PLOT_DOC"))
        .c_str());
}