#pragma once

#include "docstr.hpp"
#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <stdexcept> // For std::invalid_argument
#include <vinecopulib.hpp>

namespace nb = nanobind;
using namespace nb::literals;
using namespace vinecopulib;

// Factory function to create a Vinecop from dimensionality
inline Vinecop
vc_from_dimension(const size_t d)
{
  return Vinecop(d);
}

inline Vinecop
vc_from_structure(
  std::optional<RVineStructure> structure = std::nullopt,
  std::optional<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>> matrix =
    std::nullopt,
  const std::vector<std::vector<Bicop>>& pair_copulas = {},
  const std::vector<std::string>& var_types = {})
{
  if (structure && matrix) {
    throw std::invalid_argument(
      "Only one of 'structure' or 'matrix' can be provided, not both.");
  } else if (structure) {
    // Use the structure-based constructor
    return Vinecop(*structure, pair_copulas, var_types);
  } else if (matrix) {
    // Use the matrix-based constructor
    return Vinecop(*matrix, pair_copulas, var_types);
  } else {
    throw std::invalid_argument(
      "Either 'structure' or 'matrix' must be provided.");
  }
}

inline Vinecop
vc_from_data(
  const Eigen::MatrixXd& data,
  std::optional<RVineStructure> structure = std::nullopt,
  std::optional<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>> matrix =
    std::nullopt,
  const std::vector<std::string>& var_types = {},
  const FitControlsVinecop& controls = FitControlsVinecop())
{
  if (structure && matrix) {
    throw std::invalid_argument(
      "Only one of 'structure' or 'matrix' can be provided, not both.");
  } else if (structure) {
    // Use the structure-based constructor
    return Vinecop(data, *structure, var_types, controls);
  } else if (matrix) {
    // Use the matrix-based constructor
    return Vinecop(data, *matrix, var_types, controls);
  } else {
    // Use the default constructor
    return Vinecop(data, RVineStructure(), var_types, controls);
  }
}

// Factory function to create a Vinecop from a file
inline Vinecop
vc_from_file(const std::string& filename, bool check = true)
{
  return Vinecop(filename, check);
}

inline void
init_vinecop_class(nb::module_& module)
{

  constexpr auto& vinecop_doc = pyvinecopulib_doc.vinecopulib.Vinecop;

  const char* default_constructor_doc =
    R"""(Default constructor for the ``Vinecop`` class.

The default constructor uses ``Vinecop.from_dimension()`` to instantiate an
empty vine copula of a given dimension. It can then be used to select a model from data using ``Vinecop.select()``. Alternatives to instantiate vine copulas
are:

- ``Vinecop.from_structure()``: Instantiate a vine copula from a structure or a matrix, as well as optional pair-copulas and variable types.
- ``Vinecop.from_data()``: Instantiate a vine copula from data, as well as optional structure or matrix, pair-copulas, and variable types.

)""";

  const char* from_data_doc = R"""(
  Factory function to create a Vinecop from data.

  Parameters
  ----------
  data :
      Input data matrix.

  structure :
      RVine structure. Provide either this or `matrix`, but not both.

  matrix :
      RVine matrix. Provide either this or `structure`, but not both.

  var_types :
      Variable types for each variable (e.g., 'c' for continuous, 'd' for discrete). Defaults to all continuous.

  controls :
      Fit controls for the vinecop. Defaults to the default constructor.
  )""";

  const char* from_structure_doc = R"""(
  Factory function to create a Vinecop using either a structure or a matrix.

  Parameters
  ----------
  structure :
      Vinecop structure. Provide either this or `matrix`, but not both.

  matrix : Eigen::Matrix, optional
      Vinecop matrix. Provide either this or `structure`, but not both.

  pair_copulas : list of list of Bicop, optional
      Pairwise copulas for each edge in the vine. Defaults to an empty list.

  var_types :
      Variable types for each variable (e.g., 'c' for continuous, 'd' for discrete). Defaults to all continuous.
  )""";

  nb::class_<Vinecop>(module, "Vinecop", vinecop_doc.doc)
    .def(nb::init<const size_t>(), default_constructor_doc, "d"_a)
    .def_static(
      "from_dimension", &vc_from_dimension, "d"_a, vinecop_doc.ctor.doc_1args_d)
    .def_static("from_structure",
                &vc_from_structure,
                "structure"_a = std::nullopt,
                "matrix"_a = std::nullopt,
                "pair_copulas"_a = std::vector<std::vector<Bicop>>(),
                "var_types"_a = std::vector<std::string>(),
                from_structure_doc)
    .def_static("from_data",
                &vc_from_data,
                "data"_a,
                "structure"_a = std::nullopt,
                "matrix"_a = std::nullopt,
                "var_types"_a = std::vector<std::string>(),
                "controls"_a.sig("FitControlsVinecop()") = FitControlsVinecop(),
                from_data_doc)
    .def_static("from_file",
                &vc_from_file,
                "filename"_a,
                "check"_a = true,
                vinecop_doc.ctor.doc_2args_filename_check)
    .def("to_file", &Vinecop::to_file, "filename"_a, vinecop_doc.to_file.doc)
    .def("to_json", &Vinecop::to_json, vinecop_doc.to_json.doc)
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
    .def(
      "get_parameters",
      [](const Vinecop& self, size_t tree, size_t edge) {
        return nb::cast(self.get_parameters(tree, edge));
      },
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
      "structure", &Vinecop::get_rvine_structure, "The R-vine structure.")
    .def_prop_ro(
      "npars", &Vinecop::get_npars, "The total number of parameters.")
    .def_prop_ro(
      "matrix",
      [](const Vinecop& self) { return nb::cast(self.get_matrix()); },
      "Extracts the R-vine structure's matrix.")
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

        // Call the Python function with the C++ object and additional
        // arguments
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