#pragma once

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <algorithm>  // For std::min
#include <stdexcept>  // For std::invalid_argument
#include <vinecopulib.hpp>

#include "docstr.hpp"
#include "misc/helpers.hpp"

namespace nb = nanobind;
using namespace nb::literals;
using namespace vinecopulib;

// Convert a vinecopulib ``TriangularArray<T>`` to a nested Python list indexed
// ``[tree][edge]`` via ``TriangularArray::to_list()``. An ``Eigen::VectorXd``
// element becomes a 1-d ndarray; a ``std::vector<Eigen::MatrixXd>`` element
// becomes a list of 2-d ndarrays. Requires the GIL to be held (it builds
// Python objects).
template <typename T>
inline nb::list triangular_to_list(const TriangularArray<T>& array) {
  return nb::steal<nb::list>(nb::cast(array.to_list()).release());
}

inline void vinecop_plot_wrapper(const Vinecop& cop, nb::object tree,
                                 bool add_edge_labels,
                                 const std::string& layout,
                                 nb::object vars_names) {
  // Import the vinecop helper Python module
  auto mod = nb::module_::import_("pyvinecopulib._python_helpers.vinecop");

  // Import the Python plotting function
  auto vinecop_plot = mod.attr("vinecop_plot");

  // Call the Python function with the C++ object and additional
  // arguments
  vinecop_plot(nb::cast(cop), tree, add_edge_labels, layout, vars_names);
}

// Factory function to create a Vinecop from dimensionality
inline Vinecop vc_from_dimension(const size_t d) { return Vinecop(d); }

inline Vinecop vc_from_structure(
    std::optional<RVineStructure> structure = std::nullopt,
    std::optional<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>>
        matrix = std::nullopt,
    const std::vector<std::vector<Bicop>>& pair_copulas = {},
    const std::vector<std::string>& var_types = {}) {
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

inline Vinecop vc_from_data(
    const Eigen::MatrixXd& data,
    std::optional<RVineStructure> structure = std::nullopt,
    std::optional<Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>>
        matrix = std::nullopt,
    const std::vector<std::string>& var_types = {},
    const FitControlsVinecop& controls = FitControlsVinecop()) {
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
inline Vinecop vc_from_file(const std::string& filename, bool check = true) {
  return Vinecop(filename, check);
}

// Factory function to create a Vinecop from a JSON string
inline Vinecop vc_from_json(const std::string& json, bool check = true) {
  nlohmann::json json_obj = nlohmann::json::parse(json);
  return Vinecop(json_obj, check);
}

inline void init_vinecop_class(nb::module_& module) {
  constexpr auto& vinecop_doc = pyvinecopulib_doc.vinecopulib.Vinecop;

  const char* default_constructor_doc =
      R"""(Default constructor for the ``Vinecop`` class.

The default constructor uses ``Vinecop.from_dimension()`` to instantiate an
empty vine copula of a given dimension. It can then be used to select a model from data using ``Vinecop.select()``. Alternatives to instantiate vine copulas
are:

- ``Vinecop.from_data()``: Instantiate from data, as well as an optional ``FitControlsVinecop``, an ``RVineStructure`` or matrix, and variable types.
- ``Vinecop.from_structure()``: Instantiate from an ``RVineStructure`` or matrix, as well as optional pair-copulas and variable types.
- ``Vinecop.from_file()``: Instantiate from a file.
- ``Vinecop.from_json()``: Instantiate from a JSON string.
)""";

  const char* from_data_doc = R"""(
  Factory function to create a Vinecop from data.

  Parameters
  ----------
  data :
      Input data matrix.

  structure :
      An ``RVineStructure``. Provide either this or `matrix`, but not both.

  matrix :
      RVine matrix. Provide either this or `structure`, but not both.

  var_types :
      Variable types for each variable (e.g., 'c' for continuous, 'd' for discrete). Defaults to all continuous.

  controls :
      Fit controls for the vinecop. Defaults to the default constructor.
  )""";

  // Supplied inline rather than via the generated docstring: the C++ facade
  // returns a ``PdfWithHfuncsResult`` struct (whose inline definition trips up
  // the libclang docstring extractor), whereas the Python method returns a
  // dict. This documents the actual Python return.
  const char* pdf_full_doc =
      R"""(Evaluates the vine copula density and, optionally, the per-edge densities and h-functions.

Parameters
----------
u : ndarray, shape (n, d) or (n, 2 * d), dtype float
    Evaluation points in (0, 1); additional columns are required for discrete
    variables (see ``Vinecop.select``).
num_threads : int, default 1
    Number of threads to parallelize the row-wise evaluation over.
keep_all : bool, default True
    If True, also return the per-edge densities and h-functions.

Returns
-------
dict
    A dict with key ``"pdf"`` (an ndarray of length n, the copula density). If
    ``keep_all`` is True, it also contains ``"pdf_edges"``, ``"hfunc1"``,
    ``"hfunc2"``, ``"hfunc1_sub"`` and ``"hfunc2_sub"``, each a nested list
    indexed ``[tree][edge]`` of length-n arrays. The ``_sub`` fields hold the
    left-limit h-functions and are only populated when at least one variable is
    discrete.
)""";

  const char* from_structure_doc = R"""(
  Factory function to create a Vinecop using either a structure or a matrix.

  Parameters
  ----------
  structure :
      An ``RVineStructure``. Provide either this or `matrix`, but not both.

  matrix :
      Vinecop matrix. Provide either this or `structure`, but not both.

  pair_copulas :
      Pairwise copulas for each edge in the vine. Defaults to an empty list.

  var_types :
      Variable types for each variable (e.g., 'c' for continuous, 'd' for discrete). Defaults to all continuous.
  )""";

  nb::class_<Vinecop>(module, "Vinecop", vinecop_doc.doc)
      .def(nb::init<const size_t>(), default_constructor_doc, "d"_a,
           nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_dimension", &vc_from_dimension, "d"_a,
                  vinecop_doc.ctor.doc_1args_d,
                  nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_structure", &vc_from_structure,
                  "structure"_a = std::nullopt, "matrix"_a = std::nullopt,
                  "pair_copulas"_a = std::vector<std::vector<Bicop>>(),
                  "var_types"_a = std::vector<std::string>(),
                  from_structure_doc, nb::call_guard<nb::gil_scoped_release>())
      .def_static(
          "from_data", &vc_from_data, "data"_a, "structure"_a = std::nullopt,
          "matrix"_a = std::nullopt, "var_types"_a = std::vector<std::string>(),
          "controls"_a.sig("FitControlsVinecop()") = FitControlsVinecop(),
          from_data_doc, nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_file", &vc_from_file, "filename"_a, "check"_a = true,
                  vinecop_doc.ctor.doc_2args_filename_check,
                  nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_json", &vc_from_json, "json"_a, "check"_a = true,
                  vinecop_doc.ctor.doc_2args_input_check,
                  nb::call_guard<nb::gil_scoped_release>())
      .def("to_file", &Vinecop::to_file, "filename"_a, vinecop_doc.to_file.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def(
          "to_json",
          [](Vinecop& self) -> std::string { return self.to_json().dump(); },
          vinecop_doc.to_json.doc, nb::call_guard<nb::gil_scoped_release>())
      .def_prop_rw("var_types", &Vinecop::get_var_types,
                   &Vinecop::set_var_types, vinecop_doc.get_var_types.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      .def_prop_ro("trunc_lvl", &Vinecop::get_trunc_lvl,
                   vinecop_doc.get_trunc_lvl.doc)
      .def_prop_ro("dim", &Vinecop::get_dim, vinecop_doc.get_dim.doc)
      // The five `(tree, edge)` getters below take parameters libclang
      // can't disambiguate against the `get_all_*` bulk variants, so the
      // auto-extracted constant falls back to
      // `doc_was_unable_to_choose_unambiguous_names`. Keep short
      // hand-written literals for these — the upstream `//!` comments
      // are richer but only surface in the Doxygen-built C++ docs.
      .def("get_pair_copula", &Vinecop::get_pair_copula,
           "Gets the pair copula at the given (tree, edge) position.", "tree"_a,
           "edge"_a, nb::call_guard<nb::gil_scoped_release>())
      .def("get_family", &Vinecop::get_family,
           "Gets the family of the pair copula at the given (tree, edge).",
           "tree"_a, "edge"_a)
      .def("get_rotation", &Vinecop::get_rotation,
           "Gets the rotation (degrees) of the pair copula at the given "
           "(tree, edge).",
           "tree"_a, "edge"_a)
      .def("get_parameters", &Vinecop::get_parameters,
           "Gets the parameter matrix of the pair copula at the given "
           "(tree, edge).",
           "tree"_a, "edge"_a)
      .def("get_tau", &Vinecop::get_tau,
           "Gets Kendall's tau of the pair copula at the given (tree, edge).",
           "tree"_a, "edge"_a)
      .def_prop_ro("pair_copulas", &Vinecop::get_all_pair_copulas,
                   vinecop_doc.get_all_pair_copulas.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      .def_prop_ro("families", &Vinecop::get_all_families,
                   vinecop_doc.get_all_families.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      .def_prop_ro("rotations", &Vinecop::get_all_rotations,
                   vinecop_doc.get_all_rotations.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      .def_prop_ro("parameters", &Vinecop::get_all_parameters,
                   vinecop_doc.get_all_parameters.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      .def_prop_ro("taus", &Vinecop::get_all_taus, vinecop_doc.get_all_taus.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      .def_prop_ro("order", &Vinecop::get_order, vinecop_doc.get_order.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      .def_prop_ro("structure", &Vinecop::get_rvine_structure,
                   vinecop_doc.get_rvine_structure.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      .def_prop_ro("npars", &Vinecop::get_npars, vinecop_doc.get_npars.doc)
      .def_prop_ro("matrix", &Vinecop::get_matrix, vinecop_doc.get_matrix.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      .def(
          "get_struct_array",
          [](const Vinecop& cop, bool natural_order) -> nb::list {
            return triangular_to_list(cop.get_struct_array(natural_order));
          },
          "natural_order"_a = false,
          struct_array_list_doc(vinecop_doc.get_struct_array.doc).c_str())
      .def_prop_ro("nobs", &Vinecop::get_nobs, vinecop_doc.get_nobs.doc)
      .def_prop_ro("threshold", &Vinecop::get_threshold,
                   vinecop_doc.get_threshold.doc)
      .def("select", &Vinecop::select, "data"_a,
           "controls"_a.sig("FitControlsVinecop()") = FitControlsVinecop(),
           vinecop_doc.select.doc, nb::call_guard<nb::gil_scoped_release>())
      .def("fit", &Vinecop::fit, "data"_a,
           "controls"_a.sig("FitControlsBicop()") = FitControlsBicop(),
           "num_threads"_a = 1, vinecop_doc.fit.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def("pdf", &Vinecop::pdf, "u"_a, "num_threads"_a = 1,
           vinecop_doc.pdf.doc, nb::call_guard<nb::gil_scoped_release>())
      .def(
          "pdf_full",
          [](const Vinecop& cop, Eigen::MatrixXd u, size_t num_threads,
             bool keep_all) -> nb::dict {
            Vinecop::PdfWithHfuncsResult res;
            {
              // Release the GIL only around the C++ computation; the dict /
              // list construction below must hold it.
              nb::gil_scoped_release release;
              res = cop.pdf_full(std::move(u), num_threads, keep_all);
            }
            nb::dict out;
            out["pdf"] = nb::cast(res.pdf);
            if (keep_all) {
              out["pdf_edges"] = triangular_to_list(res.pdf_edges);
              out["hfunc1"] = triangular_to_list(res.hfunc1);
              out["hfunc2"] = triangular_to_list(res.hfunc2);
              out["hfunc1_sub"] = triangular_to_list(res.hfunc1_sub);
              out["hfunc2_sub"] = triangular_to_list(res.hfunc2_sub);
            }
            return out;
          },
          "u"_a, "num_threads"_a = 1, "keep_all"_a = true, pdf_full_doc)
      .def("cdf", &Vinecop::cdf, "u"_a, "N"_a = 10000, "num_threads"_a = 1,
           "seeds"_a = std::vector<int>(), vinecop_doc.cdf.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def("simulate", &Vinecop::simulate, "n"_a, "qrng"_a = false,
           "num_threads"_a = 1, "seeds"_a = std::vector<int>(),
           vinecop_doc.simulate.doc, nb::call_guard<nb::gil_scoped_release>())
      .def("rosenblatt", &Vinecop::rosenblatt, "u"_a, "num_threads"_a = 1,
           "randomize_discrete"_a = true, "seeds"_a = std::vector<int>(),
           vinecop_doc.rosenblatt.doc, nb::call_guard<nb::gil_scoped_release>())
      .def("inverse_rosenblatt", &Vinecop::inverse_rosenblatt, "u"_a,
           "num_threads"_a = 1, vinecop_doc.inverse_rosenblatt.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def("loglik", &Vinecop::loglik, "u"_a = Eigen::MatrixXd(),
           "num_threads"_a = 1, vinecop_doc.loglik.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def("aic", &Vinecop::aic, "u"_a = Eigen::MatrixXd(), "num_threads"_a = 1,
           vinecop_doc.aic.doc, nb::call_guard<nb::gil_scoped_release>())
      .def("bic", &Vinecop::bic, "u"_a = Eigen::MatrixXd(), "num_threads"_a = 1,
           vinecop_doc.bic.doc, nb::call_guard<nb::gil_scoped_release>())
      .def("mbicv", &Vinecop::mbicv, "u"_a = Eigen::MatrixXd(), "psi0"_a = 0.9,
           "num_threads"_a = 1, vinecop_doc.mbicv.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def("scores", &Vinecop::scores, "u"_a, "step_wise"_a = true,
           "num_threads"_a = 1, vinecop_doc.scores.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def("hessian_avg", &Vinecop::hessian_avg, "u"_a, "step_wise"_a = true,
           "num_threads"_a = 1, vinecop_doc.hessian_avg.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def("scores_cov", &Vinecop::scores_cov, "u"_a, "step_wise"_a = true,
           "num_threads"_a = 1, vinecop_doc.scores_cov.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def(
          "hessian",
          [](Vinecop& cop, Eigen::MatrixXd u, bool step_wise,
             size_t num_threads) -> nb::list {
            TriangularArray<std::vector<Eigen::MatrixXd>> hess;
            {
              // Release the GIL only around the C++ computation; the nested
              // list construction below must hold it.
              nb::gil_scoped_release release;
              hess = cop.hessian(std::move(u), step_wise, num_threads);
            }
            return triangular_to_list(hess);
          },
          "u"_a, "step_wise"_a = true, "num_threads"_a = 1,
          vinecop_doc.hessian.doc)
      .def(
          "__repr__",
          [](const Vinecop& cop) {
            return "<pyvinecopulib.core.Vinecop> " + cop.str();
          },
          vinecop_doc.str.doc)
      .def(
          "__str__",
          [](const Vinecop& cop) {
            return "<pyvinecopulib.core.Vinecop> " + cop.str();
          },
          vinecop_doc.str.doc)
      .def(
          "format",
          [](const Vinecop& cop, const std::vector<size_t>& trees = {}) {
            return cop.str(trees);
          },
          "trees"_a = std::vector<size_t>{}, vinecop_doc.str.doc)
      .def("truncate", &Vinecop::truncate, "trunc_lvl"_a,
           vinecop_doc.truncate.doc, nb::call_guard<nb::gil_scoped_release>())
      .def("plot", &vinecop_plot_wrapper, "tree"_a = nb::none(),
           "add_edge_labels"_a = true, "layout"_a = "graphviz",
           "vars_names"_a = nb::none(),
           python_doc_helper("pyvinecopulib._python_helpers.vinecop",
                             "VINECOP_PLOT_DOC",
                             "Plot the vine copula (extended doc unavailable)")
               .c_str())
      .def("__getstate__",
           [](const Vinecop& cop) {
             nb::dict state;
             state["rvine_structure"] =
                 cop.get_rvine_structure().to_json().dump();
             state["pair_copulas"] = cop.get_all_pair_copulas();
             state["var_types"] = cop.get_var_types();
             return state;
           })

      .def("__setstate__", [](Vinecop& cop, nb::dict state) {
        nlohmann::json json_obj = nlohmann::json::parse(
            nb::cast<std::string>(state["rvine_structure"]));
        RVineStructure structure(json_obj);

        new (&cop) Vinecop(
            std::move(structure),
            nb::cast<std::vector<std::vector<Bicop>>>(state["pair_copulas"]),
            nb::cast<std::vector<std::string>>(state["var_types"]));
      });
}
