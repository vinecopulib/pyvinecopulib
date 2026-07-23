#pragma once

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <tuple>
#include <vinecopulib.hpp>
#include <vinecopulib/vinecop/rvine_trees.hpp>

#include "docstr.hpp"
#include "misc/helpers.hpp"

namespace nb = nanobind;
using namespace nb::literals;
using namespace vinecopulib;

// Factory function to create RVineStructure with dimension and truncation level
inline RVineStructure rv_from_dimension(
    size_t d = static_cast<size_t>(1),
    size_t trunc_lvl = std::numeric_limits<size_t>::max()) {
  return RVineStructure(d, trunc_lvl);
}

// Factory function to create RVineStructure from a matrix
inline RVineStructure rv_from_matrix(
    const Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic>& mat,
    bool check = true) {
  return RVineStructure(mat, check);
}

// Factory function to create RVineStructure from an order and truncation level
inline RVineStructure rv_from_order(
    const std::vector<size_t>& order,
    size_t trunc_lvl = std::numeric_limits<size_t>::max(), bool check = true) {
  return RVineStructure(order, trunc_lvl, check);
}

// Factory function to create RVineStructure from an order vector and a
// structure array given as nested rows (tree i holds d - 1 - i entries)
inline RVineStructure rv_from_struct_array(
    const std::vector<size_t>& order,
    const std::vector<std::vector<size_t>>& struct_array,
    bool natural_order = false, bool check = true) {
  return RVineStructure(order, TriangularArray<size_t>(struct_array),
                        natural_order, check);
}

// Factory function to assemble an RVineStructure from a list-of-trees
// decomposition (the Dissmann finalization, reusing the upstream RVineTrees
// primitive). ``trees`` is a nested list ``[tree][edge]`` of ``(a, b,
// conditioning)`` triples with 1-based variable labels: ``a`` / ``b`` are the
// conditioned pair and ``conditioning`` the (possibly empty) conditioning set;
// tree ``t`` must hold exactly ``d - 1 - t`` edges. The upstream peel with the
// selection diagonal policy (``leaf_edges[0].back()``) chooses the diagonals,
// so the resulting matrix matches what ``Vinecop.select`` produces for the same
// trees. Pair-copulas are irrelevant here (independence placeholders); callers
// re-fit on the returned structure.
inline RVineStructure rv_from_trees(
    size_t d,
    const std::vector<
        std::vector<std::tuple<size_t, size_t, std::vector<size_t>>>>& trees) {
  if (trees.empty()) {
    return RVineStructure(d, static_cast<size_t>(0));
  }
  std::vector<RVineTrees::Tree> tree_list;
  tree_list.reserve(trees.size());
  for (const auto& tree : trees) {
    RVineTrees::Tree edges;
    edges.reserve(tree.size());
    for (const auto& e : tree) {
      edges.emplace_back(std::get<0>(e), std::get<1>(e), std::get<2>(e));
    }
    tree_list.push_back(std::move(edges));
  }
  auto dec = RVineTrees(d, std::move(tree_list))
                 .to_struct_array(
                     [](size_t, const std::vector<std::vector<size_t>>& leaf) {
                       return leaf[0].back();
                     });
  return RVineStructure(dec.order, dec.struct_array);
}

// Factory function to create RVineStructure from a file
inline RVineStructure rv_from_file(const std::string& filename,
                                   bool check = true) {
  return RVineStructure(filename, check);
}

// Factory function to create RVineStructure from a JSON string
inline RVineStructure rv_from_json(const std::string& json, bool check = true) {
  nlohmann::json json_obj = nlohmann::json::parse(json);
  return RVineStructure(json_obj, check);
}

inline void init_vinecop_rvine_structure(nb::module_& module) {
  constexpr auto& doc = pyvinecopulib_doc;
  constexpr auto& rvinestructure_doc = doc.vinecopulib.RVineStructure;
  constexpr auto& dvinestructure_doc = doc.vinecopulib.DVineStructure;
  constexpr auto& cvinestructure_doc = doc.vinecopulib.CVineStructure;

  const char* default_constructor_doc =
      R"""(Default constructor for the ``RVineStructure`` class.

The default constructor uses ``RVineStructure.from_dimension()`` to instantiate
a default structure of a given dimension and truncation level.
Alternatives to instantiate structures are:

- ``RVineStructure.from_order()``: Instantiate from an order vector.
- ``RVineStructure.from_matrix()``: Instantiate from a matrix.
- ``RVineStructure.from_struct_array()``: Instantiate from an order vector and a structure array.
- ``RVineStructure.from_file()``: Instantiate from a JSON file, or a CBOR file when the filename ends in ``.cbor``.
- ``RVineStructure.from_json()``: Instantiate from a JSON string.
)""";

  const char* from_trees_doc =
      R"""(Assemble an ``RVineStructure`` from a list-of-trees decomposition.

This is the finalization step of Dissmann-style structure selection: given the
edges chosen for each tree, it produces the corresponding R-vine matrix
(variable order + structure array), reusing the same diagonal policy as
``Vinecop.select`` so the result matches the compiled selector.

Parameters
----------
d : int
    Dimension of the vine.
trees : list of list of tuple
    Per-tree edge lists, indexed ``[tree][edge]``. Each edge is a triple
    ``(a, b, conditioning)`` of 1-based variable labels, where ``a`` and ``b``
    are the conditioned pair and ``conditioning`` is the (possibly empty)
    conditioning set. Tree ``t`` must contain exactly ``d - 1 - t`` edges. An
    empty list yields a fully truncated (independence) structure.

Returns
-------
RVineStructure
    The assembled structure. Pair-copulas are not part of the input; callers
    fit them on the returned structure.
)""";

  nb::class_<RVineStructure>(module, "RVineStructure", rvinestructure_doc.doc)
      .def(nb::init<const size_t&, const size_t&>(),
           "d"_a = static_cast<size_t>(1),
           "trunc_lvl"_a = std::numeric_limits<size_t>::max(),
           default_constructor_doc, nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_dimension", &rv_from_dimension,
                  "d"_a = static_cast<size_t>(1),
                  "trunc_lvl"_a = std::numeric_limits<size_t>::max(),
                  rvinestructure_doc.ctor.doc_2args_d_trunc_lvl,
                  nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_matrix", &rv_from_matrix, "mat"_a, "check"_a = true,
                  rvinestructure_doc.ctor.doc_2args_mat_check,
                  nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_order", &rv_from_order, "order"_a,
                  "trunc_lvl"_a = std::numeric_limits<size_t>::max(),
                  "check"_a = true,
                  rvinestructure_doc.ctor.doc_3args_order_trunc_lvl_check,
                  nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_struct_array", &rv_from_struct_array, "order"_a,
                  "struct_array"_a, "natural_order"_a = false, "check"_a = true,
                  struct_array_list_doc(
                      rvinestructure_doc.ctor
                          .doc_4args_order_struct_array_natural_order_check)
                      .c_str(),
                  nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_trees", &rv_from_trees, "d"_a, "trees"_a,
                  from_trees_doc, nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_file", &rv_from_file, "filename"_a, "check"_a = true,
                  rvinestructure_doc.ctor.doc_2args_filename_check,
                  nb::call_guard<nb::gil_scoped_release>())
      .def_static("from_json", &rv_from_json, "json"_a, "check"_a = true,
                  rvinestructure_doc.ctor.doc_2args_input_check,
                  nb::call_guard<nb::gil_scoped_release>())
      .def("to_file", &RVineStructure::to_file, "filename"_a,
           rvinestructure_doc.to_file.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def(
          "to_json",
          [](const RVineStructure& self) { return self.to_json().dump(); },
          rvinestructure_doc.to_json.doc,
          nb::call_guard<nb::gil_scoped_release>())
      .def_prop_ro("dim", &RVineStructure::get_dim,
                   rvinestructure_doc.get_dim.doc)
      .def_prop_ro("trunc_lvl", &RVineStructure::get_trunc_lvl,
                   rvinestructure_doc.get_trunc_lvl.doc)
      .def_prop_ro("order",
                   (std::vector<size_t>(RVineStructure::*)() const) &
                       RVineStructure::get_order,
                   rvinestructure_doc.get_order.doc_0args,
                   nb::call_guard<nb::gil_scoped_release>())
      .def_prop_ro("matrix", &RVineStructure::get_matrix,
                   rvinestructure_doc.get_matrix.doc,
                   nb::call_guard<nb::gil_scoped_release>())
      .def("struct_array", &RVineStructure::struct_array, "tree"_a, "edge"_a,
           "natural_order"_a = false, rvinestructure_doc.struct_array.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def(
          "get_struct_array",
          [](const RVineStructure& rvs, bool natural_order) -> nb::list {
            return triangular_to_list(rvs.get_struct_array(natural_order));
          },
          "natural_order"_a = false,
          struct_array_list_doc(rvinestructure_doc.get_struct_array.doc)
              .c_str())
      .def("min_array", &RVineStructure::min_array, "tree"_a, "edge"_a,
           rvinestructure_doc.min_array.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def("needed_hfunc1", &RVineStructure::needed_hfunc1, "tree"_a, "edge"_a,
           rvinestructure_doc.needed_hfunc1.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def("needed_hfunc2", &RVineStructure::needed_hfunc2, "tree"_a, "edge"_a,
           rvinestructure_doc.needed_hfunc2.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def("truncate", &RVineStructure::truncate, "trunc_lvl"_a,
           rvinestructure_doc.truncate.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def_static("simulate", &RVineStructure::simulate, "d"_a,
                  "natural_order"_a = false, "seeds"_a = std::vector<size_t>(),
                  rvinestructure_doc.simulate.doc,
                  nb::call_guard<nb::gil_scoped_release>())
      .def(
          "__repr__",
          [](const RVineStructure& rvs) {
            return "<pyvinecopulib.core.RVineStructure>\n" + rvs.str();
          },
          rvinestructure_doc.str.doc)
      .def(
          "__str__",
          [](const RVineStructure& rvs) {
            return "<pyvinecopulib.core.RVineStructure>\n" + rvs.str();
          },
          rvinestructure_doc.str.doc)
      .def("__getstate__",
           [](const RVineStructure& rvinestruct) {
             return rvinestruct.to_json().dump();
           })
      .def("__setstate__", [](RVineStructure& rvinestruct, std::string state) {
        nlohmann::json json_obj = nlohmann::json::parse(state);
        new (&rvinestruct) RVineStructure(json_obj);
      });

  nb::class_<DVineStructure, RVineStructure>(module, "DVineStructure",
                                             dvinestructure_doc.doc)
      .def(nb::init<const std::vector<size_t>&, size_t>(), "order"_a,
           "trunc_lvl"_a = std::numeric_limits<size_t>::max(),
           dvinestructure_doc.ctor.doc_2args,
           nb::call_guard<nb::gil_scoped_release>())
      .def("__repr__", [](const DVineStructure& rvs) {
        return "<pyvinecopulib.core.DVineStructure>\n" + rvs.str();
      });

  nb::class_<CVineStructure, RVineStructure>(module, "CVineStructure",
                                             cvinestructure_doc.doc)
      .def(nb::init<const std::vector<size_t>&, size_t>(), "order"_a,
           "trunc_lvl"_a = std::numeric_limits<size_t>::max(),
           cvinestructure_doc.ctor.doc_2args,
           nb::call_guard<nb::gil_scoped_release>())
      .def("__repr__", [](const CVineStructure& rvs) {
        return "<pyvinecopulib.core.CVineStructure>\n" + rvs.str();
      });
}
