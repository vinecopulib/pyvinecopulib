#pragma once

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/operators.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <stdexcept>
#include <string>
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
// A truncated structure stores only `trunc_lvl` rows, and upstream's per-entry
// accessors index the triangular array without a bounds check -- so reading a
// tree above the truncation crashes the interpreter. Guard the four Python
// entry points; the whole-array getters are safe because they iterate the rows
// that exist.
inline void rv_check_slot(const RVineStructure& rvs, size_t tree, size_t edge,
                          const char* name) {
  const size_t d = rvs.get_dim();
  const size_t trunc = rvs.get_trunc_lvl();
  if (tree >= trunc || edge + tree + 1 >= d) {
    throw std::runtime_error(
        std::string(name) + "(): tree " + std::to_string(tree) + ", edge " +
        std::to_string(edge) + " is outside a structure with dim " +
        std::to_string(d) + " and trunc_lvl " + std::to_string(trunc) + ".");
  }
}

inline RVineStructure rv_from_struct_array(
    const std::vector<size_t>& order,
    const std::vector<std::vector<size_t>>& struct_array,
    bool natural_order = false, bool check = true) {
  return RVineStructure(order, TriangularArray<size_t>(struct_array),
                        natural_order, check);
}

// Nested-list ``[tree][edge]`` of ``(a, b, conditioning)`` triples (1-based
// variable labels): ``a`` / ``b`` are the conditioned pair, ``conditioning``
// the (possibly empty) conditioning set.
using TreeTuples =
    std::vector<std::vector<std::tuple<size_t, size_t, std::vector<size_t>>>>;

// Normalizes both list-of-trees spellings this package returns: the
// ``(a, b, conditioning)`` triples of ``RVineStructure.get_trees()`` and the
// ``{"conditioned": (a, b), "conditioning": [...], "pair_copula": ...}``
// mappings of ``Vinecop.get_trees()``, whose pair copulas a structure has no
// place for. Requires the GIL.
inline TreeTuples to_tree_tuples(const nb::object& trees) {
  TreeTuples out;
  for (nb::handle tree : nb::cast<nb::sequence>(trees)) {
    std::vector<std::tuple<size_t, size_t, std::vector<size_t>>> edges;
    for (nb::handle edge : nb::cast<nb::sequence>(tree)) {
      if (nb::isinstance<nb::dict>(edge)) {
        const auto mapping = nb::cast<nb::dict>(edge);
        const auto conditioned =
            nb::cast<std::tuple<size_t, size_t>>(mapping["conditioned"]);
        edges.emplace_back(
            std::get<0>(conditioned), std::get<1>(conditioned),
            nb::cast<std::vector<size_t>>(mapping["conditioning"]));
      } else {
        edges.push_back(
            nb::cast<std::tuple<size_t, size_t, std::vector<size_t>>>(edge));
      }
    }
    out.push_back(std::move(edges));
  }
  return out;
}

// Faithful inverse of ``RVineStructure.get_trees()``: the RVineStructure ctor
// from RVineTrees uses the identity diagonal policy (each edge's first
// conditioned variable ``a`` on the diagonal), so
// ``from_trees(s.dim, s.get_trees())`` reproduces ``s`` exactly. This is also
// the diagonal policy ``Vinecop.select`` uses to finalize a selected structure
// (upstream), so ``VinecopBase.select`` can assemble the selected trees through
// this same factory. The diagonal choice controls the variable order (which
// variables sit at the tail — relevant for conditional sampling). Tree ``t``
// must hold exactly ``d - 1 - t`` edges.
inline RVineStructure rv_from_trees(size_t d, const TreeTuples& trees,
                                    bool check = true) {
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
  return RVineStructure(RVineTrees(d, std::move(tree_list)), check);
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

Faithful inverse of ``RVineStructure.get_trees()``: the R-vine matrix is
assembled placing each edge's first conditioned variable ``a`` on the diagonal,
so ``RVineStructure.from_trees(s.dim, s.get_trees())`` reproduces ``s`` exactly.
The diagonal choice fixes the variable order — and hence which variables sit at
its tail, which is what ``Vinecop.sample_conditional`` conditions on.

Parameters
----------
d : int
    Dimension of the vine.
trees : list of list of tuple or list of list of dict
    Per-tree edge lists, indexed ``[tree][edge]``. Each edge is either a triple
    ``(a, b, conditioning)`` of 1-based variable labels, as
    ``RVineStructure.get_trees()`` returns, or a mapping with ``"conditioned"``
    and ``"conditioning"`` keys, as ``Vinecop.get_trees()`` returns; its
    ``"pair_copula"`` entry is ignored. ``a`` and ``b`` are the conditioned pair
    (``a`` is placed on the matrix diagonal) and ``conditioning`` is the
    (possibly empty) conditioning set. Tree ``t`` must contain exactly
    ``d - 1 - t`` edges. An empty list yields a fully truncated (independence)
    structure.
check : bool, default True
    Whether to validate the assembled structure (e.g. the proximity condition).

Returns
-------
RVineStructure
    The assembled structure. Pair-copulas are not part of the input; callers
    fit them on the returned structure.

See Also
--------
RVineStructure.get_trees : The decomposition this method inverts.
Vinecop.get_trees : The same decomposition carrying the fitted pair copulas.
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
      .def_static(
          "from_trees",
          [](size_t d, const nb::object& trees, bool check) {
            const TreeTuples tuples = to_tree_tuples(trees);
            nb::gil_scoped_release release;
            return rv_from_trees(d, tuples, check);
          },
          "d"_a, "trees"_a, "check"_a = true, from_trees_doc)
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
      .def(
          "struct_array",
          [](const RVineStructure& self, size_t tree, size_t edge,
             bool natural_order) {
            rv_check_slot(self, tree, edge, "struct_array");
            return self.struct_array(tree, edge, natural_order);
          },
          "tree"_a, "edge"_a, "natural_order"_a = false,
          rvinestructure_doc.struct_array.doc,
          nb::call_guard<nb::gil_scoped_release>())
      .def(
          "get_struct_array",
          [](const RVineStructure& rvs, bool natural_order) -> nb::list {
            return triangular_to_list(rvs.get_struct_array(natural_order));
          },
          "natural_order"_a = false,
          struct_array_list_doc(rvinestructure_doc.get_struct_array.doc)
              .c_str())
      .def(
          "get_trees",
          [](const RVineStructure& self) -> nb::list {
            // Structure-only list-of-trees view: each edge is a tuple
            // (a, b, conditioning), mirroring the ``from_trees()`` input.
            // ``RVineStructure.from_trees(s.dim, s.get_trees())`` reproduces
            // ``s`` exactly (a faithful, matrix-exact inverse).
            RVineTrees rvt;
            {
              nb::gil_scoped_release release;
              rvt = self.get_trees();
            }
            nb::list out;
            for (const auto& tree : rvt.get_trees()) {
              nb::list edges;
              for (const auto& e : tree)
                edges.append(nb::make_tuple(e.a, e.b, nb::cast(e.C)));
              out.append(edges);
            }
            return out;
          },
          rvinestructure_doc.get_trees.doc)
      .def(
          "__eq__",
          [](const RVineStructure& self, const RVineStructure& other) {
            return self == other;
          },
          nb::is_operator())
      .def(
          "min_array",
          [](const RVineStructure& self, size_t tree, size_t edge) {
            rv_check_slot(self, tree, edge, "min_array");
            return self.min_array(tree, edge);
          },
          "tree"_a, "edge"_a, rvinestructure_doc.min_array.doc,
          nb::call_guard<nb::gil_scoped_release>())
      .def(
          "needed_hfunc1",
          [](const RVineStructure& self, size_t tree, size_t edge) {
            rv_check_slot(self, tree, edge, "needed_hfunc1");
            return self.needed_hfunc1(tree, edge);
          },
          "tree"_a, "edge"_a, rvinestructure_doc.needed_hfunc1.doc,
          nb::call_guard<nb::gil_scoped_release>())
      .def(
          "needed_hfunc2",
          [](const RVineStructure& self, size_t tree, size_t edge) {
            rv_check_slot(self, tree, edge, "needed_hfunc2");
            return self.needed_hfunc2(tree, edge);
          },
          "tree"_a, "edge"_a, rvinestructure_doc.needed_hfunc2.doc,
          nb::call_guard<nb::gil_scoped_release>())
      .def(
          "get_min_array",
          [](const RVineStructure& self) -> nb::list {
            return triangular_to_list(self.get_min_array());
          },
          rvinestructure_doc.get_min_array.doc)
      .def(
          "get_needed_hfunc1",
          [](const RVineStructure& self) -> nb::list {
            return triangular_to_list(self.get_needed_hfunc1());
          },
          rvinestructure_doc.get_needed_hfunc1.doc)
      .def(
          "get_needed_hfunc2",
          [](const RVineStructure& self) -> nb::list {
            return triangular_to_list(self.get_needed_hfunc2());
          },
          rvinestructure_doc.get_needed_hfunc2.doc)
      .def("truncate", &RVineStructure::truncate, "trunc_lvl"_a,
           rvinestructure_doc.truncate.doc,
           nb::call_guard<nb::gil_scoped_release>())
      .def_static("sample", &RVineStructure::simulate, "d"_a,
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
