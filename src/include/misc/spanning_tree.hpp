#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/random_spanning_tree.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/seed_seq.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <map>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace nb = nanobind;
using namespace nb::literals;

// Select a spanning tree over a candidate graph, mirroring
// vinecopulib::VinecopSelector::select_edges (boost prim / kruskal / Wilson).
// This re-exposes the exact boost routines the compiled selector uses so a
// Python-side structure selector produces the same spanning tree: for the
// deterministic MST algorithms the retained edge set is identical (unique on
// distinct weights); for the random algorithms the draw is reproducible given
// the same graph, weights and seeds.
//
// Parameters
//   n_vertices     : number of vertices (line-graph nodes = edges of the
//                    previous tree).
//   edges          : candidate edges as (v0, v1) vertex-index pairs.
//   weights        : per-candidate graph weight ``w = 1 - (crit >= threshold) *
//                    crit`` (the selector's convention; the weighted-random
//                    path draws with probability proportional to ``1 - w``).
//   tree_algorithm : "mst_prim" | "mst_kruskal" | "random_weighted" |
//                    "random_unweighted".
//   seeds          : RNG seeds for the random algorithms (ignored by the MST
//                    ones). Seeded exactly like FitControlsVinecop.
// Returns the indices (into ``edges``) of the retained spanning-tree edges.
inline std::vector<size_t> select_spanning_tree(
    size_t n_vertices, const std::vector<std::pair<size_t, size_t>>& edges,
    const std::vector<double>& weights, const std::string& tree_algorithm,
    const std::vector<int>& seeds) {
  using Graph = boost::adjacency_list<
      boost::vecS, boost::vecS, boost::undirectedS, boost::no_property,
      boost::property<boost::edge_weight_t, double,
                      boost::property<boost::edge_index_t, size_t>>>;
  using Edge = boost::graph_traits<Graph>::edge_descriptor;

  Graph g(n_vertices);
  auto w_map = boost::get(boost::edge_weight, g);
  auto idx_map = boost::get(boost::edge_index, g);
  for (size_t i = 0; i < edges.size(); ++i) {
    auto res = boost::add_edge(edges[i].first, edges[i].second, g);
    w_map[res.first] = weights[i];
    idx_map[res.first] = i;
  }

  std::vector<size_t> selected;

  // Early exit: the candidate graph is already a spanning tree.
  if (n_vertices > 0 && edges.size() == n_vertices - 1) {
    selected.resize(edges.size());
    for (size_t i = 0; i < edges.size(); ++i) {
      selected[i] = i;
    }
    return selected;
  }

  if (tree_algorithm == "mst_kruskal") {
    std::vector<Edge> tree;
    boost::kruskal_minimum_spanning_tree(g, std::back_inserter(tree));
    for (auto e : tree) {
      selected.push_back(idx_map[e]);
    }
    return selected;
  }

  // Predecessor-map algorithms: Prim and the two random (Wilson) variants.
  std::vector<size_t> pred(n_vertices);
  if (tree_algorithm == "mst_prim") {
    boost::prim_minimum_spanning_tree(g, pred.data());
  } else {
    boost::mt19937 gen;
    if (seeds.empty()) {
      std::random_device rd;
      std::vector<int> s(20);
      for (auto& x : s) {
        x = static_cast<int>(rd());
      }
      boost::random::seed_seq seq(s.begin(), s.end());
      gen.seed(seq);
    } else {
      boost::random::seed_seq seq(seeds.begin(), seeds.end());
      gen.seed(seq);
    }
    // The RNG draw order (root vertex, then the walk) mirrors the selector.
    boost::random::uniform_int_distribution<size_t> root_dist(0,
                                                              n_vertices - 1);
    size_t root = root_dist(gen);
    if (tree_algorithm == "random_unweighted") {
      boost::random_spanning_tree(
          g, gen, boost::predecessor_map(pred.data()).root_vertex(root));
    } else {
      std::map<Edge, double> inv_weights;
      for (auto e : boost::make_iterator_range(boost::edges(g))) {
        inv_weights[e] = 1.0 - w_map[e];
      }
      boost::associative_property_map<std::map<Edge, double>> inv_weight_map(
          inv_weights);
      boost::random_spanning_tree(g, gen,
                                  boost::predecessor_map(pred.data())
                                      .root_vertex(root)
                                      .weight_map(inv_weight_map));
    }
  }

  for (auto e : boost::make_iterator_range(boost::edges(g))) {
    auto s = boost::source(e, g);
    auto t = boost::target(e, g);
    if (pred[s] == t || pred[t] == s) {
      selected.push_back(idx_map[e]);
    }
  }
  return selected;
}

inline void init_spanning_tree(nb::module_& m) {
  m.def("_select_spanning_tree", &select_spanning_tree, "n_vertices"_a,
        "edges"_a, "weights"_a, "tree_algorithm"_a,
        "seeds"_a = std::vector<int>(),
        R"""(Select a spanning tree over a candidate graph.

Internal helper for Python-side vine structure selection: it re-exposes the
boost minimum-/random-spanning-tree routines used by the compiled vinecopulib
selector, so the Python selector chooses the same edges.

Parameters
----------
n_vertices : int
    Number of graph vertices (line-graph nodes).
edges : list of tuple of int
    Candidate edges as ``(v0, v1)`` vertex-index pairs.
weights : list of float
    Per-candidate graph weight ``w = 1 - (crit >= threshold) * crit``.
tree_algorithm : str
    One of ``"mst_prim"``, ``"mst_kruskal"``, ``"random_weighted"``,
    ``"random_unweighted"``.
seeds : list of int, optional
    RNG seeds for the random algorithms (ignored by the MST ones).

Returns
-------
list of int
    Indices into ``edges`` of the retained spanning-tree edges.
)""");
}
