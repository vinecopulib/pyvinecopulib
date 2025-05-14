import copy
import time

import numpy as np

import pyvinecopulib as pv


def matrix_to_trees(M: np.ndarray):
  """
  Converts a vine array (lower triangular R-vine matrix) into a list of tree graphs.

  Parameters
  ----------
  M : np.ndarray
      A (d, d) lower-triangular R-vine matrix encoding the vine structure.

  Returns
  -------
  trees : list[dict]
      A list of trees, where each tree is a dictionary with:
          - 'nodes': node ids and their associated information
          - 'edges': list of (node1, node2, (a, b, C)) triples
            where (a, b) are conditioned variables and C is the conditioning set.
  """
  d = M.shape[0]
  trees = []

  for t in range(d - 1):
    edges = []

    # In tree T0, nodes are simply the variable indices
    if t == 0:
      prev_edges = {idx: int(M[d - 1 - idx, idx]) for idx in range(d)}
    else:
      # In higher trees Tt, nodes correspond to edges of T(t-1)
      prev_edges = {idx: x[2] for idx, x in enumerate(trees[t - 1]["edges"])}

      # Build a lookup table:
      # (variable, conditioning set) -> node index
      lookup = {}
      for node_id, (a, b, cond_C) in prev_edges.items():
        lookup[(a, frozenset(cond_C.union({b})))] = node_id
        lookup[(b, frozenset(cond_C.union({a})))] = node_id

    for e in range(d - 1 - t):
      # Extract pair (a, b) and conditioning set C for each edge
      a = int(M[d - 1 - e, e])
      b = int(M[t, e])
      C = frozenset(map(int, M[:t, e]))

      edge_info = (a, b, C)

      if t == 0:
        # At tree 0, node ids are variables themselves
        node1, node2 = a, b
      else:
        # In higher trees, find corresponding node ids
        node1 = lookup.get((a, C), None)
        node2 = lookup.get((b, C), None)

      edges.append((node1, node2, edge_info))

    # Store the constructed tree
    trees.append({"nodes": prev_edges, "edges": edges})

  return trees


def matrix_to_trees2(M: np.ndarray):
  """
  Converts a vine array (lower triangular R-vine matrix) into a list of tree graphs.

  Parameters
  ----------
  M : np.ndarray
      A (d, d) lower-triangular R-vine matrix encoding the vine structure.

  Returns
  -------
  edge_infos : list[list[int, int, set[int]]]
      A list of edge_infos, where each element corresponds to a given tree,
      i.e. a list of (a, b, C) triples.
  """
  d = M.shape[0]
  edge_infos = []

  for t in range(d - 1):
    edges = []

    for e in range(d - 1 - t):
      # Extract pair (a, b) and conditioning set C for each edge
      a = int(M[d - 1 - e, e])
      b = int(M[t, e])
      C = frozenset(map(int, M[:t, e]))

      edges.append((a, b, C))

    # Store the constructed tree
    edge_infos.append(edges)

  return edge_infos


def trees_to_matrix(trees):
  """
  Converts a list of trees back into a vine array (lower triangular R-vine matrix).

  Parameters
  ----------
  trees : list[dict]
      A list of trees, as produced by `matrix_to_trees`, where each tree is a
      dictionary with 'nodes' and 'edges'.

  Returns
  -------
  M : np.ndarray
      A (d, d) lower-triangular R-vine matrix encoding the vine structure.
  """
  if len(trees) < 1:
    raise ValueError("The input trees list must contain at least one tree.")
  d = len(trees[0]["nodes"])  # Infer d from first tree nodes
  trunc_lvl = len(trees)
  M = np.zeros((d, d), dtype=np.int64)
  order = [0] * d

  # Work on a deep copy to modify edges without changing original trees
  trees = copy.deepcopy(trees)
  degrees = []
  for t in range(trunc_lvl):
    assert len(trees[t]["edges"]) == d - t - 1, (
      f"Number of edges in tree {t} is incorrect: {len(trees[t]['edges'])} != {d - t - 1}"
    )

    # Compute degree (number of connections) for each node
    degrees.append({})
    for n1, n2, _ in trees[t]["edges"]:
      degrees[t][n1] = degrees[t].get(n1, 0) + 1
      degrees[t][n2] = degrees[t].get(n2, 0) + 1
    edge_with_leaves = set()
    for idx, (n1, n2, _) in enumerate(trees[t]["edges"]):
      if degrees[t][n1] == 1 or degrees[t][n2] == 1:
        edge_with_leaves.add(idx)

  for col in range(d - 1):
    # At each column, we reconstruct one diagonal element and its upper entries
    t = max(min(trunc_lvl, d - 1 - col), 1)  # matrix above trunc_lvl is empty
    # t = d - 1 - col  # Corresponding tree level (from top down)

    tree = t - 1
    tree_edges = trees[tree]["edges"]
    tree_degrees = degrees[tree]
    for idx, (n1, n2, (a, b, C)) in enumerate(tree_edges):
      if tree_degrees[n1] == 1 or tree_degrees[n2] == 1:
        # Choose the leaf node and the "other" variable
        if tree_degrees[n1] == 1:
          leaf_var, other_var = a, b
        else:
          leaf_var, other_var = b, a

        # Set the diagonal element (leaf variable)
        order[col] = leaf_var

        # Set the first off-diagonal element
        M[t - 1, col] = other_var

        # Save the conditioning set associated with the edge
        ning_set = set(C)

        # Remove the used edge to prevent reuse
        del tree_edges[idx]
        tree_degrees[n1] -= 1
        tree_degrees[n2] -= 1
        break
    else:
      raise RuntimeError("No leaf edge found in tree.")

    # Fill the rest of the column by traversing to lower trees
    for k in range(1, t):
      check_set = set(ning_set)
      check_set.add(order[col])  # Update with current diagonal variable

      tree = t - 1 - k
      tree_edges = trees[tree]["edges"]
      tree_degrees = degrees[tree]

      # Find an edge in the lower tree matching the expanded conditioning set
      for idx, (n1, n2, (a, b, C)) in enumerate(tree_edges):
        all_indices = {a, b}.union(C)

        if all_indices == check_set:
          # Found the matching edge
          next_var = b if a == order[col] else a

          # Fill the corresponding matrix entry
          M[t - 1 - k, col] = next_var

          # Update conditioning set
          ning_set = set(C)

          # Remove the used edge
          del tree_edges[idx]
          tree_degrees[n1] -= 1
          tree_degrees[n2] -= 1
          break
      else:
        raise RuntimeError("No matching edge found in lower tree.")

  # The last diagonal element: copy from the second-to-last column
  M[0, d - 1] = M[0, d - 2]

  # Fill the remaining diagonal elements
  for col in range(d - 1):
    M[d - 1 - col, col] = order[col]

  return M


def edge_infos_to_trees(edge_infos_by_tree):
  """
  Rebuild vine trees from merged (a, b, C) triples.

  Parameters
  ----------
  edge_infos_by_tree : list[list[tuple]]
      List of (a, b, C) triples per tree level.

  Returns
  -------
  rebuilt_trees : list[dict]
      List of trees with nodes, edges, and degrees.
  """
  rebuilt_trees = []

  for t, edge_infos in enumerate(edge_infos_by_tree):
    nodes = {}
    edges = []

    if t == 0:
      # At tree 0, nodes are the variables directly
      variables = set()
      for a, b, C in edge_infos:
        variables.add(a)
        variables.add(b)

      var_list = sorted(variables)
      var_to_node = {var: idx for idx, var in enumerate(var_list)}
      nodes = {idx: var for idx, var in enumerate(var_list)}

      # Now recover edges from (a, b, C)
      for a, b, C in edge_infos:
        n1 = var_to_node[a]
        n2 = var_to_node[b]
        edges.append((n1, n2, (a, b, C)))

    else:
      # At higher trees, nodes are edges from previous tree
      prev_edges = rebuilt_trees[t - 1]["edges"]

      # Build nodes dict and a lookup table with
      # (variable, conditioning set) -> node index
      nodes = {}
      lookup = {}
      for node_id, (_, _, (a, b, C)) in enumerate(prev_edges):
        nodes[node_id] = (a, b, C)
        lookup[(a, frozenset(C.union({b})))] = node_id
        lookup[(b, frozenset(C.union({a})))] = node_id

      # Now recover edges at level t
      for a, b, C in edge_infos:
        node1 = lookup.get((a, C))
        node2 = lookup.get((b, C))
        if node1 is None or node2 is None:
          raise ValueError(
            f"No matching nodes for edge {(a, b, C)} in tree {t}"
          )
        edges.append((node1, node2, (a, b, C)))

    rebuilt_trees.append({"nodes": nodes, "edges": edges})

    # Compute degree (number of connections) for each node
    degrees = {}
    for n1, n2, _ in edges:
      degrees[n1] = degrees.get(n1, 0) + 1
      degrees[n2] = degrees.get(n2, 0) + 1
    rebuilt_trees[t]["degrees"] = degrees

  return rebuilt_trees


def trees_to_matrix2(edge_infos):
  """
  Converts a list of trees back into a vine array (lower triangular R-vine matrix).

  Parameters
  ----------
  edge_infos : list[list[tuple[int, int, set[int]]]]
      A list of edge_infos, as produced by `matrix_to_trees`, where each element corresponds to a given tree,
      i.e. a list of (a, b, C) triples.

  Returns
  -------
  M : np.ndarray
      A (d, d) lower-triangular R-vine matrix encoding the vine structure.
  """

  # Extract trees from edge_infos
  trees = edge_infos_to_trees(edge_infos)

  if len(trees) < 1:
    raise ValueError("The input trees list must contain at least one tree.")
  d = len(trees[0]["nodes"])  # Infer d from first tree nodes
  trunc_lvl = len(trees)
  M = np.zeros((d, d), dtype=np.int64)
  order = [0] * d

  for t in range(trunc_lvl):
    assert len(trees[t]["edges"]) == d - t - 1, (
      f"Number of edges in tree {t} is incorrect: {len(trees[t]['edges'])} != {d - t - 1}"
    )

  for col in range(d - 1):
    # At each column, we reconstruct one diagonal element and its upper entries
    t = max(min(trunc_lvl, d - 1 - col), 1)  # matrix above trunc_lvl is empty
    # t = d - 1 - col  # Corresponding tree level (from top down)

    tree = t - 1
    tree_edges = trees[tree]["edges"]
    tree_degrees = trees[tree]["degrees"]
    for idx, (n1, n2, (a, b, C)) in enumerate(tree_edges):
      if tree_degrees[n1] == 1 or tree_degrees[n2] == 1:
        # Choose the leaf node and the "other" variable
        if tree_degrees[n1] == 1:
          leaf_var, other_var = a, b
        else:
          leaf_var, other_var = b, a

        # Set the diagonal element (leaf variable)
        order[col] = leaf_var

        # Set the first off-diagonal element
        M[t - 1, col] = other_var

        # Save the conditioning set associated with the edge
        ning_set = set(C)

        # Remove the used edge to prevent reuse
        del tree_edges[idx]
        tree_degrees[n1] -= 1
        tree_degrees[n2] -= 1
        break
    else:
      raise RuntimeError("No leaf edge found in tree.")

    # Fill the rest of the column by traversing to lower trees
    for k in range(1, t):
      check_set = set(ning_set)
      check_set.add(order[col])  # Update with current diagonal variable

      tree = t - 1 - k
      tree_edges = trees[tree]["edges"]
      tree_degrees = trees[tree]["degrees"]

      # Find an edge in the lower tree matching the expanded conditioning set
      for idx, (n1, n2, (a, b, C)) in enumerate(tree_edges):
        all_indices = {a, b}.union(C)

        if all_indices == check_set:
          # Found the matching edge
          next_var = b if a == order[col] else a

          # Fill the corresponding matrix entry
          M[t - 1 - k, col] = next_var

          # Update conditioning set
          ning_set = set(C)

          # Remove the used edge
          del tree_edges[idx]
          tree_degrees[n1] -= 1
          tree_degrees[n2] -= 1
          break
      else:
        raise RuntimeError("No matching edge found in lower tree.")

  # The last diagonal element: copy from the second-to-last column
  M[0, d - 1] = M[0, d - 2]

  # Fill the remaining diagonal elements
  for col in range(d - 1):
    M[d - 1 - col, col] = order[col]

  return M


# Test matrix_to_trees and trees_to_matrix
d = 100
t = time.time()
rvs = pv.RVineStructure.simulate(d)
print("Simulation time:", time.time() - t)
t = time.time()
mat = rvs.matrix
print("Matrix time:", time.time() - t)
t = time.time()
trees = matrix_to_trees(mat)
print("Trees extraction time:", time.time() - t)
t = time.time()
edge_infos = matrix_to_trees2(mat)
print("Trees extraction 2 time:", time.time() - t)
t = time.time()
reconstructed_mat = trees_to_matrix(trees)
print("Reconstruction time:", time.time() - t)
assert np.array_equal(mat, reconstructed_mat)
t = time.time()
reconstructed_mat2 = trees_to_matrix2(edge_infos)
print("Reconstruction 2 time:", time.time() - t)
assert np.array_equal(mat, reconstructed_mat2)


# Test edge_infos_to_trees
# I.e., reconstruction from (a, b, C) triples
edge_infos = [
  [edge_info for _, _, edge_info in tree["edges"]] for tree in trees
]
t = time.time()
trees2 = edge_infos_to_trees(edge_infos)
reconstructed_mat2 = trees_to_matrix(trees2)
print("Trees from edge infos time:", time.time() - t)
t = time.time()
reconstructed_mat3 = trees_to_matrix2(edge_infos)
print("Trees from edge infos 2 time:", time.time() - t)
assert np.array_equal(mat, reconstructed_mat2)
assert np.array_equal(mat, reconstructed_mat3)


def truncate_matrix(mat, trunc_level):
  """
  Truncate a matrix to a given level.

  Parameters
  ----------
  mat : np.ndarray
      The matrix to truncate.
  trunc_level : int
      The level to truncate to.

  Returns
  -------
  truncated_mat : np.ndarray
      The truncated matrix.
  """
  d = mat.shape[0]
  if trunc_level >= d:
    raise ValueError("Truncation level must be less than the matrix size.")
  # Create a mask for the upper triangular part of the matrix
  mask = np.fromfunction(
    lambda i, j: (i + j < d - 1) & (i >= trunc_level),
    mat.shape,
    dtype=int,
  )
  # Set the masked elements to zero
  truncated_mat = mat.copy()
  truncated_mat[mask] = 0
  return truncated_mat


def merge_list_of_trees(list_of_trees, trunc_level=None):
  """
  Merge multiple vines structures.

  Parameters
  ----------
  trees : list[list[dict]]

  Returns
  -------
  merged_trees : list[dict]
      Merged vine structure.
  """
  if trunc_level is None:
    max_depth = max(len(trees) for trees in list_of_trees)

  edge_infos = []
  for t in range(max_depth):
    edge_infos.append([])

    for trees in list_of_trees:
      if t < len(trees):
        edge_infos[t].extend(
          [edge_info for _, _, edge_info in trees[t]["edges"]]
        )

  return edge_infos_to_trees(edge_infos)


def fill_missing_edges(trees):
  """
  Fill missing edges in a list of vine trees to make a complete vine structure.

  Parameters
  ----------
  trees : list[dict]
      List of trees with possibly missing edges.

  Returns
  -------
  completed_trees : list[dict]
      Fully connected vine trees ready for matrix conversion.
  """
  completed_trees = copy.deepcopy(trees)
  d = len(trees[0]["nodes"])

  for t, tree in enumerate(completed_trees):
    target_edges = d - 1 - t
    current_edges = len(tree["edges"])
    print(f"Tree {t}: {current_edges} edges, target: {target_edges}")

    if current_edges == target_edges:
      continue  # Already full

    if t > 0:
      completed_trees[t]["nodes"] = {
        idx: edge[2] for idx, edge in enumerate(completed_trees[t - 1]["edges"])
      }
    nodes = completed_trees[t]["nodes"]
    edges = completed_trees[t]["edges"]

    # Build lookup of node -> involved variables
    node_to_vars = {}
    for idx, info in nodes.items():
      if t == 0:
        # Node info is a variable
        node_to_vars[idx] = {info}
      else:
        # Node info is (a, b, C)
        a, b, C = info
        node_to_vars[idx] = {a, b}.union(C)

    # Add missing edges
    node_list = list(nodes.keys())
    parent = {node: node for node in node_list}

    # Union-Find structure to detect cycles
    def find(u):
      while parent[u] != u:
        parent[u] = parent[parent[u]]  # Path compression
        u = parent[u]
      return u

    def union(u, v):
      parent_u = find(u)
      parent_v = find(v)
      if parent_u == parent_v:
        return False  # Would create a cycle
      parent[parent_u] = parent_v
      return True

    # Set of already connected pairs to avoid duplicates
    connected = set()
    for n1, n2, _ in edges:
      union(n1, n2)
      connected.add((n1, n2))

    allowed_edges = []
    for i in range(len(node_list)):
      for j in range(i + 1, len(node_list)):
        n1, n2 = node_list[i], node_list[j]
        if (n1, n2) in connected:
          continue  # Already connected

        vars1 = node_to_vars[n1]
        vars2 = node_to_vars[n2]
        shared = vars1.intersection(vars2)

        if len(shared) == len(vars1) - 1:
          a, b = vars1.difference(shared).pop(), vars2.difference(shared).pop()
          C = frozenset(shared)
          allowed_edges.append((n1, n2, (a, b, C)))

    while len(edges) < target_edges and len(allowed_edges) > 0:
      n1, n2, (a, b, C) = allowed_edges.pop()
      if find(n1) == find(n2):
        continue  # Would create a cycle

      edges.append((n1, n2, (a, b, C)))
      connected.add((n1, n2))
      union(n1, n2)

    if len(edges) < target_edges:
      raise RuntimeError(f"Could not fill tree {t}.")

  return completed_trees


# # Works for any truncation level
# trunc_level = 3
# reconstructed_mat = trees_to_matrix(trees[:trunc_level])
# truncated_mat = truncate_matrix(mat, trunc_level)
# assert np.array_equal(truncated_mat, reconstructed_mat)

# # Now, let's test the merging of multiple trees
# # Create a set of local vines and a bridging vine
# n_variables = 3
# n_groups = 3
# seeds = [0, 1, 2]

# # Create local vines
# local_rvs = [
#   pv.RVineStructure.simulate(n_variables, seeds=seeds) for _ in range(n_groups)
# ]

# # Convert to matrices
# local_matrices = [
#   local_rvs[i].matrix
#   + i
#   * n_variables
#   * np.fliplr(np.triu(np.ones((n_variables, n_variables), dtype=np.uint64)))
#   for i in range(n_groups)
# ]

# # Convert to trees
# local_trees = [matrix_to_trees(local_matrices[i]) for i in range(n_groups)]

# # Create the bridging vine
# bridges = [0]
# bridges.extend([mat[0, n_variables - 1] for mat in local_matrices])
# bridges = np.array(bridges, dtype=np.uint64)
# briging_rvs = pv.RVineStructure.simulate(n_groups, seeds=seeds)
# briging_matrix = briging_rvs.matrix
# briging_matrix = bridges[briging_matrix]

# # Print the local and bridging matrices
# for i, local_matrix in enumerate(local_matrices):
#   print(f"Local matrix {i}:\n{local_matrix}")
# print("Bridging matrix:\n", briging_matrix)

# # Convert to trees
# bridging_trees = matrix_to_trees(briging_matrix)

# # Merge the two sets of trees
# merged_trees = merge_list_of_trees(local_trees + [bridging_trees])

# # Fill missing edges in the merged trees
# completed_trees = fill_missing_edges(merged_trees)

# # # Print the completed trees
# # print("Completed trees:")
# # for t, tree in enumerate(completed_trees):
# #   n_edges = len(tree["edges"])
# #   n_nodes = len(tree["nodes"])
# #   print(f"Tree {t}: {n_edges} edges, {n_nodes} nodes")
# #   print("Nodes:", tree["nodes"])
# #   print("Edges:", tree["edges"])

# # Convert completed trees back to matrix
# merged_matrix = trees_to_matrix(completed_trees)

# # Print the merged matrix
# print(merged_matrix)
