import time

import numpy as np

import pyvinecopulib as pv


def matrix_to_trees(M: np.ndarray):
  """
  Converts a vine array (lower triangular R-vine matrix) into a list of trees.

  Parameters
  ----------
  M : np.ndarray
      A (d, d) lower-triangular R-vine matrix encoding the vine structure.

  Returns
  -------
  trees : list[list[int, int, set[int]]]
      A list where each element corresponds to a given tree,
      i.e. a list of (a, b, C) triplets, where {a, b} is the
      conditioned pair and C is the conditioning set.
  """
  d = M.shape[0]
  trees = []

  for t in range(d - 1):
    edges = []

    for e in range(d - 1 - t):
      # Extract pair (a, b) and conditioning set C for each edge
      a = int(M[d - 1 - e, e])
      b = int(M[t, e])
      C = frozenset(map(int, M[:t, e]))

      edges.append((a, b, C))

    # Store the constructed tree
    trees.append(edges)

  return trees


def trees_to_augmented_trees(trees):
  """
  Reconstructs the trees from a list of edges infos and adds information about nodes and degrees in order to reconstruct the vine structure.

  Parameters
  ----------
  trees : list[list[tuple[int, int, set[int]]]]
      A list, as produced e.g. by `matrix_to_trees`,
      where each element corresponds to a given tree,
      i.e. a list of (a, b, C) triplets, where {a, b} is the
      conditioned pair and C is the conditioning set.

  Returns
  -------
  augmented_trees : list[dict]
      List of trees with nodes, edges, and degrees.
  """
  augmented_trees = []

  for t, edge_infos in enumerate(trees):
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
        node1 = var_to_node[a]
        node2 = var_to_node[b]
        edges.append((node1, node2, (a, b, C)))

    else:
      # At higher trees, nodes are edges from previous tree
      prev_edges = augmented_trees[t - 1]["edges"]

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

    augmented_trees.append({"nodes": nodes, "edges": edges})

    # Compute degree (number of connections) for each node
    degrees = {}
    for node1, node2, _ in edges:
      degrees[node1] = degrees.get(node1, 0) + 1
      degrees[node2] = degrees.get(node2, 0) + 1
    augmented_trees[t]["degrees"] = degrees

  return augmented_trees


def trees_to_matrix(trees):
  """
  Converts a list of trees back into a vine array (lower triangular R-vine matrix).

  Parameters
  ----------
  trees : list[list[tuple[int, int, set[int]]]]
      A list, as produced e.g. by `matrix_to_trees`,
      where each element corresponds to a given tree,
      i.e. a list of (a, b, C) triplets, where {a, b} is the
      conditioned pair and C is the conditioning set.

  Returns
  -------
  M : np.ndarray
      A (d, d) lower-triangular R-vine matrix encoding the vine structure.
  """

  # Extract trees from edge_infos
  trees = trees_to_augmented_trees(trees)

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
    for idx, (node1, node2, (a, b, C)) in enumerate(tree_edges):
      if tree_degrees[node1] == 1 or tree_degrees[node2] == 1:
        # Choose the leaf node and the "other" variable
        if tree_degrees[node1] == 1:
          leaf_var, other_var = a, b
        else:
          leaf_var, other_var = b, a

        # Set the diagonal element (leaf variable)
        order[col] = leaf_var

        # Set the first off-diagonal element
        M[t - 1, col] = other_var

        # Save the conditioning set associated with the edge
        ning_set = set(C)

        print(f"(a, b): {a}, {b}")
        print("check_set:", ning_set)
        print(f"node1: {node1}, node2: {node2}")
        print(f"degrees: {tree_degrees[node1]}, {tree_degrees[node2]}")

        # Remove the used edge to prevent reuse
        del tree_edges[idx]
        tree_degrees[node1] -= 1
        tree_degrees[node2] -= 1
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
      for idx, (node1, node2, (a, b, C)) in enumerate(tree_edges):
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
          tree_degrees[node1] -= 1
          tree_degrees[node2] -= 1
          break
      else:
        raise RuntimeError("No matching edge found in lower tree.")

  # The last diagonal element: copy from the second-to-last column
  M[0, d - 1] = M[0, d - 2]

  # Fill the remaining diagonal elements
  for col in range(d - 1):
    M[d - 1 - col, col] = order[col]

  return M


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
  trees : list[list[list[tuple[int, int, set[int]]]]]
      A list of vines, where each vine is represented as a list of trees,
      and each tree is represented as a list of (a, b, C) triplets,
      where {a, b} is the conditioned pair and C is the conditioning set.

  Returns
  -------
  merged_trees : list[list[tuple[int, int, set[int]]]]
      Merged vine structure.
  """
  if trunc_level is None:
    max_depth = max(len(trees) for trees in list_of_trees)

  merged_trees = []
  for t in range(max_depth):
    merged_trees.append([])

    for trees in list_of_trees:
      if t < len(trees):
        merged_trees[t].extend(trees[t])

  return merged_trees


def fill_missing_edges(trees):
  """
  Fill missing edges in a list of vine trees to make a complete vine structure.

  Parameters
  ----------
  trees : list[list[tuple[int, int, set[int]]]]
      A list of trees, where each tree is represented as a list of
      (a, b, C) triplets, where {a, b} is the conditioned pair and C is
      the conditioning set.

  Returns
  -------
  completed_trees : list[list[tuple[int, int, set[int]]]]
      The same, but fully connected trees.
  """
  completed_trees = trees_to_augmented_trees(trees)
  d = len(completed_trees[0]["nodes"])

  for t, tree in enumerate(completed_trees):
    target_edges = d - 1 - t
    current_edges = len(tree["edges"])
    # print(f"Tree {t}: {current_edges} edges, target: {target_edges}")

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
    for node1, node2, _ in edges:
      union(node1, node2)
      connected.add((node1, node2))

    allowed_edges = []
    for i in range(len(node_list)):
      for j in range(i + 1, len(node_list)):
        node1, node2 = node_list[i], node_list[j]
        if (node1, node2) in connected:
          continue  # Already connected

        vars1 = node_to_vars[node1]
        vars2 = node_to_vars[node2]
        shared = vars1.intersection(vars2)

        if len(shared) == len(vars1) - 1:
          a, b = vars1.difference(shared).pop(), vars2.difference(shared).pop()
          C = frozenset(shared)
          allowed_edges.append((node1, node2, (a, b, C)))

    while len(edges) < target_edges and len(allowed_edges) > 0:
      node1, node2, (a, b, C) = allowed_edges.pop()
      if find(node1) == find(node2):
        continue  # Would create a cycle

      edges.append((node1, node2, (a, b, C)))
      connected.add((node1, node2))
      union(node1, node2)

    if len(edges) < target_edges:
      raise RuntimeError(f"Could not fill tree {t}.")

  return [
    [edge_info for (_, _, edge_info) in tree["edges"]]
    for tree in completed_trees
  ]


# Test matrix_to_trees and trees_to_matrix
# d = 6
# t = time.time()
# rvs = pv.RVineStructure.simulate(d)
# print("Simulation time:", time.time() - t)
mat = np.array(
  [
    [5, 2, 6, 6, 6, 6, 6],
    [6, 6, 1, 2, 5, 5, 0],
    [2, 5, 2, 5, 2, 0, 0],
    [1, 1, 5, 1, 0, 0, 0],
    [3, 7, 7, 0, 0, 0, 0],
    [7, 3, 0, 0, 0, 0, 0],
    [4, 0, 0, 0, 0, 0, 0],
  ]
)
rvs = pv.RVineStructure.from_matrix(mat)
t = time.time()
mat = rvs.matrix
print("Matrix time:", time.time() - t)
t = time.time()
trees = matrix_to_trees(mat)
print("Trees extraction time:", time.time() - t)
t = time.time()
reconstructed_mat = trees_to_matrix(trees)
print("Reconstruction time:", time.time() - t)
assert np.array_equal(mat, reconstructed_mat)


# Works for any truncation level
trunc_level = 3
reconstructed_mat = trees_to_matrix(trees[:trunc_level])
truncated_mat = truncate_matrix(mat, trunc_level)
assert np.array_equal(truncated_mat, reconstructed_mat)

# # Now, let's test the merging of multiple trees
# # Create a set of local vines and a bridging vine
# n_variables = 5
# n_groups = 100
# seeds = [0, 1, 2]
# show_matrices = False

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
# if show_matrices:
#   print("Local matrices:")
#   for i, local_matrix in enumerate(local_matrices):
#     print(f"Local matrix {i}:\n{local_matrix}")
#   print("Bridging matrix:\n", briging_matrix)

# # Convert to trees
# bridging_trees = matrix_to_trees(briging_matrix)

# # Merge the two sets of trees
# t = time.time()
# merged_trees = merge_list_of_trees(local_trees + [bridging_trees])
# print("Merging time:", time.time() - t)

# # Fill missing edges in the merged trees
# t = time.time()
# completed_trees = fill_missing_edges(merged_trees)
# print("Filling missing edges time:", time.time() - t)

# # Convert completed trees back to matrix
# t = time.time()
# merged_matrix = trees_to_matrix(completed_trees)
# print("Reconstruction time:", time.time() - t)
# assert pv.RVineStructure.from_matrix(merged_matrix) is not None

# # Print the merged matrix
# if show_matrices:
#   print("Merged matrix:")
#   print(merged_matrix)
