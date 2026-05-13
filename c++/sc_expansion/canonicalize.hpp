#pragma once
#include <cstdint>
#include <vector>

namespace sc_expansion {

  struct CanonicalResult {
    std::vector<uint8_t> canonical_matrix;
    std::vector<unsigned int> canonical_permutation;
    uint64_t automorphism_count;
  };

  // Compute the canonical form of a directed multigraph using bliss.
  //
  // adj            : V*V row-major multiplicity matrix (uint8 entries are edge
  //                  multiplicities; matrix is treated as directed).
  // V              : number of vertices.
  // initial_color  : size-V vector of vertex colors used to seed bliss's
  //                  partition refinement. Vertices in different colors are
  //                  never identified by the canonical labeling. Pass all-zero
  //                  for an uncolored graph; pass {0,...,1,...} to mark a
  //                  subset (e.g. for rooted graphs).
  //
  // Returns the canonically relabeled adjacency, the permutation that maps the
  // input labeling to the canonical one (canonical_permutation[i] = image of
  // old vertex i), and the size of the vertex-automorphism group of the
  // multidigraph (does NOT include multi-edge factorial corrections).
  CanonicalResult canonicalize(std::vector<uint8_t> const &adj, int V,
                               std::vector<int> const &initial_color);

} // namespace sc_expansion
