#include "canonicalize.hpp"

#include <cassert>
#include <cmath>
#include <stdexcept>

#include "digraph.hh"
#include "stats.hh"

namespace sc_expansion {

  CanonicalResult canonicalize(std::vector<uint8_t> const &adj, int V,
                               std::vector<int> const &initial_color) {

    if (V < 0 || adj.size() != static_cast<size_t>(V) * V)
      throw std::invalid_argument("canonicalize: adj size does not match V*V");
    if (initial_color.size() != static_cast<size_t>(V))
      throw std::invalid_argument("canonicalize: initial_color size != V");

    // Build the bliss multidigraph. Multi-edges are encoded as repeated
    // add_edge calls; bliss::Digraph handles them natively.
    bliss::Digraph g(V);
    for (int v = 0; v < V; ++v) g.change_color(v, initial_color[v]);
    for (int i = 0; i < V; ++i) {
      for (int j = 0; j < V; ++j) {
        uint8_t mult = adj[i * V + j];
        for (uint8_t k = 0; k < mult; ++k) g.add_edge(i, j);
      }
    }

    bliss::Stats stats;
    // canonical_form returns a pointer to a permutation P of length V such
    // that P[i] is the canonical label of input vertex i.
    unsigned int const *perm = g.canonical_form(stats);

    CanonicalResult result;
    result.canonical_permutation.assign(perm, perm + V);

    // Apply the permutation to produce the canonical adjacency. With
    // P[i] = canonical label of i, the canonical entry (P[i], P[j]) equals
    // the original entry (i, j) — i.e. canonical[P[i]*V + P[j]] = adj[i*V + j].
    result.canonical_matrix.assign(static_cast<size_t>(V) * V, 0);
    for (int i = 0; i < V; ++i) {
      for (int j = 0; j < V; ++j) {
        result.canonical_matrix[perm[i] * V + perm[j]] = adj[i * V + j];
      }
    }

    // Group size. Without GMP, bliss exposes only the approximate
    // long-double value. For V up to ~14 (|Aut| <= 14!) this round-trips
    // exactly through long double; we round to nearest integer and verify
    // it landed within tolerance.
    long double g_approx = stats.get_group_size_approx();
    long double rounded = std::round(g_approx);
    if (std::fabs(g_approx - rounded) > 1e-6L)
      throw std::runtime_error(
          "canonicalize: bliss group size not integral within tolerance — "
          "graph may be too large for long-double precision (build bliss with GMP)");
    result.automorphism_count = static_cast<uint64_t>(rounded);

    return result;
  }

} // namespace sc_expansion
