#include "canonicalize.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "digraph.hh"
#include "stats.hh"

namespace sc_expansion {

  CanonicalResult canonicalize(std::vector<uint8_t> const &adj, int V,
                               std::vector<int> const &initial_color) {

    if (V < 0 || adj.size() != static_cast<size_t>(V) * V)
      throw std::invalid_argument("canonicalize: adj size does not match V*V");
    if (initial_color.size() != static_cast<size_t>(V))
      throw std::invalid_argument("canonicalize: initial_color size != V");

    // Encode the directed multigraph as a colored simple bliss::Digraph by
    // subdividing every edge. bliss::Digraph::add_edge silently drops
    // duplicates (see digraph.hh), so repeated add_edge calls do NOT produce
    // parallel edges. Instead, for each entry adj[i,j] = m, insert m
    // auxiliary "edge" vertices, each connected as i -> aux -> j. Aux
    // vertices share a color distinct from any original vertex color so the
    // canonicalization never mixes the two classes. The encoded graph's
    // automorphism group equals |Aut(multigraph)| * prod_e (m_e!), since the
    // m_e auxiliaries on edge e can be permuted freely among themselves.
    int aux_color = 0;
    for (int c : initial_color)
      if (c + 1 > aux_color) aux_color = c + 1;

    size_t total_aux = 0;
    for (auto m : adj) total_aux += m;

    size_t N = static_cast<size_t>(V) + total_aux;
    bliss::Digraph g(static_cast<unsigned int>(N));
    for (int v = 0; v < V; ++v) g.change_color(v, initial_color[v]);
    for (size_t a = V; a < N; ++a) g.change_color(static_cast<unsigned int>(a), aux_color);

    unsigned int next_aux        = static_cast<unsigned int>(V);
    long double factorial_product = 1.0L;
    for (int i = 0; i < V; ++i) {
      for (int j = 0; j < V; ++j) {
        uint8_t mult = adj[i * V + j];
        for (uint8_t k = 0; k < mult; ++k) {
          g.add_edge(static_cast<unsigned int>(i), next_aux);
          g.add_edge(next_aux, static_cast<unsigned int>(j));
          ++next_aux;
        }
        long double f = 1.0L;
        for (uint8_t k = 2; k <= mult; ++k) f *= k;
        factorial_product *= f;
      }
    }

    bliss::Stats stats;
    // perm has length N: perm[v] = canonical label of input vertex v in the
    // encoded graph. We only care about its restriction to the original V.
    unsigned int const *perm = g.canonical_form(stats);

    // Rank the originals by their encoded-canonical labels. orig_rank[i] is
    // the position of original vertex i in the canonical ordering of just
    // the originals (0..V-1). Aux vertices are dropped from the output.
    std::vector<std::pair<unsigned int, int>> label_idx;
    label_idx.reserve(V);
    for (int i = 0; i < V; ++i) label_idx.emplace_back(perm[i], i);
    std::sort(label_idx.begin(), label_idx.end());
    std::vector<int> orig_rank(V);
    for (int r = 0; r < V; ++r) orig_rank[label_idx[r].second] = r;

    CanonicalResult result;
    result.canonical_permutation.assign(V, 0u);
    for (int i = 0; i < V; ++i)
      result.canonical_permutation[i] = static_cast<unsigned int>(orig_rank[i]);

    // Project the canonical relabeling back onto a V*V multiplicity matrix.
    result.canonical_matrix.assign(static_cast<size_t>(V) * V, 0);
    for (int i = 0; i < V; ++i)
      for (int j = 0; j < V; ++j)
        result.canonical_matrix[orig_rank[i] * V + orig_rank[j]] = adj[i * V + j];

    // Recover |Aut(multigraph)| = |Aut(encoded)| / prod_e (m_e!).
    long double g_approx = stats.get_group_size_approx() / factorial_product;
    long double rounded  = std::round(g_approx);
    if (std::fabs(g_approx - rounded) > 1e-6L)
      throw std::runtime_error(
          "canonicalize: bliss group size not integral within tolerance — "
          "graph may be too large for long-double precision (build bliss with GMP)");
    result.automorphism_count = static_cast<uint64_t>(rounded);

    return result;
  }

} // namespace sc_expansion
