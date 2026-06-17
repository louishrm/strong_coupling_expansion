#pragma once

#include <map>
#include <vector>
#include <algorithm>
#include "graph.hpp"
#include "fock_space.hpp"

namespace sc_expansion {

  // ---------------------------------------------------------------------------
  // POD structs shared by atomic/ and dimer/ diagrams.
  // ---------------------------------------------------------------------------

  struct HoppingLine {
    int from_vertex;
    int to_vertex;
  };

  struct Lines {
    std::vector<HoppingLine> lines;
  };

  // Per-leg info at a vertex: which hopping line it belongs to and its role.
  struct LegInfo {
    int line_index;
    bool is_source; // true = annihilation (from_vertex), false = creation (to_vertex)
  };

  // A symmetry-reduced global configuration with its accumulated weight.
  // `config` stores the full op_id vector: legs at v0, then v1, ..., in
  // hopping-line iteration order.
  struct ValidGlobalConfig {
    std::vector<uint8_t> config;
    double weight;
  };

  // ---------------------------------------------------------------------------
  // Templated shared helpers.
  // ---------------------------------------------------------------------------

  template <int N_sites> struct GlobalConfiguration {
    std::vector<uint8_t> config; // size = total legs in the diagram

    bool operator<(GlobalConfiguration const &other) const { return this->config < other.config; }
    bool operator==(GlobalConfiguration const &other) const { return this->config == other.config; }
  };

  // Cumulant-symmetry orbit under spin flip (and, for N_sites=2, a no-op
  // reflection placeholder). The dimer's site-swap symmetry is handled inside
  // its spatial-embedding canonicalisation, not here, so get_orbit() considers
  // spin flip only — including reflection here would double-count embeddings.
  template <int N_sites, typename T> struct SymmetryGroup {
    using Config = GlobalConfiguration<N_sites>;

    static Config apply_spin_flip(Config const &c) {
      Config res = c;
      for (auto &op_id : res.config) { op_id = FermionOperator<N_sites, T>(op_id).apply_spin_flip().op; }
      return res;
    }

    static std::vector<Config> get_orbit(Config const &c) {
      std::vector<Config> orbit;
      orbit.push_back(c);
      orbit.push_back(apply_spin_flip(c));
      std::sort(orbit.begin(), orbit.end());
      orbit.erase(std::unique(orbit.begin(), orbit.end()), orbit.end());
      return orbit;
    }

    static Config get_canonical(Config const &c) {
      auto orbit = get_orbit(c);
      return orbit[0];
    }
  };

  // ---------------------------------------------------------------------------
  // Free helpers shared across diagram implementations.
  // ---------------------------------------------------------------------------

  // Square-lattice (bipartite) or triangular-lattice (non-bipartite) embedding
  // count for the single-site case. Lives here rather than on Graph because
  // it presupposes a *lattice*, which is a diagram-level concept.
  int compute_lattice_free_multiplicity(Graph const &graph);

  // Shell-indexed embedding count for rooted graphs. Anchors mark0_vertex at
  // the origin and recurses over the lattice; the returned map buckets final
  // positions of mark1_vertex by squared distance d² = Δx²+Δy². Pass
  // mark1_vertex = -1 for the single-mark case; the returned map then has a
  // single entry at d² = 0 holding the total embedding count.
  std::map<int, int> compute_rooted_shell_multiplicity(Graph const &graph, int mark0_vertex, int mark1_vertex);

  // Enumerate directed hopping lines from the adjacency matrix.
  Lines compute_hopping_lines(Graph const &graph);

  // For each vertex, list its legs (line index + source/dest role).
  std::vector<std::vector<LegInfo>> compute_legs_per_vertex(Graph const &graph, Lines const &hopping_lines);

  // Fermion-loop count → sign (+1 for even loop count, -1 for odd).
  int compute_diagram_sign(int V, Lines const &hopping_lines, std::vector<std::vector<LegInfo>> const &legs_per_vertex);

} // namespace sc_expansion
