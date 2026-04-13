#pragma once

#include <vector>
#include <map>
#include <set>
#include <numeric>   // For std::accumulate
#include <algorithm> // For std::next_permutation
#include <queue>
#include "./hubbard_solver.hpp"
#include "./cumulant.hpp"
#include "./graph.hpp"
#include "./dual.hpp"
#include "./fock_space.hpp"
#include "./vertex.hpp"

namespace sc_expansion {

  struct HoppingLine {
    int from_vertex;
    int to_vertex;
  };

  struct Lines {
    std::vector<HoppingLine> lines;
  };

  // Helper to store valid physics states during generation
  template <int N_sites> struct GlobalConfiguration {
    std::vector<uint8_t> config; // Size = Total number of legs in the diagram

    bool operator<(const GlobalConfiguration &other) const { return this->config < other.config; }
    bool operator==(const GlobalConfiguration &other) const { return this->config == other.config; }
  };

  // Helper for symmetry operations (Spin Flip and Lattice Reflection)
  template <int N_sites, typename T> struct SymmetryGroup {
    using Config = GlobalConfiguration<N_sites>;

    static Config apply_spin_flip(Config const &c) {
      Config res = c;
      for (auto &op_id : res.config) { op_id = FermionOperator<N_sites, T>(op_id).apply_spin_flip().op; }
      return res;
    }

    static Config apply_reflection(Config const &c) {
      Config res = c;
      if constexpr (N_sites == 2) {
        for (auto &op_id : res.config) { op_id = FermionOperator<N_sites, T>(op_id).apply_reflection().op; }
      }
      return res;
    }

    // Computes all unique configurations in the orbit of the given configuration.
    // Only SpinFlip is used here — for N_sites=2, the Reflect (dimer site-swap)
    // symmetry is equivalent to lattice inversion, which is already applied in
    // canonicalize_directions when merging spatial configurations. Including
    // Reflect here would double-count those embeddings.
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
      return orbit[0]; // Lexicographically smallest
    }
  };

  // Per-leg info at a vertex: which hopping line it belongs to and its role
  struct LegInfo {
    int line_index; // Index into hopping_lines
    bool is_source; // true = annihilation (from_vertex), false = creation (to_vertex)
  };

  // A symmetry-reduced global configuration with its total weight.
  // config stores the full op_id vector: legs at v0, then v1, ..., in hopping-line iteration order.
  // weight = spatial_weight * orbit_size / automorphism_count
  struct ValidGlobalConfig {
    std::vector<uint8_t> config;
    double weight;
  };

  // Stores a unique spatial embedding pattern for the dimer (N_sites=2) case.
  // Each entry in `directions` is a bond label (0-3) for the corresponding hopping line:
  //   0 = horizontal rightward (source site 1, dest site 0)
  //   1 = horizontal leftward  (source site 0, dest site 1)
  //   2 = vertical, site-0 bond (source site 0, dest site 0)
  //   3 = vertical, site-1 bond (source site 1, dest site 1)
  // `weight` is the total number of lattice embeddings producing this pattern.
  struct SpatialConfiguration {
    std::vector<uint8_t> directions; // Per hopping line: bond label 0-3
    double weight;
  };

  // Lattice embedding count for a graph on the square (bipartite) or triangular (non-bipartite) lattice.
  // Moved here from Graph — this is the N_sites=1 (single-site) free multiplicity.
  int compute_lattice_free_multiplicity(Graph const &graph);

  template <int N_sites, typename T> class Diagram {

    public:
    explicit Diagram(Graph const &graph, std::vector<VertexType<N_sites, T> *> const &vertex_types);

    // Cluster-restricted embedding: spatial configs computed from only the given positions.
    // Weights are divided by n_cluster_sites to give per-site (per-dimer) values.
    Diagram(Graph const &graph, std::vector<VertexType<N_sites, T> *> const &vertex_types,
            std::vector<std::pair<int, int>> const &cluster_positions, int n_cluster_sites);

    T evaluate(std::vector<double> const &taus, HubbardSolver<N_sites, T> const &solver, bool infinite_U);

    // Mark all VertexInstances that depend on tau_index as dirty
    void mark_tau_dirty(int tau_index);

    // Mark all VertexInstances as dirty (e.g. after full tau replacement)
    void mark_all_dirty();

    const std::vector<SpatialConfiguration> &get_spatial_configurations() const { return this->spatial_configurations; }

    double get_free_multiplicity() const {
      double total = 0.0;
      for (auto const &sc : this->spatial_configurations) { total += sc.weight; }
      return total;
    }

    const std::vector<ValidGlobalConfig> &get_valid_configurations() const { return this->valid_configurations; }
    double get_diagram_sign() const { return (double)this->diagram_sign; }
    Graph const &get_graph() const { return this->graph; }

    private:
    Graph const &graph;
    std::vector<VertexType<N_sites, T> *> vertex_type_ptrs; // Per-vertex VertexType* for caching (may be nullptr)
    std::vector<ValidGlobalConfig> valid_configurations;
    std::vector<SpatialConfiguration> spatial_configurations;
    std::vector<std::vector<LegInfo>> legs_per_vertex;
    Lines hopping_lines;
    int diagram_sign = 1;

    // VertexInstances: [config_idx][vertex_idx] — local cache per (config, vertex) pair
    std::vector<std::vector<VertexInstance<N_sites, T>>> vertex_instances;

    // Precomputed inverse map: tau_index → list of vertex indices that depend on it.
    // All configs share the same vertex-to-tau mapping, so we only store vertex indices.
    std::vector<std::vector<int>> tau_to_vertices;

    // Factored evaluation data (N_sites=2 only)
    // Per-vertex: the distinct op_id tuples that appear at this vertex across all global configs.
    std::vector<std::vector<std::vector<uint8_t>>> local_states;    // [vertex][state_idx] -> op_ids
    std::vector<std::vector<T>> local_values;                       // [vertex][state_idx] -> cached value
    std::vector<std::vector<T>> local_values_infinite;              // [vertex][state_idx] -> cached value (inf-U)
    std::vector<bool> vertex_dirty_finite;                          // [vertex]
    std::vector<bool> vertex_dirty_infinite;                        // [vertex]
    std::vector<std::vector<int>> config_to_local;                  // [gc_idx][vertex] -> state_idx

    void compute_hopping_lines();
    void setup_vertices(std::vector<VertexType<N_sites, T> *> const &vertex_types);
    void compute_spatial_configurations();
    void compute_valid_configurations();
    void compute_diagram_sign();
    void build_vertex_instances();
    void build_local_state_tables();
    T evaluate_factored(std::vector<double> const &taus, HubbardSolver<N_sites, T> const &solver, bool infinite_U);

    // Helpers for dimer spatial embedding on the rectangular (columnar) superlattice
    void solve_dimer_embedding(int placed_count, std::vector<bool> &placed, std::vector<std::pair<int, int>> &coords,
                               std::map<std::vector<uint8_t>, int> &config_counts) const;

    // Cluster-restricted embedding: only place vertices at positions in cluster_positions
    void solve_cluster_embedding(int placed_count, std::vector<bool> &placed, std::vector<std::pair<int, int>> &coords,
                                 std::map<std::vector<uint8_t>, int> &config_counts,
                                 std::vector<std::pair<int, int>> const &cluster_positions) const;

    void compute_spatial_configurations_cluster(std::vector<std::pair<int, int>> const &cluster_positions, int n_cluster_sites);

    std::vector<uint8_t> canonicalize_directions(std::vector<uint8_t> const &dirs,
                                                 std::vector<std::vector<int>> const &automorphisms) const;

    std::vector<uint8_t> apply_automorphism_to_directions(std::vector<uint8_t> const &dirs, std::vector<int> const &perm) const;
  };
}; // namespace sc_expansion
