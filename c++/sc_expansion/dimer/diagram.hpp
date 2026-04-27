#pragma once

#include <vector>
#include <map>
#include <utility>
#include "../graph.hpp"
#include "../hubbard_solver.hpp"
#include "../cumulant.hpp"
#include "../diagram_common.hpp"
#include "vertex.hpp"

namespace sc_expansion::dimer {

  // A unique spatial embedding pattern on the staggered (triangular) dimer
  // superlattice. Dimer at (u,v) covers physical sites (2u + v%2, v) and
  // (2u + v%2 + 1, v); each dimer has 6 NN dimers, all of which cross the A/B
  // sublattice (src ≠ dst). Since the local cumulant value depends only on
  // (source_site, dest_site) — not on which of the three directions in each
  // class — `directions[k]` is a binary bond label per hopping line k:
  //   0 = src=1, dst=0   (dimer offsets (+1, 0), (0, +1), (0, -1))
  //   1 = src=0, dst=1   (dimer offsets (-1, 0), (-1, +1), (-1, -1))
  // `weight` is the total number of lattice embeddings producing this pattern.
  struct SpatialConfiguration {
    std::vector<uint8_t> directions;
    double weight;
  };

  template <typename T> class Diagram {
    public:
    explicit Diagram(Graph const &graph, std::vector<VertexType<T> *> const &vertex_types);

    // Cluster-restricted embedding: spatial configs computed from only the given
    // positions. Weights are divided by n_cluster_sites for per-dimer values.
    Diagram(Graph const &graph, std::vector<VertexType<T> *> const &vertex_types,
            std::vector<std::pair<int, int>> const &cluster_positions, int n_cluster_sites);

    T evaluate(std::vector<double> const &taus, HubbardSolver<2, T> const &solver);

    // Diagnostic: per-config signed contributions (diagram_sign * w_c * prod_v C_v),
    // WITHOUT the -1/beta prefactor. Forces full recomputation.
    std::vector<T> evaluate_per_config(std::vector<double> const &taus, HubbardSolver<2, T> const &solver);

    void mark_tau_dirty(int tau_index);
    void mark_all_dirty();

    std::vector<SpatialConfiguration> const &get_spatial_configurations() const { return this->spatial_configurations; }

    double get_free_multiplicity() const {
      double total = 0.0;
      for (auto const &sc : this->spatial_configurations) total += sc.weight;
      return total;
    }

    std::vector<ValidGlobalConfig> const &get_valid_configurations() const { return this->valid_configurations; }
    double get_diagram_sign() const { return (double)this->diagram_sign; }
    Graph const &get_graph() const { return this->graph; }
    std::pair<long, long> get_local_cache_stats() const { return {this->local_cache_hits, this->local_cache_misses}; }
    std::vector<int> get_local_state_counts() const {
      std::vector<int> counts;
      for (auto const &ls : this->local_states) counts.push_back((int)ls.size());
      return counts;
    }
    std::vector<std::vector<std::vector<uint8_t>>> const &get_local_states() const { return this->local_states; }

    // Cumulative time (seconds) in Phase 1 (cumulant eval) and Phase 2 (config sum).
    double get_phase1_time() const { return this->phase1_seconds; }
    double get_phase2_time() const { return this->phase2_seconds; }

    private:
    Graph const &graph;
    std::vector<VertexType<T> *> vertex_type_ptrs;
    std::vector<ValidGlobalConfig> valid_configurations;
    std::vector<SpatialConfiguration> spatial_configurations;
    std::vector<std::vector<LegInfo>> legs_per_vertex;
    Lines hopping_lines;
    int diagram_sign = 1;

    // [config_idx][vertex_idx] — fallback path when local_states isn't built.
    std::vector<std::vector<VertexInstance<T>>> vertex_instances;

    // tau_index → list of vertex indices that depend on it.
    std::vector<std::vector<int>> tau_to_vertices;

    // Factored evaluation tables.
    std::vector<std::vector<std::vector<uint8_t>>> local_states; // [vertex][state_idx] -> op_ids
    std::vector<std::vector<T>> local_values;                    // [vertex][state_idx]
    std::vector<std::vector<CumulantPlan>> local_plans_finite;
    bool local_plans_built = false;
    std::vector<bool> vertex_dirty_finite;
    mutable long local_cache_hits   = 0;
    mutable long local_cache_misses = 0;
    mutable double phase1_seconds   = 0.0;
    mutable double phase2_seconds   = 0.0;
    std::vector<std::vector<int>> config_to_local; // [gc_idx][vertex] -> state_idx

    void setup_vertices(std::vector<VertexType<T> *> const &vertex_types);
    void build_local_plans(HubbardSolver<2, T> const &solver);
    void compute_spatial_configurations();
    void compute_spatial_configurations_cluster(std::vector<std::pair<int, int>> const &cluster_positions, int n_cluster_sites);
    void compute_valid_configurations();
    void build_vertex_instances();
    void build_local_state_tables();
    T evaluate_factored(std::vector<double> const &taus, HubbardSolver<2, T> const &solver);

    void solve_dimer_embedding(int placed_count, std::vector<bool> &placed, std::vector<std::pair<int, int>> &coords,
                               std::map<std::vector<uint8_t>, int> &config_counts) const;

    void solve_cluster_embedding(int placed_count, std::vector<bool> &placed, std::vector<std::pair<int, int>> &coords,
                                 std::map<std::vector<uint8_t>, int> &config_counts,
                                 std::vector<std::pair<int, int>> const &cluster_positions) const;

    std::vector<uint8_t> canonicalize_directions(std::vector<uint8_t> const &dirs, std::vector<std::vector<int>> const &automorphisms) const;
    std::vector<uint8_t> apply_automorphism_to_directions(std::vector<uint8_t> const &dirs, std::vector<int> const &perm) const;
  };

} // namespace sc_expansion::dimer
