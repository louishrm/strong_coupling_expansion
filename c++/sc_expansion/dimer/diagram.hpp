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

  // Density-insertion encoding for a marked (rooted) vertex. The dimer
  // correlator uses only StaticDensity: the n_σ(0) density is attached to the
  // per-vertex cumulant leaf as an external decoration via add_static_density,
  // contributing nothing to op_ids and never entering the partition lattice
  // (see dimer_density_density_correlator-01.md, task 1). Unlike atomic there
  // is no BlockConstraint path — the dimer never folds a real (c†, c) pair into
  // op_ids — so the enum carries a single value, kept as a named type for
  // signature parity with atomic::MarkEncoding.
  enum class MarkEncoding { StaticDensity };

  template <typename T> class Diagram {
    public:
    explicit Diagram(Graph const &graph, std::vector<VertexType<T> *> const &vertex_types);

    // Cluster-restricted embedding: spatial configs computed from only the given
    // positions. Weights are divided by n_cluster_sites for per-dimer values.
    Diagram(Graph const &graph, std::vector<VertexType<T> *> const &vertex_types,
            std::vector<std::pair<int, int>> const &cluster_positions, int n_cluster_sites);

    // Rooted (density-density correlator) constructor. `marks` are the marked-
    // vertex indices in `graph`'s (canonical) labeling; `sites` is the within-
    // dimer site (0 or 1) per mark; `mark_spins` carries σ ∈ {0=↓, 1=↑} per
    // mark; `r` is the target physical displacement (x, y). For coincident
    // marks both entries of `marks` are the same vertex index.
    //
    // The mark-constrained embedding (compute_spatial_configurations_rooted) pins
    // the two marks at physical (0,0)/r on the infinite staggered superlattice and
    // sums both anchorings via the ±r enumeration; the marked vertices' cumulants
    // carry a static density n_σ(0) via add_static_density, decorated in
    // build_local_plans. is_rooted is set true.
    Diagram(Graph const &graph, std::vector<VertexType<T> *> const &vertex_types, std::vector<int> marks, std::vector<int> sites,
            std::vector<int> mark_spins, std::vector<int> r, MarkEncoding mark_encoding = MarkEncoding::StaticDensity);

    // Cluster-restricted rooted constructor (finite-cluster density-density). Same
    // mark pinning and density decoration as the rooted ctor above, but every dimer
    // — mark1's forced dimer and all interior vertices — must land on a
    // `cluster_positions` superlattice cell.
    //
    // pin_origin selects the normalisation convention:
    //   false (default): mark0's home dimer is SWEPT over all cluster positions and
    //     the summed embedding weight divided by n_cluster_sites — the per-dimer
    //     translation average, mirroring the EXTENSIVE vacuum free-energy cluster
    //     ctor. Correct for per-dimer free-energy-like quantities, but on an
    //     inhomogeneous open cluster it averages over INEQUIVALENT reference sites.
    //   true: mark0's home dimer is PINNED at cluster_positions[0] (no sweep, no
    //     ÷n_cluster_sites). This is the INTENSIVE/local convention for a
    //     ⟨n(r)n(0)⟩ correlator anchored at one reference site — it reproduces a
    //     finite-cluster ED measurement at that single site (e.g. the pendant
    //     site (0,0)) instead of smearing over all cells.
    Diagram(Graph const &graph, std::vector<VertexType<T> *> const &vertex_types, std::vector<int> marks, std::vector<int> sites,
            std::vector<int> mark_spins, std::vector<int> r, std::vector<std::pair<int, int>> const &cluster_positions, int n_cluster_sites,
            MarkEncoding mark_encoding = MarkEncoding::StaticDensity, bool pin_origin = false);

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

    // Rooted (density-density) accessors. is_rooted_diagram() is false for the
    // vacuum/cluster ctors; in rooted mode the mark/site/spin state and target
    // displacement are exposed for diagnostics and the task-3 embedding.
    bool is_rooted_diagram() const { return this->is_rooted; }
    std::vector<int> const &get_marks() const { return this->marks; }
    std::vector<int> const &get_sites() const { return this->sites; }
    std::vector<int> const &get_mark_spins() const { return this->mark_spins; }
    std::vector<int> const &get_target_r() const { return this->target_r; }
    // Diagnostic: rooted weight normaliser (|mark-fixing automorphisms| ×
    // multi-edge factorial), used in place of the vacuum symmetry factor by the
    // rooted branch of compute_valid_configurations. Exposed for the on-site
    // same-spin weight-bookkeeping cross-check.
    double get_rooted_sym_factor() const { return this->rooted_sym_factor; }
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

    // Rooted (density-density) state. is_rooted == false ⇒ vacuum/cluster mode.
    bool is_rooted             = false;
    MarkEncoding mark_encoding = MarkEncoding::StaticDensity;
    std::vector<int> marks;      // marked-vertex indices (canonical labeling)
    std::vector<int> sites;      // within-dimer site (0/1) per mark
    std::vector<int> mark_spins; // σ ∈ {0=↓, 1=↑} per mark
    std::vector<int> target_r;   // target physical displacement (x, y)

    // Rooted weight normalisation: |pointwise mark-stabiliser automorphisms| ×
    // multi-edge factorial. Computed in compute_spatial_configurations_rooted
    // and used (in place of the vacuum graph.get_symmetry_factor()) by the
    // rooted branch of compute_valid_configurations.
    double rooted_sym_factor = 1.0;

    // True iff some graph automorphism swaps the two marks (perm[m0]=m1,
    // perm[m1]=m0). Set alongside rooted_sym_factor in compute_rooted_automorphisms.
    // Drives the r=0 label-swap multiplicity (see compute_valid_configurations):
    // a distinct-mark graph WITHOUT such an automorphism represents two physically
    // identical but distinctly-labelled mark assignments, only one of which the
    // mark-fixing canonicalisation keeps, so its weight is doubled at r=0.
    bool mark_swap_exists = false;

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
    void compute_spatial_configurations_rooted();
    void compute_spatial_configurations_rooted_cluster(std::vector<std::pair<int, int>> const &cluster_positions, int n_cluster_sites,
                                                        bool pin_origin);

    // Graph automorphisms that fix every marked vertex pointwise; also sets
    // rooted_sym_factor (= |these automorphisms| × multi-edge factorial). Shared
    // by the infinite and cluster-restricted rooted spatial enumerations.
    std::vector<std::vector<int>> compute_rooted_automorphisms();

    void compute_valid_configurations();
    void build_vertex_instances();
    void build_local_state_tables();
    T evaluate_factored(std::vector<double> const &taus, HubbardSolver<2, T> const &solver);

    void solve_dimer_embedding(int placed_count, std::vector<bool> &placed, std::vector<std::pair<int, int>> &coords,
                               std::map<std::vector<uint8_t>, int> &config_counts) const;

    void solve_cluster_embedding(int placed_count, std::vector<bool> &placed, std::vector<std::pair<int, int>> &coords,
                                 std::map<std::vector<uint8_t>, int> &config_counts,
                                 std::vector<std::pair<int, int>> const &cluster_positions) const;

    std::vector<uint8_t> canonicalize_directions(std::vector<uint8_t> const &dirs, std::vector<std::vector<int>> const &automorphisms,
                                                 bool include_inversion = true) const;
    std::vector<uint8_t> apply_automorphism_to_directions(std::vector<uint8_t> const &dirs, std::vector<int> const &perm) const;
  };

} // namespace sc_expansion::dimer
