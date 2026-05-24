#pragma once

#include <vector>
#include "../graph.hpp"
#include "../hubbard_solver.hpp"
#include "../diagram_common.hpp"
#include "vertex.hpp"

namespace sc_expansion::atomic {

  // Two equivalent encodings for the density insertion on a marked vertex:
  //
  //   BlockConstraint (default, existing): the mark contributes a real
  //     (c†_σ, c_σ) operator pair pinned into op_ids at τ=0; the cumulant
  //     plan keeps them in one partition block via add_block_constraint.
  //   StaticDensity (new): the mark contributes nothing to op_ids; the
  //     density is attached to the per-vertex cumulant leaf as an external
  //     decoration via add_static_density. The (c†, c) pair never enters
  //     the partition lattice; each n_σ(0) is one bosonic atom routed
  //     independently across blocks (see static_density_correlator.md).
  //
  // The two paths produce identical diagram values for the same (graph,
  // marks, mark_spins) — cross-validated in the test suite.
  enum class MarkEncoding { BlockConstraint, StaticDensity };

  template <typename T> class Diagram {
    public:
    // Free-energy (vacuum) constructor.
    explicit Diagram(Graph const &graph, std::vector<VertexType<T> *> const &vertex_types);

    // Rooted (correlator) constructor. `marks` are the marked-vertex indices
    // in `graph`'s labeling (the rooted graph's canonical labeling, by
    // convention); `mark_spins` carries one entry per mark with σ∈{0=↓,1=↑}.
    // Caller must pre-set graph's symmetry_factor to the rooted symmetry
    // factor and free_multiplicity = 1 via the override Graph constructor;
    // the shell-indexed lattice multiplier is the measurement layer's job.
    // `evaluate(taus)` expects taus of size n_lines+1 with taus.back() = 0;
    // for MarkEncoding::StaticDensity the pinned slot is unused but still
    // expected (keeps the callers' tau-layout convention uniform).
    // `flip_mark_order` is honored only for BlockConstraint encoding; it is
    // a no-op for StaticDensity (no mark pair to flip).
    Diagram(Graph const &graph, std::vector<VertexType<T> *> const &vertex_types,
            std::vector<int> marks, std::vector<int> mark_spins,
            bool flip_mark_order = false,
            MarkEncoding mark_encoding = MarkEncoding::BlockConstraint);

    T evaluate(std::vector<double> const &taus, HubbardSolver<1, T> const &solver, bool infinite_U);

    // Diagnostic: per-config (sign × weight × Π_v vertex_value) for the rooted path.
    // Returns one entry per valid configuration, same order as get_valid_configurations().
    std::vector<T> evaluate_per_config(std::vector<double> const &taus, HubbardSolver<1, T> const &solver, bool infinite_U);

    void mark_tau_dirty(int tau_index);
    void mark_all_dirty();

    double get_free_multiplicity() const { return this->graph.get_free_multiplicity(); }
    std::vector<ValidGlobalConfig> const &get_valid_configurations() const { return this->valid_configurations; }
    double get_diagram_sign() const { return (double)this->diagram_sign; }
    Graph const &get_graph() const { return this->graph; }
    bool is_rooted_diagram() const { return this->is_rooted; }
    std::vector<int> const &get_marks() const { return this->marks; }
    std::vector<int> const &get_mark_spins() const { return this->mark_spins; }
    int get_pinned_tau_index() const { return (int)this->hopping_lines.lines.size(); }

    private:
    Graph const &graph;
    std::vector<VertexType<T> *> vertex_type_ptrs;
    std::vector<ValidGlobalConfig> valid_configurations;
    std::vector<std::vector<LegInfo>> legs_per_vertex;
    Lines hopping_lines;
    int diagram_sign = 1;

    bool is_rooted = false;
    bool flip_mark_order = false; // if true, push (c_σ, c†_σ) instead of (c†_σ, c_σ)
    MarkEncoding mark_encoding = MarkEncoding::BlockConstraint;
    std::vector<int> marks;       // vertex indices of marked vertices
    std::vector<int> mark_spins;  // σ ∈ {0=↓, 1=↑} for each mark

    // [config_idx][vertex_idx] — local cache per (config, vertex) pair
    std::vector<std::vector<VertexInstance<T>>> vertex_instances;

    // tau_index → vertex indices that depend on it.
    std::vector<std::vector<int>> tau_to_vertices;

    void setup_vertices(std::vector<VertexType<T> *> const &vertex_types);
    void compute_valid_configurations();
    void build_vertex_instances();
  };

} // namespace sc_expansion::atomic
