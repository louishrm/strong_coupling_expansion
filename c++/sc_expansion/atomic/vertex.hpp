#pragma once

#include <vector>
#include "../args.hpp"
#include "../cumulant.hpp"
#include "../hubbard_solver.hpp"

namespace sc_expansion::atomic {

  // Single-site vertex. Templated on scalar type only; the Hilbert size is
  // fixed to 1 inside the namespace. Mirrors the original sc_expansion::VertexType
  // / VertexInstance interface so call sites only need a namespace qualifier.
  template <typename T> class VertexType {
    public:
    explicit VertexType(int n_legs);

    std::pair<long, long> get_local_cache_stats() const { return {this->local_cache_hits, this->local_cache_misses}; }

    void record_local_hit() const { ++this->local_cache_hits; }
    void record_local_miss() const { ++this->local_cache_misses; }

    private:
    int n_legs;
    mutable long local_cache_hits   = 0;
    mutable long local_cache_misses = 0;
  };

  template <typename T> class VertexInstance {
    public:
    VertexInstance(VertexType<T> *type_, std::vector<int> tau_indices_, std::vector<uint8_t> op_ids_);

    T get_value(std::vector<double> const &global_taus, HubbardSolver<1, T> const &solver, bool infinite_U) const;

    void mark_dirty() {
      this->is_dirty_finite   = true;
      this->is_dirty_infinite = true;
    }

    // Mark this vertex as carrying a density insertion n_σ(0) = c†_σ·c_σ
    // — the cumulant plan keeps the (c†, c) pair in the same partition
    // block, producing the partial cumulant κ_partial rather than the full
    // κ. Call once per mark on this vertex (a vertex with two coincident
    // marks gets two block constraints — one per density insertion).
    // u_input_idx / p_input_idx are positions in the unprimed / primed
    // operator lists *as produced by Args::split_from_raw on op_ids*.
    void add_block_constraint(int u_input_idx, int p_input_idx) {
      this->block_u_inputs.push_back(u_input_idx);
      this->block_p_inputs.push_back(p_input_idx);
    }

    // Declare a coincidence group of marks sitting at the same (τ=0, σ=
    // orbital). block_indices reference the order of prior
    // add_block_constraint calls on this vertex.
    void add_coincidence_group(int orbital, std::vector<int> block_indices) {
      this->coincidence_orbitals.push_back(orbital);
      this->coincidence_block_indices.push_back(std::move(block_indices));
    }

    // Register one static density n_σ(0) decoration on this vertex (for the
    // static-density correlator path). Unlike add_block_constraint, this
    // mechanism does NOT require the mark's (c†_σ, c_σ) pair to live in
    // op_ids — the density is attached as an external decoration of the
    // per-vertex cumulant leaf via CumulantSolver::add_static_density.
    // The two mechanisms are mutually exclusive per mark in practice.
    void add_static_density(int orbital) { this->static_density_orbitals.push_back(orbital); }

    std::vector<int> const &get_tau_indices() const { return this->tau_indices; }

    private:
    VertexType<T> *type;
    std::vector<int> tau_indices;
    std::vector<uint8_t> op_ids;
    std::vector<int> block_u_inputs;
    std::vector<int> block_p_inputs;
    std::vector<int> coincidence_orbitals;
    std::vector<std::vector<int>> coincidence_block_indices;
    std::vector<int> static_density_orbitals;

    mutable T local_cache_finite;
    mutable T local_cache_infinite;
    mutable bool is_dirty_finite   = true;
    mutable bool is_dirty_infinite = true;

    mutable CumulantPlan plan;
    mutable bool plan_built = false;
  };

} // namespace sc_expansion::atomic
