#pragma once
#include <vector>
#include <deque>
#include "../graph.hpp"
#include "../hubbard_solver.hpp"
#include "diagram.hpp"
#include "vertex.hpp"

namespace sc_expansion::dimer {

  template <typename T> class FreeEnergyCalculator {
    public:
    FreeEnergyCalculator(Parameters<T> const &params, int order, int override_fm = -1);

    // Construct from pre-built graphs (avoids redundant diagram generation across MPI ranks).
    FreeEnergyCalculator(Parameters<T> const &params, int order, std::vector<Graph> const &prebuilt_graphs, int override_fm = -1);

    // Construct with finite-cluster embedding: spatial configurations are restricted
    // to embeddings that fit on the given cluster (in superlattice (u,v) coordinates).
    // Per-config weights are normalised by n_cluster_sites.
    FreeEnergyCalculator(Parameters<T> const &params, int order, std::vector<std::pair<int, int>> const &cluster_positions, int n_cluster_sites);

    T compute_sum_diagrams(std::vector<double> const &taus) const;

    void mark_tau_dirty(int tau_index);
    void mark_all_dirty();

    // Drop every cached value so the next compute_sum_diagrams call recomputes
    // from scratch. Used by Configuration::compute_omega to guarantee a clean
    // state between MC steps independent of the move's per-tau dirty marking.
    void clear_all_caches() { this->mark_all_dirty(); }

    std::deque<Graph> const &get_graphs() const { return this->graphs; }
    std::deque<VertexType<T>> const &get_vertex_types() const { return this->vertex_types; }
    std::deque<Diagram<T>> const &get_diagrams() const { return this->diagrams; }
    HubbardSolver<2, T> const &get_solver() const { return this->solver; }

    int get_n_diagrams() const { return (int)this->diagrams.size(); }

    private:
    void init_from_graphs(std::vector<Graph> const &source_graphs, int override_fm);

    Parameters<T> const &params;
    int order;
    std::deque<Graph> graphs;
    std::deque<Diagram<T>> diagrams;
    std::deque<VertexType<T>> vertex_types; // index = cumulant_order - 1
    HubbardSolver<2, T> solver;
  };

} // namespace sc_expansion::dimer
