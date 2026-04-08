#pragma once
#include <vector>
#include <deque>
#include "diagram.hpp"
#include "hubbard_solver.hpp"

namespace sc_expansion {

  template <int N_sites, typename T>
  class FreeEnergyCalculator {
    public:
    FreeEnergyCalculator(Parameters<T> const &params, int order, int override_fm = -1);

    // Construct from pre-built graphs (avoids redundant diagram generation across MPI ranks)
    FreeEnergyCalculator(Parameters<T> const &params, int order, std::vector<Graph> const &prebuilt_graphs);

    // Cluster-restricted embedding: spatial configs computed from the given cluster positions.
    FreeEnergyCalculator(Parameters<T> const &params, int order,
                         std::vector<std::pair<int, int>> const &cluster_positions, int n_cluster_sites);

    T compute_sum_diagrams(std::vector<double> const &taus, bool infinite_U) const;
    T compute_sum_diagrams_dimer(std::vector<double> const &taus, bool infinite_U) const;

    std::pair<double, double> compute_infinite_U_coefficient(bool dimer = false) const;

    void clear_all_caches() const;
    void mark_tau_dirty(int tau_index);
    void mark_all_dirty();

    std::deque<Graph> const &get_graphs() const { return this->graphs; }
    std::deque<VertexType<N_sites, T>> const &get_vertex_types() const { return this->vertex_types; }

    private:
    Parameters<T> const &params;
    int order;
    std::deque<Graph> graphs;
    std::deque<Diagram<N_sites, T>> diagrams;
    std::deque<VertexType<N_sites, T>> vertex_types; // One per cumulant order (index = cumulant_order - 1)
    HubbardSolver<N_sites, T> solver;
  };

} // namespace sc_expansion
