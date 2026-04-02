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

    // Cluster-restricted embedding: spatial configs computed from the given cluster positions.
    FreeEnergyCalculator(Parameters<T> const &params, int order,
                         std::vector<std::pair<int, int>> const &cluster_positions, int n_cluster_sites);

    T compute_sum_diagrams(std::vector<double> const &taus, bool infinite_U, bool use_cache) const;
    T compute_sum_diagrams_dimer(std::vector<double> const &taus, bool infinite_U, bool use_cache) const;

    std::pair<double, double> compute_infinite_U_coefficient(bool dimer = false) const;

    private:
    Parameters<T> const &params;
    int order;
    std::deque<Graph> graphs;
    std::deque<Diagram<N_sites, T>> diagrams;
    HubbardSolver<N_sites, T> solver;
  };

} // namespace sc_expansion
