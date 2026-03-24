#pragma once
#include <vector>
#include <deque>
#include "diagram.hpp"
#include "hubbard_solver.hpp"

namespace sc_expansion {

  template <int N_sites, typename T>
  class FreeEnergyCalculator {
    public:
    FreeEnergyCalculator(Parameters<T> const &params, int order);

    T compute_sum_diagrams(std::vector<double> const &taus, bool infinite_U, bool use_cache) const;
    T compute_sum_diagrams_dimer(std::vector<double> const &taus, bool infinite_U, bool use_cache) const;

    std::pair<double, double> compute_infinite_U_coefficient(bool dimer = false) const;

    private:
    Parameters<T> const &params;
    int order;
    std::deque<Diagram> diagrams;
    std::vector<DiagramEvaluator<N_sites, T>> evaluators;
  };

} // namespace sc_expansion
