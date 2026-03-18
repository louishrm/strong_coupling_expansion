#pragma once
#include <vector>
#include <deque>
#include "diagram.hpp"
#include "hubbard_solver.hpp"

namespace sc_expansion {

  template <typename T>
  class FreeEnergyCalculator {
    public:
    FreeEnergyCalculator(Parameters<T> const &params, int order);

    T compute_sum_diagrams(std::vector<double> const &taus, bool infinite_U, bool use_cache) const;
    T compute_sum_diagrams_dimer(std::vector<double> const &taus, bool infinite_U, bool use_cache) const;

    private:
    Parameters<T> const &params;
    int order;
    std::deque<Diagram> diagrams;
    std::vector<DiagramEvaluator<T>> evaluators;
  };

} // namespace sc_expansion
