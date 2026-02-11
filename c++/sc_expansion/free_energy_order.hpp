#pragma once
#include <vector>
#include <deque>
#include "diagram.hpp"
#include "hubbard_atom.hpp"

namespace sc_expansion {

  class FreeEnergyCalculator {
    public:
    FreeEnergyCalculator(Parameters const &params, int order);

    double compute_sum_diagrams(std::vector<double> const &taus, bool infinite_U) const;

    private:
    Parameters params;
    int order;
    std::deque<Diagram> diagrams;
    std::vector<DiagramEvaluator> evaluators;
  };

} // namespace sc_expansion
