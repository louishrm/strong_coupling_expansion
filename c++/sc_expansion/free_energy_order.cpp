#include "free_energy_order.hpp"
#include "generate_diagrams.hpp"
#include <cmath>

namespace sc_expansion {

  FreeEnergyCalculator::FreeEnergyCalculator(Parameters const &params_, int order_) : params(params_), order(order_) {
    VacuumDiagramGenerator gen(this->order);
    gen.generate();
    const auto &unique_graphs = gen.get_unique_graphs();

    for (auto const &g_vec : unique_graphs) {
      int V = std::sqrt(g_vec.size());
      this->diagrams.emplace_back(Graph(g_vec, V));
    }

    for (auto const &diag : this->diagrams) { this->evaluators.emplace_back(diag, this->params); }
  }

  double FreeEnergyCalculator::compute_sum_diagrams(std::vector<double> const &taus, bool infinite_U) const {
    double sum = 0.0;
    for (auto const &evaluator : this->evaluators) { sum += evaluator.evaluate_at_taus(taus, infinite_U); }
    return sum;
  }

} // namespace sc_expansion
