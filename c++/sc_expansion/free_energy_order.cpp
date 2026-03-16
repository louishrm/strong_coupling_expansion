#include "free_energy_order.hpp"
#include "generate_diagrams.hpp"
#include <cmath>

#include "dual.hpp"

namespace sc_expansion {

  template <typename T>
  FreeEnergyCalculator<T>::FreeEnergyCalculator(Parameters<T> const &params_, int order_) : params(params_), order(order_) {
    VacuumDiagramGenerator gen(this->order, params_.bipartite);
    gen.generate();
    const auto &unique_graphs = gen.get_unique_graphs();

    for (auto const &g : unique_graphs) { this->diagrams.emplace_back(g); }

    for (auto const &diag : this->diagrams) { this->evaluators.emplace_back(diag, this->params); }
  }

  template <typename T>
  T FreeEnergyCalculator<T>::compute_sum_diagrams(std::vector<double> const &taus, bool infinite_U, bool use_cache) const {
    T sum = T(0.0);
    for (auto const &evaluator : this->evaluators) { sum = sum + evaluator.evaluate_at_taus(taus, infinite_U, use_cache); }
    return sum;
  }

  template <typename T>
  T FreeEnergyCalculator<T>::compute_sum_diagrams_dimer(std::vector<double> const &taus, bool infinite_U, bool use_cache) const {
    T sum = T(0.0);
    for (auto const &evaluator : this->evaluators) { sum = sum + evaluator.evaluate_at_taus_dimer(taus, infinite_U, use_cache); }
    return sum;
  }

  template class FreeEnergyCalculator<double>;
  template class FreeEnergyCalculator<Dual>;

} // namespace sc_expansion
