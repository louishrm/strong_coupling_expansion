#include "free_energy_calculator.hpp"
#include "generate_diagrams.hpp"
#include "combinatorics.hpp"
#include <cmath>
#include <numeric>
#include <algorithm>

#include "dual.hpp"

namespace sc_expansion {

  template <int N_sites, typename T>
  FreeEnergyCalculator<N_sites, T>::FreeEnergyCalculator(Parameters<T> const &params_, int order_) : params(params_), order(order_) {
    VacuumDiagramGenerator gen(this->order, params.bipartite);
    gen.generate();
    const auto &unique_graphs = gen.get_unique_graphs();

    for (auto const &g : unique_graphs) { this->diagrams.emplace_back(g); }

    for (auto const &diag : this->diagrams) { this->evaluators.emplace_back(diag, this->params); }
  }

  template <int N_sites, typename T>
  T FreeEnergyCalculator<N_sites, T>::compute_sum_diagrams(std::vector<double> const &taus, bool infinite_U, bool use_cache) const {
    T sum = T(0.0);
    for (auto const &evaluator : this->evaluators) { sum = sum + evaluator.evaluate_at_taus(taus, infinite_U, use_cache); }
    return sum;
  }

  template <int N_sites, typename T>
  T FreeEnergyCalculator<N_sites, T>::compute_sum_diagrams_dimer(std::vector<double> const &taus, bool infinite_U, bool use_cache) const {
    T sum = T(0.0);
    for (auto const &evaluator : this->evaluators) { sum = sum + evaluator.evaluate_at_taus_dimer(taus, infinite_U, use_cache); }
    return sum;
  }

  template <int N_sites, typename T>
  std::pair<double, double> FreeEnergyCalculator<N_sites, T>::compute_infinite_U_coefficient(bool dimer) const {
    int n = this->order;
    std::vector<double> taus(n);
    std::iota(taus.begin(), taus.end(), 0.0);

    double sum_abs    = 0.0;
    double sum_signed = 0.0;
    do {
      T val_T = dimer ? this->compute_sum_diagrams_dimer(taus, true, false) : this->compute_sum_diagrams(taus, true, false);
      double val;
      if constexpr (std::is_same_v<T, Dual>) {
        val = val_T.derivative;
      } else {
        val = (double)val_T;
      }
      sum_abs += std::abs(val);
      sum_signed += val;
    } while (std::next_permutation(taus.begin(), taus.end()));

    double fact = 1.0;
    for (int i = 1; i <= n; ++i) fact *= i;

    double beta_val;
    if constexpr (std::is_same_v<T, Dual>) {
      beta_val = this->params.beta.value;
    } else {
      beta_val = (double)this->params.beta;
    }

    std::pair<double, double> result;
    result.first  = (std::pow(beta_val, n) / fact) * sum_abs;
    result.second = (std::pow(beta_val, n) / fact) * sum_signed;
    return result;
  }

  template class FreeEnergyCalculator<1, double>;
  template class FreeEnergyCalculator<1, Dual>;
  template class FreeEnergyCalculator<2, double>;
  template class FreeEnergyCalculator<2, Dual>;

} // namespace sc_expansion
