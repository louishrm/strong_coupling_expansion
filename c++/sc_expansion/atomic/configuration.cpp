#include "configuration.hpp"
#include "../dual.hpp"
#include <triqs/mc_tools/random_generator.hpp>

namespace sc_expansion::atomic {

  template <typename T>
  Configuration<T>::Configuration(Parameters<T> const &params_, int order_, double alpha_, FreeEnergyCalculator<T> &calculator_)
     : ConfigurationBase<T>(params_, order_), alpha(alpha_), integrand(0.0), reference_integrand(0.0), proposed_integrand(0.0),
       proposed_reference_integrand(0.0), calculator(calculator_) {

    this->state.resize(this->order);
    triqs::mc_tools::random_generator RNG("mt19937", 23432);
    for (int i = 0; i < this->order; ++i) this->state[i] = RNG(this->beta);

    auto [finite_U, infinite_U] = this->compute_integrands();
    this->integrand              = finite_U;
    this->reference_integrand    = infinite_U;
    this->metropolis_weight      = this->compute_weight(finite_U, infinite_U);
  }

  template <typename T> std::pair<double, double> Configuration<T>::compute_integrands() const {
    double finite_U   = 0.0;
    double infinite_U = 0.0;

    T finite_U_T   = this->calculator.compute_sum_diagrams(this->state, false);
    T infinite_U_T = this->calculator.compute_sum_diagrams(this->state, true);

    if constexpr (std::is_same_v<T, Dual>) {
      finite_U   = finite_U_T.derivative;
      infinite_U = infinite_U_T.derivative;
    } else {
      finite_U   = finite_U_T;
      infinite_U = infinite_U_T;
    }
    return {finite_U, infinite_U};
  }

  template <typename T> double Configuration<T>::compute_weight(double finite_U, double infinite_U) const {
    return std::abs(finite_U - infinite_U) + this->alpha * std::abs(infinite_U);
  }

  template <typename T> double Configuration<T>::evaluate_proposed() {
    auto [finite_U, infinite_U]        = this->compute_integrands();
    this->proposed_integrand           = finite_U;
    this->proposed_reference_integrand = infinite_U;
    return this->compute_weight(finite_U, infinite_U);
  }

  template <typename T> void Configuration<T>::commit_proposal() {
    this->integrand           = this->proposed_integrand;
    this->reference_integrand = this->proposed_reference_integrand;
    this->metropolis_weight   = this->compute_weight(this->integrand, this->reference_integrand);
  }

  template <typename T> double Configuration<T>::get_integrand() const { return this->integrand; }

  template <typename T> double Configuration<T>::get_reference_integrand() const { return this->reference_integrand; }

  template <typename T> void Configuration<T>::mark_tau_dirty(int tau_index) { this->calculator.mark_tau_dirty(tau_index); }

  template class Configuration<double>;
  template class Configuration<Dual>;

} // namespace sc_expansion::atomic
