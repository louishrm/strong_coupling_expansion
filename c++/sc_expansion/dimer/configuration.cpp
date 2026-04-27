#include "configuration.hpp"
#include "../dual.hpp"
#include <triqs/mc_tools/random_generator.hpp>
#include <cmath>

namespace sc_expansion::dimer {

  template <typename T>
  Configuration<T>::Configuration(Parameters<T> const &params_, int order_, FreeEnergyCalculator<T> &calculator_, double alpha_)
     : ConfigurationBase<T>(params_, order_), omega(0.0), proposed_omega(0.0), alpha(alpha_), calculator(calculator_) {

    this->state.resize(this->order);
    triqs::mc_tools::random_generator RNG("mt19937", 23432);
    for (int i = 0; i < this->order; ++i) this->state[i] = RNG(this->beta);

    this->omega             = this->compute_omega();
    this->metropolis_weight = std::abs(this->omega + this->alpha);
  }

  template <typename T> double Configuration<T>::compute_omega() const {
    T val = this->calculator.compute_sum_diagrams(this->state);
    this->calculator.clear_all_caches();
    if constexpr (std::is_same_v<T, Dual>) {
      return val.derivative;
    } else {
      return val;
    }
  }

  template <typename T> double Configuration<T>::evaluate_proposed() {
    this->proposed_omega = this->compute_omega();
    return std::abs(this->proposed_omega + this->alpha);
  }

  template <typename T> void Configuration<T>::commit_proposal() {
    this->omega             = this->proposed_omega;
    this->metropolis_weight = std::abs(this->omega + this->alpha);
  }

  template <typename T> double Configuration<T>::get_integrand() const { return this->omega; }

  template <typename T> double Configuration<T>::get_reference_integrand() const { return 0.0; }

  template <typename T> void Configuration<T>::mark_tau_dirty(int tau_index) { this->calculator.mark_tau_dirty(tau_index); }

  template class Configuration<double>;
  template class Configuration<Dual>;

} // namespace sc_expansion::dimer
