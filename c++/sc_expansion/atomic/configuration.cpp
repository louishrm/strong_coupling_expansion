#include "configuration.hpp"
#include "../dual.hpp"
#include <triqs/mc_tools/random_generator.hpp>
#include <cmath>

namespace sc_expansion::atomic {

  template <typename T>
  Configuration<T>::Configuration(Parameters<T> const &params_, int order_, double alpha_, SumDiagrams<T> &calculator_)
     : ConfigurationBase<T>(params_, order_), alpha(alpha_), integrand(0.0), proposed_integrand(0.0), calculator(calculator_),
       target_d_sq(calculator_.is_density_density_mode() ? calculator_.get_target_d_sq() : -1) {

    this->state.resize(this->order);
    triqs::mc_tools::random_generator RNG("mt19937", 23432);
    for (int i = 0; i < this->order; ++i) this->state[i] = RNG(this->beta);

    this->integrand         = this->compute_integrand();
    this->metropolis_weight = std::abs(this->integrand + this->alpha);
  }

  template <typename T> double Configuration<T>::compute_integrand() const {
    T val_T;
    if (this->target_d_sq >= 0) {
      val_T = this->calculator.density_density_single(this->state, false);
    } else {
      val_T = this->calculator.free_energy(this->state, false);
    }

    if constexpr (std::is_same_v<T, Dual>) {
      return val_T.derivative;
    } else {
      return val_T;
    }
  }

  template <typename T> double Configuration<T>::evaluate_proposed() {
    this->proposed_integrand = this->compute_integrand();
    return std::abs(this->proposed_integrand + this->alpha);
  }

  template <typename T> void Configuration<T>::commit_proposal() {
    this->integrand         = this->proposed_integrand;
    this->metropolis_weight = std::abs(this->integrand + this->alpha);
  }

  template <typename T> double Configuration<T>::get_integrand() const { return this->integrand; }

  template <typename T> double Configuration<T>::get_reference_integrand() const { return 0.0; }

  template <typename T> void Configuration<T>::mark_tau_dirty(int tau_index) { this->calculator.mark_tau_dirty(tau_index); }

  template class Configuration<double>;
  template class Configuration<Dual>;

} // namespace sc_expansion::atomic
