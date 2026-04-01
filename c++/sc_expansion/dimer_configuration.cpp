#include "dimer_configuration.hpp"
#include "dual.hpp"
#include <triqs/mc_tools/random_generator.hpp>

template <typename T>
DimerConfiguration<T>::DimerConfiguration(sc_expansion::Parameters<T> const &params_, int order_)
   : ConfigurationBase<T>(params_, order_), omega(0.0), proposed_omega(0.0), calculator(params_, order_) {

  this->state.resize(this->order);
  triqs::mc_tools::random_generator RNG("mt19937", 23432);
  for (int i = 0; i < this->order; i++) { this->state[i] = RNG(this->beta); }

  this->omega             = this->compute_omega();
  this->metropolis_weight = std::abs(this->omega);
}

template <typename T> double DimerConfiguration<T>::compute_omega() const {
  T val = this->calculator.compute_sum_diagrams_dimer(this->state, false, true);
  if constexpr (std::is_same_v<T, Dual>) {
    return val.derivative;
  } else {
    return val;
  }
}

template <typename T> double DimerConfiguration<T>::evaluate_proposed() {
  this->proposed_omega = this->compute_omega();
  return std::abs(this->proposed_omega);
}

template <typename T> void DimerConfiguration<T>::commit_proposal() {
  this->omega             = this->proposed_omega;
  this->metropolis_weight = std::abs(this->omega);
}

template <typename T> double DimerConfiguration<T>::get_integrand() const { return this->omega; }

template <typename T> double DimerConfiguration<T>::get_reference_integrand() const { return 0.0; }

template class DimerConfiguration<double>;
template class DimerConfiguration<Dual>;
