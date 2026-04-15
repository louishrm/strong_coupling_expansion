#include "diagmc_configuration.hpp"
#include "configuration_base.hpp"
#include "dual.hpp"
#include <triqs/mc_tools/random_generator.hpp>
#include <cmath>

template <typename T>
DiagMCConfiguration<T>::DiagMCConfiguration(sc_expansion::Parameters<T> const &params_, int order_,
                                             sc_expansion::FreeEnergyCalculator<2, T> &calculator_, double alpha_)
   : params(params_), beta(sc_expansion::get_val(params_.beta)), order(order_), n_diagrams(calculator_.get_n_diagrams()), alpha(alpha_),
     defensive_per_diagram(alpha_ / calculator_.get_n_diagrams()), calculator(calculator_) {

  // Initialize tau values randomly in [0, beta)
  this->state.resize(this->order);
  triqs::mc_tools::random_generator RNG("mt19937", 23432);
  for (int i = 0; i < this->order; i++) { this->state[i] = RNG(this->beta); }

  // Start on diagram 1, evaluate it to get the initial weight
  this->current_diagram = 1;
  this->mark_diagram_all_dirty(1);
  double val            = this->evaluate_diagram(1);
  this->current_signed  = val;
  this->current_abs     = std::abs(val) + this->defensive_per_diagram;
}

template <typename T> double DiagMCConfiguration<T>::evaluate_diagram(int d) {
  T val = this->calculator.evaluate_single_diagram(d - 1, this->state, false);
  this->calculator.clear_all_caches();
  if constexpr (std::is_same_v<T, Dual>) {
    return val.derivative;
  } else {
    return val;
  }
}

template <typename T> void DiagMCConfiguration<T>::mark_diagram_tau_dirty(int d, int tau_index) {
  this->calculator.mark_single_diagram_dirty(d - 1, tau_index);
}

template <typename T> void DiagMCConfiguration<T>::mark_diagram_all_dirty(int d) {
  this->calculator.mark_single_diagram_all_dirty(d - 1);
}

template <typename T> double DiagMCConfiguration<T>::get_alpha() const { return this->alpha; }

template <typename T> double DiagMCConfiguration<T>::get_defensive_per_diagram() const { return this->defensive_per_diagram; }

template <typename T> int DiagMCConfiguration<T>::get_n_diagrams() const { return this->n_diagrams; }

template <typename T> int DiagMCConfiguration<T>::get_order() const { return this->order; }

template <typename T> double DiagMCConfiguration<T>::get_U() const { return sc_expansion::get_val(this->params.U); }

template <typename T>
sc_expansion::FreeEnergyCalculator<2, T> const &DiagMCConfiguration<T>::get_calculator() const {
  return this->calculator;
}

template class DiagMCConfiguration<double>;
template class DiagMCConfiguration<Dual>;
