#pragma once

#include "../configuration_base.hpp"
#include "free_energy_calculator.hpp"

namespace sc_expansion::dimer {

  // Two-site MCMC configuration. Reference weight is uniform: weight = |omega| + alpha,
  // and the reference integrand exposed to the measure is the constant `alpha`.
  // Contrast with sc_expansion::atomic::Configuration, which uses the infinite-U
  // integrand as the reference.
  template <typename T> class Configuration : public ConfigurationBase<T> {
    public:
    Configuration(Parameters<T> const &params, int order, FreeEnergyCalculator<T> &calculator, double alpha = 0.001);

    double evaluate_proposed() override;
    void commit_proposal() override;

    double get_integrand() const override;
    double get_reference_integrand() const override;

    void mark_tau_dirty(int tau_index) override;

    double get_alpha() const { return this->alpha; }
    FreeEnergyCalculator<T> const &get_calculator() const { return this->calculator; }

    private:
    double omega;
    double proposed_omega;
    double alpha;
    FreeEnergyCalculator<T> &calculator;

    double compute_omega() const;
  };

} // namespace sc_expansion::dimer
