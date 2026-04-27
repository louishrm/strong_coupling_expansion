#pragma once

#include "../configuration_base.hpp"
#include "free_energy_calculator.hpp"

namespace sc_expansion::atomic {

  // Single-site MCMC configuration. Reference weight is the infinite-U integrand
  // (defensive importance sampling). The dimer namespace defines its own
  // Configuration with a uniform reference weight.
  template <typename T> class Configuration : public ConfigurationBase<T> {
    public:
    Configuration(Parameters<T> const &params, int order, double alpha, FreeEnergyCalculator<T> &calculator);

    double evaluate_proposed() override;
    void commit_proposal() override;

    double get_integrand() const override;
    double get_reference_integrand() const override;

    void mark_tau_dirty(int tau_index) override;

    FreeEnergyCalculator<T> const &get_calculator() const { return this->calculator; }

    private:
    double alpha;

    double integrand;
    double reference_integrand;

    double proposed_integrand;
    double proposed_reference_integrand;

    FreeEnergyCalculator<T> &calculator;

    double compute_weight(double finite_U, double infinite_U) const;
    std::pair<double, double> compute_integrands() const;
  };

} // namespace sc_expansion::atomic
