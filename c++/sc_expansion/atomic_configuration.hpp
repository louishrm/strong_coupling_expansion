#pragma once

#include "configuration_base.hpp"
#include "free_energy_calculator.hpp"

template <typename T> class AtomicConfiguration : public ConfigurationBase<T> {

  public:
  AtomicConfiguration(sc_expansion::Parameters<T> const &params, int order, double alpha, int override_fm = -1);

  double evaluate_proposed() override;
  void commit_proposal() override;

  double get_integrand() const override;
  double get_reference_integrand() const override;

  sc_expansion::FreeEnergyCalculator<1, T> const &get_calculator() const { return this->calculator; }

  private:
  double alpha;

  // Committed cache
  double integrand;
  double reference_integrand;

  // Tentative cache (filled by evaluate_proposed, consumed by commit_proposal)
  double proposed_integrand;
  double proposed_reference_integrand;

  // Diagram evaluation
  sc_expansion::FreeEnergyCalculator<1, T> calculator;

  double compute_weight(double finite_U, double infinite_U) const;
  std::pair<double, double> compute_integrands() const;
};
