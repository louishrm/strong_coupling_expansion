#pragma once

#include "configuration_base.hpp"
#include "free_energy_calculator.hpp"

template <typename T> class DimerConfiguration : public ConfigurationBase<T> {

  public:
  DimerConfiguration(sc_expansion::Parameters<T> const &params, int order,
                     sc_expansion::FreeEnergyCalculator<2, T> &calculator);

  double evaluate_proposed() override;
  void commit_proposal() override;

  double get_integrand() const override;
  double get_reference_integrand() const override;

  void mark_tau_dirty(int tau_index) override;

  double evaluate_at(std::vector<double> const &taus) const;

  private:
  // Committed cache
  double omega;

  // Tentative cache
  double proposed_omega;

  // Diagram evaluation (non-owning reference to externally-owned calculator)
  sc_expansion::FreeEnergyCalculator<2, T> &calculator;

  double compute_omega() const;
};
