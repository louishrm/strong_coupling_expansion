#pragma once

#include "../configuration_base.hpp"
#include "free_energy_calculator.hpp"

namespace sc_expansion::dimer {

  // Two-site MCMC configuration. Reference weight is uniform: weight = |omega| + alpha,
  // and the reference integrand exposed to the measure is the constant `alpha`.
  // Contrast with sc_expansion::atomic::Configuration, which uses the infinite-U
  // integrand as the reference.
  //
  // `Calculator` is the diagram-series engine driven each MC step. It defaults to
  // FreeEnergyCalculator<T> (the free-energy / Omega series). Any type exposing the
  // same surface works — in particular SumDiagrams<T>, which carries the rooted
  // ⟨n(r)n(0)⟩ density-density series; for that calculator the integrand `omega`
  // is the correlator value at the sampled imaginary times rather than the free
  // energy. Required surface:
  //   T    compute_sum_diagrams(std::vector<double> const &taus) const;
  //   void clear_all_caches();
  //   void mark_tau_dirty(int tau_index);
  template <typename T, typename Calculator = FreeEnergyCalculator<T>> class Configuration : public ConfigurationBase<T> {
    public:
    Configuration(Parameters<T> const &params, int order, Calculator &calculator, double alpha = 0.001);

    double evaluate_proposed() override;
    void commit_proposal() override;

    double get_integrand() const override;
    double get_reference_integrand() const override;

    void mark_tau_dirty(int tau_index) override;

    double get_alpha() const { return this->alpha; }
    Calculator const &get_calculator() const { return this->calculator; }

    private:
    double omega;
    double proposed_omega;
    double alpha;
    Calculator &calculator;

    double compute_omega() const;
  };

} // namespace sc_expansion::dimer
