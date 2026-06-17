#pragma once

#include "../configuration_base.hpp"
#include "sum_diagrams.hpp"

namespace sc_expansion::atomic {

  // Single-site MCMC configuration with uniform-reference defensive weight:
  //   W = |f(τ) + α|
  // where f(τ) is either the free-energy or rooted density-density integrand
  // (selected by the SumDiagrams calculator's mode). Mirrors the scheme used
  // by sc_expansion::dimer::Configuration so the same measure_dimer estimator
  // serves both expansions.
  template <typename T> class Configuration : public ConfigurationBase<T> {
    public:
    Configuration(Parameters<T> const &params, int order, double alpha, SumDiagrams<T> &calculator);

    double evaluate_proposed() override;
    void commit_proposal() override;

    double get_integrand() const override;
    double get_reference_integrand() const override;

    void mark_tau_dirty(int tau_index) override;

    double get_alpha() const { return this->alpha; }
    SumDiagrams<T> const &get_calculator() const { return this->calculator; }

    private:
    double alpha;

    double integrand;
    double proposed_integrand;

    SumDiagrams<T> &calculator;
    int target_d_sq;

    double compute_integrand() const;
  };

} // namespace sc_expansion::atomic
