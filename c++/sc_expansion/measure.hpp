#pragma once
// DEPRECATED: these estimators assume the old atomic sampling weight
//   W = |f - f_ref| + alpha * |f_ref|
// (with f_ref = the U=infinity reference integrand). The atomic Configuration
// now uses the uniform-reference weight W = |f + alpha| (matching the dimer
// scheme), and the corresponding estimator lives in dimer/measure_dimer.hpp.
// Nothing currently includes this file. Kept as a reference for the old scheme.
#include "configuration_base.hpp"
#include <triqs/arrays.hpp>
#include <triqs/stat/accumulator.hpp>
#include "myjackknife.hpp"
#include <cmath>
#include <iostream>
#include <chrono>
#include <memory>

struct MeasureResult {
  double mean  = 0.0;
  double error = 0.0;
};

// Per-step inputs handed to an estimator. The finite_U / infinite_U values are
// the ones the configuration already cached when it computed the metropolis
// weight, so no diagram re-evaluation happens in the estimator.
struct EstimatorInputs {
  double finite_U;
  double infinite_U;
  double weight;
};

struct EstimatorSample {
  double integrand;   // pushed to acc_integrand
  double denominator; // pushed to acc_denominator
};

// Estimator for the free-energy / density-mode series. Reference is the
// atomic-limit integrand, which has a non-vanishing absolute integral R_a; the
// denominator <|f_ref|/W> calibrates the ratio against R_a as the multiplier.
//
//   I = R_s + R_a · <(f - f_ref) / W> / <|f_ref| / W>
struct free_energy_estimator {
  double reference_integral;        // R_a
  double signed_reference_integral; // R_s

  EstimatorSample sample(EstimatorInputs const &x) const {
    return {(x.finite_U - x.infinite_U) / x.weight, std::abs(x.infinite_U) / x.weight};
  }

  double combine(double avg_int, double avg_den) const {
    if (std::abs(avg_den) < 1e-300) return 0.0;
    return reference_integral * (avg_int / avg_den) + signed_reference_integral;
  }

  void print_header(std::ostream &os) const {
    os << "--- Measurement Results (free_energy_estimator, I = R_s + R_a·<(f-f_ref)/W>/<|f_ref|/W>) ---" << std::endl;
    os << "Reference Integral (abs):    " << reference_integral << std::endl;
    os << "Reference Integral (signed): " << signed_reference_integral << std::endl;
  }
};

// Estimator for the density-density correlator series. At half-filling off-site
// R_a → 0, so the free-energy estimator's denominator collapses; here we
// replace it with <1/W> → V/Z_W where V is the MC sampling domain volume.
// The atomic move samples each τ independently on [0, β], so V = β^order
// (hypercube), NOT β^order / order!.
//
//   I = R_s + V · <(f - f_ref) / W> / <1 / W>
struct density_density_estimator {
  double domain_volume;             // V = β^order
  double signed_reference_integral; // R_s
  double reference_integral;        // R_a, kept for diagnostic printout only

  EstimatorSample sample(EstimatorInputs const &x) const {
    return {(x.finite_U - x.infinite_U) / x.weight, 1.0 / x.weight};
  }

  double combine(double avg_int, double avg_den) const {
    if (std::abs(avg_den) < 1e-300) return 0.0;
    return domain_volume * (avg_int / avg_den) + signed_reference_integral;
  }

  void print_header(std::ostream &os) const {
    os << "--- Measurement Results (density_density_estimator, I = R_s + V·<(f-f_ref)/W>/<1/W>) ---" << std::endl;
    os << "Reference Integral (abs):    " << reference_integral << std::endl;
    os << "Reference Integral (signed): " << signed_reference_integral << std::endl;
    os << "Domain Volume:               " << domain_volume << std::endl;
  }
};

template <typename T, typename Estimator> struct measure {

  ConfigurationBase<T> *config;
  Estimator estimator;

  triqs::stat::accumulator<double> acc_integrand;
  triqs::stat::accumulator<double> acc_denominator;

  // Shared result struct — survives copy into mc_generic internals
  std::shared_ptr<MeasureResult> result;

  // Progress tracking
  long step_count  = 0;
  int report_every = 5000;
  int verbosity    = 0;
  std::chrono::high_resolution_clock::time_point last_report;

  measure(ConfigurationBase<T> *config_, Estimator estimator_, int n_bins, int block_size, int verbosity_ = 0)
     : config(config_),
       estimator(estimator_),
       acc_integrand(0.0, 0, n_bins, block_size + 100),
       acc_denominator(0.0, 0, n_bins, block_size + 100),
       result(std::make_shared<MeasureResult>()),
       verbosity(verbosity_),
       last_report(std::chrono::high_resolution_clock::now()) {}

  void accumulate(double) {
    double W = config->metropolis_weight;

    if (W > 0.0) {
      EstimatorInputs x{config->get_integrand(), config->get_reference_integrand(), W};
      EstimatorSample s = estimator.sample(x);
      acc_integrand << s.integrand;
      acc_denominator << s.denominator;
    }

    this->step_count++;
    if (this->verbosity > 0 && this->step_count % this->report_every == 0) {
      auto now             = std::chrono::high_resolution_clock::now();
      double dt            = std::chrono::duration<double>(now - this->last_report).count();
      double steps_per_sec = this->report_every / dt;
      std::cout << "[measure] step " << this->step_count << " | " << steps_per_sec << " steps/s" << std::endl;
      this->last_report = now;
    }
  }

  void collect_results(mpi::communicator c) {

    auto estimator_func = [this](double avg_int, double avg_den) { return this->estimator.combine(avg_int, avg_den); };

    auto jk = triqs::stat::local::jackknife_mpi(c, estimator_func, acc_integrand, acc_denominator);

    this->result->mean  = std::get<0>(jk);
    this->result->error = std::get<1>(jk);

    if (c.rank() == 0) {
      estimator.print_header(std::cout);
      std::cout << "Jackknife Mean:              " << this->result->mean << std::endl;
      std::cout << "Jackknife Error:             " << this->result->error << std::endl;
    }
  }
};
