#pragma once
#include "configuration_base.hpp"
#include <triqs/arrays.hpp>
#include <triqs/stat/accumulator.hpp>
#include "myjackknife.hpp"
#include <iostream>
#include <chrono>
#include <memory>

struct MeasureResult {
  double mean  = 0.0;
  double error = 0.0;
};

template <typename T> struct measure {

  ConfigurationBase<T> *config;

  // Accumulators for defensive importance sampling
  // We want to estimate I = I_ref * <integrand/W> / <reference_integrand/W>
  triqs::stat::accumulator<double> acc_integrand;
  triqs::stat::accumulator<double> acc_reference;

  double reference_integral;
  double signed_reference_integral;
  double mu;

  // Shared result struct — survives copy into mc_generic internals
  std::shared_ptr<MeasureResult> result;

  // Progress tracking
  long step_count  = 0;
  int report_every = 5000;
  int verbosity    = 0;
  std::chrono::high_resolution_clock::time_point last_report;

  measure(ConfigurationBase<T> *config_, double reference_integral_, double signed_reference_integral_, int n_bins, int block_size, double mu_,
          int verbosity_ = 0)
     : config(config_),
       acc_integrand(0.0, 0, n_bins, block_size + 100),
       acc_reference(0.0, 0, n_bins, block_size + 100),
       reference_integral(reference_integral_),
       signed_reference_integral(signed_reference_integral_),
       mu(mu_),
       result(std::make_shared<MeasureResult>()),
       verbosity(verbosity_),
       last_report(std::chrono::high_resolution_clock::now()) {}

  void accumulate(double) {
    double W = config->metropolis_weight;

    // Safety check for W=0, though MC should not visit such states
    if (W > 0.0) {
      acc_integrand << ((config->get_integrand() - config->get_reference_integrand()) / W);
      acc_reference << (std::abs(config->get_reference_integrand()) / W);
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

    // The ratio estimator: I = I_ref * (avg(integrand/W) / avg(ref_integrand/W))
    auto ratio_func = [this](double avg_int, double avg_ref) {
      if (std::abs(avg_ref) < 1e-18) return 0.0;
      return (avg_int / avg_ref) * this->reference_integral + this->signed_reference_integral;
    };

    // Perform Jackknife on the ratio of the two accumulators
    auto jk = triqs::stat::local::jackknife_mpi(c, ratio_func, acc_integrand, acc_reference);

    this->result->mean  = std::get<0>(jk);
    this->result->error = std::get<1>(jk);

    if (c.rank() == 0) {
      std::cout << "--- Measurement Results (Defensive Importance Sampling) ---" << std::endl;
      std::cout << "Reference Integral: " << reference_integral << std::endl;
      std::cout << "Jackknife Mean:     " << this->result->mean << std::endl;
      std::cout << "Jackknife Error:    " << this->result->error << std::endl;
    }
  }
};
