#pragma once
#include "configuration.hpp"
#include <triqs/stat/accumulator.hpp>
#include "../myjackknife.hpp"
#include <iostream>
#include <chrono>
#include <memory>
#include <cmath>

struct DimerMeasureResult {
  double coeff          = 0.0;
  double error          = 0.0;
  double mean_sign      = 0.0;
  double sign_error     = 0.0;
  double mean_abs       = 0.0;
  double abs_error      = 0.0;
  double mean_omega_abs = 0.0; // <|Omega|/W>, used for alpha auto-tuning
};

template <typename T> struct measure_dimer {

  sc_expansion::dimer::Configuration<T> *config;

  // Accumulators for defensive importance sampling ratio estimator:
  //   acc_integrand: accumulates Omega / W
  //   acc_reference: accumulates alpha / W
  // where W = |Omega + alpha| is the Metropolis weight.
  triqs::stat::accumulator<double> acc_integrand;
  triqs::stat::accumulator<double> acc_reference;
  triqs::stat::accumulator<double> acc_abs_integrand;

  double alpha;
  double mu;

  // Shared result struct — survives copy into mc_generic internals
  std::shared_ptr<DimerMeasureResult> result;

  // Progress tracking
  long step_count  = 0;
  int report_every = 5000;
  int verbosity    = 0;
  std::chrono::high_resolution_clock::time_point last_report;

  measure_dimer(sc_expansion::dimer::Configuration<T> *config_, int n_bins, int block_size, double mu_, int verbosity_ = 0)
     : config(config_),
       acc_integrand(0.0, 0, n_bins, block_size + 100),
       acc_reference(0.0, 0, n_bins, block_size + 100),
       acc_abs_integrand(0.0, 0, n_bins, block_size + 100),
       alpha(config_->get_alpha()),
       mu(mu_),
       result(std::make_shared<DimerMeasureResult>()),
       verbosity(verbosity_),
       last_report(std::chrono::high_resolution_clock::now()) {}

  void accumulate(double) {
    double W = config->metropolis_weight;

    if (W > 0.0) {
      double omega = config->get_integrand();
      acc_integrand << (omega / W);
      acc_reference << (this->alpha / W);
      acc_abs_integrand << (std::abs(omega) / W);
    }

    this->step_count++;
    if (this->verbosity > 0 && this->step_count % this->report_every == 0) {
      auto now             = std::chrono::high_resolution_clock::now();
      double dt            = std::chrono::duration<double>(now - this->last_report).count();
      double steps_per_sec = this->report_every / dt;
      std::cout << "[measure_dimer] step " << this->step_count << " | " << steps_per_sec << " steps/s" << std::endl;
      this->last_report = now;
    }
  }

  void collect_results(mpi::communicator c) {

    // Ratio estimator: coeff = alpha * beta^n * <Omega/W> / <alpha/W>
    int order   = config->get_order();
    double norm = this->alpha * std::pow(config->beta, order);

    auto ratio_func = [norm](double avg_int, double avg_ref) {
      if (std::abs(avg_ref) < 1e-18) return 0.0;
      return (avg_int / avg_ref) * norm;
    };

    auto jk = triqs::stat::local::jackknife_mpi(c, ratio_func, acc_integrand, acc_reference);

    this->result->coeff = std::get<0>(jk);
    this->result->error = std::get<1>(jk);

    auto pick_first  = [](double x, double) { return x; };
    auto pick_second = [](double, double y) { return y; };
    auto int_jk      = triqs::stat::local::jackknife_mpi(c, pick_first, acc_integrand, acc_reference);
    auto ref_jk      = triqs::stat::local::jackknife_mpi(c, pick_second, acc_integrand, acc_reference);

    this->result->mean_sign  = std::get<0>(int_jk);
    this->result->sign_error = std::get<1>(int_jk);
    this->result->mean_abs   = std::get<0>(ref_jk);
    this->result->abs_error  = std::get<1>(ref_jk);

    auto abs_jk                  = triqs::stat::local::jackknife_mpi(c, pick_first, acc_abs_integrand, acc_reference);
    this->result->mean_omega_abs = std::get<0>(abs_jk);

    if (c.rank() == 0) {
      std::cout << "--- Measurement Results (Dimer, defensive ratio estimator) ---" << std::endl;
      std::cout << "Alpha:               " << this->alpha << std::endl;
      std::cout << "Mean Omega/W:        " << this->result->mean_sign << std::endl;
      std::cout << "Mean alpha/W:        " << this->result->mean_abs << std::endl;
      std::cout << "Jackknife Coeff:     " << this->result->coeff << std::endl;
      std::cout << "Jackknife Error:     " << this->result->error << std::endl;
    }
  }
};
