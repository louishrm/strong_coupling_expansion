#pragma once
#include "dimer_configuration.hpp"
#include <triqs/stat/accumulator.hpp>
#include "myjackknife.hpp"
#include <iostream>
#include <random>
#include <chrono>
#include <memory>

struct DimerMeasureResult {
  double coeff      = 0.0;
  double error      = 0.0;
  double mean_sign  = 0.0;
  double sign_error = 0.0;
  double mean_abs   = 0.0;
  double abs_error  = 0.0;
};

template <typename T> struct measure_dimer {

  DimerConfiguration<T> *config;

  // Accumulator for the sign of Omega (from |Omega| Metropolis sampling)
  triqs::stat::accumulator<double> acc_sign;

  // Accumulator for |Omega| from uniform sampling (normalization estimate)
  triqs::stat::accumulator<double> acc_norm;

  double mu;

  // Shared result struct — survives copy into mc_generic internals
  std::shared_ptr<DimerMeasureResult> result;

  // RNG for drawing uniform tau samples (independent of the Metropolis chain)
  std::mt19937 rng;
  std::uniform_real_distribution<double> tau_dist;

  // Progress tracking
  long step_count  = 0;
  int report_every = 5000;
  int verbosity    = 0;
  std::chrono::high_resolution_clock::time_point last_report;

  measure_dimer(DimerConfiguration<T> *config_, int n_bins, int block_size, double mu_, int random_seed, int verbosity_ = 0)
     : config(config_), acc_sign(0.0, 0, n_bins, block_size + 100), acc_norm(0.0, 0, n_bins, block_size + 100), mu(mu_),
       result(std::make_shared<DimerMeasureResult>()),
       rng(random_seed), tau_dist(0.0, config_->beta), verbosity(verbosity_),
       last_report(std::chrono::high_resolution_clock::now()) {}

  void accumulate(double) {
    double W = config->metropolis_weight;

    // Sign estimation from |Omega| Metropolis sampling
    if (W > 0.0) { acc_sign << (config->get_integrand() >= 0.0 ? 1.0 : -1.0); }

    // Normalization estimation: evaluate |Omega| at uniform random taus
    int order = config->get_order();
    std::vector<double> uniform_taus(order);
    for (int i = 0; i < order; i++) { uniform_taus[i] = tau_dist(rng); }
    acc_norm << std::abs(config->evaluate_at(uniform_taus));

    this->step_count++;
    if (this->verbosity > 0 && this->step_count % this->report_every == 0) {
      auto now     = std::chrono::high_resolution_clock::now();
      double dt    = std::chrono::duration<double>(now - this->last_report).count();
      double steps_per_sec = this->report_every / dt;
      std::cout << "[measure_dimer] step " << this->step_count << " | " << steps_per_sec << " steps/s" << std::endl;
      this->last_report = now;
    }
  }

  void collect_results(mpi::communicator c) {

    auto identity_func = [](double x) { return x; };

    auto sign_jk = triqs::stat::local::jackknife_mpi(c, identity_func, acc_sign);
    auto norm_jk = triqs::stat::local::jackknife_mpi(c, identity_func, acc_norm);

    int order   = config->get_order();
    double beta_n = std::pow(config->beta, order);

    this->result->mean_sign  = std::get<0>(sign_jk);
    this->result->sign_error = std::get<1>(sign_jk);
    this->result->mean_abs   = std::get<0>(norm_jk);
    this->result->abs_error  = std::get<1>(norm_jk);

    // Combined coefficient: beta^n * <|Omega|>_uniform * <sign>
    this->result->coeff = beta_n * this->result->mean_abs * this->result->mean_sign;
    // Error propagation (independent estimators):
    // delta(coeff) = beta^n * sqrt( (<|Omega|> * delta_sign)^2 + (<sign> * delta_|Omega|)^2 )
    this->result->error = beta_n * std::sqrt(this->result->mean_abs * this->result->mean_abs * this->result->sign_error * this->result->sign_error
                                             + this->result->mean_sign * this->result->mean_sign * this->result->abs_error * this->result->abs_error);

    if (c.rank() == 0) {
      std::cout << "--- Measurement Results (Dimer, uniform reference) ---" << std::endl;
      std::cout << "Mean sign:           " << this->result->mean_sign << std::endl;
      std::cout << "Sign error:          " << this->result->sign_error << std::endl;
      std::cout << "Mean |Omega| (unif): " << this->result->mean_abs << std::endl;
      std::cout << "|Omega| error:       " << this->result->abs_error << std::endl;
      std::cout << "Jackknife Mean:      " << this->result->coeff << std::endl;
      std::cout << "Jackknife Error:     " << this->result->error << std::endl;
    }
  }
};
