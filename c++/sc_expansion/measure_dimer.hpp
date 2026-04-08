#pragma once
#include "dimer_configuration.hpp"
#include <triqs/stat/accumulator.hpp>
#include "myjackknife.hpp"
#include <iostream>
#include <h5/h5.hpp>
#include <filesystem>
#include <random>
#include <chrono>

template <typename T> struct measure_dimer {

  DimerConfiguration<T> *config;

  // Accumulator for the sign of Omega (from |Omega| Metropolis sampling)
  triqs::stat::accumulator<double> acc_sign;

  // Accumulator for |Omega| from uniform sampling (normalization estimate)
  triqs::stat::accumulator<double> acc_norm;

  double mu;

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

    auto sign_result = triqs::stat::local::jackknife_mpi(c, identity_func, acc_sign);
    auto norm_result = triqs::stat::local::jackknife_mpi(c, identity_func, acc_norm);

    if (c.rank() == 0) {
      double mean_sign      = std::get<0>(sign_result);
      double sign_error     = std::get<1>(sign_result);
      double mean_abs       = std::get<0>(norm_result);
      double abs_error      = std::get<1>(norm_result);
      int order             = config->get_order();
      double beta_n         = std::pow(config->beta, order);

      // Combined coefficient: beta^n * <|Omega|>_uniform * <sign>
      double coeff          = beta_n * mean_abs * mean_sign;
      // Error propagation (independent estimators):
      // delta(coeff) = beta^n * sqrt( (<|Omega|> * delta_sign)^2 + (<sign> * delta_|Omega|)^2 )
      double coeff_error    = beta_n * std::sqrt(mean_abs * mean_abs * sign_error * sign_error
                                                 + mean_sign * mean_sign * abs_error * abs_error);

      std::cout << "--- Measurement Results (Dimer, uniform reference) ---" << std::endl;
      std::cout << "Mean sign:           " << mean_sign << std::endl;
      std::cout << "Sign error:          " << sign_error << std::endl;
      std::cout << "Mean |Omega| (unif): " << mean_abs << std::endl;
      std::cout << "|Omega| error:       " << abs_error << std::endl;
      std::cout << "Jackknife Mean:      " << coeff << std::endl;
      std::cout << "Jackknife Error:     " << coeff_error << std::endl;

      if (std::filesystem::is_directory("./results")) {
        std::string filename = "./results/dimer_data_order_" + std::to_string(order) + "_U_" + std::to_string(config->get_U()) + "_beta_"
           + std::to_string(config->beta) + "_mu_" + std::to_string(mu) + ".h5";
        h5::file file(filename, 'w');
        h5_write(file, "mean", coeff);
        h5_write(file, "error", coeff_error);
        h5_write(file, "mean_sign", mean_sign);
        h5_write(file, "sign_error", sign_error);
        h5_write(file, "mean_abs_integrand", mean_abs);
        h5_write(file, "abs_integrand_error", abs_error);
        h5_write(file, "mu", mu);
      }
    }
  }
};
