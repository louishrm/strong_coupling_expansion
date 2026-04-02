#pragma once
#include "dimer_configuration.hpp"
#include <triqs/stat/accumulator.hpp>
#include "myjackknife.hpp"
#include <iostream>
#include <h5/h5.hpp>
#include <filesystem>
#include <random>

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

  measure_dimer(DimerConfiguration<T> *config_, int n_bins, int block_size, double mu_, int random_seed)
     : config(config_), acc_sign(0.0, 0, n_bins, block_size + 100), acc_norm(0.0, 0, n_bins, block_size + 100), mu(mu_),
       rng(random_seed), tau_dist(0.0, config_->beta) {}

  void accumulate(double) {
    double W = config->metropolis_weight;

    // Sign estimation from |Omega| Metropolis sampling
    if (W > 0.0) { acc_sign << (config->get_integrand() >= 0.0 ? 1.0 : -1.0); }

    // Normalization estimation: evaluate |Omega| at uniform random taus
    int order = config->get_order();
    std::vector<double> uniform_taus(order);
    for (int i = 0; i < order; i++) { uniform_taus[i] = tau_dist(rng); }
    acc_norm << std::abs(config->evaluate_at(uniform_taus));
  }

  void collect_results(mpi::communicator c) {

    auto identity_func = [](double x) { return x; };

    auto sign_result = triqs::stat::local::jackknife_mpi(c, identity_func, acc_sign);
    auto norm_result = triqs::stat::local::jackknife_mpi(c, identity_func, acc_norm);

    if (c.rank() == 0) {
      std::cout << "--- Measurement Results (Dimer, uniform reference) ---" << std::endl;
      std::cout << "Mean sign:           " << std::get<0>(sign_result) << std::endl;
      std::cout << "Sign error:          " << std::get<1>(sign_result) << std::endl;
      std::cout << "Mean |Omega| (unif): " << std::get<0>(norm_result) << std::endl;
      std::cout << "|Omega| error:       " << std::get<1>(norm_result) << std::endl;

      if (std::filesystem::is_directory("./results")) {
        std::string filename = "./results/dimer_data_order_" + std::to_string(config->get_order()) + "_U_" + std::to_string(config->get_U()) + "_beta_"
           + std::to_string(config->beta) + "_mu_" + std::to_string(mu) + ".h5";
        h5::file file(filename, 'w');
        h5_write(file, "mean_sign", std::get<0>(sign_result));
        h5_write(file, "sign_error", std::get<1>(sign_result));
        h5_write(file, "mean_abs_integrand", std::get<0>(norm_result));
        h5_write(file, "abs_integrand_error", std::get<1>(norm_result));
        h5_write(file, "mu", mu);
      }
    }
  }
};
