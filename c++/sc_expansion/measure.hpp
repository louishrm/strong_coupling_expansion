#pragma once
#include "configuration.hpp"
#include <triqs/arrays.hpp>
#include <triqs/stat/accumulator.hpp>
#include "myjackknife.hpp"
#include <iostream>
#include <h5/h5.hpp>
#include <fstream>
#include <cmath>

struct measure {

  Configuration *config;

  // Accumulators for defensive importance sampling
  // We want to estimate I = I_ref * <integrand/W> / <reference_integrand/W>
  triqs::stat::accumulator<double> acc_integrand;
  triqs::stat::accumulator<double> acc_reference;

  double reference_integral;
  double mu;

  measure(Configuration *config_, double reference_integral_, int n_bins, int block_size, double mu_)
     : config(config_),
       acc_integrand(0.0, 0, n_bins, block_size + 100),
       acc_reference(0.0, 0, n_bins, block_size + 100),
       reference_integral(reference_integral_),
       mu(mu_) {}

  void accumulate(double) {
    double W = config->metropolis_weight;

    // Safety check for W=0, though MC should not visit such states
    if (W > 0.0) {
      acc_integrand << (config->integrand / W);
      acc_reference << (config->reference_integrand / W);
    }
  }

  void collect_results(mpi::communicator c) {

    // The ratio estimator: I = I_ref * (avg(integrand/W) / avg(ref_integrand/W))
    auto ratio_func = [this](double avg_int, double avg_ref) {
      if (std::abs(avg_ref) < 1e-18) return 0.0;
      return (avg_int / avg_ref) * this->reference_integral;
    };

    // Perform Jackknife on the ratio of the two accumulators
    auto result = triqs::stat::local::jackknife_mpi(c, ratio_func, acc_integrand, acc_reference);

    if (c.rank() == 0) {
      std::cout << "--- Measurement Results (Defensive Importance Sampling) ---" << std::endl;
      std::cout << "Reference Integral: " << reference_integral << std::endl;
      std::cout << "Jackknife Mean:     " << std::get<0>(result) << std::endl;
      std::cout << "Jackknife Error:    " << std::get<1>(result) << std::endl;

      // Ensure directory exists (basic check, usually handled by build system or user)
      // For now, we assume 'results' directory exists or we write to current dir if needed.
      try {
        std::string filename = "./results/full_lattice_data_mu_" + std::to_string(mu) + ".h5";
        h5::file file(filename, 'w');
        h5_write(file, "mean", std::get<0>(result));
        h5_write(file, "error", std::get<1>(result));
        h5_write(file, "mu", mu);
        h5_write(file, "reference_integral", reference_integral);
      } catch (...) { std::cerr << "Warning: Could not write HDF5 results file." << std::endl; }
    }
  }
};
