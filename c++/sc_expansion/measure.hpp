#pragma once
#include "configuration.hpp"
#include <triqs/arrays.hpp>
#include <triqs/stat/accumulator.hpp>
//#include <triqs/stat.hpp>
#include "myjackknife.hpp"
#include <iostream>
#include <h5/h5.hpp>
#include <fstream>

struct measure {

  Configuration *config;

  // TRIQS Accumulators
  // These replace the huge vectors. They store the binned statistics
  // required for the Jackknife automatically.
  triqs::stat::accumulator<double> acc_sign;
  triqs::stat::accumulator<double> acc_ref;
  double reference_integral;
  double mu;

  // Constructor
  // We need n_bins (stat quality) and block_size (capacity per bin).
  // Typically n_bins ~ 100-500 is sufficient.
  // block_size should be total_mc_steps / n_bins.
  measure(Configuration *config_, double reference_integral_, int n_bins, int block_size, double mu_)
     : config(config_),
       acc_sign(0.0, 0, n_bins, block_size + 100), // +100 buffer for safety
       acc_ref(0.0, 0, n_bins, block_size + 100),
       reference_integral(reference_integral_),
       mu(mu_) {}

  // Accumulate is now O(1) memory and very fast
  void accumulate(double) {
    // Calculate the ratio quantity
    double current_ref_val = config->get_reference_weight() / config->weight;

    // Push directly to accumulators
    acc_sign << config->sign;
    acc_ref << current_ref_val;
  }

  // Perform the distributed Jackknife
  // We pass 'reference_integral' here as it is needed for the final calculation
  void collect_results(mpi::communicator c) {

    // Define the function f(avg_sign, avg_ref) -> physical_result
    auto ratio_func = [this](double avg_sign, double avg_ref) { return avg_sign * this->reference_integral / avg_ref; };

    // triqs::stat::jackknife_mpi automatically:
    // 1. Reduces the accumulators across all MPI nodes
    // 2. Performs the Jackknife resampling
    // 3. Returns the Mean and Error
    auto result = triqs::stat::local::jackknife_mpi(c, ratio_func, acc_sign, acc_ref);

    if (c.rank() == 0) {
      std::cout << "--- Measurement Results ---" << std::endl;

      // Optional: Check total samples used
      // We reduce the count just for display purposes
      //long total_samples = mpi::reduce(acc_sign.count(), c);
      //std::cout << "Total Samples:   " << total_samples << std::endl;

      std::cout << "Jackknife Mean:  " << std::get<0>(result) << std::endl;
      std::cout << "Jackknife Error: " << std::get<1>(result) << std::endl;

      {
        std::string filename = "results/data_mu_" + std::to_string(mu) + ".h5";
        h5::file file(filename, 'w'); // 'w' is safe now, we are the only owner of this file!
        h5_write(file, "mean", std::get<0>(result));
        h5_write(file, "error", std::get<1>(result));
        h5_write(file, "mu", mu);
      }
    }
  }
};