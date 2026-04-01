#pragma once
#include "configuration_base.hpp"
#include <triqs/stat/accumulator.hpp>
#include "myjackknife.hpp"
#include <iostream>
#include <h5/h5.hpp>
#include <filesystem>

template <typename T> struct measure_dimer {

  ConfigurationBase<T> *config;

  // Accumulator for the sign of Omega
  triqs::stat::accumulator<double> acc_sign;

  double mu;

  measure_dimer(ConfigurationBase<T> *config_, int n_bins, int block_size, double mu_)
     : config(config_), acc_sign(0.0, 0, n_bins, block_size + 100), mu(mu_) {}

  void accumulate(double) {
    double W = config->metropolis_weight;

    // W = |Omega|; only accumulate if nonzero
    if (W > 0.0) { acc_sign << (config->get_integrand() >= 0.0 ? 1.0 : -1.0); }
  }

  void collect_results(mpi::communicator c) {

    auto identity_func = [](double avg_sign) { return avg_sign; };

    auto result = triqs::stat::local::jackknife_mpi(c, identity_func, acc_sign);

    if (c.rank() == 0) {
      std::cout << "--- Measurement Results (Dimer, |Omega| sampling) ---" << std::endl;
      std::cout << "Mean sign:     " << std::get<0>(result) << std::endl;
      std::cout << "Sign error:    " << std::get<1>(result) << std::endl;

      if (std::filesystem::is_directory("./results")) {
        std::string filename = "./results/dimer_data_order_" + std::to_string(config->get_order()) + "_U_" + std::to_string(config->get_U()) + "_beta_"
           + std::to_string(config->beta) + "_mu_" + std::to_string(mu) + ".h5";
        h5::file file(filename, 'w');
        h5_write(file, "mean_sign", std::get<0>(result));
        h5_write(file, "sign_error", std::get<1>(result));
        h5_write(file, "mu", mu);
      }
    }
  }
};
