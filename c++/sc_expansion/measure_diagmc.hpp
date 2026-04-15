#pragma once
#include "diagmc_configuration.hpp"
#include "measure_dimer.hpp"
#include <triqs/stat/accumulator.hpp>
#include "myjackknife.hpp"
#include <iostream>
#include <chrono>
#include <memory>
#include <vector>

template <typename T> struct MeasureDiagMC {

  DiagMCConfiguration<T> *config;

  // Defensive ratio estimator accumulators:
  //   acc_numerator:   D_d(tau) / W    where W = |D_d| + c
  //   acc_denominator: c / W
  // Result: coeff = alpha * beta^n * <num> / <den>
  triqs::stat::accumulator<double> acc_numerator;
  triqs::stat::accumulator<double> acc_denominator;

  double alpha;
  double c; // per-diagram defensive constant = alpha / N

  // Shared result struct (reuses DimerMeasureResult) — survives copy into mc_generic internals
  std::shared_ptr<DimerMeasureResult> result;

  // Diagram visit statistics — shared_ptr so copies inside mc_generic update the same counters
  std::shared_ptr<std::vector<long>> visit_counts;

  // Progress tracking
  long step_count  = 0;
  int report_every = 5000;
  int verbosity    = 0;
  std::chrono::high_resolution_clock::time_point last_report;

  MeasureDiagMC(DiagMCConfiguration<T> *config_, int n_bins, int block_size, int verbosity_ = 0)
     : config(config_), acc_numerator(0.0, 0, n_bins, block_size + 100), acc_denominator(0.0, 0, n_bins, block_size + 100),
       alpha(config_->get_alpha()), c(config_->get_defensive_per_diagram()), result(std::make_shared<DimerMeasureResult>()),
       visit_counts(std::make_shared<std::vector<long>>(config_->get_n_diagrams() + 1, 0)), verbosity(verbosity_),
       last_report(std::chrono::high_resolution_clock::now()) {}

  void accumulate(double) {
    int d    = config->current_diagram;
    double W = config->current_abs; // |D_d| + c, always > 0

    (*this->visit_counts)[d]++;

    // Continuous defensive ratio estimator
    this->acc_numerator << (config->current_signed / W);
    this->acc_denominator << (this->c / W);

    this->step_count++;
    if (this->verbosity > 1 && this->step_count % this->report_every == 0) {
      auto now             = std::chrono::high_resolution_clock::now();
      double dt            = std::chrono::duration<double>(now - this->last_report).count();
      double steps_per_sec = this->report_every / dt;
      std::cout << "[MeasureDiagMC rank " << mpi::communicator().rank() << "] step " << this->step_count << " | " << steps_per_sec
                << " steps/s | d=" << d << std::endl;
      this->last_report = now;
    }
  }

  void collect_results(mpi::communicator comm) {
    int order   = config->get_order();
    double norm = this->alpha * std::pow(config->beta, order);

    auto ratio_func = [norm](double avg_num, double avg_den) {
      if (std::abs(avg_den) < 1e-18) return 0.0;
      return (avg_num / avg_den) * norm;
    };

    auto jk = triqs::stat::local::jackknife_mpi(comm, ratio_func, acc_numerator, acc_denominator);

    this->result->coeff = std::get<0>(jk);
    this->result->error = std::get<1>(jk);

    // Diagnostic quantities from individual accumulators
    auto identity_func = [](double x) { return x; };
    auto num_jk        = triqs::stat::local::jackknife_mpi(comm, identity_func, acc_numerator);
    auto den_jk        = triqs::stat::local::jackknife_mpi(comm, identity_func, acc_denominator);

    this->result->mean_sign  = std::get<0>(num_jk);
    this->result->sign_error = std::get<1>(num_jk);
    this->result->mean_abs   = std::get<0>(den_jk);
    this->result->abs_error  = std::get<1>(den_jk);

    if (comm.rank() == 0) {
      std::cout << "--- DiagMC Results (defensive ratio estimator) ---" << std::endl;
      std::cout << "Alpha:               " << this->alpha << std::endl;
      std::cout << "c (per-diagram):     " << this->c << std::endl;
      std::cout << "Mean D_d/W:          " << this->result->mean_sign << std::endl;
      std::cout << "Mean c/W:            " << this->result->mean_abs << std::endl;
      std::cout << "Jackknife Coeff:     " << this->result->coeff << std::endl;
      std::cout << "Jackknife Error:     " << this->result->error << std::endl;

      // Visit statistics
      auto const &vc = *this->visit_counts;
      long total_visits = 0;
      for (auto v : vc) total_visits += v;
      std::cout << "\nTotal number of measures: " << total_visits << std::endl;
      std::cout << "\n--- Diagram visit statistics ---" << std::endl;
      auto const &graphs = config->get_calculator().get_graphs();
      for (int di = 1; di <= config->get_n_diagrams(); di++) {
        std::cout << "d=" << di << " (V=" << graphs[di - 1].get_V() << "): " << vc[di] << " (" << 100.0 * vc[di] / total_visits << "%)"
                  << std::endl;
      }
    }
  }
};
