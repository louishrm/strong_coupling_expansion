/*
 * Diagnostic tests for dimer MC sampling at U=8, beta=2 at orders 2, 4, 6.
 *
 * Designed to be CHEAP: no long MCMC runs at high orders.
 *
 * Test 1: Integrand structure — evaluate Omega at O(500) uniform tau samples.
 *         Reports peakedness (CV), dynamic range, effective sample size, sign fraction.
 *         This is the key diagnostic: if the integrand is very peaked, uniform
 *         normalization is doomed regardless of MCMC quality.
 *
 * Test 2: Short MCMC mixing — run ~10k cycles just to measure acceptance rate
 *         and check if the chain explores well.
 *
 * Usage:  mpirun -np 4 ./test_dimer_diagnostics
 */

#include <gtest/gtest.h>
#include "sc_expansion/dimer_configuration.hpp"
#include "sc_expansion/free_energy_calculator.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/measure_dimer.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <mpi/mpi.hpp>
#include <cmath>
#include <memory>
#include <random>
#include <filesystem>
#include <algorithm>

using namespace sc_expansion;

static const std::vector<std::pair<int, int>> triangle_cluster = {{0, 0}, {1, 0}, {0, 1}};
static const int n_cluster_sites = 3;

// =====================================================================
// Test 1: Integrand structure from uniform tau samples
//
// No MCMC — just evaluate Omega at uniform random tau points.
// Reports distribution statistics that reveal whether uniform
// normalization can work.
// =====================================================================

void run_integrand_structure(int order, int n_samples) {
  mpi::communicator world;
  if (world.rank() != 0) return;

  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  bool bipartite = (order == 2);
  Parameters<double> params{U, beta, mu, t_intra, bipartite};

  std::unique_ptr<FreeEnergyCalculator<2, double>> calc;
  if (order == 2) {
    calc = std::make_unique<FreeEnergyCalculator<2, double>>(params, order);
  } else {
    calc = std::make_unique<FreeEnergyCalculator<2, double>>(params, order, triangle_cluster, n_cluster_sites);
  }

  std::mt19937 rng(54321 + order);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);

  std::vector<double> abs_values;
  abs_values.reserve(n_samples);
  int n_positive = 0, n_negative = 0, n_zero = 0;
  double sum_val = 0, sum_val2 = 0;

  for (int s = 0; s < n_samples; s++) {
    std::vector<double> taus(order);
    for (int i = 0; i < order; i++) taus[i] = tau_dist(rng);

    double omega = calc->compute_sum_diagrams_dimer(taus, false);
    double abs_omega = std::abs(omega);
    abs_values.push_back(abs_omega);
    sum_val += abs_omega;
    sum_val2 += abs_omega * abs_omega;

    if (omega > 1e-15) n_positive++;
    else if (omega < -1e-15) n_negative++;
    else n_zero++;
  }

  std::sort(abs_values.begin(), abs_values.end());

  double mean = sum_val / n_samples;
  double mean_sq = sum_val2 / n_samples;
  double var  = mean_sq - mean * mean;
  double std_dev = std::sqrt(std::max(0.0, var));

  double p10  = abs_values[(int)(0.10 * n_samples)];
  double p25  = abs_values[(int)(0.25 * n_samples)];
  double p50  = abs_values[(int)(0.50 * n_samples)];
  double p75  = abs_values[(int)(0.75 * n_samples)];
  double p90  = abs_values[(int)(0.90 * n_samples)];
  double p99  = abs_values[(int)(0.99 * n_samples)];
  double pmax = abs_values.back();

  // n_eff = n * (<|Omega|>)^2 / <|Omega|^2>
  double n_eff = (mean * mean) / mean_sq * n_samples;

  std::cout << "\n========================================" << std::endl;
  std::cout << "  INTEGRAND STRUCTURE — Order " << order << std::endl;
  std::cout << "  U=8 beta=2, " << n_samples << " uniform samples" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "  Mean |Omega|     = " << mean << std::endl;
  std::cout << "  Std dev          = " << std_dev << std::endl;
  std::cout << "  Coeff of var     = " << (mean > 0 ? std_dev / mean : 0) << std::endl;
  std::cout << std::endl;
  std::cout << "  Percentiles of |Omega|:" << std::endl;
  std::cout << "    10%:  " << p10 << std::endl;
  std::cout << "    25%:  " << p25 << std::endl;
  std::cout << "    50%:  " << p50 << " (median)" << std::endl;
  std::cout << "    75%:  " << p75 << std::endl;
  std::cout << "    90%:  " << p90 << std::endl;
  std::cout << "    99%:  " << p99 << std::endl;
  std::cout << "    max:  " << pmax << std::endl;
  std::cout << std::endl;
  std::cout << "  Dynamic range (max/median) = " << (p50 > 0 ? pmax / p50 : 999) << std::endl;
  std::cout << "  Dynamic range (p99/p50)    = " << (p50 > 0 ? p99 / p50 : 999) << std::endl;
  std::cout << "  Fraction |Omega|~0         = " << (double)n_zero / n_samples << std::endl;
  std::cout << std::endl;
  std::cout << "  Effective sample size      = " << n_eff << " / " << n_samples << std::endl;
  std::cout << "  Efficiency                 = " << n_eff / n_samples * 100 << "%" << std::endl;
  std::cout << "  (100%=flat, low=peaked → uniform sampling wastes work)" << std::endl;
  std::cout << std::endl;
  std::cout << "  Sign: +" << n_positive << "  -" << n_negative << "  zero=" << n_zero << std::endl;
  std::cout << "  <sign>_uniform = " << (double)(n_positive - n_negative) / (n_positive + n_negative + 1e-30) << std::endl;
  std::cout << "========================================\n" << std::endl;
}

TEST(DimerDiagnostics, IntegrandStructure_Order2) { run_integrand_structure(2, 2000); }
TEST(DimerDiagnostics, IntegrandStructure_Order4) { run_integrand_structure(4, 500); }
TEST(DimerDiagnostics, IntegrandStructure_Order6) { run_integrand_structure(6, 200); }

// =====================================================================
// Test 2: Short MCMC mixing diagnostic
//
// Run a brief MCMC (10k cycles) and report acceptance rate + whether
// the chain visits both signs.  Uses the standard |Omega| weight.
// =====================================================================

void run_mixing_diagnostic(int order, int n_cycles) {
  mpi::communicator world;

  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  bool bipartite = (order == 2);
  Parameters<double> params{U, beta, mu, t_intra, bipartite};

  std::unique_ptr<FreeEnergyCalculator<2, double>> calculator;
  if (order == 2) {
    calculator = std::make_unique<FreeEnergyCalculator<2, double>>(params, order);
  } else {
    calculator = std::make_unique<FreeEnergyCalculator<2, double>>(params, order, triangle_cluster, n_cluster_sites);
  }
  auto config = std::make_unique<DimerConfiguration<double>>(params, order, *calculator);

  int n_warmup     = 2000;
  int length_cycle = 1;
  int random_seed  = 55186222 + world.rank() * 786512 + order * 1000;
  int verbosity    = (world.rank() == 0 ? 2 : 0);

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 20;
  int block_size = (n_cycles / n_bins) + 1;

  if (world.rank() == 0) { std::filesystem::create_directory("./results"); }

  int measure_seed = 88871234 + world.rank() * 271828 + order * 999;
  measure_dimer<double> meas(config.get(), n_bins, block_size, mu, measure_seed);
  mc.add_move(move<double>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "dimer_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  if (world.rank() == 0) {
    std::string filename = "./results/dimer_data_order_" + std::to_string(order) + "_U_" + std::to_string(config->get_U()) + "_beta_"
       + std::to_string(config->beta) + "_mu_" + std::to_string(mu) + ".h5";

    double mean_sign = 0, sign_error = 0, mean_abs = 0, abs_error = 0;
    {
      h5::file file(filename, 'r');
      h5_read(file, "mean_sign", mean_sign);
      h5_read(file, "sign_error", sign_error);
      h5_read(file, "mean_abs_integrand", mean_abs);
      h5_read(file, "abs_integrand_error", abs_error);
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "  MCMC MIXING — Order " << order << " (" << n_cycles << " cycles x4 ranks)" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "  <sign>         = " << mean_sign << std::endl;
    std::cout << "  sign error     = " << sign_error << std::endl;
    std::cout << "  <|Omega|>_unif = " << mean_abs << std::endl;
    std::cout << "  |Omega| error  = " << abs_error << std::endl;
    std::cout << "  (check MC output above for acceptance rate)" << std::endl;
    std::cout << "========================================\n" << std::endl;

    std::filesystem::remove_all("./results");
  }
}

TEST(DimerDiagnostics, Mixing_Order2) { run_mixing_diagnostic(2, 50000); }
TEST(DimerDiagnostics, Mixing_Order4) { run_mixing_diagnostic(4, 10000); }
TEST(DimerDiagnostics, Mixing_Order6) { run_mixing_diagnostic(6, 5000); }

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
