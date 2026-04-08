/*
 * Dimer-atomic consistency test at t0 = 0.
 *
 * At t0 = 0 the intra-dimer hopping vanishes, so each dimer is two
 * decoupled atomic sites.  The dimer expansion on the triangular
 * superlattice should then reduce to the atomic expansion:
 *
 *   a_n^dimer(t0=0) / 2  ==  a_n^atom   (per site)
 *
 * Usage:  mpirun -np <N> ./test_vandermonde_dimer
 */

#include <gtest/gtest.h>
#include "sc_expansion/dimer_configuration.hpp"
#include "sc_expansion/atomic_configuration.hpp"
#include "sc_expansion/free_energy_calculator.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/measure_dimer.hpp"
#include "sc_expansion/measure.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <mpi/mpi.hpp>
#include <h5/h5.hpp>
#include <cmath>
#include <memory>
#include <filesystem>
#include <iostream>
#include <iomanip>

using namespace sc_expansion;

// =====================================================================
// Helper: run dimer MCMC at a given order with t0, return coefficient
// =====================================================================
struct DimerMCResult {
  double coeff;
  double error;
  double mean_sign;
  double sign_error;
  double mean_abs;
  double abs_error;
};

static DimerMCResult run_dimer_mc(mpi::communicator &world, double U, double beta, double mu,
                                  double t0, int order, int n_cycles) {

  Parameters<double> params{U, beta, mu, t0, false};

  int n_warmup     = 5000;
  int length_cycle = 1;

  FreeEnergyCalculator<2, double> dimer_calculator(params, order);
  auto config = std::make_unique<DimerConfiguration<double>>(params, order, dimer_calculator);

  int random_seed = 52186000 + world.rank() * 786512;
  int verbosity   = 0;

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  if (world.rank() == 0) { std::filesystem::create_directory("./results"); }

  int measure_seed = 88871000 + world.rank() * 271828;
  measure_dimer<double> meas(config.get(), n_bins, block_size, mu, measure_seed);
  mc.add_move(move<double>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "dimer_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  DimerMCResult result{};

  if (world.rank() == 0) {
    std::string filename = "./results/dimer_data_order_" + std::to_string(order) + "_U_"
       + std::to_string(config->get_U()) + "_beta_" + std::to_string(config->beta)
       + "_mu_" + std::to_string(mu) + ".h5";

    {
      h5::file file(filename, 'r');
      h5_read(file, "mean_sign", result.mean_sign);
      h5_read(file, "sign_error", result.sign_error);
      h5_read(file, "mean_abs_integrand", result.mean_abs);
      h5_read(file, "abs_integrand_error", result.abs_error);
    }

    double abs_integral = std::pow(beta, order) * result.mean_abs;
    result.coeff        = abs_integral * result.mean_sign;
    result.error        = std::pow(beta, order)
       * std::sqrt(result.mean_abs * result.mean_abs * result.sign_error * result.sign_error
                   + result.mean_sign * result.mean_sign * result.abs_error * result.abs_error);

    std::filesystem::remove(filename);
  }

  mpi::broadcast(result.coeff, world);
  mpi::broadcast(result.error, world);
  mpi::broadcast(result.mean_sign, world);
  mpi::broadcast(result.sign_error, world);
  mpi::broadcast(result.mean_abs, world);
  mpi::broadcast(result.abs_error, world);

  return result;
}

// =====================================================================
// Helper: run atomic MCMC at a given order, return coefficient
// =====================================================================
struct AtomicMCResult {
  double coeff;
  double error;
};

static AtomicMCResult run_atomic_mc(mpi::communicator &world, double U, double beta, double mu,
                                    int order, int n_cycles) {

  Parameters<double> params{U, beta, mu, 0.0, true};

  double alpha     = 0.01;
  int override_fm  = -1;
  int n_warmup     = 2000;
  int length_cycle = 1;

  FreeEnergyCalculator<1, double> calculator(params, order, override_fm);

  double reference_integral        = 0.0;
  double signed_reference_integral = 0.0;

  if (world.rank() == 0) {
    auto [ref_abs, ref_signed] = calculator.compute_infinite_U_coefficient(false);
    reference_integral         = ref_abs;
    signed_reference_integral  = ref_signed;
  }

  mpi::broadcast(reference_integral, world);
  mpi::broadcast(signed_reference_integral, world);

  auto config = std::make_unique<AtomicConfiguration<double>>(params, order, alpha, calculator);

  int random_seed = 72186222 + world.rank() * 786512;
  int verbosity   = 0;

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  if (world.rank() == 0) { std::filesystem::create_directory("./results"); }

  measure<double> meas(config.get(), reference_integral, signed_reference_integral, n_bins, block_size, mu);
  mc.add_move(move<double>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "defensive_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  AtomicMCResult result{};

  if (world.rank() == 0) {
    std::string filename = "./results/full_lattice_data_order_" + std::to_string(order)
       + "_U_" + std::to_string(U) + "_beta_" + std::to_string(beta)
       + "_mu_" + std::to_string(mu) + "_bipartite.h5";

    {
      h5::file file(filename, 'r');
      h5_read(file, "mean", result.coeff);
      h5_read(file, "error", result.error);
    }

    std::filesystem::remove(filename);
  }

  mpi::broadcast(result.coeff, world);
  mpi::broadcast(result.error, world);

  return result;
}

// =====================================================================
// Test: dimer at t0=0 vs atomic expansion
// =====================================================================

TEST(DimerAtomicConsistency, Order6_t0_zero) {
  mpi::communicator world;

  double U    = 8.0;
  double beta = 2.0;
  double mu   = 3.0;
  int order   = 6;

  int n_cycles = 100000;

  // --- Dimer MCMC at t0 = 0 ---
  if (world.rank() == 0) {
    std::cout << "\n=== Dimer MC: order " << order << ", t0 = 0 ===" << std::endl;
  }
  DimerMCResult dimer = run_dimer_mc(world, U, beta, mu, 0.0, order, n_cycles);

  // --- Atomic MCMC ---
  if (world.rank() == 0) {
    std::cout << "\n=== Atomic MC: order " << order << " ===" << std::endl;
  }
  AtomicMCResult atom = run_atomic_mc(world, U, beta, mu, order, n_cycles);

  // --- Comparison ---
  double dimer_per_site       = dimer.coeff / 2.0;
  double dimer_error_per_site = dimer.error / 2.0;

  if (world.rank() == 0) {
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "\n========================================" << std::endl;
    std::cout << "  Order " << order << " consistency check (t0 = 0)" << std::endl;
    std::cout << "  U = " << U << ", beta = " << beta << ", mu = " << mu << std::endl;
    std::cout << "  MC cycles per rank: " << n_cycles << std::endl;
    std::cout << "  MPI ranks: " << world.size() << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "  a_" << order << "^dimer(t0=0) / 2  = " << dimer_per_site
              << "  +/- " << dimer_error_per_site << std::endl;
    std::cout << "  a_" << order << "^atom             = " << atom.coeff
              << "  +/- " << atom.error << std::endl;

    double diff       = std::abs(dimer_per_site - atom.coeff);
    double combined_e = std::sqrt(dimer_error_per_site * dimer_error_per_site + atom.error * atom.error);
    double n_sigma    = (combined_e > 0) ? diff / combined_e : 0.0;

    std::cout << "  |diff| = " << diff << ",  combined error = " << combined_e
              << ",  n_sigma = " << n_sigma << std::endl;
    std::cout << "========================================\n" << std::endl;

    EXPECT_LT(n_sigma, 5.0) << "Dimer (t0=0) and atomic results differ by more than 5 sigma";

    std::filesystem::remove_all("./results");
  }
}

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
