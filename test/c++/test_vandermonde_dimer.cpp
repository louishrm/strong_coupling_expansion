/*
 * Vandermonde consistency test: dimer order-6 coefficient vs atomic expansion.
 *
 * Computes a_6^dimer(t0) via MCMC for t0 = {0, 0.2, 0.4, 0.6, 0.8, 1.0},
 * then uses Vandermonde interpolation to extract polynomial coefficients c_k
 * such that
 *   a_6^dimer(t0) ≈ c_0 + c_1 * t0 + ... + c_5 * t0^5
 *
 * Also runs the atomic (N_sites=1) MCMC at order 6 for comparison.
 * The constant term c_0 = a_6^dimer(0) corresponds to the decoupled-dimer
 * limit and should relate to the atomic expansion on the triangular
 * superlattice formed by the staggered dimer tiling.
 *
 * Designed for cluster execution (high cycle counts).
 *
 * Usage:  mpirun -np 4 ./test_vandermonde_dimer
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
#include <vector>
#include <iostream>
#include <iomanip>

using namespace sc_expansion;

// =====================================================================
// Gaussian elimination for a dense N x N system  (column-major safe)
// =====================================================================
static bool solve_linear_system(std::vector<std::vector<double>> A,
                                std::vector<double> b,
                                std::vector<double> &x) {
  int n = static_cast<int>(b.size());
  x.resize(n);

  // Forward elimination with partial pivoting
  for (int col = 0; col < n; col++) {
    // Find pivot
    int pivot = col;
    double max_val = std::abs(A[col][col]);
    for (int row = col + 1; row < n; row++) {
      if (std::abs(A[row][col]) > max_val) {
        max_val = std::abs(A[row][col]);
        pivot   = row;
      }
    }
    if (max_val < 1e-14) return false; // singular

    std::swap(A[col], A[pivot]);
    std::swap(b[col], b[pivot]);

    for (int row = col + 1; row < n; row++) {
      double factor = A[row][col] / A[col][col];
      for (int j = col; j < n; j++) { A[row][j] -= factor * A[col][j]; }
      b[row] -= factor * b[col];
    }
  }

  // Back substitution
  for (int row = n - 1; row >= 0; row--) {
    x[row] = b[row];
    for (int j = row + 1; j < n; j++) { x[row] -= A[row][j] * x[j]; }
    x[row] /= A[row][row];
  }
  return true;
}

// =====================================================================
// Helper: run dimer MCMC at order 4 for a given t0, return coefficient
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
                                  double t0, int order, int n_cycles, int seed_offset) {

  Parameters<double> params{U, beta, mu, t0, false};

  int n_warmup     = 5000;
  int length_cycle = 1;

  FreeEnergyCalculator<2, double> dimer_calculator(params, order);
  auto config = std::make_unique<DimerConfiguration<double>>(params, order, dimer_calculator);

  int random_seed = 52186000 + seed_offset + world.rank() * 786512;
  int verbosity   = 0; // quiet for batch runs

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  if (world.rank() == 0) { std::filesystem::create_directory("./results"); }

  int measure_seed = 88871000 + seed_offset + world.rank() * 271828;
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

    // Clean up h5 so next run doesn't collide
    std::filesystem::remove(filename);
  }

  // Broadcast result to all ranks
  mpi::broadcast(result.coeff, world);
  mpi::broadcast(result.error, world);
  mpi::broadcast(result.mean_sign, world);
  mpi::broadcast(result.sign_error, world);
  mpi::broadcast(result.mean_abs, world);
  mpi::broadcast(result.abs_error, world);

  return result;
}

// =====================================================================
// Helper: run atomic MCMC at order 4, return coefficient
// =====================================================================
struct AtomicMCResult {
  double coeff;
  double error;
};

static AtomicMCResult run_atomic_mc(mpi::communicator &world, double U, double beta, double mu,
                                    int order, int n_cycles) {

  Parameters<double> params{U, beta, mu, 0.0, true};

  double alpha    = 0.01;
  int override_fm = -1; // full lattice multiplicities
  int n_warmup    = 2000;
  int length_cycle = 1;

  // Construct calculator once, compute reference integral on master, broadcast
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
// Main test
// =====================================================================

TEST(VandermondeDimer, Order6PolynomialDecomposition) {
  mpi::communicator world;

  // Physics parameters
  double U    = 8.0;
  double beta = 2.0;
  double mu   = 3.0;
  int order   = 6;

  // t0 values for Vandermonde interpolation (6 points for degree-5 polynomial)
  std::vector<double> t0_values = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
  int n_points = static_cast<int>(t0_values.size());

  int n_cycles_dimer = 100000;
  int n_cycles_atom  = 100000;

  // ---------------------------------------------------------------
  // Step 1: Run dimer MCMC for each t0
  // ---------------------------------------------------------------
  std::vector<DimerMCResult> dimer_results(n_points);

  for (int i = 0; i < n_points; i++) {
    if (world.rank() == 0) {
      std::cout << "\n=== Dimer MC: t0 = " << t0_values[i]
                << " (run " << (i + 1) << "/" << n_points << ") ===" << std::endl;
    }
    dimer_results[i] = run_dimer_mc(world, U, beta, mu, t0_values[i], order,
                                    n_cycles_dimer, i * 1000);
  }

  // ---------------------------------------------------------------
  // Step 2: Run atomic MCMC for reference
  // ---------------------------------------------------------------
  if (world.rank() == 0) {
    std::cout << "\n=== Atomic MC: order " << order << " ===" << std::endl;
  }
  AtomicMCResult atom_result = run_atomic_mc(world, U, beta, mu, order, n_cycles_atom);

  // ---------------------------------------------------------------
  // Step 3: Vandermonde interpolation (rank 0 only for output)
  // ---------------------------------------------------------------
  if (world.rank() == 0) {

    // Build Vandermonde matrix: V[i][j] = t0[i]^j
    std::vector<std::vector<double>> V(n_points, std::vector<double>(n_points));
    std::vector<double> a(n_points);

    for (int i = 0; i < n_points; i++) {
      double t0_pow = 1.0;
      for (int j = 0; j < n_points; j++) {
        V[i][j] = t0_pow;
        t0_pow *= t0_values[i];
      }
      a[i] = dimer_results[i].coeff;
    }

    std::vector<double> c; // polynomial coefficients
    bool ok = solve_linear_system(V, a, c);
    ASSERT_TRUE(ok) << "Vandermonde system is singular";

    // ---------------------------------------------------------------
    // Step 4: Print results
    // ---------------------------------------------------------------
    std::cout << std::fixed << std::setprecision(10);

    std::cout << "\n========================================" << std::endl;
    std::cout << "  Dimer order-6 coefficients a_6(t0)" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "  U = " << U << ", beta = " << beta << ", mu = " << mu << std::endl;
    std::cout << "  MC cycles (dimer): " << n_cycles_dimer << std::endl;
    std::cout << "  MC cycles (atom):  " << n_cycles_atom << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    for (int i = 0; i < n_points; i++) {
      std::cout << "  t0 = " << std::setw(6) << t0_values[i]
                << "  =>  a_4 = " << std::setw(16) << dimer_results[i].coeff
                << "  +/- " << std::setw(12) << dimer_results[i].error
                << "  (sign = " << dimer_results[i].mean_sign << ")" << std::endl;
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "  Vandermonde polynomial coefficients" << std::endl;
    std::cout << "  a_6(t0) = c_0 + c_1*t0 + ... + c_5*t0^5" << std::endl;
    std::cout << "========================================" << std::endl;

    for (int k = 0; k < n_points; k++) {
      std::cout << "  c_" << k << " = " << std::setw(16) << c[k] << std::endl;
    }

    // Verify interpolation: evaluate polynomial at each t0
    std::cout << "\n  Interpolation check:" << std::endl;
    for (int i = 0; i < n_points; i++) {
      double poly_val = 0.0;
      double t0_pow   = 1.0;
      for (int k = 0; k < n_points; k++) {
        poly_val += c[k] * t0_pow;
        t0_pow *= t0_values[i];
      }
      std::cout << "  t0 = " << std::setw(6) << t0_values[i]
                << "  poly = " << std::setw(16) << poly_val
                << "  MC = " << std::setw(16) << dimer_results[i].coeff
                << "  diff = " << std::setw(12) << (poly_val - dimer_results[i].coeff) << std::endl;
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "  Atomic expansion comparison" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "  a_6^atom (square lattice, per site) = " << atom_result.coeff
              << "  +/- " << atom_result.error << std::endl;
    std::cout << "  c_0 = a_6^dimer(t0=0) (triang. superlattice, per dimer) = " << c[0] << std::endl;
    std::cout << "  c_0 / 2 (per site)  = " << c[0] / 2.0 << std::endl;
    std::cout << "========================================\n" << std::endl;

    // Cleanup results directory
    std::filesystem::remove_all("./results");
  }
}

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
