/*
 * Vandermonde consistency test: dimer vs atomic expansion.
 *
 * The physical free energy per site F(t) can be computed two ways:
 *
 *   Atomic:  F(t) = sum_{m=0}^{N} a_m^atom * t^m
 *
 *   Dimer:   F(t) = G(t) + (1/2) * sum_{n=1}^{N} a_n^dimer(t) * t^n
 *
 * where G(t) = -ln(Z_2(t)) / (2*beta) is the isolated-dimer free energy
 * per site and a_n^dimer(t_0) depends on the intra-dimer hopping t_0.
 * At the physical point t_0 = t both must agree.
 *
 * Expanding G(t) = sum_m G_m * t^m and a_n(t_0) = sum_k p_{n,k} * t_0^k,
 * then setting t_0 = t and matching powers of t^m:
 *
 *   a_m^atom = G_m + (1/2) * sum_{n=1}^{m} p_{n, m-n}
 *
 * This test extracts G_m and p_{n,k} via Vandermonde interpolation at
 * 7 values of t_0, then checks the identity at even orders m = 2, 4, 6.
 * (Odd orders vanish on the bipartite square lattice.)
 *
 * Usage:  mpirun -np <N> ./test_compare_atom_dimer
 */

#include <gtest/gtest.h>
#include "sc_expansion/dimer_configuration.hpp"
#include "sc_expansion/atomic_configuration.hpp"
#include "sc_expansion/free_energy_calculator.hpp"
#include "sc_expansion/hubbard_solver.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/measure_dimer.hpp"
#include "sc_expansion/measure.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <mpi/mpi.hpp>
#include <cmath>
#include <memory>
#include <vector>
#include <iostream>
#include <iomanip>

using namespace sc_expansion;

// =====================================================================
// Gaussian elimination for a dense N x N system
// =====================================================================
static bool solve_linear_system(std::vector<std::vector<double>> A,
                                std::vector<double> b,
                                std::vector<double> &x) {
  int n = static_cast<int>(b.size());
  x.resize(n);

  for (int col = 0; col < n; col++) {
    int pivot = col;
    double max_val = std::abs(A[col][col]);
    for (int row = col + 1; row < n; row++) {
      if (std::abs(A[row][col]) > max_val) {
        max_val = std::abs(A[row][col]);
        pivot   = row;
      }
    }
    if (max_val < 1e-14) return false;

    std::swap(A[col], A[pivot]);
    std::swap(b[col], b[pivot]);

    for (int row = col + 1; row < n; row++) {
      double factor = A[row][col] / A[col][col];
      for (int j = col; j < n; j++) { A[row][j] -= factor * A[col][j]; }
      b[row] -= factor * b[col];
    }
  }

  for (int row = n - 1; row >= 0; row--) {
    x[row] = b[row];
    for (int j = row + 1; j < n; j++) { x[row] -= A[row][j] * x[j]; }
    x[row] /= A[row][row];
  }
  return true;
}

// =====================================================================
// Build Vandermonde matrix and solve for polynomial coefficients.
// Given values y_i = f(t_i), finds c_k such that f(t) ~ sum_k c_k t^k.
// =====================================================================
static bool vandermonde_fit(const std::vector<double> &t_values,
                            const std::vector<double> &y_values,
                            std::vector<double> &coeffs) {
  int K = static_cast<int>(t_values.size());
  std::vector<std::vector<double>> V(K, std::vector<double>(K));
  for (int i = 0; i < K; i++) {
    double t_pow = 1.0;
    for (int j = 0; j < K; j++) {
      V[i][j] = t_pow;
      t_pow *= t_values[i];
    }
  }
  return solve_linear_system(V, y_values, coeffs);
}

// =====================================================================
// Compute isolated-dimer free energy per site: G(t0) = -ln(Z_2(t0))/(2*beta)
// =====================================================================
static double dimer_free_energy_per_site(double U, double beta, double mu, double t0) {
  Parameters<double> params{U, beta, mu, t0, false};
  HubbardSolver<2, double> solver(params);
  return -std::log(solver.get_Z()) / (2.0 * beta);
}

// =====================================================================
// Helper: run dimer MCMC at a given order and t0, return coefficient
// =====================================================================
static DimerMeasureResult run_dimer_mc(mpi::communicator &world, double U, double beta, double mu,
                                       double t0, int order, int n_cycles, int seed_offset) {

  Parameters<double> params{U, beta, mu, t0, false};

  int n_warmup     = 5000;
  int length_cycle = 1;

  FreeEnergyCalculator<2, double> calculator(params, order);
  auto config = std::make_unique<DimerConfiguration<double>>(params, order, calculator);

  int random_seed  = 52186000 + seed_offset + world.rank() * 786512;
  int measure_seed = 88871000 + seed_offset + world.rank() * 271828;
  int verbosity    = 0;

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  measure_dimer<double> meas(config.get(), n_bins, block_size, mu, measure_seed);
  mc.add_move(move<double>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "dimer_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  return *meas.result;
}

// =====================================================================
// Helper: run atomic MCMC at a given order, return coefficient
// =====================================================================
static MeasureResult run_atomic_mc(mpi::communicator &world, double U, double beta, double mu,
                                   int order, int n_cycles, int seed_offset) {

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

  int random_seed = 72186222 + seed_offset + world.rank() * 786512;
  int verbosity   = 0;

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  measure<double> meas(config.get(), reference_integral, signed_reference_integral, n_bins, block_size, mu);
  mc.add_move(move<double>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "defensive_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  return *meas.result;
}

// =====================================================================
// Main test: Vandermonde consistency check up to order 6
// =====================================================================

TEST(VandermondeConsistency, DimerAtomicTaylorMatching) {
  mpi::communicator world;

  // Physics parameters
  double U    = 8.0;
  double beta = 2.0;
  double mu   = 3.0;

  int max_order = 6;
  int n_cycles  = 100000;

  // 7 interpolation points for degree-6 polynomial (Chebyshev nodes on [0, 1])
  int K = max_order + 1;
  std::vector<double> t0_values(K);
  for (int i = 0; i < K; i++) {
    t0_values[i] = 0.5 * (1.0 + std::cos(M_PI * (2 * i + 1) / (2.0 * K)));
  }

  // =================================================================
  // Step 1: Dimer free energy G(t0) Taylor coefficients (analytic)
  // =================================================================
  std::vector<double> G_coeffs;

  if (world.rank() == 0) {
    std::vector<double> G_values(K);
    for (int i = 0; i < K; i++) {
      G_values[i] = dimer_free_energy_per_site(U, beta, mu, t0_values[i]);
    }
    vandermonde_fit(t0_values, G_values, G_coeffs);

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "\n========================================" << std::endl;
    std::cout << "  G(t0) Taylor coefficients (per site)" << std::endl;
    std::cout << "========================================" << std::endl;
    for (int m = 0; m <= max_order; m++) {
      std::cout << "  G_" << m << " = " << G_coeffs[m] << std::endl;
    }
  }

  // Broadcast G_coeffs to all ranks
  if (world.rank() != 0) { G_coeffs.resize(K); }
  MPI_Bcast(G_coeffs.data(), K, MPI_DOUBLE, 0, world.get());

  // =================================================================
  // Step 2: Dimer expansion polynomial coefficients p_{n,k}
  //         For each even order n=2,4,...,max_order, run at K t0 values.
  // =================================================================

  // p_coeffs[n][k] = coefficient of t0^k in a_n^dimer(t0) / 2 (per site)
  // Odd orders vanish (bipartite square lattice), so only n = 2, 4, 6.
  std::vector<std::vector<double>> p_coeffs(max_order + 1, std::vector<double>(K, 0.0));

  for (int n = 2; n <= max_order; n += 2) {
    if (world.rank() == 0) {
      std::cout << "\n--- Dimer order " << n << ": running " << K << " t0 values ---" << std::endl;
    }

    std::vector<double> dimer_values(K);
    for (int i = 0; i < K; i++) {
      if (world.rank() == 0) {
        std::cout << "  t0 = " << t0_values[i] << " (run " << (i + 1) << "/" << K << ")" << std::endl;
      }
      int seed_offset = n * 10000 + i * 1000;
      auto result = run_dimer_mc(world, U, beta, mu, t0_values[i], n, n_cycles, seed_offset);
      // Store per-site coefficient
      dimer_values[i] = result.coeff / 2.0;
    }

    // Vandermonde fit (all ranks have the same values via mpi broadcast in collect_results)
    vandermonde_fit(t0_values, dimer_values, p_coeffs[n]);

    if (world.rank() == 0) {
      std::cout << "  Polynomial coefficients p_{" << n << ",k}:" << std::endl;
      for (int k = 0; k < K; k++) {
        std::cout << "    p_{" << n << "," << k << "} = " << p_coeffs[n][k] << std::endl;
      }
    }
  }

  // =================================================================
  // Step 3: Atomic expansion coefficients a_m^atom
  //         Only even orders (odd orders vanish on bipartite lattice).
  // =================================================================
  std::vector<double> a_atom(max_order + 1, 0.0);

  for (int m = 2; m <= max_order; m += 2) {
    if (world.rank() == 0) {
      std::cout << "\n--- Atomic order " << m << " ---" << std::endl;
    }
    int seed_offset = m * 20000;
    auto result = run_atomic_mc(world, U, beta, mu, m, n_cycles, seed_offset);
    a_atom[m] = result.mean;
  }

  // =================================================================
  // Step 4: Check consistency identity (even orders only)
  //   a_m^atom = G_m + sum_{n=2,4,...,m} p_{n, m-n}
  // =================================================================
  if (world.rank() == 0) {
    std::cout << "\n========================================================" << std::endl;
    std::cout << "  Consistency check: a_m^atom = G_m + sum p_{n, m-n}" << std::endl;
    std::cout << "  U = " << U << ", beta = " << beta << ", mu = " << mu << std::endl;
    std::cout << "  MC cycles per rank: " << n_cycles << std::endl;
    std::cout << "  MPI ranks: " << world.size() << std::endl;
    std::cout << "========================================================" << std::endl;

    for (int m = 2; m <= max_order; m += 2) {
      double dimer_sum = G_coeffs[m];
      for (int n = 2; n <= m; n += 2) { dimer_sum += p_coeffs[n][m - n]; }

      std::cout << "  m = " << m << ":  dimer_side = " << std::setw(16) << dimer_sum
                << "  atom_side = " << std::setw(16) << a_atom[m]
                << "  diff = " << std::setw(14) << (dimer_sum - a_atom[m]) << std::endl;

      double diff  = std::abs(dimer_sum - a_atom[m]);
      double scale = std::max(std::abs(a_atom[m]), 1e-10);

      EXPECT_LT(diff / scale, 0.20)
         << "Order " << m << " mismatch: dimer_side=" << dimer_sum
         << " atom_side=" << a_atom[m] << " rel_diff=" << diff / scale;
    }

    std::cout << "========================================================\n" << std::endl;
  }
}

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
