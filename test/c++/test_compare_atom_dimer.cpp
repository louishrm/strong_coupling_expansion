/*
 * Consistency test: dimer vs atomic expansion.
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
 * Rather than extracting individual Taylor coefficients via Vandermonde
 * interpolation (which is extremely ill-conditioned for the monomial
 * basis on [0,1]), this test compares the truncated free energy directly
 * at a specific small t value.  The identity
 *
 *   a_m^atom = G_m + (1/2) * sum_{n=1}^{m} p_{n, m-n}
 *
 * holding at every order implies the function-value identity
 *
 *   sum_{m=0}^{M} a_m^atom * t^m = G(t) + (1/2) * sum_{n} a_n^dimer(t) * t^n + O(t^{M+2})
 *
 * For t = 0.3 and M = 6, the truncation error is O(t^8) ~ 7e-5.
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
// Compute isolated-atom free energy per site: -ln(Z_1)/(beta)
// =====================================================================
static double atom_free_energy_per_site(double U, double beta, double mu) {
  Parameters<double> params{U, beta, mu, 0.0, true};
  HubbardSolver<1, double> solver(params);
  return -std::log(solver.get_Z()) / beta;
}

// =====================================================================
// Main test: function-value consistency check up to order 6
// =====================================================================

TEST(ConsistencyCheck, DimerAtomicFunctionValue) {
  mpi::communicator world;

  // Physics parameters
  double U    = 8.0;
  double beta = 2.0;
  double mu   = 3.0;
  double t    = 0.3;   // Physical hopping for function-value comparison

  int max_order = 6;
  int n_cycles  = 100000;

  // =================================================================
  // Step 1: Analytic zeroth-order contributions
  // =================================================================
  double a0_atom = atom_free_energy_per_site(U, beta, mu);             // F_atom(0)
  double G_at_t  = dimer_free_energy_per_site(U, beta, mu, t);         // G(t)
  double G_at_0  = dimer_free_energy_per_site(U, beta, mu, 0.0);       // G(0) = a0_atom

  if (world.rank() == 0) {
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "\n========================================" << std::endl;
    std::cout << "  Analytic zeroth-order (per site)" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "  a0_atom  = " << a0_atom << std::endl;
    std::cout << "  G(0)     = " << G_at_0 << std::endl;
    std::cout << "  G(t=" << t << ") = " << G_at_t << std::endl;
  }

  // =================================================================
  // Step 2: Atomic expansion coefficients a_m^atom (even m = 2, 4, 6)
  // =================================================================
  double atom_sum = a0_atom;
  std::vector<double> a_atom(max_order + 1, 0.0);

  for (int m = 2; m <= max_order; m += 2) {
    if (world.rank() == 0) {
      std::cout << "\n--- Atomic order " << m << " ---" << std::endl;
    }
    int seed_offset = m * 20000;
    auto result = run_atomic_mc(world, U, beta, mu, m, n_cycles, seed_offset);
    a_atom[m] = result.mean;
    atom_sum += a_atom[m] * std::pow(t, m);
  }

  // =================================================================
  // Step 3: Dimer expansion coefficients a_n^dimer(t) at t0 = t
  //         Only one MC run per order (no Vandermonde interpolation).
  // =================================================================
  double dimer_sum = G_at_t;
  std::vector<double> c_dimer(max_order + 1, 0.0);

  for (int n = 2; n <= max_order; n += 2) {
    if (world.rank() == 0) {
      std::cout << "\n--- Dimer order " << n << " at t0 = " << t << " ---" << std::endl;
    }
    int seed_offset = n * 10000;
    auto result = run_dimer_mc(world, U, beta, mu, t, n, n_cycles, seed_offset);
    // Per-site contribution: (1/2) * a_n^dimer(t) * t^n
    c_dimer[n] = result.coeff / 2.0;
    dimer_sum += c_dimer[n] * std::pow(t, n);
  }

  // =================================================================
  // Step 4: Compare truncated free energies
  //   atom_sum  = a0 + sum_{m=2,4,6} a_m^atom * t^m
  //   dimer_sum = G(t) + (1/2) * sum_{n=2,4,6} a_n^dimer(t) * t^n
  //   These agree up to O(t^8) truncation error.
  // =================================================================
  if (world.rank() == 0) {
    std::cout << "\n========================================================" << std::endl;
    std::cout << "  Consistency check: F_atom(t) vs F_dimer(t)" << std::endl;
    std::cout << "  U = " << U << ", beta = " << beta << ", mu = " << mu << ", t = " << t << std::endl;
    std::cout << "  MC cycles per rank: " << n_cycles << std::endl;
    std::cout << "  MPI ranks: " << world.size() << std::endl;
    std::cout << "========================================================" << std::endl;

    std::cout << "\n  Atomic coefficients:" << std::endl;
    for (int m = 2; m <= max_order; m += 2) {
      std::cout << "    a_" << m << " = " << a_atom[m] << std::endl;
    }

    std::cout << "\n  Dimer coefficients (per site, at t0 = " << t << "):" << std::endl;
    for (int n = 2; n <= max_order; n += 2) {
      std::cout << "    c_" << n << "/2 = " << c_dimer[n] << std::endl;
    }

    // Compare the perturbative corrections (subtract the common zeroth order)
    double delta_atom  = atom_sum - a0_atom;
    double delta_dimer = dimer_sum - G_at_0;

    std::cout << "\n  atom_sum   = " << atom_sum << std::endl;
    std::cout << "  dimer_sum  = " << dimer_sum << std::endl;
    std::cout << "  delta_atom  (atom_sum - a0)  = " << delta_atom << std::endl;
    std::cout << "  delta_dimer (dimer_sum - G0) = " << delta_dimer << std::endl;
    std::cout << "  diff = " << (delta_atom - delta_dimer) << std::endl;

    double diff  = std::abs(delta_atom - delta_dimer);
    double scale = std::max(std::abs(delta_atom), 1e-10);

    EXPECT_LT(diff / scale, 0.20)
       << "Function-value mismatch at t=" << t
       << ": delta_atom=" << delta_atom
       << " delta_dimer=" << delta_dimer
       << " rel_diff=" << diff / scale;

    std::cout << "  relative diff = " << diff / scale << std::endl;
    std::cout << "========================================================\n" << std::endl;
  }
}

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
