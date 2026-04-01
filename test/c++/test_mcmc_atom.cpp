/*
 * MCMC test for AtomicConfiguration with the TRIQS MC framework.
 *
 * Uses AtomicConfiguration<double> with override_fm=1 (matching the 2-site
 * ED benchmark from analytical/benchmark_atomic_expansion.py) and the TRIQS
 * mc_generic loop with move/measure classes (as in apps/mcmc_atom.cpp).
 *
 * Checks that the MC estimate for the order-4 coefficient converges to
 * the exact value -0.019604213442906773.
 *
 * Usage:  mpirun -np 4 ./test_mcmc_atom
 *         (also works with 1 rank: ./test_mcmc_atom)
 */

#include "sc_expansion/atomic_configuration.hpp"
#include "sc_expansion/free_energy_calculator.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/measure.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <mpi/mpi.hpp>
#include <h5/h5.hpp>
#include <iostream>
#include <cmath>
#include <memory>
#include <filesystem>

int main(int argc, char **argv) {

  mpi::environment env(argc, argv);
  mpi::communicator world;

  // =====================================================================
  // Parameters (matching analytical/benchmark_atomic_expansion.py)
  // =====================================================================
  double U    = 8.0;
  double beta = 2.0;
  double mu   = 3.0;
  int order   = 4;

  double alpha     = 0.01;
  int override_fm  = 1;
  int n_cycles     = 100000;
  int n_warmup     = 2000;
  int length_cycle = 1;

  // Exact order-4 coefficient from the Python ED benchmark
  double exact_coeff = -0.019604213442906773;

  sc_expansion::Parameters<double> params{U, beta, mu, 0.0, true};

  // =====================================================================
  // Compute exact infinite-U reference integral on master rank, broadcast
  // =====================================================================
  double reference_integral        = 0.0;
  double signed_reference_integral = 0.0;

  if (world.rank() == 0) {
    sc_expansion::FreeEnergyCalculator<1, double> calculator(params, order, override_fm);
    auto [ref_abs, ref_signed] = calculator.compute_infinite_U_coefficient(false);
    reference_integral         = ref_abs;
    signed_reference_integral  = ref_signed;
  }

  mpi::broadcast(reference_integral, world);
  mpi::broadcast(signed_reference_integral, world);

  // =====================================================================
  // Construct AtomicConfiguration with override_fm=1
  // =====================================================================
  auto config = std::make_unique<AtomicConfiguration<double>>(params, order, alpha, override_fm);

  // =====================================================================
  // Set up TRIQS MC loop (same pattern as apps/mcmc_atom.cpp)
  // =====================================================================
  int random_seed = 32186222 + world.rank() * 786512;
  int verbosity   = (world.rank() == 0 ? 2 : 0);

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  // Create results directory so that measure writes h5 output
  if (world.rank() == 0) { std::filesystem::create_directory("./results"); }

  measure<double> meas(config.get(), reference_integral, signed_reference_integral, n_bins, block_size, mu);
  mc.add_move(move<double>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "defensive_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  // =====================================================================
  // Read back h5 result on master and verify against exact coefficient
  // =====================================================================
  if (world.rank() == 0) {
    std::string filename = "./results/full_lattice_data_order_" + std::to_string(order) + "_U_" + std::to_string(U) + "_beta_" + std::to_string(beta)
       + "_mu_" + std::to_string(mu) + "_bipartite.h5";

    double mc_mean  = 0.0;
    double mc_error = 0.0;

    {
      h5::file file(filename, 'r');
      h5_read(file, "mean", mc_mean);
      h5_read(file, "error", mc_error);
    }

    double rel_err = std::abs(mc_mean - exact_coeff) / std::abs(exact_coeff);

    std::cout << "=== MCMC AtomicConfiguration Test (order " << order << ", fm=" << override_fm << ") ===" << std::endl;
    std::cout << "MPI ranks:          " << world.size() << std::endl;
    std::cout << "Samples per rank:   " << n_cycles << std::endl;
    std::cout << "Exact (Python ED):  " << exact_coeff << std::endl;
    std::cout << "MC estimate:        " << mc_mean << std::endl;
    std::cout << "MC error:           " << mc_error << std::endl;
    std::cout << "Relative error:     " << rel_err << std::endl;

    // Cleanup
    std::filesystem::remove_all("./results");

    if (rel_err < 0.15) {
      std::cout << "PASS" << std::endl;
    } else {
      std::cerr << "FAIL: relative error " << rel_err << " exceeds 15% tolerance" << std::endl;
      return 1;
    }
  }

  return 0;
}
