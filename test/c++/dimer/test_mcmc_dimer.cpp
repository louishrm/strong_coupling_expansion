/*
 * MCMC test for sc_expansion::dimer::Configuration on a finite 3-dimer cluster.
 *
 * 3-dimer cluster (6 sites) — triangle on the staggered superlattice.
 * Dimers in superlattice (u,v) coords:  A=(0,0), B=(1,0), C=(0,1).
 * Inter-dimer bonds:  A<->B, A<->C, B<->C.
 *
 * Coefficients are compared against exact ED reference values supplied in the
 * task description.
 *
 * Usage:  mpirun -np 4 ./test_mcmc_dimer
 */

#include <gtest/gtest.h>
#include "sc_expansion/dimer/configuration.hpp"
#include "sc_expansion/dimer/free_energy_calculator.hpp"
#include "sc_expansion/dimer/measure_dimer.hpp"
#include "sc_expansion/hubbard_solver.hpp"
#include "sc_expansion/move.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <mpi/mpi.hpp>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

static void run_mcmc_dimer_check(int order, double exact_coeff, double rel_tol, int n_cycles) {
  mpi::communicator world;

  double U    = 8.0;
  double beta = 2.0;
  double mu   = 3.0;
  double t    = 1.0; // intra-dimer hopping

  // 3-dimer triangular cluster on the staggered superlattice. The dimer
  // superlattice is non-bipartite (each dimer has 6 NN), so bipartite=false
  // is required to include odd-order diagrams. n_cluster_sites is the
  // per-dimer normaliser (number of dimers in the cluster), not physical sites.
  std::vector<std::pair<int, int>> cluster_positions = {{0, 0}, {1, 0}, {0, 1}};
  int n_cluster_sites                                = 3;

  double alpha     = 0.001;
  int n_warmup     = 2000;
  int length_cycle = 1;

  sc_expansion::Parameters<double> params{U, beta, mu, t, /*bipartite=*/false};
  sc_expansion::dimer::FreeEnergyCalculator<double> calculator(params, order, cluster_positions, n_cluster_sites);

  auto config = std::make_unique<sc_expansion::dimer::Configuration<double>>(params, order, calculator, alpha);

  int random_seed = 32186222 + world.rank() * 786512;
  int verbosity   = (world.rank() == 0 ? 2 : 0);

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  measure_dimer<double> meas(config.get(), n_bins, block_size, mu);
  mc.add_move(move<double>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "dimer_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  if (world.rank() == 0) {
    double mc_mean  = meas.result->coeff;
    double mc_error = meas.result->error;
    double rel_err  = std::abs(mc_mean - exact_coeff) / std::abs(exact_coeff);

    std::cout << "Exact (ED):     " << exact_coeff << std::endl;
    std::cout << "MC estimate:    " << mc_mean << std::endl;
    std::cout << "MC error:       " << mc_error << std::endl;
    std::cout << "Relative error: " << rel_err << std::endl;

    EXPECT_LT(rel_err, rel_tol) << "MC estimate " << mc_mean << " deviates from exact " << exact_coeff << " by " << rel_err * 100 << "%";
  }
}

TEST(McmcDimerCluster, Order2Coefficient) {
  run_mcmc_dimer_check(/*order=*/2, /*exact_coeff=*/-0.132935639042, /*rel_tol=*/0.10, /*n_cycles=*/100000);
}

TEST(McmcDimerCluster, Order3Coefficient) {
  run_mcmc_dimer_check(/*order=*/3, /*exact_coeff=*/0.001298542799, /*rel_tol=*/0.15, /*n_cycles=*/300000);
}

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
