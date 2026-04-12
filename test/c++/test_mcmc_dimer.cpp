/*
 * Tests for the columnar dimer (N_sites=2) expansion at orders 2 and 4.
 *
 * Cluster geometry: 3 dimers on the rectangular (columnar) superlattice.
 *   A=(0,0), B=(0,1), C=(1,0)
 *   A--B vertical (2 bonds), A--C horizontal (1 bond)
 *
 * Tests:
 *   1. Order 2 quadrature: 1D trapezoidal rule on the 3-dimer cluster.
 *   2. Order 2 MCMC: cluster-restricted sign-based sampling.
 *   3. Order 4 MCMC: cluster-restricted sign-based sampling.
 *
 * ED references from analytical/benchmark_columnar_dimer_expansion.py:
 *   Order 2 coefficient on the 3-dimer cluster = -0.130006787492
 *   Order 4 coefficient on the 3-dimer cluster = -0.037169097044
 *
 * Usage:  mpirun -np 4 ./test_mcmc_dimer
 */

#include <gtest/gtest.h>
#include "sc_expansion/diagram.hpp"
#include "sc_expansion/graph.hpp"
#include "sc_expansion/dimer_configuration.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/measure_dimer.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <mpi/mpi.hpp>
#include <cmath>
#include <memory>

using namespace sc_expansion;

// =====================================================================
// MCMC test: order 2, cluster-restricted sign-based sampling
//
// Uses FreeEnergyCalculator with cluster-restricted embedding on the
// 3-dimer L-shaped cluster. The integrand already has per-dimer weights.
//
// coefficient = beta^n * <|Omega|>_uniform * <sign>
// =====================================================================

TEST(DimerExpansion, Order2MCMC) {
  mpi::communicator world;

  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order = 2;
  Parameters<double> params{U, beta, mu, t_intra, true};

  // 3-dimer L-shaped cluster on the rectangular superlattice
  std::vector<std::pair<int, int>> cluster_positions = {{0, 0}, {0, 1}, {1, 0}};
  int n_cluster_sites                                = 3;

  // --- MCMC: cluster-restricted ---
  int n_cycles     = 100000;
  int n_warmup     = 2000;
  int length_cycle = 1;

  FreeEnergyCalculator<2, double> calculator(params, order, cluster_positions, n_cluster_sites);
  auto config = std::make_unique<DimerConfiguration<double>>(params, order, calculator);

  int random_seed = 32186222 + world.rank() * 786512;
  int verbosity   = (world.rank() == 0 ? 2 : 0);

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  int measure_seed = 99871234 + world.rank() * 314159;
  measure_dimer<double> meas(config.get(), n_bins, block_size, mu, measure_seed);
  mc.add_move(move<double>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "dimer_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  if (world.rank() == 0) {
    // coefficient = beta^n * <|Omega|>_uniform * <sign>  (no fm division — cluster weights)
    double abs_integral = std::pow(beta, order) * meas.result->mean_abs;
    double mc_coeff     = abs_integral * meas.result->mean_sign;
    double exact        = -0.130006787492;
    double rel_err      = std::abs(mc_coeff - exact) / std::abs(exact);

    std::cout << "Exact (Python ED):       " << exact << std::endl;
    std::cout << "MC coefficient:          " << mc_coeff << std::endl;
    std::cout << "|Omega| integral (MC):   " << abs_integral << std::endl;
    std::cout << "Mean |Omega| (uniform):  " << meas.result->mean_abs << std::endl;
    std::cout << "Mean sign:               " << meas.result->mean_sign << std::endl;
    std::cout << "Sign error:              " << meas.result->sign_error << std::endl;
    std::cout << "Relative error:          " << rel_err << std::endl;

    EXPECT_LT(rel_err, 0.025) << "MC estimate " << mc_coeff << " deviates from exact " << exact << " by " << rel_err * 100 << "%";
  }
}

// =====================================================================
// MCMC test: order 4 on the 3-dimer L-shaped cluster
//
// Cluster: A=(0,0), B=(0,1), C=(1,0) on the rectangular superlattice.
// Inter-dimer bonds:
//   A--B vertical: 2 bonds (site0-site0, site1-site1)
//   A--C horizontal: 1 bond (A.site1 <-> C.site0)
//
// DimerConfiguration uses cluster-restricted spatial embeddings, so
// the integrand already has per-dimer weights — no fm division needed.
//
// ED reference from analytical/benchmark_columnar_dimer_expansion.py:
//   Order 4 coefficient on the 3-dimer cluster = -0.037169097044
// =====================================================================

TEST(DimerExpansion, Order4MCMC_3DimerCluster) {
  mpi::communicator world;

  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order = 4;
  Parameters<double> params{U, beta, mu, t_intra, true}; // bipartite (rectangular superlattice)

  // 3-dimer L-shaped cluster positions on the rectangular superlattice
  std::vector<std::pair<int, int>> cluster_positions = {{0, 0}, {0, 1}, {1, 0}};
  int n_cluster_sites                                = 3;

  // --- MCMC: cluster-restricted DimerConfiguration ---
  int n_cycles     = 200000;
  int n_warmup     = 5000;
  int length_cycle = 1;

  FreeEnergyCalculator<2, double> calculator(params, order, cluster_positions, n_cluster_sites);
  auto config = std::make_unique<DimerConfiguration<double>>(params, order, calculator);

  int random_seed = 42186333 + world.rank() * 786512;
  int verbosity   = (world.rank() == 0 ? 2 : 0);

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  int measure_seed = 77871234 + world.rank() * 271828;
  measure_dimer<double> meas(config.get(), n_bins, block_size, mu, measure_seed);
  mc.add_move(move<double>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "dimer_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  if (world.rank() == 0) {
    // Cluster per-dimer weights are already in the integrand — no fm division
    double abs_integral = std::pow(beta, order) * meas.result->mean_abs;
    double mc_coeff     = abs_integral * meas.result->mean_sign;
    double exact        = -0.037169097044;
    double rel_err      = std::abs(mc_coeff - exact) / std::abs(exact);

    std::cout << "\n=== Order 4, 3-dimer L-shaped cluster (columnar) ===" << std::endl;
    std::cout << "Exact (Python ED):       " << exact << std::endl;
    std::cout << "MC coefficient:          " << mc_coeff << std::endl;
    std::cout << "|Omega| integral (MC):   " << abs_integral << std::endl;
    std::cout << "Mean |Omega| (uniform):  " << meas.result->mean_abs << std::endl;
    std::cout << "Mean sign:               " << meas.result->mean_sign << std::endl;
    std::cout << "Sign error:              " << meas.result->sign_error << std::endl;
    std::cout << "Relative error:          " << rel_err << std::endl;

    EXPECT_LT(rel_err, 0.05) << "MC estimate " << mc_coeff << " deviates from exact " << exact << " by " << rel_err * 100 << "%";
  }
}

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
