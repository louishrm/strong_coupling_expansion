/*
 * Tests for Diagrammatic Monte Carlo (DiagMC) on the dimer expansion.
 *
 * Cluster geometry: 3 dimers on the rectangular (columnar) superlattice.
 *   A=(0,0), B=(0,1), C=(1,0)
 *   A--B vertical (2 bonds), A--C horizontal (1 bond)
 *
 * Tests:
 *   1. Order 2 DiagMC: validate against exact ED benchmark.
 *   2. Order 4 DiagMC: validate against exact ED benchmark.
 *   3. Diagram visit statistics: verify all diagrams are visited.
 *
 * ED references from analytical/benchmark_columnar_dimer_expansion.py:
 *   Order 2 coefficient on the 3-dimer cluster = -0.130006787492
 *   Order 4 coefficient on the 3-dimer cluster = -0.037169097044
 *
 * Usage:  mpirun -np 4 ./test_diagmc
 */

#include <gtest/gtest.h>
#include "sc_expansion/diagmc_configuration.hpp"
#include "sc_expansion/diagmc_move.hpp"
#include "sc_expansion/measure_diagmc.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <mpi/mpi.hpp>
#include <cmath>
#include <memory>

using namespace sc_expansion;

// =====================================================================
// DiagMC test: order 2, 3-dimer L-shaped cluster
//
// At order 2 there is only 1 physical diagram, so DiagMC alternates
// between the reference (d=0) and the single physical diagram (d=1).
// Validates that the DiagMC estimator produces the correct coefficient.
// =====================================================================

TEST(DiagMC, Order2) {
  mpi::communicator world;

  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order    = 2;
  double alpha = 0.01;
  Parameters<double> params{U, beta, mu, t_intra, true};

  // 3-dimer L-shaped cluster on the rectangular superlattice
  std::vector<std::pair<int, int>> cluster_positions = {{0, 0}, {0, 1}, {1, 0}};
  int n_cluster_sites                                = 3;

  int n_cycles     = 200000;
  int n_warmup     = 5000;
  int length_cycle = 1;

  FreeEnergyCalculator<2, double> calculator(params, order, cluster_positions, n_cluster_sites);
  auto config = std::make_unique<DiagMCConfiguration<double>>(params, order, calculator, alpha);

  if (world.rank() == 0) {
    std::cout << "\n=== DiagMC Order 2 ===" << std::endl;
    std::cout << "N_diagrams = " << config->get_n_diagrams() << std::endl;
  }

  int random_seed = 32186222 + world.rank() * 786512;
  int verbosity   = (world.rank() == 0 ? 2 : 0);

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  MeasureDiagMC<double> meas(config.get(), n_bins, block_size, verbosity);
  mc.add_move(DiagMCMove<double>(config.get(), mc.get_rng()), "diagmc_move");
  mc.add_measure(meas, "diagmc_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  if (world.rank() == 0) {
    double mc_coeff = meas.result->coeff;
    double exact    = -0.130006787492;
    double rel_err  = std::abs(mc_coeff - exact) / std::abs(exact);

    std::cout << "Exact (Python ED):       " << exact << std::endl;
    std::cout << "DiagMC coefficient:      " << mc_coeff << std::endl;
    std::cout << "DiagMC error:            " << meas.result->error << std::endl;
    std::cout << "Relative error:          " << rel_err << std::endl;

    // DiagMC sign estimator has higher variance than defensive scheme on small clusters
    EXPECT_LT(rel_err, 0.10) << "DiagMC estimate " << mc_coeff << " deviates from exact " << exact << " by " << rel_err * 100 << "%";
  }
}

// =====================================================================
// DiagMC test: order 4, 3-dimer L-shaped cluster
//
// Order 4 has 3 physical diagrams (V=2, V=3, V=4). The walker should
// spend time proportional to each diagram's integrand magnitude.
// =====================================================================

TEST(DiagMC, Order4) {
  mpi::communicator world;

  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order    = 4;
  double alpha = 0.001;
  Parameters<double> params{U, beta, mu, t_intra, true};

  std::vector<std::pair<int, int>> cluster_positions = {{0, 0}, {0, 1}, {1, 0}};
  int n_cluster_sites                                = 3;

  int n_cycles     = 500000;
  int n_warmup     = 5000;
  int length_cycle = 1;

  FreeEnergyCalculator<2, double> calculator(params, order, cluster_positions, n_cluster_sites);
  auto config = std::make_unique<DiagMCConfiguration<double>>(params, order, calculator, alpha);

  if (world.rank() == 0) {
    std::cout << "\n=== DiagMC Order 4 ===" << std::endl;
    std::cout << "N_diagrams = " << config->get_n_diagrams() << std::endl;
  }

  int random_seed = 42186333 + world.rank() * 786512;
  int verbosity   = (world.rank() == 0 ? 2 : 0);

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  MeasureDiagMC<double> meas(config.get(), n_bins, block_size, verbosity);
  mc.add_move(DiagMCMove<double>(config.get(), mc.get_rng()), "diagmc_move");
  mc.add_measure(meas, "diagmc_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  if (world.rank() == 0) {
    double mc_coeff = meas.result->coeff;
    double exact    = -0.037169097044;
    double rel_err  = std::abs(mc_coeff - exact) / std::abs(exact);

    std::cout << "\n=== Order 4, 3-dimer L-shaped cluster (DiagMC) ===" << std::endl;
    std::cout << "Exact (Python ED):       " << exact << std::endl;
    std::cout << "DiagMC coefficient:      " << mc_coeff << std::endl;
    std::cout << "DiagMC error:            " << meas.result->error << std::endl;
    std::cout << "Relative error:          " << rel_err << std::endl;

    // DiagMC sign estimator has higher variance than defensive scheme on small clusters
    EXPECT_LT(rel_err, 0.10) << "DiagMC estimate " << mc_coeff << " deviates from exact " << exact << " by " << rel_err * 100 << "%";
  }
}

// =====================================================================
// Diagnostic: diagram visit statistics at order 4
//
// Verifies that all diagrams (including reference) are visited,
// and prints the visit distribution.
// =====================================================================

TEST(DiagMC, DiagramVisitStats) {
  mpi::communicator world;
  if (world.rank() != 0) return; // single-rank test

  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order    = 4;
  double alpha = 0.001;
  Parameters<double> params{U, beta, mu, t_intra, true};

  std::vector<std::pair<int, int>> cluster_positions = {{0, 0}, {0, 1}, {1, 0}};
  int n_cluster_sites                                = 3;

  int n_cycles     = 100000;
  int n_warmup     = 2000;
  int length_cycle = 1;

  FreeEnergyCalculator<2, double> calculator(params, order, cluster_positions, n_cluster_sites);
  auto config = std::make_unique<DiagMCConfiguration<double>>(params, order, calculator, alpha);

  triqs::mc_tools::mc_generic<double> mc("", 99999, 0);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  MeasureDiagMC<double> meas(config.get(), n_bins, block_size, 0);
  mc.add_move(DiagMCMove<double>(config.get(), mc.get_rng()), "diagmc_move");
  mc.add_measure(meas, "diagmc_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  std::cout << "\n=== Diagram Visit Statistics (Order 4, DiagMC) ===" << std::endl;
  auto const &vc = *meas.visit_counts;
  long total_visits = 0;
  for (auto v : vc) total_visits += v;

  // Verify all physical diagrams were visited at least once
  for (int d = 1; d <= config->get_n_diagrams(); d++) {
    EXPECT_GT(vc[d], 0) << "Diagram d=" << d << " was never visited";
    double frac = 100.0 * (double)vc[d] / (double)total_visits;
    std::cout << "d=" << d << ": " << vc[d] << " visits (" << frac << "%)" << std::endl;
  }
  std::cout << "Total visits: " << total_visits << std::endl;
}

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
