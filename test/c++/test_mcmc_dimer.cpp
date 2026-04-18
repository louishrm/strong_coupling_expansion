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
#include <random>

using namespace sc_expansion;

// =====================================================================
// MCMC test: order 2, defensive ratio estimator
//
// Uses FreeEnergyCalculator with cluster-restricted embedding on the
// 3-dimer L-shaped cluster. The integrand already has per-dimer weights.
//
// Defensive scheme: W = |Omega + alpha|, ratio estimator with known
// normalization alpha * beta^n.
// =====================================================================

TEST(DimerExpansion, Order2MCMC) {
  mpi::communicator world;

  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order    = 2;
  double alpha = 0.01;
  Parameters<double> params{U, beta, mu, t_intra, true};

  // 3-dimer L-shaped cluster on the rectangular superlattice
  std::vector<std::pair<int, int>> cluster_positions = {{0, 0}, {0, 1}, {1, 0}};
  int n_cluster_sites                                = 3;

  // --- MCMC: cluster-restricted ---
  int n_cycles     = 100000;
  int n_warmup     = 2000;
  int length_cycle = 1;

  FreeEnergyCalculator<2, double> calculator(params, order, cluster_positions, n_cluster_sites);
  auto config = std::make_unique<DimerConfiguration<double>>(params, order, calculator, alpha);

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
    double mc_coeff = meas.result->coeff;
    double exact    = -0.130006787492;
    double rel_err  = std::abs(mc_coeff - exact) / std::abs(exact);

    std::cout << "Exact (Python ED):       " << exact << std::endl;
    std::cout << "MC coefficient:          " << mc_coeff << std::endl;
    std::cout << "MC error:                " << meas.result->error << std::endl;
    std::cout << "Relative error:          " << rel_err << std::endl;

    // Profiling and cache statistics
    std::cout << "\n--- Profile & Cache (Order 2) ---" << std::endl;
    auto const &vt2 = config->get_calculator().get_vertex_types();
    for (size_t i = 0; i < vt2.size(); i++) {
      auto [hits, misses] = vt2[i].get_cache_stats();
      long total          = hits + misses;
      double hit_rate     = total > 0 ? 100.0 * hits / total : 0.0;
      std::cout << "Global cache (cumulant order " << (i + 1) << "): hits=" << hits << " misses=" << misses << " hit_rate=" << hit_rate << "%"
                << std::endl;
    }
    auto const &diags2 = config->get_calculator().get_diagrams();
    double total_p1 = 0, total_p2 = 0;
    long total_local_hits = 0, total_local_misses = 0;
    for (size_t d = 0; d < diags2.size(); d++) {
      auto [lh, lm] = diags2[d].get_local_cache_stats();
      total_local_hits += lh;
      total_local_misses += lm;
      total_p1 += diags2[d].get_phase1_time();
      total_p2 += diags2[d].get_phase2_time();
    }
    long total_local      = total_local_hits + total_local_misses;
    double local_hit_rate = total_local > 0 ? 100.0 * total_local_hits / total_local : 0.0;
    std::cout << "Local cache: hits=" << total_local_hits << " misses=" << total_local_misses << " hit_rate=" << local_hit_rate << "%" << std::endl;
    std::cout << "Phase 1 (cumulants): " << total_p1 << " s" << std::endl;
    std::cout << "Phase 2 (contraction): " << total_p2 << " s" << std::endl;
    std::cout << "Phase 1 fraction: " << (total_p1 + total_p2 > 0 ? 100.0 * total_p1 / (total_p1 + total_p2) : 0) << "%" << std::endl;

    EXPECT_LT(rel_err, 0.04) << "MC estimate " << mc_coeff << " deviates from exact " << exact << " by " << rel_err * 100 << "%";
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
// Defensive ratio estimator: W = |Omega + alpha|.
//
// ED reference from analytical/benchmark_columnar_dimer_expansion.py:
//   Order 4 coefficient on the 3-dimer cluster = -0.037169097044
// =====================================================================

TEST(DimerExpansion, Order4MCMC_3DimerCluster) {
  mpi::communicator world;

  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order    = 4;
  double alpha = 0.001;
  Parameters<double> params{U, beta, mu, t_intra, true}; // bipartite (rectangular superlattice)

  // 3-dimer L-shaped cluster positions on the rectangular superlattice
  std::vector<std::pair<int, int>> cluster_positions = {{0, 0}, {0, 1}, {1, 0}};
  int n_cluster_sites                                = 3;

  // --- MCMC: cluster-restricted DimerConfiguration ---
  int n_cycles     = 200000;
  int n_warmup     = 5000;
  int length_cycle = 1;

  FreeEnergyCalculator<2, double> calculator(params, order, cluster_positions, n_cluster_sites);
  auto config = std::make_unique<DimerConfiguration<double>>(params, order, calculator, alpha);

  int random_seed = 42186333 + world.rank() * 786512;
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
    double mc_coeff = meas.result->coeff;
    double exact    = -0.037169097044;
    double rel_err  = std::abs(mc_coeff - exact) / std::abs(exact);

    std::cout << "\n=== Order 4, 3-dimer L-shaped cluster (columnar) ===" << std::endl;
    std::cout << "Exact (Python ED):       " << exact << std::endl;
    std::cout << "MC coefficient:          " << mc_coeff << std::endl;
    std::cout << "MC error:                " << meas.result->error << std::endl;
    std::cout << "Relative error:          " << rel_err << std::endl;

    // Profiling and cache statistics
    std::cout << "\n--- Profile & Cache (Order 4) ---" << std::endl;
    auto const &vt4 = config->get_calculator().get_vertex_types();
    for (size_t i = 0; i < vt4.size(); i++) {
      auto [hits, misses] = vt4[i].get_cache_stats();
      long total          = hits + misses;
      double hit_rate     = total > 0 ? 100.0 * hits / total : 0.0;
      std::cout << "Global cache (cumulant order " << (i + 1) << "): hits=" << hits << " misses=" << misses << " hit_rate=" << hit_rate << "%"
                << std::endl;
    }
    auto const &diags4 = config->get_calculator().get_diagrams();
    double total_p1 = 0, total_p2 = 0;
    long total_local_hits = 0, total_local_misses = 0;
    for (size_t d = 0; d < diags4.size(); d++) {
      auto [lh, lm] = diags4[d].get_local_cache_stats();
      total_local_hits += lh;
      total_local_misses += lm;
      total_p1 += diags4[d].get_phase1_time();
      total_p2 += diags4[d].get_phase2_time();
    }
    long total_local      = total_local_hits + total_local_misses;
    double local_hit_rate = total_local > 0 ? 100.0 * total_local_hits / total_local : 0.0;
    std::cout << "Local cache: hits=" << total_local_hits << " misses=" << total_local_misses << " hit_rate=" << local_hit_rate << "%" << std::endl;
    std::cout << "Phase 1 (cumulants): " << total_p1 << " s" << std::endl;
    std::cout << "Phase 2 (contraction): " << total_p2 << " s" << std::endl;
    std::cout << "Phase 1 fraction: " << (total_p1 + total_p2 > 0 ? 100.0 * total_p1 / (total_p1 + total_p2) : 0) << "%" << std::endl;

    EXPECT_LT(rel_err, 0.05) << "MC estimate " << mc_coeff << " deviates from exact " << exact << " by " << rel_err * 100 << "%";
  }
}

// =====================================================================
// Diagnostic: per-diagram contribution at order 4
//
// Evaluates each diagram independently at several random tau points
// to measure the relative magnitude of each diagram's integrand.
// Shows how much the V=2 diagram contributes vs the total.
// =====================================================================

TEST(DimerExpansion, Order4PerDiagramContribution) {
  mpi::communicator world;
  if (world.rank() != 0) return; // single-rank test

  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order = 4;
  Parameters<double> params{U, beta, mu, t_intra, true};

  std::vector<std::pair<int, int>> cluster_positions = {{0, 0}, {0, 1}, {1, 0}};
  int n_cluster_sites                                = 3;

  FreeEnergyCalculator<2, double> calculator(params, order, cluster_positions, n_cluster_sites);
  HubbardSolver<2, double> solver(params);

  auto const &diags  = calculator.get_diagrams();
  auto const &graphs = calculator.get_graphs();
  int n_diags        = (int)diags.size();

  // Evaluate each diagram at many random tau points and accumulate |D_d|
  std::mt19937 rng(12345);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);
  int n_samples = 5000;

  std::vector<double> avg_abs(n_diags, 0.0);
  double avg_abs_total = 0.0;

  for (int s = 0; s < n_samples; s++) {
    std::vector<double> taus(order);
    for (int i = 0; i < order; i++) taus[i] = tau_dist(rng);

    double total = 0.0;
    for (int d = 0; d < n_diags; d++) {
      auto &diagram = const_cast<Diagram<2, double> &>(diags[d]);
      diagram.mark_all_dirty();
      double val = diagram.evaluate(taus, solver, false);
      avg_abs[d] += std::abs(val);
      total += val;
    }
    avg_abs_total += std::abs(total);
  }

  std::cout << "\n=== Per-diagram contribution (Order 4, " << n_samples << " samples) ===" << std::endl;
  std::cout << "Diag\tV\tSymF\t<|D_d|>\t\t\tfraction of <|Omega|>" << std::endl;
  std::cout << std::string(70, '-') << std::endl;

  double sum_avg_abs = 0;
  for (int d = 0; d < n_diags; d++) sum_avg_abs += avg_abs[d];

  for (int d = 0; d < n_diags; d++) {
    avg_abs[d] /= n_samples;
    std::cout << d << "\t" << graphs[d].get_V() << "\t" << graphs[d].get_symmetry_factor() << "\t" << avg_abs[d] << "\t\t"
              << 100.0 * avg_abs[d] * n_samples / sum_avg_abs << "%" << std::endl;
  }
  avg_abs_total /= n_samples;
  std::cout << "Total <|Omega|> = " << avg_abs_total << std::endl;

  // Show V=2 contribution specifically
  double v2_total = 0;
  for (int d = 0; d < n_diags; d++) {
    if (graphs[d].get_V() <= 2) v2_total += avg_abs[d];
  }
  std::cout << "\nV<=2 fraction of sum(|D_d|): " << 100.0 * v2_total * n_samples / sum_avg_abs << "%" << std::endl;
  std::cout << "======================================================\n" << std::endl;
}

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
