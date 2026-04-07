/*
 * Tests for the dimer (N_sites=2) expansion at order 2.
 *
 * At order 2 there is a single diagram (2-cycle, D2a) with fm=6 on the
 * triangular superlattice and 1 spatial config.  The 2-dimer cluster
 * coefficient is the infinite-lattice result / fm.
 *
 * Tests:
 *   1. Quadrature: 1D trapezoidal rule (fix tau_2=0, integrate tau_1).
 *   2. MCMC: sign-based sampling with |Omega| weight via DimerConfiguration,
 *      combined with a deterministic |Omega| integral.
 *
 * ED reference from analytical/benchmark_staggered_dimer_expansion.py:
 *   Order 2 coefficient on the 2-dimer cluster = -0.066467819521
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
#include <h5/h5.hpp>
#include <cmath>
#include <memory>
#include <filesystem>

using namespace sc_expansion;

// =====================================================================
// Deterministic test: 1D trapezoidal rule
//
// Time-translation invariance: the integrand depends only on (tau_1 - tau_2).
// Fix tau_2 = 0, integrate tau_1 over [0, beta], multiply by beta.
// Divide by fm to get the 2-dimer cluster coefficient.
// =====================================================================

TEST(DimerExpansion, Order2Quadrature) {
  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  Parameters<double> params{U, beta, mu, t_intra, true};

  Graph graph({0, 1, 1, 0}, 2);
  std::vector<VertexType<2, double> *> vt;
  Diagram<2, double> diagram(graph, vt);
  HubbardSolver<2, double> solver(params);

  double fm = diagram.get_free_multiplicity();
  EXPECT_DOUBLE_EQ(fm, 6.0);

  int N           = 5000;
  double h        = beta / N;
  double integral = 0.0;

  for (int i = 0; i <= N; i++) {
    double wi                = (i == 0 || i == N) ? 0.5 : 1.0;
    std::vector<double> taus = {i * h, 0.0};
    integral += wi * diagram.evaluate(taus, solver, false);
  }
  integral *= h * beta;

  double cluster_coeff = integral / fm;

  std::cout << "Quadrature coeff:    " << cluster_coeff << std::endl;
  std::cout << "Exact (Python ED):   " << -0.066467819521 << std::endl;

  // ED result from analytical/benchmark_staggered_dimer_expansion.py
  EXPECT_NEAR(cluster_coeff, -0.066467819521, 1e-4);
}

// =====================================================================
// MCMC test: sign-based sampling with |Omega| weight
//
// The DimerConfiguration samples taus according to |Omega(tau)|.
// measure_dimer accumulates <sign(Omega)> and estimates the
// normalization integral(|Omega|) = beta^n * <|Omega|>_uniform
// via independent uniform tau samples.
//
// The coefficient is:  beta^n * <|Omega|>_uniform * <sign> / fm
// =====================================================================

TEST(DimerExpansion, Order2MCMC) {
  mpi::communicator world;

  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order = 2;
  Parameters<double> params{U, beta, mu, t_intra, true};

  Graph graph({0, 1, 1, 0}, 2);
  std::vector<VertexType<2, double> *> vt;
  Diagram<2, double> diagram(graph, vt);
  HubbardSolver<2, double> solver(params);
  double fm = diagram.get_free_multiplicity();

  // --- MCMC: estimate <sign> from |Omega| sampling, <|Omega|> from uniform sampling ---
  int n_cycles     = 100000;
  int n_warmup     = 2000;
  int length_cycle = 1;

  FreeEnergyCalculator<2, double> calculator(params, order);
  auto config = std::make_unique<DimerConfiguration<double>>(params, order, calculator);

  int random_seed = 32186222 + world.rank() * 786512;
  int verbosity   = (world.rank() == 0 ? 2 : 0);

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  if (world.rank() == 0) { std::filesystem::create_directory("./results"); }

  int measure_seed = 99871234 + world.rank() * 314159;
  measure_dimer<double> meas(config.get(), n_bins, block_size, mu, measure_seed);
  mc.add_move(move<double>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "dimer_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  // --- Read results and compute coefficient ---
  if (world.rank() == 0) {
    std::string filename = "./results/dimer_data_order_" + std::to_string(order) + "_U_" + std::to_string(config->get_U()) + "_beta_"
       + std::to_string(config->beta) + "_mu_" + std::to_string(mu) + ".h5";

    double mean_sign          = 0.0;
    double sign_error         = 0.0;
    double mean_abs_integrand = 0.0;
    double abs_integrand_error = 0.0;

    {
      h5::file file(filename, 'r');
      h5_read(file, "mean_sign", mean_sign);
      h5_read(file, "sign_error", sign_error);
      h5_read(file, "mean_abs_integrand", mean_abs_integrand);
      h5_read(file, "abs_integrand_error", abs_integrand_error);
    }

    // coefficient = beta^n * <|Omega|>_uniform * <sign> / fm
    double abs_integral = std::pow(beta, order) * mean_abs_integrand;
    double mc_coeff     = abs_integral * mean_sign / fm;
    double exact        = -0.066467819521;
    double rel_err      = std::abs(mc_coeff - exact) / std::abs(exact);

    std::cout << "Exact (Python ED):       " << exact << std::endl;
    std::cout << "MC coefficient:          " << mc_coeff << std::endl;
    std::cout << "|Omega| integral (MC):   " << abs_integral << std::endl;
    std::cout << "Mean |Omega| (uniform):  " << mean_abs_integrand << std::endl;
    std::cout << "Mean sign:               " << mean_sign << std::endl;
    std::cout << "Sign error:              " << sign_error << std::endl;
    std::cout << "Free multiplicity:       " << fm << std::endl;
    std::cout << "Relative error:          " << rel_err << std::endl;

    EXPECT_LT(rel_err, 0.10) << "MC estimate " << mc_coeff << " deviates from exact " << exact << " by " << rel_err * 100 << "%";

    std::filesystem::remove_all("./results");
  }
}

// =====================================================================
// MCMC test: order 4 on the 3-dimer triangle cluster
//
// Cluster: A=(0,0), B=(1,0), C=(0,1) on the triangular superlattice.
// Inter-dimer directions: A-B RIGHT, A-C RIGHT, B-C LEFT.
//
// DimerConfiguration uses cluster-restricted spatial embeddings, so
// the integrand already has per-dimer weights — no fm division needed.
//
// ED reference from analytical/benchmark_staggered_dimer_expansion.py:
//   Order 4 coefficient on the 3-dimer cluster = -0.037479608143
// =====================================================================

TEST(DimerExpansion, Order4MCMC_3DimerCluster) {
  mpi::communicator world;

  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order = 4;
  Parameters<double> params{U, beta, mu, t_intra, false}; // non-bipartite (triangular superlattice)

  // 3-dimer triangle cluster positions on the triangular superlattice
  std::vector<std::pair<int, int>> cluster_positions = {{0, 0}, {1, 0}, {0, 1}};
  int n_cluster_sites = 3;

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

  if (world.rank() == 0) { std::filesystem::create_directory("./results"); }

  int measure_seed = 77871234 + world.rank() * 271828;
  measure_dimer<double> meas(config.get(), n_bins, block_size, mu, measure_seed);
  mc.add_move(move<double>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "dimer_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  // --- Read results and compute coefficient ---
  if (world.rank() == 0) {
    std::string filename = "./results/dimer_data_order_" + std::to_string(order) + "_U_" + std::to_string(config->get_U()) + "_beta_"
       + std::to_string(config->beta) + "_mu_" + std::to_string(mu) + ".h5";

    double mean_sign           = 0.0;
    double sign_error          = 0.0;
    double mean_abs_integrand  = 0.0;
    double abs_integrand_error = 0.0;

    {
      h5::file file(filename, 'r');
      h5_read(file, "mean_sign", mean_sign);
      h5_read(file, "sign_error", sign_error);
      h5_read(file, "mean_abs_integrand", mean_abs_integrand);
      h5_read(file, "abs_integrand_error", abs_integrand_error);
    }

    // Cluster per-dimer weights are already in the integrand — no fm division
    double abs_integral = std::pow(beta, order) * mean_abs_integrand;
    double mc_coeff     = abs_integral * mean_sign;
    double exact        = -0.037479608143;
    double rel_err      = std::abs(mc_coeff - exact) / std::abs(exact);

    std::cout << "\n=== Order 4, 3-dimer triangle cluster ===" << std::endl;
    std::cout << "Exact (Python ED):       " << exact << std::endl;
    std::cout << "MC coefficient:          " << mc_coeff << std::endl;
    std::cout << "|Omega| integral (MC):   " << abs_integral << std::endl;
    std::cout << "Mean |Omega| (uniform):  " << mean_abs_integrand << std::endl;
    std::cout << "Mean sign:               " << mean_sign << std::endl;
    std::cout << "Sign error:              " << sign_error << std::endl;
    std::cout << "Relative error:          " << rel_err << std::endl;

    EXPECT_LT(rel_err, 0.15) << "MC estimate " << mc_coeff << " deviates from exact " << exact << " by " << rel_err * 100 << "%";

    std::filesystem::remove_all("./results");
  }
}

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
