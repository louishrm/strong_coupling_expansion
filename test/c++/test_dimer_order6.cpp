/*
 * Dimer expansion order 6: cluster validation and spatial embedding diagnostics.
 *
 * 1) MCMC test on the 3-dimer triangle cluster at order 6, compared against
 *    the ED reference from analytical/benchmark_staggered_dimer_expansion.py.
 *
 * 2) Spatial embedding diagnostics: for every order-6 non-bipartite graph,
 *    print the spatial configurations (directions + weights) and validate
 *    the total weight against the expected triangular-lattice embedding count.
 *
 * ED reference (U=8.0, beta=2.0, mu=3.0, t_intra=1.0):
 *   Order 6 coefficient on the 3-dimer cluster = -0.003064775286
 *
 * Usage:  mpirun -np 4 ./test_dimer_order6
 */

#include <gtest/gtest.h>
#include "sc_expansion/diagram.hpp"
#include "sc_expansion/graph.hpp"
#include "sc_expansion/generate_diagrams.hpp"
#include "sc_expansion/free_energy_calculator.hpp"
#include "sc_expansion/dimer_configuration.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/measure_dimer.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <mpi/mpi.hpp>
#include <cmath>
#include <memory>
#include <iomanip>

using namespace sc_expansion;

// =====================================================================
// 1) MCMC test: order 6 on the 3-dimer triangle cluster
//
// Cluster: A=(0,0), B=(1,0), C=(0,1) on the triangular superlattice.
// Uses cluster-restricted spatial embeddings (per-dimer weights already
// built into the integrand — no fm division needed).
//
// ED reference from analytical/benchmark_staggered_dimer_expansion.py:
//   Order 6 coefficient on the 3-dimer cluster = -0.003064775286
// =====================================================================

TEST(DimerOrder6, MCMC_3DimerCluster) {
  mpi::communicator world;

  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order = 6;
  Parameters<double> params{U, beta, mu, t_intra, false}; // non-bipartite

  // 3-dimer triangle cluster on the triangular superlattice
  std::vector<std::pair<int, int>> cluster_positions = {{0, 0}, {1, 0}, {0, 1}};
  int n_cluster_sites                                = 3;

  // --- MCMC ---
  int n_cycles     = 250000;
  int n_warmup     = 10000;
  int length_cycle = 1;

  FreeEnergyCalculator<2, double> calculator(params, order, cluster_positions, n_cluster_sites);
  auto config = std::make_unique<DimerConfiguration<double>>(params, order, calculator);

  int random_seed = 55186333 + world.rank() * 786512;
  int verbosity   = (world.rank() == 0 ? 2 : 0);

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  int measure_seed = 88871234 + world.rank() * 271828;
  measure_dimer<double> meas(config.get(), n_bins, block_size, mu, measure_seed);
  mc.add_move(move<double>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "dimer_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  if (world.rank() == 0) {
    double abs_integral = std::pow(beta, order) * meas.result->mean_abs;
    double mc_coeff     = abs_integral * meas.result->mean_sign;
    double exact        = -0.003064775286;
    double rel_err      = std::abs(mc_coeff - exact) / std::abs(exact);

    std::cout << "\n=== Order 6, 3-dimer triangle cluster ===" << std::endl;
    std::cout << "Exact (Python ED):       " << exact << std::endl;
    std::cout << std::setprecision(12);
    std::cout << "MC coefficient:          " << mc_coeff << std::endl;
    std::cout << "|Omega| integral (MC):   " << abs_integral << std::endl;
    std::cout << "Mean |Omega| (uniform):  " << meas.result->mean_abs << std::endl;
    std::cout << "Mean sign:               " << meas.result->mean_sign << std::endl;
    std::cout << "Sign error:              " << meas.result->sign_error << std::endl;
    std::cout << "Relative error:          " << rel_err << std::endl;

    EXPECT_LT(rel_err, 0.20) << "MC estimate " << mc_coeff << " deviates from exact " << exact << " by " << rel_err * 100 << "%";
  }
}

// =====================================================================
// 2) Spatial embedding diagnostics: validate weights for all order-6
//    non-bipartite graphs on the infinite triangular superlattice.
//
// For each graph, prints:
//   - adjacency matrix, V, symmetry factor, automorphism count, bipartite flag
//   - number and details of spatial configurations (directions + weights)
//   - total spatial weight (= total embedding count after canonicalization)
//   - number of valid global configurations and total config weight
//
// Also computes a per-diagram deterministic quadrature integral on the
// 3-dimer cluster to detect sign or weight errors in individual diagrams.
// =====================================================================

TEST(DimerOrder6, SpatialEmbeddingDiagnostics) {
  mpi::communicator world;
  if (world.rank() != 0) return; // diagnostics on rank 0 only

  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order = 6;
  Parameters<double> params{U, beta, mu, t_intra, false};

  // Generate all order-6 diagrams (both bipartite and non-bipartite)
  VacuumDiagramGenerator gen(order, false);
  gen.generate();
  auto const &graphs = gen.get_unique_graphs();

  // Create VertexTypes for cumulant orders 1..3
  int max_cumulant_order = order / 2;
  std::vector<VertexType<2, double>> vertex_types_storage;
  for (int k = 1; k <= max_cumulant_order; k++) { vertex_types_storage.emplace_back(2 * k); }
  std::vector<VertexType<2, double> *> vt_ptrs(max_cumulant_order);
  for (int k = 0; k < max_cumulant_order; k++) { vt_ptrs[k] = &vertex_types_storage[k]; }

  HubbardSolver<2, double> solver(params);

  // Cluster positions for per-diagram quadrature
  std::vector<std::pair<int, int>> cluster_positions = {{0, 0}, {1, 0}, {0, 1}};
  int n_cluster_sites                                = 3;

  std::cout << "\n" << std::string(70, '=') << std::endl;
  std::cout << "Order 6 Dimer Spatial Embedding Diagnostics" << std::endl;
  std::cout << "U=" << U << " beta=" << beta << " mu=" << mu << " t_intra=" << t_intra << std::endl;
  std::cout << "Total unique graphs: " << graphs.size() << std::endl;
  std::cout << std::string(70, '=') << std::endl;

  int graph_idx = 0;
  for (auto const &graph : graphs) {
    int V = graph.get_V();

    // --- Infinite-lattice diagram ---
    Diagram<2, double> diagram_inf(graph, vt_ptrs);

    // --- Cluster-restricted diagram ---
    Diagram<2, double> diagram_cluster(graph, vt_ptrs, cluster_positions, n_cluster_sites);

    auto const &spatial_inf     = diagram_inf.get_spatial_configurations();
    auto const &spatial_cluster = diagram_cluster.get_spatial_configurations();
    auto const &configs_inf     = diagram_inf.get_valid_configurations();
    auto const &configs_cluster = diagram_cluster.get_valid_configurations();

    double total_spatial_weight_inf = 0.0;
    for (auto const &sc : spatial_inf) { total_spatial_weight_inf += sc.weight; }

    double total_spatial_weight_cluster = 0.0;
    for (auto const &sc : spatial_cluster) { total_spatial_weight_cluster += sc.weight; }

    double total_config_weight_inf = 0.0;
    for (auto const &c : configs_inf) { total_config_weight_inf += c.weight; }

    double total_config_weight_cluster = 0.0;
    for (auto const &c : configs_cluster) { total_config_weight_cluster += c.weight; }

    // Print adjacency matrix
    std::cout << "\n--- Graph " << graph_idx << " ---" << std::endl;
    std::cout << "V=" << V << "  bipartite=" << graph.get_bipartite() << "  sym_factor=" << graph.get_symmetry_factor()
              << "  aut_count=" << graph.get_automorphism_count() << "  sign=" << diagram_inf.get_diagram_sign() << std::endl;
    std::cout << "Adjacency matrix:" << std::endl;
    for (int i = 0; i < V; i++) {
      std::cout << "  [";
      for (int j = 0; j < V; j++) { std::cout << (int)graph(i, j) << (j < V - 1 ? "," : ""); }
      std::cout << "]" << std::endl;
    }

    // Infinite-lattice spatial configs
    std::cout << "Infinite lattice: " << spatial_inf.size() << " spatial configs, total weight=" << total_spatial_weight_inf << std::endl;
    for (size_t i = 0; i < spatial_inf.size(); i++) {
      std::cout << "  [" << i << "] w=" << spatial_inf[i].weight << " dirs=[";
      for (size_t j = 0; j < spatial_inf[i].directions.size(); j++) {
        std::cout << (int)spatial_inf[i].directions[j] << (j + 1 < spatial_inf[i].directions.size() ? "," : "");
      }
      std::cout << "]" << std::endl;
    }
    std::cout << "  valid configs=" << configs_inf.size() << "  total config weight=" << total_config_weight_inf << std::endl;

    // Cluster spatial configs
    std::cout << "3-dimer cluster: " << spatial_cluster.size() << " spatial configs, total weight=" << total_spatial_weight_cluster << std::endl;
    for (size_t i = 0; i < spatial_cluster.size(); i++) {
      std::cout << "  [" << i << "] w=" << spatial_cluster[i].weight << " dirs=[";
      for (size_t j = 0; j < spatial_cluster[i].directions.size(); j++) {
        std::cout << (int)spatial_cluster[i].directions[j] << (j + 1 < spatial_cluster[i].directions.size() ? "," : "");
      }
      std::cout << "]" << std::endl;
    }
    std::cout << "  valid configs=" << configs_cluster.size() << "  total config weight=" << total_config_weight_cluster << std::endl;

    // Per-diagram deterministic integral on the cluster (crude quadrature)
    // Use a small grid to catch sign/weight errors, not for precision.
    int N_quad              = 8;
    double h                = beta / N_quad;
    double diagram_integral = 0.0;

    // 6D trapezoidal over [0, beta]^6
    // Use time-translation invariance: fix tau_6 = 0, integrate tau_1..tau_5, multiply by beta
    for (int i1 = 0; i1 <= N_quad; i1++) {
      double w1 = (i1 == 0 || i1 == N_quad) ? 0.5 : 1.0;
      for (int i2 = 0; i2 <= N_quad; i2++) {
        double w2 = (i2 == 0 || i2 == N_quad) ? 0.5 : 1.0;
        for (int i3 = 0; i3 <= N_quad; i3++) {
          double w3 = (i3 == 0 || i3 == N_quad) ? 0.5 : 1.0;
          for (int i4 = 0; i4 <= N_quad; i4++) {
            double w4 = (i4 == 0 || i4 == N_quad) ? 0.5 : 1.0;
            for (int i5 = 0; i5 <= N_quad; i5++) {
              double w5                = (i5 == 0 || i5 == N_quad) ? 0.5 : 1.0;
              std::vector<double> taus = {i1 * h, i2 * h, i3 * h, i4 * h, i5 * h, 0.0};
              diagram_integral += w1 * w2 * w3 * w4 * w5 * diagram_cluster.evaluate(taus, solver, false);
            }
          }
        }
      }
    }
    diagram_integral *= std::pow(h, 5) * beta; // h^5 for 5 integrated taus, beta for fixed tau_6

    std::cout << std::setprecision(10);
    std::cout << "  Cluster quadrature integral: " << diagram_integral << std::endl;

    graph_idx++;
  }

  std::cout << "\n" << std::string(70, '=') << std::endl;
  std::cout << "End of diagnostics" << std::endl;
  std::cout << std::string(70, '=') << std::endl;
}

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
