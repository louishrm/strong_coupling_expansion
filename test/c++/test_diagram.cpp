#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <gtest/gtest.h>
#include "sc_expansion/diagram.hpp"
#include "sc_expansion/graph.hpp"
#include "sc_expansion/free_energy_calculator.hpp"

using namespace sc_expansion;

// =====================================================================
// Helper
// =====================================================================

static double single_site_free_multiplicity(std::vector<uint8_t> const &adjmat, int V, bool bipartite_only = true) {
  Graph graph(adjmat, V, bipartite_only);
  std::vector<VertexType<double> *> vt;
  Diagram<double> diagram(graph, vt);
  return diagram.get_free_multiplicity();
}

// =====================================================================
// Free multiplicity tests (single-site expansion)
// =====================================================================

TEST(DiagramFreeMultiplicity, D2a) { EXPECT_DOUBLE_EQ(single_site_free_multiplicity({0, 1, 1, 0}, 2), 4); }

TEST(DiagramFreeMultiplicity, D3a) { EXPECT_DOUBLE_EQ(single_site_free_multiplicity({0, 1, 0, 0, 0, 1, 1, 0, 0}, 3, false), 12); }

TEST(DiagramFreeMultiplicity, D4a) { EXPECT_DOUBLE_EQ(single_site_free_multiplicity({0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0}, 4), 36); }

TEST(DiagramFreeMultiplicity, D4b) { EXPECT_DOUBLE_EQ(single_site_free_multiplicity({0, 1, 1, 1, 0, 0, 1, 0, 0}, 3), 16); }

TEST(DiagramFreeMultiplicity, D4c) { EXPECT_DOUBLE_EQ(single_site_free_multiplicity({0, 2, 2, 0}, 2), 4); }

TEST(DiagramFreeMultiplicity, D6_1) { EXPECT_DOUBLE_EQ(single_site_free_multiplicity({0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0}, 4), 64); }

TEST(DiagramFreeMultiplicity, D6_2) { EXPECT_DOUBLE_EQ(single_site_free_multiplicity({0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0}, 4), 64); }

TEST(DiagramFreeMultiplicity, D6a) {
  EXPECT_DOUBLE_EQ(
     single_site_free_multiplicity({0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0}, 6),
     400);
}

TEST(DiagramFreeMultiplicity, D6b) { EXPECT_DOUBLE_EQ(single_site_free_multiplicity({0, 3, 3, 0}, 2), 4); }

TEST(DiagramFreeMultiplicity, D6c) { EXPECT_DOUBLE_EQ(single_site_free_multiplicity({0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0}, 4), 64); }

TEST(DiagramFreeMultiplicity, D6d) {
  EXPECT_DOUBLE_EQ(single_site_free_multiplicity({0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0}, 5), 144);
}

TEST(DiagramFreeMultiplicity, D6e) { EXPECT_DOUBLE_EQ(single_site_free_multiplicity({0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0}, 4), 64); }

TEST(DiagramFreeMultiplicity, D6f) { EXPECT_DOUBLE_EQ(single_site_free_multiplicity({0, 2, 1, 2, 0, 0, 1, 0, 0}, 3), 16); }

TEST(DiagramFreeMultiplicity, D6g) { EXPECT_DOUBLE_EQ(single_site_free_multiplicity({0, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0}, 4), 36); }

TEST(DiagramFreeMultiplicity, D8a) {
  std::vector<uint8_t> D8a(64, 0);
  for (int i = 0; i < 8; ++i) {
    D8a[i * 8 + (i + 1) % 8]   = 1;
    D8a[((i + 1) % 8) * 8 + i] = 1;
  }
  EXPECT_DOUBLE_EQ(single_site_free_multiplicity(D8a, 8), 4900);
}

TEST(DiagramFreeMultiplicity, D8b) { EXPECT_DOUBLE_EQ(single_site_free_multiplicity({0, 4, 4, 0}, 2), 4); }

// =====================================================================
// Global configuration tests
// =====================================================================

TEST(DiagramGlobalConfigs, D2aAtom) {
  Graph graph({0, 1, 1, 0}, 2);
  std::vector<VertexType<double> *> vt;
  Diagram<double> diagram(graph, vt);
  auto const &configs = diagram.get_valid_configurations();

  EXPECT_EQ(configs.size(), 1u);
  EXPECT_DOUBLE_EQ(configs[0].weight, 4.0);
}

TEST(DiagramGlobalConfigs, D4bAtom) {
  Graph graph({0, 1, 1, 1, 0, 0, 1, 0, 0}, 3);
  std::vector<VertexType<double> *> vt;
  Diagram<double> diagram(graph, vt);
  auto const &configs = diagram.get_valid_configurations();

  EXPECT_EQ(configs.size(), 2u);
  double total = 0;
  for (auto const &c : configs) total += c.weight;
  EXPECT_DOUBLE_EQ(total, 32.0);
}

// =====================================================================
// Numerical evaluation tests
// =====================================================================

TEST(DiagramEvaluation, D2aAtom) {
  double U = 8.0, beta = 1.0, mu = 2.0;
  Parameters<double> params{U, beta, mu, true};
  std::vector<double> taus = {0.5, 0.0};

  Graph graph({0, 1, 1, 0}, 2);
  VertexType<double> vt1(2);
  std::vector<VertexType<double> *> vt = {&vt1};
  Diagram<double> diagram(graph, vt);
  HubbardSolver<double> solver(params);

  double val = diagram.evaluate(taus, solver, false);
  EXPECT_TRUE(std::isfinite(val));
  EXPECT_NE(val, 0.0);

  double val2 = diagram.evaluate(taus, solver, false);
  EXPECT_DOUBLE_EQ(val, val2);
}

// =====================================================================
// Order-1 self-insertion graph: a single vertex carrying one density
// insertion c^dag_sigma c_sigma (the {1} adjacency). Pen-and-paper
// prediction:
//
//   compute_sum_diagrams([tau], false) with delta = 1
//       = (-delta)^1 * evaluate(tau)   (per the (-t)^{n-k} delta^k factor)
//       = -2 * delta * n_sigma / beta
//
// which is tau-independent (single vertex -> only one time variable, and
// the atomic density is a static expectation value). The overall minus is
// the implicit (-t) that the evaluator carries for every line; with one
// self-insertion and zero hops at order 1 we pick up one factor of -1.
// =====================================================================
TEST(DiagramEvaluation, Order1SelfLoopGivesSingleSpinDensity) {
  double U     = 8.0;
  double beta  = 1.0;
  double mu    = 2.0;
  double delta = 1.0;
  Parameters<double> params{U, beta, mu, /*bipartite=*/true, delta};

  FreeEnergyCalculator<double> calculator(params, /*order=*/1, /*override_fm=*/-1, /*allow_self_loops=*/true);
  ASSERT_EQ(calculator.get_n_diagrams(), 1);

  // Pick a nonzero random time in (0, beta) to confirm tau-independence isn't masked by tau=0.
  std::vector<double> taus = {0.37 * beta};
  double val               = calculator.compute_sum_diagrams(taus, false);

  // Single-spin density on the Hubbard atom:
  //   Z        = 1 + 2 e^{beta mu} + e^{beta (2 mu - U)}
  //   <n_sigma> = (e^{beta mu} + e^{beta (2 mu - U)}) / Z
  double exp_mu  = std::exp(beta * mu);
  double exp_dbl = std::exp(beta * (2.0 * mu - U));
  double Z       = 1.0 + 2.0 * exp_mu + exp_dbl;
  double n_sigma = (exp_mu + exp_dbl) / Z;

  double expected = 2.0 * delta * n_sigma / beta;
  EXPECT_NEAR(val, expected, 1e-12);
}

TEST(DiagramDiagnostics, Order8AtomGlobalConfigCounts) {
  double U = 8.0, beta = 2.0, mu = 3.0;
  Parameters<double> params{U, beta, mu, true};

  FreeEnergyCalculator<double> calculator(params, 8);

  auto const &graphs   = calculator.get_graphs();
  auto const &diagrams = calculator.get_diagrams();

  std::cout << "\n===== Order-8 Atomic: Global Config Counts =====" << std::endl;
  std::cout << std::left << std::setw(8) << "Graph" << std::setw(6) << "V" << std::setw(12) << "SymFactor" << std::setw(12) << "FreeMult"
            << std::setw(12) << "GlobalCfgs" << std::setw(14) << "TotalWeight" << std::endl;
  std::cout << std::string(64, '-') << std::endl;

  double grand_total_weight = 0;
  for (size_t i = 0; i < diagrams.size(); i++) {
    auto const &diagram = diagrams[i];
    auto const &graph   = graphs[i];

    auto const &configs = diagram.get_valid_configurations();
    double free_mult    = diagram.get_free_multiplicity();

    double total_weight = 0;
    for (auto const &c : configs) total_weight += c.weight;
    grand_total_weight += total_weight;

    std::cout << std::left << std::setw(8) << i << std::setw(6) << graph.get_V() << std::setw(12) << graph.get_symmetry_factor() << std::setw(12)
              << free_mult << std::setw(12) << configs.size() << std::setw(14) << total_weight << std::endl;
  }

  std::cout << std::string(64, '-') << std::endl;
  std::cout << "Total graphs: " << diagrams.size() << "  |  Grand total weight: " << grand_total_weight << std::endl;
  std::cout << "=========================================================\n" << std::endl;
}
