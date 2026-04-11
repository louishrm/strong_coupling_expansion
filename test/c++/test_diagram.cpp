#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
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
  std::vector<VertexType<1, double> *> vt;
  Diagram<1, double> diagram(graph, vt);
  return diagram.get_free_multiplicity();
}

// =====================================================================
// Free multiplicity tests (N_sites=1, single-site expansion)
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
// Dimer spatial configuration tests (N_sites=2)
// =====================================================================

TEST(DiagramSpatialConfigs, D4b) {
  Graph graph({0, 1, 1, 1, 0, 0, 1, 0, 0}, 3);
  std::vector<VertexType<2, double> *> vt;
  Diagram<2, double> diagram(graph, vt);
  auto const &spatial = diagram.get_spatial_configurations();

  EXPECT_EQ(spatial.size(), 2u);
  double total = 0;
  for (auto const &sc : spatial) total += sc.weight;
  EXPECT_DOUBLE_EQ(total, 36.0);
}

TEST(DiagramSpatialConfigs, D2a) {
  Graph graph({0, 1, 1, 0}, 2);
  std::vector<VertexType<2, double> *> vt;
  Diagram<2, double> diagram(graph, vt);
  auto const &spatial = diagram.get_spatial_configurations();

  EXPECT_EQ(spatial.size(), 1u);
  EXPECT_DOUBLE_EQ(spatial[0].weight, 6.0);
}

TEST(DiagramSpatialConfigs, D4c) {
  Graph graph({0, 2, 2, 0}, 2);
  std::vector<VertexType<2, double> *> vt;
  Diagram<2, double> diagram(graph, vt);
  auto const &spatial = diagram.get_spatial_configurations();

  EXPECT_EQ(spatial.size(), 1u);
  EXPECT_DOUBLE_EQ(spatial[0].weight, 6.0);
}

TEST(DiagramSpatialConfigs, D6c) {
  Graph graph({0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0}, 4);
  std::vector<VertexType<2, double> *> vt;
  Diagram<2, double> diagram(graph, vt);
  auto const &spatial = diagram.get_spatial_configurations();

  EXPECT_EQ(spatial.size(), 2u);
  EXPECT_DOUBLE_EQ(spatial[0].weight, 54.0);
  EXPECT_DOUBLE_EQ(spatial[1].weight, 162.0);
}

// =====================================================================
// Global configuration tests
// =====================================================================

TEST(DiagramGlobalConfigs, D2aAtom) {
  Graph graph({0, 1, 1, 0}, 2);
  std::vector<VertexType<1, double> *> vt;
  Diagram<1, double> diagram(graph, vt);
  auto const &configs = diagram.get_valid_configurations();

  EXPECT_EQ(configs.size(), 1u);
  EXPECT_DOUBLE_EQ(configs[0].weight, 4.0);
}

TEST(DiagramGlobalConfigs, D2aDimer) {
  Graph graph({0, 1, 1, 0}, 2);
  std::vector<VertexType<2, double> *> vt;
  Diagram<2, double> diagram(graph, vt);
  auto const &configs = diagram.get_valid_configurations();

  EXPECT_EQ(configs.size(), 1u);
  // weight = spatial_weight * orbit_size / sym_factor = 6 * 2 / 2 = 6
  // orbit_size=2 from SpinFlip only (Reflect is in spatial canonicalization)
  EXPECT_DOUBLE_EQ(configs[0].weight, 6.0);
}

TEST(DiagramGlobalConfigs, D4bAtom) {
  Graph graph({0, 1, 1, 1, 0, 0, 1, 0, 0}, 3);
  std::vector<VertexType<1, double> *> vt;
  Diagram<1, double> diagram(graph, vt);
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
  Parameters<double> params{U, beta, mu, 0.0, true};
  std::vector<double> taus = {0.5, 0.0};

  Graph graph({0, 1, 1, 0}, 2);
  VertexType<1, double> vt1(2);
  std::vector<VertexType<1, double> *> vt = {&vt1};
  Diagram<1, double> diagram(graph, vt);
  HubbardSolver<1, double> solver(params);

  double val = diagram.evaluate(taus, solver, false);
  EXPECT_TRUE(std::isfinite(val));
  EXPECT_NE(val, 0.0);

  double val2 = diagram.evaluate(taus, solver, false);
  EXPECT_DOUBLE_EQ(val, val2);
}

TEST(DiagramEvaluation, D2aDimer) {
  double U = 8.0, beta = 1.0, mu = 2.0;
  Parameters<double> params{U, beta, mu, 1.0, true};
  std::vector<double> taus = {0.5, 0.0};

  Graph graph({0, 1, 1, 0}, 2);
  VertexType<2, double> vt1(2);
  std::vector<VertexType<2, double> *> vt = {&vt1};
  Diagram<2, double> diagram(graph, vt);
  HubbardSolver<2, double> solver(params);

  double val = diagram.evaluate(taus, solver, false);
  EXPECT_TRUE(std::isfinite(val));
  EXPECT_NE(val, 0.0);
}

// =====================================================================
// Diagnostics: [2,2,2] V=3 non-bipartite graphs (dimer, N_sites=2)
// =====================================================================

TEST(DiagramDiagnostics, V3_222_K3) {
  // K3: complete directed graph on 3 vertices (one edge each direction between all pairs)
  // adj = {0,1,1, 1,0,1, 1,1,0}
  std::vector<uint8_t> adjmat = {0, 1, 1, 1, 0, 1, 1, 1, 0};
  Graph graph(adjmat, 3, false); // bipartite_only = false

  double U = 8.0, beta = 2.0, mu = 3.0, t_hop = 1.0;
  Parameters<double> params{U, beta, mu, t_hop, false};

  VertexType<2, double> vt2(4); // n_legs = 4 for cumulant_order = 2
  std::vector<VertexType<2, double> *> vt = {nullptr, &vt2};
  Diagram<2, double> diagram(graph, vt);
  HubbardSolver<2, double> solver(params);

  auto const &spatial = diagram.get_spatial_configurations();
  auto const &configs = diagram.get_valid_configurations();

  std::cout << "\n========== [2,2,2] V=3: K3 (complete directed graph) ==========" << std::endl;
  std::cout << "Symmetry factor:        " << graph.get_symmetry_factor() << std::endl;
  std::cout << "Automorphism count:     " << graph.get_automorphism_count() << std::endl;
  std::cout << "Diagram sign:           " << diagram.get_diagram_sign() << std::endl;
  std::cout << "Num spatial configs:    " << spatial.size() << std::endl;

  double total_spatial_weight = 0;
  for (size_t i = 0; i < spatial.size(); i++) {
    std::cout << "  spatial[" << i << "]: weight=" << spatial[i].weight << " dirs=[";
    for (size_t j = 0; j < spatial[i].directions.size(); j++) {
      std::cout << (int)spatial[i].directions[j] << (j + 1 < spatial[i].directions.size() ? "," : "");
    }
    std::cout << "]" << std::endl;
    total_spatial_weight += spatial[i].weight;
  }
  std::cout << "Total spatial weight:   " << total_spatial_weight << std::endl;

  std::cout << "Num valid configs:      " << configs.size() << std::endl;
  double total_config_weight = 0;
  for (auto const &c : configs) total_config_weight += c.weight;
  std::cout << "Total config weight:    " << total_config_weight << std::endl;

  std::vector<double> taus = {0.1, 0.4, 0.7, 1.2, 1.5, 1.8};
  double val = diagram.evaluate(taus, solver, false);
  std::cout << "evaluate(taus):         " << val << std::endl;
  std::cout << "=========================================================\n" << std::endl;

  EXPECT_TRUE(std::isfinite(val));
}

TEST(DiagramDiagnostics, V3_222_Doubled3Cycle) {
  // Doubled 3-cycle: 0→1 (×2), 1→2 (×2), 2→0 (×2)
  // adj = {0,2,0, 0,0,2, 2,0,0}
  std::vector<uint8_t> adjmat = {0, 2, 0, 0, 0, 2, 2, 0, 0};
  Graph graph(adjmat, 3, false); // bipartite_only = false

  double U = 8.0, beta = 2.0, mu = 3.0, t_hop = 1.0;
  Parameters<double> params{U, beta, mu, t_hop, false};

  VertexType<2, double> vt2(4); // n_legs = 4 for cumulant_order = 2
  std::vector<VertexType<2, double> *> vt = {nullptr, &vt2};
  Diagram<2, double> diagram(graph, vt);
  HubbardSolver<2, double> solver(params);

  auto const &spatial = diagram.get_spatial_configurations();
  auto const &configs = diagram.get_valid_configurations();

  std::cout << "\n========== [2,2,2] V=3: Doubled 3-cycle ==========" << std::endl;
  std::cout << "Symmetry factor:        " << graph.get_symmetry_factor() << std::endl;
  std::cout << "Automorphism count:     " << graph.get_automorphism_count() << std::endl;
  std::cout << "Diagram sign:           " << diagram.get_diagram_sign() << std::endl;
  std::cout << "Num spatial configs:    " << spatial.size() << std::endl;

  double total_spatial_weight = 0;
  for (size_t i = 0; i < spatial.size(); i++) {
    std::cout << "  spatial[" << i << "]: weight=" << spatial[i].weight << " dirs=[";
    for (size_t j = 0; j < spatial[i].directions.size(); j++) {
      std::cout << (int)spatial[i].directions[j] << (j + 1 < spatial[i].directions.size() ? "," : "");
    }
    std::cout << "]" << std::endl;
    total_spatial_weight += spatial[i].weight;
  }
  std::cout << "Total spatial weight:   " << total_spatial_weight << std::endl;

  std::cout << "Num valid configs:      " << configs.size() << std::endl;
  double total_config_weight = 0;
  for (auto const &c : configs) total_config_weight += c.weight;
  std::cout << "Total config weight:    " << total_config_weight << std::endl;

  std::vector<double> taus = {0.1, 0.4, 0.7, 1.2, 1.5, 1.8};
  double val = diagram.evaluate(taus, solver, false);
  std::cout << "evaluate(taus):         " << val << std::endl;
  std::cout << "=========================================================\n" << std::endl;

  EXPECT_TRUE(std::isfinite(val));
}

// =====================================================================
// Benchmark: atomic expansion on the dimer (exact diagonalization)
//
// On a 2-site dimer every graph has free multiplicity = 1.
// Using FreeEnergyCalculator<1, double> with override_fm=1,
// we sum all order-4 diagrams over the full simplex and compare
// the infinite-U coefficient against the ED result from
// analytical/benchmark_atomic_expansion.py.
// =====================================================================

TEST(DiagramBenchmark, AtomicExpansionDimerOrder4InfiniteU) {
  Parameters<double> params{8.0, 2.0, 3.0, 0.0, true};
  FreeEnergyCalculator<1, double> calculator(params, 4, /*override_fm=*/1);
  auto [abs_coeff, signed_coeff] = calculator.compute_infinite_U_coefficient();

  EXPECT_NEAR(signed_coeff, -4.0904630472238777e-04, 1e-12);
}
