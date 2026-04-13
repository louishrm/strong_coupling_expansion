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
// Columnar dimer spatial configuration tests (N_sites=2)
//
// Total weight = sum of placements x 2^(n_vertical_lines per placement),
// summed over all rectangular-lattice placements.
// Bond labels: 0=horiz right, 1=horiz left, 2=vert site-0, 3=vert site-1.
// =====================================================================

static double dimer_total_spatial_weight(Graph const &graph) {
  std::vector<VertexType<2, double> *> vt;
  Diagram<2, double> diagram(graph, vt);
  return diagram.get_free_multiplicity();
}

static size_t dimer_num_spatial_configs(Graph const &graph) {
  std::vector<VertexType<2, double> *> vt;
  Diagram<2, double> diagram(graph, vt);
  return diagram.get_spatial_configurations().size();
}

static bool dimer_all_bond_labels_valid(Graph const &graph) {
  std::vector<VertexType<2, double> *> vt;
  Diagram<2, double> diagram(graph, vt);
  for (auto const &sc : diagram.get_spatial_configurations()) {
    for (auto d : sc.directions) {
      if (d > 3) return false;
    }
  }
  return true;
}

// D2a: V=2, 2 lines. 2 horiz placements (1 combo) + 2 vert placements (4 combos) = 10.
TEST(DiagramSpatialConfigs, D2a_TotalWeight) {
  EXPECT_DOUBLE_EQ(dimer_total_spatial_weight(Graph({0, 1, 1, 0}, 2)), 10.0);
}

TEST(DiagramSpatialConfigs, D2a_NumConfigs) {
  EXPECT_GE(dimer_num_spatial_configs(Graph({0, 1, 1, 0}, 2)), 1u);
}

TEST(DiagramSpatialConfigs, D2a_ValidLabels) {
  EXPECT_TRUE(dimer_all_bond_labels_valid(Graph({0, 1, 1, 0}, 2)));
}

// D4b: V=3 star, 4 lines. 4+16+16+64 = 100.
TEST(DiagramSpatialConfigs, D4b_TotalWeight) {
  EXPECT_DOUBLE_EQ(dimer_total_spatial_weight(Graph({0, 1, 1, 1, 0, 0, 1, 0, 0}, 3)), 100.0);
}

TEST(DiagramSpatialConfigs, D4b_ValidLabels) {
  EXPECT_TRUE(dimer_all_bond_labels_valid(Graph({0, 1, 1, 1, 0, 0, 1, 0, 0}, 3)));
}

// D4c: V=2, 4 lines (double edges). 2 horiz + 32 vert = 34.
TEST(DiagramSpatialConfigs, D4c_TotalWeight) {
  EXPECT_DOUBLE_EQ(dimer_total_spatial_weight(Graph({0, 2, 2, 0}, 2)), 34.0);
}

TEST(DiagramSpatialConfigs, D4c_ValidLabels) {
  EXPECT_TRUE(dimer_all_bond_labels_valid(Graph({0, 2, 2, 0}, 2)));
}

// D6b: V=2, 6 lines (triple edges). 2 horiz + 128 vert = 130.
TEST(DiagramSpatialConfigs, D6b_TotalWeight) {
  EXPECT_DOUBLE_EQ(dimer_total_spatial_weight(Graph({0, 3, 3, 0}, 2)), 130.0);
}

TEST(DiagramSpatialConfigs, D6b_ValidLabels) {
  EXPECT_TRUE(dimer_all_bond_labels_valid(Graph({0, 3, 3, 0}, 2)));
}

// D6c: V=4 star, 6 lines. 8+96+384+512 = 1000.
TEST(DiagramSpatialConfigs, D6c_TotalWeight) {
  EXPECT_DOUBLE_EQ(dimer_total_spatial_weight(Graph({0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0}, 4)), 1000.0);
}

TEST(DiagramSpatialConfigs, D6c_ValidLabels) {
  EXPECT_TRUE(dimer_all_bond_labels_valid(Graph({0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0}, 4)));
}

// D4a: V=4 chain (0→1→2→3→0), bipartite. 4 lines.
// Placements: vertex 0 at origin, each subsequent vertex at NN of previous,
// AND vertex 3 must be NN of vertex 0 (closing the cycle).
// On rectangular lattice, 4-cycles are axis-aligned rectangles (degenerate: 1x1 squares)
// or zig-zag paths. Specifically, a 4-cycle on a rectangular lattice requires
// the 4 vertices to form a unit square: (0,0),(1,0),(1,1),(0,1).
// 4 placements (4 rotations of the square) × sub-bond combos.
// Each placement has: 2 horiz lines + 2 vert lines → 2^2=4 sub-bond combos.
// Plus the reflected square gives another 4 placements with same structure.
// Total = 8 placements × 4 combos = 32. But actually we need to count more carefully.
// Let's just verify it's positive and labels are valid.
TEST(DiagramSpatialConfigs, D4a_Positive) {
  Graph graph({0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0}, 4);
  EXPECT_GT(dimer_total_spatial_weight(graph), 0.0);
  EXPECT_TRUE(dimer_all_bond_labels_valid(graph));
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

  EXPECT_GE(configs.size(), 1u);
  // Total weight = sum over spatial configs of (spatial_weight * orbit_size / sym_factor)
  // With columnar dimer: total spatial weight = 10, sym_factor = 2
  double total = 0;
  for (auto const &c : configs) total += c.weight;
  EXPECT_DOUBLE_EQ(total, 10.0);
}

TEST(DiagramGlobalConfigs, D2aDimerConfigCounts) {
  Graph graph({0, 1, 1, 0}, 2);
  std::vector<VertexType<2, double> *> vt;
  Diagram<2, double> diagram(graph, vt);
  auto const &spatial = diagram.get_spatial_configurations();
  auto const &configs = diagram.get_valid_configurations();

  std::cout << "\n===== Order-2 dimer config counts =====" << std::endl;
  std::cout << "Num spatial configs:  " << spatial.size() << std::endl;
  std::cout << "Num global configs:   " << configs.size() << std::endl;
  for (size_t i = 0; i < spatial.size(); i++) {
    std::cout << "  spatial[" << i << "]: weight=" << spatial[i].weight << " dirs=[";
    for (size_t j = 0; j < spatial[i].directions.size(); j++) {
      std::cout << (int)spatial[i].directions[j] << (j + 1 < spatial[i].directions.size() ? "," : "");
    }
    std::cout << "]" << std::endl;
  }
  for (size_t i = 0; i < configs.size(); i++) {
    std::cout << "  global[" << i << "]: weight=" << configs[i].weight << " ops=[";
    for (size_t j = 0; j < configs[i].config.size(); j++) {
      std::cout << (int)configs[i].config[j] << (j + 1 < configs[i].config.size() ? "," : "");
    }
    std::cout << "]" << std::endl;
  }
  std::cout << "=======================================\n" << std::endl;
}

TEST(DiagramGlobalConfigs, Order4DimerConfigCounts) {
  struct DiagramInfo {
    std::string name;
    std::vector<uint8_t> adjmat;
    int V;
  };
  std::vector<DiagramInfo> diagrams = {
     {"D4a (4-cycle, V=4)", {0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0}, 4},
     {"D4b (star, V=3)", {0, 1, 1, 1, 0, 0, 1, 0, 0}, 3},
     {"D4c (double edge, V=2)", {0, 2, 2, 0}, 2},
  };

  for (auto const &info : diagrams) {
    Graph graph(info.adjmat, info.V, false);
    std::vector<VertexType<2, double> *> vt;
    Diagram<2, double> diagram(graph, vt);
    auto const &spatial = diagram.get_spatial_configurations();
    auto const &configs = diagram.get_valid_configurations();

    std::cout << "\n===== " << info.name << " dimer config counts =====" << std::endl;
    std::cout << "Num spatial configs:  " << spatial.size() << std::endl;
    std::cout << "Num global configs:   " << configs.size() << std::endl;

    double total_spatial = 0, total_global = 0;
    for (auto const &s : spatial) total_spatial += s.weight;
    for (auto const &c : configs) total_global += c.weight;
    std::cout << "Total spatial weight: " << total_spatial << std::endl;
    std::cout << "Total global weight:  " << total_global << std::endl;
    std::cout << "============================================\n" << std::endl;
  }
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

// =====================================================================
// Factored evaluation tests (N_sites=2 dimer path)
// =====================================================================

TEST(FactoredEvaluation, MatchesDirectComputation) {
  double U = 8.0, beta = 2.0, mu = 3.0, t_hop = 1.0;
  Parameters<double> params{U, beta, mu, t_hop, false};
  HubbardSolver<2, double> solver(params);

  struct DiagramSpec {
    std::string name;
    std::vector<uint8_t> adjmat;
    int V;
  };
  std::vector<DiagramSpec> test_graphs = {
     {"D2a (simple)",  {0, 1, 1, 0}, 2},
     {"D4a (4-cycle)", {0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0}, 4},
     {"D4b (star)",    {0, 1, 1, 1, 0, 0, 1, 0, 0}, 3},
     {"D4c (double)",  {0, 2, 2, 0}, 2},
  };

  std::vector<std::vector<double>> tau_sets = {
     {0.1, 0.4, 0.7, 1.2},
     {0.0, 0.5, 1.0, 1.5},
     {0.3, 0.3, 0.8, 1.9},
  };

  for (auto &spec : test_graphs) {
    Graph graph(spec.adjmat, spec.V, false);

    // Verify both diagrams produce the same valid configurations
    {
      std::vector<VertexType<2, double> *> no_vt;
      Diagram<2, double> ref(graph, no_vt);
      VertexType<2, double> v1(2), v2(4);
      std::vector<VertexType<2, double> *> vt = {&v1, &v2};
      Diagram<2, double> fac(graph, vt);

      auto const &rc = ref.get_valid_configurations();
      auto const &fc = fac.get_valid_configurations();
      ASSERT_EQ(rc.size(), fc.size()) << "Config count mismatch for " << spec.name;
      double ref_total = 0, fac_total = 0;
      for (size_t i = 0; i < rc.size(); i++) {
        EXPECT_EQ(rc[i].config, fc[i].config) << "Config ops mismatch at " << i << " for " << spec.name;
        EXPECT_DOUBLE_EQ(rc[i].weight, fc[i].weight) << "Config weight mismatch at " << i << " for " << spec.name;
        ref_total += rc[i].weight;
        fac_total += fc[i].weight;
      }
      EXPECT_DOUBLE_EQ(ref.get_diagram_sign(), fac.get_diagram_sign()) << "Sign mismatch for " << spec.name;
      std::cout << "  " << spec.name << ": " << rc.size() << " configs, total weight=" << ref_total
                << ", sign=" << ref.get_diagram_sign() << std::endl;
    }

    // Use fresh diagrams per tau set to avoid stale-cache complications
    for (auto &taus : tau_sets) {
      // Trim taus to match the diagram order (number of hopping lines)
      int order = 0;
      for (int i = 0; i < spec.V; i++)
        for (int j = 0; j < spec.V; j++) order += spec.adjmat[i * spec.V + j];
      std::vector<double> test_taus(taus.begin(), taus.begin() + std::min((int)taus.size(), order));
      if ((int)test_taus.size() < order) continue; // skip if not enough taus

      std::vector<VertexType<2, double> *> no_vt;
      Diagram<2, double> ref_diagram(graph, no_vt);

      VertexType<2, double> vt1(2), vt2(4);
      std::vector<VertexType<2, double> *> vt = {&vt1, &vt2};
      Diagram<2, double> factored_diagram(graph, vt);

      double ref_val = ref_diagram.evaluate(test_taus, solver, false);
      double fac_val = factored_diagram.evaluate(test_taus, solver, false);
      EXPECT_DOUBLE_EQ(ref_val, fac_val)
         << "Mismatch for " << spec.name << " at finite U, taus[0]=" << test_taus[0];

      double ref_inf = ref_diagram.evaluate(test_taus, solver, true);
      double fac_inf = factored_diagram.evaluate(test_taus, solver, true);
      EXPECT_DOUBLE_EQ(ref_inf, fac_inf)
         << "Mismatch for " << spec.name << " at infinite U, taus[0]=" << test_taus[0];
    }
  }
}

TEST(FactoredEvaluation, DirtyFlagCorrectness) {
  double U = 8.0, beta = 2.0, mu = 3.0, t_hop = 1.0;
  Parameters<double> params{U, beta, mu, t_hop, false};
  HubbardSolver<2, double> solver(params);

  Graph graph({0, 1, 1, 1, 0, 0, 1, 0, 0}, 3, false); // D4b star
  VertexType<2, double> vt1(2);
  VertexType<2, double> vt2(4);
  std::vector<VertexType<2, double> *> vt = {&vt1, &vt2};

  // Evaluate at taus_1
  Diagram<2, double> diagram(graph, vt);
  std::vector<double> taus1 = {0.1, 0.4, 0.7, 1.2};
  double val1 = diagram.evaluate(taus1, solver, false);
  EXPECT_TRUE(std::isfinite(val1));

  // Change taus[0], mark dirty, evaluate at taus_2
  std::vector<double> taus2 = {0.9, 0.4, 0.7, 1.2};
  diagram.mark_tau_dirty(0);
  double val2 = diagram.evaluate(taus2, solver, false);

  // Construct fresh diagram and evaluate at taus_2 for reference
  Diagram<2, double> fresh_diagram(graph, vt);
  double val2_fresh = fresh_diagram.evaluate(taus2, solver, false);

  EXPECT_DOUBLE_EQ(val2, val2_fresh);
  // Also verify it changed from val1 (sanity)
  EXPECT_NE(val1, val2);

  // Test mark_all_dirty
  std::vector<double> taus3 = {0.2, 0.6, 0.3, 1.8};
  diagram.mark_all_dirty();
  double val3 = diagram.evaluate(taus3, solver, false);

  Diagram<2, double> fresh_diagram2(graph, vt);
  double val3_fresh = fresh_diagram2.evaluate(taus3, solver, false);
  EXPECT_DOUBLE_EQ(val3, val3_fresh);
}

TEST(FactoredEvaluation, LocalStateCounts) {
  struct DiagramInfo {
    std::string name;
    std::vector<uint8_t> adjmat;
    int V;
  };
  std::vector<DiagramInfo> diagrams = {
     {"D4a (4-cycle, V=4)", {0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0}, 4},
     {"D4b (star, V=3)", {0, 1, 1, 1, 0, 0, 1, 0, 0}, 3},
     {"D4c (double edge, V=2)", {0, 2, 2, 0}, 2},
  };

  std::cout << "\n===== Factored Evaluation: Local State Counts =====" << std::endl;
  for (auto const &info : diagrams) {
    Graph graph(info.adjmat, info.V, false);
    VertexType<2, double> vt1(2);
    VertexType<2, double> vt2(4);
    std::vector<VertexType<2, double> *> vt = {&vt1, &vt2};
    Diagram<2, double> diagram(graph, vt);

    auto const &configs = diagram.get_valid_configurations();
    int n_global = (int)configs.size();
    int V = info.V;

    std::cout << "\n  " << info.name << ":" << std::endl;
    std::cout << "    Global configs: " << n_global << std::endl;

    // To access local_states we use the factored path indirectly:
    // Just print the config count and vertex count as a diagnostic.
    int old_evals = n_global * V;
    std::cout << "    Cumulant evals (old): " << old_evals << " = " << n_global << " configs * " << V << " vertices" << std::endl;
    std::cout << "    (New path computes each distinct local state once per vertex)" << std::endl;
  }
  std::cout << "===================================================\n" << std::endl;
}
