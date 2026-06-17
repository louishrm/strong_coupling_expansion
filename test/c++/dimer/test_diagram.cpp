#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "../c++/sc_expansion/graph.hpp"
#include "../c++/sc_expansion/hubbard_solver.hpp"
#include "../c++/sc_expansion/dimer/diagram.hpp"
#include "../c++/sc_expansion/dimer/vertex.hpp"

using namespace sc_expansion;
using namespace sc_expansion::dimer;

// ============================================================================
//  Spatial-configuration counts for representative low-order diagrams on the
//  staggered (triangular) dimer superlattice.
// ============================================================================

TEST(DiagramSpatialConfigs, D4b) {
  Graph graph({0, 1, 1, 1, 0, 0, 1, 0, 0}, 3, false);
  std::vector<VertexType<double> *> vt;
  Diagram<double> diagram(graph, vt);
  auto const &spatial = diagram.get_spatial_configurations();

  EXPECT_EQ(spatial.size(), 2u);
  double total = 0;
  for (auto const &sc : spatial) total += sc.weight;
  EXPECT_DOUBLE_EQ(total, 36.0);
}

TEST(DiagramSpatialConfigs, D2a) {
  Graph graph({0, 1, 1, 0}, 2);
  std::vector<VertexType<double> *> vt;
  Diagram<double> diagram(graph, vt);
  auto const &spatial = diagram.get_spatial_configurations();

  EXPECT_EQ(spatial.size(), 1u);
  EXPECT_DOUBLE_EQ(spatial[0].weight, 6.0);
}

TEST(DiagramSpatialConfigs, D4c) {
  Graph graph({0, 2, 2, 0}, 2);
  std::vector<VertexType<double> *> vt;
  Diagram<double> diagram(graph, vt);
  auto const &spatial = diagram.get_spatial_configurations();

  EXPECT_EQ(spatial.size(), 1u);
  EXPECT_DOUBLE_EQ(spatial[0].weight, 6.0);
}

TEST(DiagramSpatialConfigs, D6c) {
  Graph graph({0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0}, 4);
  std::vector<VertexType<double> *> vt;
  Diagram<double> diagram(graph, vt);
  auto const &spatial = diagram.get_spatial_configurations();

  EXPECT_EQ(spatial.size(), 2u);
  EXPECT_DOUBLE_EQ(spatial[0].weight, 54.0);
  EXPECT_DOUBLE_EQ(spatial[1].weight, 162.0);
}

// ============================================================================
//  Global (spin-resolved) configurations.
// ============================================================================

TEST(DiagramGlobalConfigs, D2aDimer) {
  Graph graph({0, 1, 1, 0}, 2);
  std::vector<VertexType<double> *> vt;
  Diagram<double> diagram(graph, vt);
  auto const &configs = diagram.get_valid_configurations();

  EXPECT_EQ(configs.size(), 1u);
  // weight = spatial_weight * orbit_size / sym_factor = 6 * 2 / 2 = 6
  // orbit_size = 2 from SpinFlip only (lattice symmetries already absorbed
  // into spatial-config canonicalization).
  EXPECT_DOUBLE_EQ(configs[0].weight, 6.0);
}

// ============================================================================
//  End-to-end evaluation: a finite, non-zero diagram value at a sample (τ, U).
// ============================================================================

TEST(DiagramEvaluation, D2aDimer) {
  double U = 8.0, beta = 1.0, mu = 2.0;
  Parameters<double> params{U, beta, mu, 1.0, true};
  std::vector<double> taus = {0.5, 0.0};

  Graph graph({0, 1, 1, 0}, 2);
  VertexType<double> vt1(2);
  std::vector<VertexType<double> *> vt = {&vt1};
  Diagram<double> diagram(graph, vt);
  HubbardSolver<2, double> solver(params);

  double val = diagram.evaluate(taus, solver);
  EXPECT_TRUE(std::isfinite(val));
  EXPECT_NE(val, 0.0);
}
