#include <gtest/gtest.h>
#include <algorithm>
#include <numeric>
#include "sc_expansion/diagram2.hpp"
#include "sc_expansion/graph.hpp"

using namespace sc_expansion;

class Diagram2Test : public ::testing::Test {};

// Helper: create a Diagram2<1, double> and return its free multiplicity.
static double single_site_free_multiplicity(std::vector<uint8_t> const &adjmat, int V, bool bipartite_only = true) {
  Graph graph(adjmat, V, bipartite_only);
  std::vector<VertexType<1, double> *> vt;
  Diagram2<1, double> diagram(graph, vt);
  return diagram.get_free_multiplicity();
}

// =====================================================================
// Free multiplicity tests (N_sites=1, single-site expansion)
// These were moved from test_graph.cpp into Diagram2.
// =====================================================================

TEST_F(Diagram2Test, FreeMultiplicityIsCorrect) {

  // D2a (digon)
  std::vector<uint8_t> D2a = {0, 1, 1, 0};
  EXPECT_EQ((int)single_site_free_multiplicity(D2a, 2), 4);

  // D3a (3-cycle, non-bipartite)
  std::vector<uint8_t> D3a = {0, 1, 0, 0, 0, 1, 1, 0, 0};
  EXPECT_EQ((int)single_site_free_multiplicity(D3a, 3, false), 12);

  // D4a (4-cycle)
  std::vector<uint8_t> D4a = {0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0};
  EXPECT_EQ((int)single_site_free_multiplicity(D4a, 4), 36);

  // D4b
  std::vector<uint8_t> D4b = {0, 1, 1, 1, 0, 0, 1, 0, 0};
  EXPECT_EQ((int)single_site_free_multiplicity(D4b, 3), 16);

  // D4c (double digon)
  std::vector<uint8_t> D4c = {0, 2, 2, 0};
  EXPECT_EQ((int)single_site_free_multiplicity(D4c, 2), 4);

  // D6 variants
  std::vector<uint8_t> D6_1 = {0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0};
  EXPECT_EQ((int)single_site_free_multiplicity(D6_1, 4), 64);

  std::vector<uint8_t> D6_2 = {0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  EXPECT_EQ((int)single_site_free_multiplicity(D6_2, 4), 64);

  std::vector<uint8_t> D6a = {0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0};
  EXPECT_EQ((int)single_site_free_multiplicity(D6a, 6), 400);

  std::vector<uint8_t> D6b = {0, 3, 3, 0};
  EXPECT_EQ((int)single_site_free_multiplicity(D6b, 2), 4);

  std::vector<uint8_t> D6c = {0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  EXPECT_EQ((int)single_site_free_multiplicity(D6c, 4), 64);

  std::vector<uint8_t> D6d = {0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0};
  EXPECT_EQ((int)single_site_free_multiplicity(D6d, 5), 144);

  std::vector<uint8_t> D6e = {0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0};
  EXPECT_EQ((int)single_site_free_multiplicity(D6e, 4), 64);

  std::vector<uint8_t> D6f = {0, 2, 1, 2, 0, 0, 1, 0, 0};
  EXPECT_EQ((int)single_site_free_multiplicity(D6f, 3), 16);

  std::vector<uint8_t> D6g = {0, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0};
  EXPECT_EQ((int)single_site_free_multiplicity(D6g, 4), 36);

  // D8a (8-cycle)
  std::vector<uint8_t> D8a(64, 0);
  for (int i = 0; i < 8; ++i) {
    D8a[i * 8 + (i + 1) % 8]   = 1;
    D8a[((i + 1) % 8) * 8 + i] = 1;
  }
  EXPECT_EQ((int)single_site_free_multiplicity(D8a, 8), 4900);

  // D8b (quadruple digon)
  std::vector<uint8_t> D8b = {0, 4, 4, 0};
  EXPECT_EQ((int)single_site_free_multiplicity(D8b, 2), 4);
}

// =====================================================================
// Dimer spatial configuration tests (N_sites=2)
// =====================================================================

// D4b: 3 vertices, adjacency matrix {0,1,1, 1,0,0, 1,0,0}
// Vertex 0 (degree 4) connected to vertex 1 and vertex 2 (both degree 2).
//
// On the staggered-dimer triangular superlattice (N_sites=2):
//   - 6 NN directions, split 3 right / 3 left
//   - Vertex 0 at origin; v1 and v2 each at any of 6 positions (independent) -> 36 embeddings
//   - Direction pattern per embedding: (dir_to_v1, dir_to_v2) in {L,R}^2
//   - Grouping under graph automorphism (swap v1<->v2) + lattice inversion (L<->R):
//       {LL, RR} merge -> weight 9+9 = 18
//       {LR, RL} merge -> weight 9+9 = 18
//   => 2 distinct spatial configurations, weight 18 each.
TEST_F(Diagram2Test, D4bDimerSpatialConfigs) {

  std::vector<uint8_t> D4b = {0, 1, 1, 1, 0, 0, 1, 0, 0};
  Graph graph(D4b, 3);

  std::vector<VertexType<2, double> *> vertex_types; // empty, not needed for spatial config test
  Diagram2<2, double> diagram(graph, vertex_types);

  auto const &spatial = diagram.get_spatial_configurations();

  // Expect exactly 2 distinct spatial configurations
  ASSERT_EQ(spatial.size(), 2u);

  // Both should have weight 18
  std::vector<double> weights;
  for (auto const &sc : spatial) { weights.push_back(sc.weight); }
  std::sort(weights.begin(), weights.end());

  EXPECT_DOUBLE_EQ(weights[0], 18.0);
  EXPECT_DOUBLE_EQ(weights[1], 18.0);

  // Total weight should equal 36 (6 * 6 triangular lattice embeddings)
  double total = std::accumulate(weights.begin(), weights.end(), 0.0);
  EXPECT_DOUBLE_EQ(total, 36.0);
}

// D2a (digon): 2 vertices, adjacency {0,1, 1,0}. Single edge each way.
// On triangular superlattice: 6 embeddings (v1 at any of 6 neighbours of v0).
// 3 right + 3 left -> merged by lattice inversion -> 1 spatial config, weight 6.
TEST_F(Diagram2Test, D2aDimerSpatialConfigs) {

  std::vector<uint8_t> D2a = {0, 1, 1, 0};
  Graph graph(D2a, 2);

  std::vector<VertexType<2, double> *> vertex_types;
  Diagram2<2, double> diagram(graph, vertex_types);

  auto const &spatial = diagram.get_spatial_configurations();

  ASSERT_EQ(spatial.size(), 1u);
  EXPECT_DOUBLE_EQ(spatial[0].weight, 6.0);
}

// D4c (double digon / watermelon): 2 vertices, adjacency {0,2, 2,0}.
// On triangular superlattice: 6 embeddings, all have both lines in the same direction.
// -> 1 spatial config, weight 6.
TEST_F(Diagram2Test, D4cDimerSpatialConfigs) {

  std::vector<uint8_t> D4c = {0, 2, 2, 0};
  Graph graph(D4c, 2);

  std::vector<VertexType<2, double> *> vertex_types;
  Diagram2<2, double> diagram(graph, vertex_types);

  auto const &spatial = diagram.get_spatial_configurations();

  ASSERT_EQ(spatial.size(), 1u);
  EXPECT_DOUBLE_EQ(spatial[0].weight, 6.0);
}

TEST_F(Diagram2Test, D6cDimerSpatialConfigs) {

  std::vector<uint8_t> D6c = {0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  Graph graph(D6c, 4);

  std::vector<VertexType<2, double> *> vertex_types;
  Diagram2<2, double> diagram(graph, vertex_types);

  auto const &spatial = diagram.get_spatial_configurations();

  ASSERT_EQ(spatial.size(), 2u);
  EXPECT_DOUBLE_EQ(spatial[0].weight, 54.0);
  EXPECT_DOUBLE_EQ(spatial[1].weight, 162.0);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
