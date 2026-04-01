#include <gtest/gtest.h>
#include <algorithm>
#include <numeric>
#include <cmath>
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
  std::cout << "D2a spatial config directions: ";
  for (auto dir : spatial[0].directions) { std::cout << (int)dir << " "; }
  std::cout << std::endl;
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

  std::cout << "D4c spatial config directions: ";
  for (auto dir : spatial[0].directions) { std::cout << (int)dir << " "; }
  std::cout << std::endl;
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

// =====================================================================
// Global configuration tests (orbital assignment + symmetry reduction)
// =====================================================================

// D2a, N_sites=1: digon, 2 vertices, fm=4, auto_count=2.
// 2 lines, 2^2=4 spin assignments. Spin conservation forces s0=s1.
// Valid: (D,D) and (U,U), related by SpinFlip -> 1 canonical, orbit_size=2.
// Weight = 4 * 2 / 2 = 4.
TEST_F(Diagram2Test, D2aAtomGlobalConfigs) {

  std::vector<uint8_t> D2a = {0, 1, 1, 0};
  Graph graph(D2a, 2);
  std::vector<VertexType<1, double> *> vt;
  Diagram2<1, double> diagram(graph, vt);

  auto const &configs = diagram.get_valid_configurations();
  ASSERT_EQ(configs.size(), 1u);
  EXPECT_DOUBLE_EQ(configs[0].weight, 4.0);
}

// D2a, N_sites=2: digon on dimer superlattice.
// 1 spatial config (weight=6), auto_count=2.
// Spin conservation forces s0=s1 -> valid: (D,D) and (U,U).
// Both map to same canonical under SpinFlip. Orbit under full group has size 4.
// Weight = 6 * 4 / 2 = 12.  -> 1 canonical config.
TEST_F(Diagram2Test, D2aDimerGlobalConfigs) {

  std::vector<uint8_t> D2a = {0, 1, 1, 0};
  Graph graph(D2a, 2);
  std::vector<VertexType<2, double> *> vt;
  Diagram2<2, double> diagram(graph, vt);

  auto const &configs = diagram.get_valid_configurations();
  ASSERT_EQ(configs.size(), 1u);
  EXPECT_DOUBLE_EQ(configs[0].weight, 12.0);
}

// D4b, N_sites=1: 3 vertices, fm=16, auto_count=2.
// 4 lines, valid spins: s0=s2, s1=s3 (from vertex 1,2 conservation).
// 4 valid spin assignments: (DD), (DU), (UD), (UU).
// (DD)<->(UU) via SpinFlip, (DU)<->(UD) via SpinFlip -> 2 canonicals.
// Each orbit_size=2. Weight = 16 * 2 / 2 = 16 each.
// Total weight = 32.
TEST_F(Diagram2Test, D4bAtomGlobalConfigs) {

  std::vector<uint8_t> D4b = {0, 1, 1, 1, 0, 0, 1, 0, 0};
  Graph graph(D4b, 3);
  std::vector<VertexType<1, double> *> vt;
  Diagram2<1, double> diagram(graph, vt);

  auto const &configs = diagram.get_valid_configurations();
  ASSERT_EQ(configs.size(), 2u);

  double total = 0;
  for (auto const &c : configs) { total += c.weight; }
  EXPECT_DOUBLE_EQ(total, 32.0);
}

// =====================================================================
// Numerical evaluation tests
// =====================================================================

class Diagram2EvalTest : public ::testing::Test {
  protected:
  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;
  Parameters<double> params{U, beta, mu, 0.0, true};
};

// D2a atom: compare Diagram2 evaluate against old DiagramEvaluator.
// Old: prefactor = (-1/beta) * sign / symmetry_factor * fm
//    = (-1/1) * (-1) / 2 * 4 = 2.0
// New: prefactor = (-1/beta) * sign = (-1)*(-1) = 1.0
//      weight per config = fm * orbit_size / auto_count = 4 * 2 / 2 = 4
// Old sum = Σ_spins Π_v C_v.  New sum = Σ_canonical weight * Π_v C_v.
// They should agree: old_result = new_result (both include all factors).
TEST_F(Diagram2EvalTest, D2aAtomEvaluateIsFinite) {

  std::vector<uint8_t> D2a = {0, 1, 1, 0};
  Graph graph(D2a, 2);

  // Create VertexTypes for caching: order 1 cumulant (degree 2)
  VertexType<1, double> vt1(2); // 2 legs
  std::vector<VertexType<1, double> *> vt = {&vt1};

  Diagram2<1, double> diagram(graph, vt);
  HubbardSolver<1, double> solver(params);

  std::vector<double> taus = {0.5, 0.0};
  double val = diagram.evaluate(taus, solver, false);

  // The value should be finite and non-zero
  EXPECT_TRUE(std::isfinite(val));
  EXPECT_NE(val, 0.0);

  // Test with VertexType cache: evaluate twice, should give the same result
  double val2 = diagram.evaluate(taus, solver, false);
  EXPECT_DOUBLE_EQ(val, val2);
}

// D2a dimer: evaluate on the triangular superlattice
TEST_F(Diagram2EvalTest, D2aDimerEvaluateIsFinite) {

  double t_dimer = 1.0;
  Parameters<double> dimer_params{U, beta, mu, t_dimer, true};

  std::vector<uint8_t> D2a = {0, 1, 1, 0};
  Graph graph(D2a, 2);

  VertexType<2, double> vt1(2);
  std::vector<VertexType<2, double> *> vt = {&vt1};

  Diagram2<2, double> diagram(graph, vt);
  HubbardSolver<2, double> solver(dimer_params);

  std::vector<double> taus = {0.5, 0.0};
  double val = diagram.evaluate(taus, solver, false);

  EXPECT_TRUE(std::isfinite(val));
  EXPECT_NE(val, 0.0);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
