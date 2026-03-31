#include <gtest/gtest.h>
#include <cmath>
#include <numeric>
#include <memory>
#include "../../c++/sc_expansion/hubbard_solver.hpp"
#include "../../c++/sc_expansion/diagram.hpp"
#include "../../c++/sc_expansion/diagram2.hpp"
#include "../../c++/sc_expansion/graph.hpp"

using namespace sc_expansion;

class DiagramDimerTest : public ::testing::Test {
  protected:
  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;
  double t    = 1.0;
  Parameters<double> params{U, beta, mu, t, true};
};

TEST_F(DiagramDimerTest, GeometricalOrbits) {
  // C1C2C1 graph on triangular lattice (not bipartite)
  std::vector<uint8_t> D4b = {0, 1, 1, 1, 0, 0, 1, 0, 0};
  Graph g(D4b, 3, false);
  std::vector<VertexType<2, double>*> vtypes;
  Diagram2<2, double> diag(g, vtypes);

  auto embeddings = diag.find_spatial_embeddings();
  auto orbits = diag.group_spatial_embeddings(embeddings);

  EXPECT_EQ(embeddings.size(), 36);
  EXPECT_EQ(orbits.size(), 2);

  for (const auto& o : orbits) {
      EXPECT_EQ(o.spatial_weight, 18);
  }
}

TEST_F(DiagramDimerTest, DimerSumRule) {
  // D2a on triangular lattice
  std::vector<uint8_t> D2a = {0, 1, 1, 0};
  Graph g(D2a, 2, false);
  std::vector<VertexType<2, double>*> vtypes;
  Diagram2<2, double> diag(g, vtypes);

  double total_weight = 0;
  for (auto const& conf : diag.get_valid_configs()) {
      total_weight += conf.weight;
  }
  
  // Free multiplicity of D2a on triangular lattice:
  // Spatial mult = 6
  // Spin assignments = 2 (up, down)
  // Sym factor = 2
  // Weight = 6 * 2 / 2 = 6
  EXPECT_EQ(total_weight, 6.0);
}

TEST_F(DiagramDimerTest, AtomicSumRule) {
  // D2a on square lattice (atomic)
  std::vector<uint8_t> D2a = {0, 1, 1, 0};
  Graph g(D2a, 2, true);
  std::vector<VertexType<1, double>*> vtypes;
  Diagram2<1, double> diag(g, vtypes);
  
  double total_weight = 0;
  for (auto const& conf : diag.get_valid_configs()) {
      total_weight += conf.weight;
  }
  
  // Free multiplicity = 4 (square lattice)
  // Spin assignments = 2 (up, down)
  // Sym factor = 2
  // Weight = 4 * 2 / 2 = 4
  EXPECT_EQ(total_weight, 4.0);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
