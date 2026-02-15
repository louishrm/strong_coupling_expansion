#include <gtest/gtest.h>
#include <cmath>
#include "sc_expansion/hubbard_atom.hpp"
#include "sc_expansion/cumulant.hpp"
#include "sc_expansion/diagram.hpp"
#include "sc_expansion/graph.hpp"

using namespace sc_expansion;

// The Test Fixture: Sets up the data common to all tests
class DiagramTest : public ::testing::Test {
  protected:
  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;
  Parameters params{U, beta, mu};

  std::unique_ptr<HubbardAtom> atom;

  void SetUp() override { atom = std::make_unique<HubbardAtom>(U, beta, mu); }
};

TEST_F(DiagramTest, DiagramSignIsCorrect) {
  // D2a = {0, 1, 1, 0} is a 2-cycle (0->1, 1->0).
  // It has 1 loop, so sign should be -1.
  std::vector<uint8_t> D2a = {0, 1, 1, 0};
  Graph g(D2a, 2);
  Diagram diagram(g);
  EXPECT_EQ(diagram.get_diagram_sign(), -1.0);

  std::vector<uint8_t> D4a = {0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0}; // 4-cycle
  std::vector<uint8_t> D4b = {0, 1, 1, 1, 0, 0, 1, 0, 0};                      // 3-cycle with double lines (0->1, 0->2, 1->0, 2->0)
  std::vector<uint8_t> D4c = {0, 2, 2, 0};                                     // 2-cycle with double lines

  Graph g_a(D4a, 4);
  Graph g_b(D4b, 3);
  Graph g_c(D4c, 2);

  Diagram diagram_a(g_a);
  Diagram diagram_b(g_b);
  Diagram diagram_c(g_c);

  EXPECT_EQ(diagram_a.get_diagram_sign(), -1.0); // 4-cycle has 2 loops, so sign should be
  EXPECT_EQ(diagram_b.get_diagram_sign(), 1.0);  // 3-cycle with double lines has 3 loops, so sign should be
  EXPECT_EQ(diagram_c.get_diagram_sign(), 1.0);  // 2-cycle with double lines has 4 loops, so sign should be
}

TEST_F(DiagramTest, Order6DiagramSignIsCorrect) {
  // D6a: 6-cycle, 1 loop -> sign -1
  std::vector<uint8_t> D6a = {0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0};
  // D6b: watermelon triple, 3 loops -> sign -1
  std::vector<uint8_t> D6b = {0, 3, 3, 0};
  // D6c: petal with 4 vertices, 3 loops -> sign -1
  std::vector<uint8_t> D6c = {0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  // D6d: square + digon, 2 loops -> sign 1
  std::vector<uint8_t> D6d = {0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0};
  // D6e: crab diagram, 3 loops -> sign -1
  std::vector<uint8_t> D6e = {0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0};
  // D6f: watermelon double + digon, 3 loops -> sign -1
  std::vector<uint8_t> D6f = {0, 2, 1, 2, 0, 0, 1, 0, 0};
  // D6g: square with one double line, 2 loops -> sign 1
  std::vector<uint8_t> D6g = {0, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0};

  EXPECT_EQ(Diagram(Graph(D6a, 6)).get_diagram_sign(), -1.0);
  EXPECT_EQ(Diagram(Graph(D6b, 2)).get_diagram_sign(), -1.0);
  EXPECT_EQ(Diagram(Graph(D6c, 4)).get_diagram_sign(), -1.0);
  EXPECT_EQ(Diagram(Graph(D6d, 5)).get_diagram_sign(), 1.0);
  EXPECT_EQ(Diagram(Graph(D6e, 4)).get_diagram_sign(), -1.0);
  EXPECT_EQ(Diagram(Graph(D6f, 3)).get_diagram_sign(), -1.0);
  EXPECT_EQ(Diagram(Graph(D6g, 4)).get_diagram_sign(), 1.0);
}

TEST_F(DiagramTest, SymmetryFactorOfDiagramIsCorrect) {

  std::vector<uint8_t> D4a = {0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0}; // 4-cycle
  std::vector<uint8_t> D4b = {0, 1, 1, 1, 0, 0, 1, 0, 0};                      // 3-cycle with double lines (0->1, 0->2, 1->0, 2->0)
  std::vector<uint8_t> D4c = {0, 2, 2, 0};                                     // 2-cycle with double lines

  Graph g_a(D4a, 4);
  Graph g_b(D4b, 3);
  Graph g_c(D4c, 2);

  Diagram diagram_a(g_a);
  Diagram diagram_b(g_b);
  Diagram diagram_c(g_c);

  EXPECT_EQ(diagram_a.get_graph().get_symmetry_factor(), 4);
  EXPECT_EQ(diagram_b.get_graph().get_symmetry_factor(), 2);
  EXPECT_EQ(diagram_c.get_graph().get_symmetry_factor(), 8);
}

TEST_F(DiagramTest, InfiniteUDiagramConstantInSimplex) {

  std::vector<uint8_t> D4a = {0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0}; // 4-cycle
  Graph g4a(D4a, 4);
  Diagram d4a(g4a);
  DiagramEvaluator eval4a(d4a, params);

  std::vector<double> taus_41 = {0.1, 0.2, 0.3, 0.4};
  std::vector<double> taus_42 = {0.15, 0.23, 0.31, 0.76};

  double val_41 = eval4a.evaluate_at_taus(taus_41, true);
  double val_42 = eval4a.evaluate_at_taus(taus_42, true);

  std::vector<uint8_t> D6c = {0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  Graph g6c(D6c, 4);
  Diagram d6c(g6c);
  DiagramEvaluator eval6c(d6c, params);

  std::vector<double> taus_61 = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
  std::vector<double> taus_62 = {0.15, 0.23, 0.31, 0.45, 0.52, 0.78};

  double val_61 = eval6c.evaluate_at_taus(taus_61, true);
  double val_62 = eval6c.evaluate_at_taus(taus_62, true);

  EXPECT_DOUBLE_EQ(val_61, val_62);
  EXPECT_DOUBLE_EQ(val_41, val_42);
}
