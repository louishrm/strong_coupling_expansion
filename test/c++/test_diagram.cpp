#include <gtest/gtest.h>
#include <cmath>
#include "sc_expansion/hubbard_atom.hpp"
#include "sc_expansion/cumulant.hpp"
#include "sc_expansion/diagram.hpp"

using namespace sc_expansion;

// The Test Fixture: Sets up the data common to all tests
class DiagramTest : public ::testing::Test {
  protected:
  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;

  std::unique_ptr<HubbardAtom> atom;

  void SetUp() override { atom = std::make_unique<HubbardAtom>(U, beta, mu); }
};

TEST_F(DiagramTest, SymmetryFactorOfDiagramIsCorrect) {

  adjmat D4a = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}; //4-cycle
  adjmat D4b = {{0, 1, 1}, {1, 0, 0}, {1, 0, 0}};                        //3-cycle with double lines
  adjmat D4c = {{0, 2}, {2, 0}};                                         //2-cycle with double lines
  Diagram diagram_a(D4a, U, beta, mu);
  Diagram diagram_b(D4b, U, beta, mu);
  Diagram diagram_c(D4c, U, beta, mu);

  auto lines = diagram_b.get_hopping_lines();
  for (const auto &line : lines) { std::cout << "Line from vertex " << line.from_vertex << " to vertex " << line.to_vertex << std::endl; }

  EXPECT_EQ(diagram_a.get_symmetry_factor(), 4);
  EXPECT_EQ(diagram_b.get_symmetry_factor(), 2);
  EXPECT_EQ(diagram_c.get_symmetry_factor(), 8);
}

TEST_F(DiagramTest, InfiniteUDiagramConstantInSimplex) {

  adjmat D4a = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}; //4-cycle
  Diagram diagram_4a(D4a, U, beta, mu);

  std::vector<double> taus_41 = {0.1, 0.2, 0.3, 0.4};
  std::vector<double> taus_42 = {0.15, 0.23, 0.31, 0.76};

  double val_41 = diagram_4a.evaluate_at_taus(taus_41, true);
  double val_42 = diagram_4a.evaluate_at_taus(taus_42, true);

  adjmat D6c = {{0, 1, 1, 1}, {1, 0, 0, 0}, {1, 0, 0, 0}, {1, 0, 0, 0}};
  Diagram diagram_6c(D6c, U, beta, mu);

  std::vector<double> taus_61 = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
  std::vector<double> taus_62 = {0.15, 0.23, 0.31, 0.45, 0.52, 0.78};

  double val_61 = diagram_6c.evaluate_at_taus(taus_61, true);
  double val_62 = diagram_6c.evaluate_at_taus(taus_62, true);
  EXPECT_DOUBLE_EQ(val_61, val_62);
  EXPECT_DOUBLE_EQ(val_41, val_42);
}