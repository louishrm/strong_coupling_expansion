#include <gtest/gtest.h>
#include <cmath>
#include <numeric>
#include <memory>
#include "../c++/sc_expansion/hubbard_solver.hpp"
#include "../c++/sc_expansion/diagram.hpp"
#include "../c++/sc_expansion/graph.hpp"

using namespace sc_expansion;

class DiagramDimerTest : public ::testing::Test {
  protected:
  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;
  double t    = 1.0;
  Parameters<double> params{U, beta, mu, t, true};
};

TEST_F(DiagramDimerTest, DimerDiagramEvaluation) {
  // A simple diagram: 2nd order (2-cycle)
  // D2a = {0, 1, 1, 0}
  std::vector<uint8_t> D2a = {0, 1, 1, 0};
  Graph g(D2a, 2);
  Diagram d(g);
  
  // Evaluator for N_sites=2 (dimer)
  DiagramEvaluator<2, double> eval(d, params);

  // For a 2nd order diagram, we have 2 taus.
  // evaluate_at_taus computes the diagram value.
  // The user requested that we treat a vertex as a dimer.
  // In the evaluator, we fixed site=0 for the operators.
  
  std::vector<double> taus = {0.5, 0.0};
  double val = eval.evaluate_at_taus(taus, false, false);
  
  // This value should be finite and reproducible.
  EXPECT_TRUE(std::isfinite(val));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
