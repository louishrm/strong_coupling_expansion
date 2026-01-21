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

  EXPECT_EQ(diagram_a.get_symmetry_factor(), 4);
  EXPECT_EQ(diagram_b.get_symmetry_factor(), 2);
  EXPECT_EQ(diagram_c.get_symmetry_factor(), 8);
}

TEST_F(DiagramTest, DiagramFreeMultiplicityIsCorrect) {

  adjmat D4a = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}; //4-cycle
  adjmat D4b = {{0, 1, 1}, {1, 0, 0}, {1, 0, 0}};                        //3-cycle with double lines
  adjmat D4c = {{0, 2}, {2, 0}};                                         //2-cycle with double lines

  int free_multiplicity_a = compute_free_multiplicity(D4a, 4);
  int free_multiplicity_b = compute_free_multiplicity(D4b, 4);
  int free_multiplicity_c = compute_free_multiplicity(D4c, 4);

  EXPECT_EQ(free_multiplicity_a, 36);
  EXPECT_EQ(free_multiplicity_b, 16);
  EXPECT_EQ(free_multiplicity_c, 4);
}