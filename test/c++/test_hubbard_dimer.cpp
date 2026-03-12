#include <gtest/gtest.h>
#include <iostream>
#include "../../c++/sc_expansion/hubbard_dimer.hpp"

using namespace sc_expansion;

TEST(HubbardDimerTest, PartitionFunctionComputation) {
  // Model parameters
  double t    = 1.0;
  double U    = 8.0;
  double beta = 1.0;
  double mu   = 4.0;

  // Instantiate the class, which prints Z to std::cout
  HubbardDimer<double> dimer(t, U, beta, mu);

  // Basic check for physical value: Z should be positive
  EXPECT_GT(dimer.Z, 0.0);
  
  // The test essentially verifies that instantiation is successful
  // and that Z is computed (it is printed by the constructor itself).
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
