#include <gtest/gtest.h>
#include "sc_expansion/generate_diagrams.hpp"
#include <vector>
#include <iostream>

using namespace sc_expansion;

TEST(GenerateDiagramsTest, Order2Diagrams) {
  VacuumDiagramGenerator gen(2);
  gen.generate();

  const auto &unique_graphs = gen.get_unique_graphs();
  EXPECT_EQ(unique_graphs.size(), 1);
}

TEST(GenerateDiagramsTest, Order4Diagrams) {
  VacuumDiagramGenerator gen(4);
  gen.generate();

  const auto &unique_graphs = gen.get_unique_graphs();
  EXPECT_EQ(unique_graphs.size(), 3);
}

TEST(GenerateDiagramsTest, Order6Diagrams) {
  VacuumDiagramGenerator gen(6);
  gen.generate();

  const auto &unique_graphs = gen.get_unique_graphs();
  EXPECT_EQ(unique_graphs.size(), 7);
}

TEST(GenerateDiagramsTest, Order8Diagrams) {
  VacuumDiagramGenerator gen(8);
  gen.generate();

  const auto &unique_graphs = gen.get_unique_graphs();
  std::cout << "Detected " << unique_graphs.size() << " unique bipartite topologies of order 8." << std::endl;
}
