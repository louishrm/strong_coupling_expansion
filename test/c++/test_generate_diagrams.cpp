#include <gtest/gtest.h>
#include "sc_expansion/generate_diagrams.hpp"
#include <vector>
#include <iostream>
#include <cmath>

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

  std::vector<std::vector<uint8_t>> expected_mats = {
      {0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0}, // D4a
      {0, 2, 2, 0},                                     // D4c
      {0, 1, 1, 1, 0, 0, 1, 0, 0}                       // D4b
  };

  for (auto const &m : expected_mats) {
    int V                          = std::sqrt(m.size());
    std::vector<uint8_t> canonical = Graph(m, V).get_canonical_form();
    EXPECT_TRUE(unique_graphs.count(canonical)) << "Missing order 4 diagram with V=" << V;
  }
}

TEST(GenerateDiagramsTest, Order6Diagrams) {
  VacuumDiagramGenerator gen(6);
  gen.generate();

  const auto &unique_graphs = gen.get_unique_graphs();
  EXPECT_EQ(unique_graphs.size(), 7);

  std::vector<std::vector<uint8_t>> expected_mats = {
      {0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0}, // D6a
      {0, 3, 3, 0},                                                                                                 // D6b
      {0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0},                                                             // D6c
      {0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},                                  // D6d
      {0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0},                                                             // D6e
      {0, 2, 1, 2, 0, 0, 1, 0, 0},                                                                                  // D6f
      {0, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0}                                                              // D6g
  };

  for (auto const &m : expected_mats) {
    int V                          = std::sqrt(m.size());
    std::vector<uint8_t> canonical = Graph(m, V).get_canonical_form();
    EXPECT_TRUE(unique_graphs.count(canonical)) << "Missing order 6 diagram with V=" << V;
  }
}

TEST(GenerateDiagramsTest, Order8Diagrams) {
  VacuumDiagramGenerator gen(8);
  gen.generate();

  const auto &unique_graphs = gen.get_unique_graphs();
  std::cout << "Detected " << unique_graphs.size() << " unique bipartite topologies of order 8." << std::endl;
}
