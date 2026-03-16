#include <gtest/gtest.h>
#include "sc_expansion/generate_diagrams.hpp"
#include <vector>
#include <iostream>
#include <cmath>

using namespace sc_expansion;

TEST(GenerateDiagramsTest, NCycleAdjacencyMatrix) {
  // Test for n=2
  std::vector<uint8_t> adjmat2   = generate_n_cycle_adjacency_matrix(2);
  std::vector<uint8_t> expected2 = {0, 1, 1, 0};
  EXPECT_EQ(adjmat2, expected2);

  // Test for n=4
  std::vector<uint8_t> adjmat4 = generate_n_cycle_adjacency_matrix(4);
  Graph g4(adjmat4, 4, false);
  EXPECT_EQ(g4.get_canonical_form(), adjmat4);

  //test for n=6
  std::vector<uint8_t> adjmat6 = generate_n_cycle_adjacency_matrix(6);
  Graph g6(adjmat6, 6, false);
  EXPECT_EQ(g6.get_canonical_form(), adjmat6);

  //test for n=3
  std::vector<uint8_t> adjmat3 = generate_n_cycle_adjacency_matrix(3);
  Graph g3(adjmat3, 3, false);
  EXPECT_EQ(g3.get_canonical_form(), adjmat3);

  //test for n=5
  std::vector<uint8_t> adjmat5 = generate_n_cycle_adjacency_matrix(5);
  Graph g5(adjmat5, 5, false);
  EXPECT_EQ(g5.get_canonical_form(), adjmat5);
}

TEST(GenerateDiagramsTest, NCycleFreeMultiplicityIsCorrect) {

  //n=2
  int fm2               = calculate_n_cycle_free_multiplicity(2, true);
  int fm2_non_bipartite = calculate_n_cycle_free_multiplicity(2, false);

  EXPECT_EQ(fm2, 4);
  EXPECT_EQ(fm2_non_bipartite, 6);

  //n=3
  int fm3_non_bipartite = calculate_n_cycle_free_multiplicity(3, false);
  EXPECT_EQ(fm3_non_bipartite, 12);

  //n=4
  int fm4               = calculate_n_cycle_free_multiplicity(4, true);
  int fm4_non_bipartite = calculate_n_cycle_free_multiplicity(4, false);

  EXPECT_EQ(fm4, 36);
  EXPECT_EQ(fm4_non_bipartite, 90);
}

TEST(GenerateDiagramsTest, Order2Diagrams) {
  VacuumDiagramGenerator gen(2);
  gen.generate();

  const auto &graphs = gen.get_unique_graphs();
  EXPECT_EQ(graphs.size(), 1);

  // Check n-cycle optimization for n=2
  EXPECT_EQ(graphs[0].get_symmetry_factor(), 2);
  uint64_t nCn2 = binomial_coefficient(2, 1);
  EXPECT_EQ(graphs[0].get_free_multiplicity(), nCn2 * nCn2);
}

TEST(GenerateDiagramsTest, Order3DiagramsNonBipartite) {
  VacuumDiagramGenerator gen(3, false); // Allow non-bipartite diagrams
  gen.generate();

  const auto &graphs = gen.get_unique_graphs();
  EXPECT_EQ(graphs.size(), 1);
}

TEST(GenerateDiagramsTest, Order4DiagramsNonBipartite) {
  VacuumDiagramGenerator gen(4, false); // Allow non-bipartite diagrams
  gen.generate();

  const auto &graphs = gen.get_unique_graphs();
  EXPECT_EQ(graphs.size(), 3);

  auto has_canonical = [&](const std::vector<uint8_t> &canonical) {
    for (const auto &g : graphs) {
      if (g.get_canonical_form() == canonical) return true;
    }
    return false;
  };

  std::vector<std::vector<uint8_t>> expected_mats = {
     {0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0}, // D4a
     {0, 2, 2, 0},                                     // D4c
     {0, 1, 1, 1, 0, 0, 1, 0, 0}                       // D4b
  };

  for (auto const &m : expected_mats) {
    int V                          = std::sqrt(m.size());
    std::vector<uint8_t> canonical = Graph(m, V).get_canonical_form();
    EXPECT_TRUE(has_canonical(canonical)) << "Missing order 4 diagram with V=" << V;
  }

  // Check n-cycle optimization for n=4 (D4a is the 4-cycle)
  bool found_n_cycle = false;
  for (const auto &g : graphs) {
    if (g.get_V() == 4) {
      // EXPECT_EQ(g.get_symmetry_factor(), 4);
      // uint64_t nCn2 = binomial_coefficient(4, 2);
      // EXPECT_EQ(g.get_free_multiplicity(), nCn2 * nCn2);
      found_n_cycle = true;
    }
  }
  EXPECT_TRUE(found_n_cycle);
}

TEST(GenerateDiagramsTest, Order5DiagramsNonBipartite) {

  VacuumDiagramGenerator gen(5, false); // Allow non-bipartite diagrams
  gen.generate();

  const auto &graphs = gen.get_unique_graphs();
  EXPECT_EQ(graphs.size(), 3);

  std::vector<std::vector<uint8_t>> expected_mats = {
     generate_n_cycle_adjacency_matrix(5),
     {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0}, // D5b
     {0, 0, 1, 1, 0, 1, 0, 2, 0}                       // D5c
  };

  for (auto const &m : expected_mats) {
    int V                          = std::sqrt(m.size());
    std::vector<uint8_t> canonical = Graph(m, V, false).get_canonical_form();
    EXPECT_TRUE(std::any_of(graphs.begin(), graphs.end(), [&](const Graph &g) { return g.get_canonical_form() == canonical; }))
       << "Missing order 5 diagram with V=" << V;
  }
}

TEST(GenerateDiagramsTest, Order6DiagramsBipartite) {
  VacuumDiagramGenerator gen(6);
  gen.generate();

  const auto &graphs = gen.get_unique_graphs();
  EXPECT_EQ(graphs.size(), 7);

  auto has_canonical = [&](const std::vector<uint8_t> &canonical) {
    for (const auto &g : graphs) {
      if (g.get_canonical_form() == canonical) return true;
    }
    return false;
  };

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
    EXPECT_TRUE(has_canonical(canonical)) << "Missing order 6 diagram with V=" << V;
  }
}

TEST(GenerateDiagramsTest, Order6DiagramsNoNBipartite) {
  VacuumDiagramGenerator gen(6, false);
  gen.generate();

  const auto &graphs = gen.get_unique_graphs();
  EXPECT_EQ(graphs.size(), 12);
}

TEST(GenerateDiagramsTest, Order8DiagramsBipartite) {
  VacuumDiagramGenerator gen(8);
  gen.generate();

  const auto &graphs = gen.get_unique_graphs();
  EXPECT_EQ(graphs.size(), 32);
}
