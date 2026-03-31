#include <gtest/gtest.h>
#include <iostream>
#include "sc_expansion/graph.hpp"

using namespace sc_expansion;

class GraphTest : public ::testing::Test {};

TEST_F(GraphTest, GraphIsConnected) {

  std::vector<uint8_t> adjacency_matrix = {0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0}; //4 vertices, edges: 0->1, 1->2, 2->3, 3->0
  Graph graph(adjacency_matrix, 4);
  EXPECT_TRUE(graph.get_connectivity());

  std::vector<uint8_t> adjacency_matrix2 = {0, 2, 1, 2, 0, 0, 1, 0, 1};
  Graph graph2(adjacency_matrix2, 3);
  EXPECT_TRUE(graph2.get_connectivity()); //3 vertices, 0->1 (twice) 0->2

  std::vector<uint8_t> adjacency_matrix3 = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  Graph graph3(adjacency_matrix3, 3);
  EXPECT_FALSE(graph3.get_connectivity()); //3 vertices, no edges

  std::vector<uint8_t> adjacency_matrix4 = {0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0}; //Two disconnected digons
  Graph graph4(adjacency_matrix4, 4);
  EXPECT_FALSE(graph4.get_connectivity());
}

TEST_F(GraphTest, GraphIsBipartite) {
  // 4-cycle is bipartite
  std::vector<uint8_t> adjacency_matrix = {0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0};
  Graph graph(adjacency_matrix, 4);
  EXPECT_TRUE(graph.get_bipartite());

  // 3-cycle is not bipartite
  std::vector<uint8_t> adjacency_matrix2 = {0, 1, 1, 1, 0, 1, 1, 1, 0};
  Graph graph2(adjacency_matrix2, 3);
  EXPECT_FALSE(graph2.get_bipartite());

  // Digon is bipartite
  std::vector<uint8_t> adjacency_matrix3 = {0, 2, 2, 0};
  Graph graph3(adjacency_matrix3, 2);
  EXPECT_TRUE(graph3.get_bipartite());
}

TEST_F(GraphTest, GraphSymmetryFactorIsCorrect) {

  std::vector<uint8_t> adjacency_matrix = {0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0}; //4 vertices, edges: 0->1, 1->2, 2->3, 3->0
  Graph graph(adjacency_matrix, 4);
  EXPECT_EQ(graph.get_symmetry_factor(), 4);

  std::vector<uint8_t> adjacency_matrix2 = {0, 1, 1, 1, 0, 0, 1, 0, 0};
  Graph graph2(adjacency_matrix2, 3);
  EXPECT_EQ(graph2.get_symmetry_factor(), 2); //3 vertices, 0->1, 0->2, 1->0

  std::vector<uint8_t> adjacency_matrix3 = {0, 2, 2, 0};
  Graph graph3(adjacency_matrix3, 2);
  EXPECT_EQ(graph3.get_symmetry_factor(), 8);
}

TEST_F(GraphTest, GraphCanonicalFormIsCorrect) {

  std::vector<uint8_t> adjacency_matrix = {0, 1, 1, 1, 0, 0, 1, 0, 0}; //3 vertices, 0->1, 0->2, 1->0
  Graph graph(adjacency_matrix, 3);
  EXPECT_EQ(graph.get_canonical_form(), std::vector<uint8_t>({0, 0, 1, 0, 0, 1, 1, 1, 0}));
}

TEST_F(GraphTest, GraphConstructorOverride) {
  std::vector<uint8_t> adjacency_matrix = {0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0}; // D6e crab diagram
  int V                                 = 4;
  int automorphism_count                = 2;
  int symmetry_factor                   = 2;
  int free_multiplicity                 = 64;

  Graph graph(adjacency_matrix, V, automorphism_count, symmetry_factor, free_multiplicity);

  EXPECT_EQ(graph.get_V(), V);
  EXPECT_EQ(graph.get_order(), 6);
  EXPECT_EQ(graph.get_symmetry_factor(), symmetry_factor);
  EXPECT_EQ(graph.get_free_multiplicity(), free_multiplicity);
  EXPECT_TRUE(graph.get_connectivity());
  EXPECT_TRUE(graph.get_bipartite());
  EXPECT_EQ(graph.get_canonical_form(), adjacency_matrix);
}
