#include <gtest/gtest.h>
#include <iostream>
#include "sc_expansion/graph.hpp"

using namespace sc_expansion;

class GraphTest : public ::testing::Test {};

// void print_all_order6_multiplicities() {
//   // Diagrams D6a to D6g from c++/sc_expansion/dimer_order_4.hpp
//   std::vector<uint8_t> D6a = {0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0};
//   std::vector<uint8_t> D6b = {0, 3, 3, 0};
//   std::vector<uint8_t> D6c = {0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
//   std::vector<uint8_t> D6d = {0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0};
//   std::vector<uint8_t> D6e = {0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0};
//   std::vector<uint8_t> D6f = {0, 2, 1, 2, 0, 0, 1, 0, 0};
//   std::vector<uint8_t> D6g = {0, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0};

//   std::cout << "D6a (6-cycle) Free Multiplicity: " << Graph(D6a, 6).get_free_multiplicity() << std::endl;
//   std::cout << "D6b (watermelon triple) Free Multiplicity: " << Graph(D6b, 2).get_free_multiplicity() << std::endl;
//   std::cout << "D6c (petal with 4 vertices) Free Multiplicity: " << Graph(D6c, 4).get_free_multiplicity() << std::endl;
//   std::cout << "D6d (square + digon) Free Multiplicity: " << Graph(D6d, 5).get_free_multiplicity() << std::endl;
//   std::cout << "D6e (crab diagram) Free Multiplicity: " << Graph(D6e, 4).get_free_multiplicity() << std::endl;
//   std::cout << "D6f (watermelon double + digon) Free Multiplicity: " << Graph(D6f, 3).get_free_multiplicity() << std::endl;
//   std::cout << "D6g (square with one double line) Free Multiplicity: " << Graph(D6g, 4).get_free_multiplicity() << std::endl;
// }

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

TEST_F(GraphTest, GraphFreeMultiplicityIsCorrect) {

  // D2a = {{0, 1}, {1, 0}};
  std::vector<uint8_t> D2a = {0, 1, 1, 0};
  Graph graph_2a(D2a, 2);
  EXPECT_EQ((int)graph_2a.get_free_multiplicity(), 4);

  std::vector<uint8_t> D3a = {0, 1, 0, 0, 0, 1, 1, 0, 0};
  Graph graph_3a(D3a, 3, false);
  EXPECT_EQ((int)graph_3a.get_free_multiplicity(), 12);

  std::vector<uint8_t> D4a = {0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0};
  Graph graph_4a(D4a, 4);
  EXPECT_EQ((int)graph_4a.get_free_multiplicity(), 36);

  std::vector<uint8_t> D4b = {0, 1, 1, 1, 0, 0, 1, 0, 0};
  Graph graph_4b(D4b, 3);
  EXPECT_EQ((int)graph_4b.get_free_multiplicity(), 16);

  std::vector<uint8_t> D4c = {0, 2, 2, 0};
  Graph graph_4c(D4c, 2);
  EXPECT_EQ((int)graph_4c.get_free_multiplicity(), 4);

  //D6
  std::vector<uint8_t> D6_1 = {0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0};
  Graph graph_6_1(D6_1, 4);
  EXPECT_EQ((int)graph_6_1.get_free_multiplicity(), 64);

  std::vector<uint8_t> D6_2 = {0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  Graph graph_6_2(D6_2, 4);
  EXPECT_EQ((int)graph_6_2.get_free_multiplicity(), 64);

  std::vector<uint8_t> D6a = {0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0};
  EXPECT_EQ((int)Graph(D6a, 6).get_free_multiplicity(), 400);

  std::vector<uint8_t> D6b = {0, 3, 3, 0};
  EXPECT_EQ((int)Graph(D6b, 2).get_free_multiplicity(), 4);

  std::vector<uint8_t> D6c = {0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  EXPECT_EQ((int)Graph(D6c, 4).get_free_multiplicity(), 64);

  std::vector<uint8_t> D6d = {0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0};
  EXPECT_EQ((int)Graph(D6d, 5).get_free_multiplicity(), 144);

  std::vector<uint8_t> D6e = {0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0};
  EXPECT_EQ((int)Graph(D6e, 4).get_free_multiplicity(), 64);

  std::vector<uint8_t> D6f = {0, 2, 1, 2, 0, 0, 1, 0, 0};
  EXPECT_EQ((int)Graph(D6f, 3).get_free_multiplicity(), 16);

  std::vector<uint8_t> D6g = {0, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0};
  EXPECT_EQ((int)Graph(D6g, 4).get_free_multiplicity(), 36);

  std::vector<uint8_t> D8a(64, 0);
  for (int i = 0; i < 8; ++i) {
    D8a[i * 8 + (i + 1) % 8]   = 1;
    D8a[((i + 1) % 8) * 8 + i] = 1;
  }
  EXPECT_EQ((int)Graph(D8a, 8).get_free_multiplicity(), 4900);

  std::vector<uint8_t> D8b = {0, 4, 4, 0};
  Graph graph_8b(D8b, 2);
  EXPECT_EQ((int)graph_8b.get_free_multiplicity(), 4);
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

//TEST_F(GraphTest, PrintOtherFreeMultiplicities) { print_all_order6_multiplicities(); }
