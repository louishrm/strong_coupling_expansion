#include <gtest/gtest.h>
#include "sc_expansion/graph.hpp"
#include "sc_expansion/generate_diagrams.hpp"

using namespace sc_expansion;

class GraphTest : public ::testing::Test {};

// =====================================================================
// Disconnected-graph handling (order-agnostic sanity checks).
// =====================================================================
TEST_F(GraphTest, DisconnectedGraphs) {
  // No edges → disconnected.
  {
    std::vector<uint8_t> adj(9, 0);
    Graph g(adj, 3);
    EXPECT_FALSE(g.get_connectivity());
  }

  // Two disconnected digons.
  {
    std::vector<uint8_t> adj = {0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0};
    Graph g(adj, 4);
    EXPECT_FALSE(g.get_connectivity());
  }
}

// =====================================================================
// Order 2: pure-hop digon (the only order-2 graph in the no-self-loop branch).
// =====================================================================
TEST_F(GraphTest, Order2) {
  // (a) Pure-hop digon, bipartite_only=true.
  {
    std::vector<uint8_t> adj = {0, 1, 1, 0};
    Graph g(adj, 2);
    EXPECT_EQ(g.get_order(), 2);
    EXPECT_TRUE(g.get_connectivity());
    EXPECT_TRUE(g.get_bipartite());
    EXPECT_EQ(g.get_symmetry_factor(), 2);
  }

  // (b) Same digon built with bipartite_only=false: graph is still bipartite.
  {
    std::vector<uint8_t> adj = {0, 1, 1, 0};
    Graph g(adj, 2, /*bipartite_only=*/false);
    EXPECT_TRUE(g.get_bipartite());
  }
}

// =====================================================================
// Order 3: the 3-cycle (non-bipartite). With self-loops removed there are
// no other connected order-3 graphs at this branch.
// =====================================================================
TEST_F(GraphTest, Order3) {
  std::vector<uint8_t> adj = {0, 1, 0, 0, 0, 1, 1, 0, 0}; // 0->1, 1->2, 2->0
  Graph g(adj, 3, /*bipartite_only=*/false);
  EXPECT_EQ(g.get_order(), 3);
  EXPECT_TRUE(g.get_connectivity());
  EXPECT_FALSE(g.get_bipartite());
}

// =====================================================================
// Order 4: 4-cycle and the V=3 order-4 graph (canonical-form witness).
// =====================================================================
TEST_F(GraphTest, Order4) {
  // 4-cycle: connected, bipartite, symmetry_factor = 4 (cyclic rotations).
  {
    std::vector<uint8_t> adj = {0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0};
    Graph g(adj, 4);
    EXPECT_EQ(g.get_order(), 4);
    EXPECT_TRUE(g.get_connectivity());
    EXPECT_EQ(g.get_symmetry_factor(), 4);
  }

  // Double-edge digon {0,2,2,0}: symmetry_factor = 2*(2!)^2 = 8.
  {
    std::vector<uint8_t> adj = {0, 2, 2, 0};
    Graph g(adj, 2);
    EXPECT_EQ(g.get_order(), 4);
    EXPECT_EQ(g.get_symmetry_factor(), 8);
  }

  // Canonical-form witness on the V=3 order-4 graph.
  {
    std::vector<uint8_t> adj = {0, 1, 1, 1, 0, 0, 1, 0, 0};
    Graph g(adj, 3);
    EXPECT_EQ(g.get_order(), 4);
    EXPECT_EQ(g.get_symmetry_factor(), 2);
    EXPECT_EQ(g.get_canonical_form(), std::vector<uint8_t>({0, 0, 1, 0, 0, 1, 1, 1, 0}));
  }
}

// =====================================================================
// Constructor that accepts precomputed properties (used by the disk
// cache and the MPI broadcast path).
// =====================================================================
TEST_F(GraphTest, ConstructorOverride) {
  std::vector<uint8_t> adj = {0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0}; // D6e "crab"
  int V                    = 4;
  int automorphism_count   = 2;
  int symmetry_factor      = 2;
  int free_multiplicity    = 64;

  Graph g(adj, V, automorphism_count, symmetry_factor, free_multiplicity);

  EXPECT_EQ(g.get_V(), V);
  EXPECT_EQ(g.get_order(), 6);
  EXPECT_EQ(g.get_symmetry_factor(), symmetry_factor);
  EXPECT_EQ(g.get_free_multiplicity(), free_multiplicity);
  EXPECT_TRUE(g.get_connectivity());
  EXPECT_TRUE(g.get_bipartite());
  EXPECT_EQ(g.get_canonical_form(), adj);
}
