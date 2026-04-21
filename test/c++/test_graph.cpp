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
// Order 1: the lone self-insertion graph {1}.
// Only exists when self-loops are allowed; a single density operator
// c^dag_sigma c_sigma at one vertex.
// =====================================================================
TEST_F(GraphTest, Order1) {
  std::vector<uint8_t> adj = {1};
  Graph g(adj, 1);

  EXPECT_EQ(g.get_order(), 1);
  EXPECT_EQ(g.get_n_self_loops(), 1);
  EXPECT_TRUE(g.get_connectivity());
  EXPECT_TRUE(g.get_bipartite()); // self-loop ignored by the bipartite check
  EXPECT_EQ(g.get_symmetry_factor(), 1);
  EXPECT_EQ(g.get_free_multiplicity(), 1);
  EXPECT_EQ(g.get_canonical_form(), adj);
}

// =====================================================================
// Order 2: cover the four flag-combinations (bipartite/non-bipartite,
// self-inserting/non-self-inserting). At order 2 no connected graph is
// non-bipartite — both available graphs are bipartite — so the
// "non-bipartite" slot is vacuous and noted.
// =====================================================================
TEST_F(GraphTest, Order2) {
  // (a) Non-self-inserting, bipartite: the pure-hop digon.
  {
    std::vector<uint8_t> adj = {0, 1, 1, 0};
    Graph g(adj, 2);
    EXPECT_EQ(g.get_order(), 2);
    EXPECT_EQ(g.get_n_self_loops(), 0);
    EXPECT_TRUE(g.get_connectivity());
    EXPECT_TRUE(g.get_bipartite());
    EXPECT_EQ(g.get_symmetry_factor(), 2);
  }

  // (b) Same digon built with bipartite_only=false: the bipartite property
  // depends on the graph, not the flag — it is still bipartite.
  {
    std::vector<uint8_t> adj = {0, 1, 1, 0};
    Graph g(adj, 2, /*bipartite_only=*/false);
    EXPECT_TRUE(g.get_bipartite());
  }

  // (c) Self-inserting, bipartite: two density operators at one vertex.
  //     Symmetry factor = 2! from the two interchangeable self-loops.
  {
    std::vector<uint8_t> adj = {2};
    Graph g(adj, 1);
    EXPECT_EQ(g.get_order(), 2);
    EXPECT_EQ(g.get_n_self_loops(), 2);
    EXPECT_TRUE(g.get_connectivity());
    EXPECT_TRUE(g.get_bipartite());
    EXPECT_EQ(g.get_symmetry_factor(), 2);
    EXPECT_EQ(g.get_free_multiplicity(), 1);
  }
}

// =====================================================================
// Order 3: the 3-cycle (non-bipartite) plus the two self-insertion
// graphs that survive the bipartite filter when self-loops are enabled.
// =====================================================================
TEST_F(GraphTest, Order3) {
  // (a) 3-cycle (triangle) — connected, non-bipartite.
  {
    std::vector<uint8_t> adj = {0, 1, 0, 0, 0, 1, 1, 0, 0}; // 0->1, 1->2, 2->0
    Graph g(adj, 3, /*bipartite_only=*/false);
    EXPECT_EQ(g.get_order(), 3);
    EXPECT_EQ(g.get_n_self_loops(), 0);
    EXPECT_TRUE(g.get_connectivity());
    EXPECT_FALSE(g.get_bipartite());
  }

  // (b) With self-insertions enabled and bipartite_only=true, the order-3
  // generator must return exactly two graphs:
  //     {3}       — three density insertions at one vertex;
  //     {1,1,1,0} — one density insertion plus a digon (V=2).
  // The 3-cycle is discarded because it is non-bipartite.
  {
    VacuumDiagramGenerator gen(3, /*bipartite_only=*/true, /*allow_self_loops=*/true);
    gen.generate();
    auto const &graphs = gen.get_unique_graphs();
    ASSERT_EQ(graphs.size(), 2);

    auto canonical_of         = [](std::vector<uint8_t> a, int V) { return Graph(a, V).get_canonical_form(); };
    auto target_triple_self   = canonical_of({3}, 1);
    auto target_hop_plus_self = canonical_of({1, 1, 1, 0}, 2);

    bool found_triple_self   = false;
    bool found_hop_plus_self = false;
    for (auto const &g : graphs) {
      if (g.get_canonical_form() == target_triple_self) {
        found_triple_self = true;
        EXPECT_EQ(g.get_n_self_loops(), 3);
        EXPECT_EQ(g.get_V(), 1);
      }
      if (g.get_canonical_form() == target_hop_plus_self) {
        found_hop_plus_self = true;
        EXPECT_EQ(g.get_n_self_loops(), 1);
        EXPECT_EQ(g.get_V(), 2);
      }
    }
    EXPECT_TRUE(found_triple_self);
    EXPECT_TRUE(found_hop_plus_self);
  }
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
