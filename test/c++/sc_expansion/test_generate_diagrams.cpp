#include <gtest/gtest.h>
#include "sc_expansion/generate_diagrams.hpp"
#include <vector>
#include <iostream>
#include <cmath>

using namespace sc_expansion;

/*----- N-cycle adjmat------*/
// The raw n-cycle adjacency happened to be the lex-min canonical form under
// the old implementation; with bliss it is a different (equally valid)
// representative. We test invariants here: same number of edges, and any two
// labelings of the same n-cycle map to the same canonical form.
TEST(GenerateDiagramsTest, NCycleAdjacencyMatrix) {
  // n=2 is special: the digon is uniquely labeled.
  std::vector<uint8_t> adjmat2 = generate_n_cycle_adjacency_matrix(2);
  EXPECT_EQ(adjmat2, std::vector<uint8_t>({0, 1, 1, 0}));

  for (int n : {3, 4, 5, 6}) {
    std::vector<uint8_t> adj = generate_n_cycle_adjacency_matrix(n);
    Graph g(adj, n, /*bipartite_only=*/false);
    EXPECT_EQ(g.get_order(), n);
    // A relabeling of the same n-cycle must produce the same canonical form.
    std::vector<uint8_t> shifted(n * n, 0);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) shifted[((i + 1) % n) * n + ((j + 1) % n)] = adj[i * n + j];
    Graph g_shift(shifted, n, /*bipartite_only=*/false);
    EXPECT_EQ(g.get_canonical_form(), g_shift.get_canonical_form()) << "n=" << n;
  }
}

/*----- N-cycle free multiplicity------*/
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

/*----- Order 2 Diagrams ------*/
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

/*----- Order 3 Diagrams ------*/
TEST(GenerateDiagramsTest, Order3DiagramsNonBipartite) {
  VacuumDiagramGenerator gen(3, false); // Allow non-bipartite diagrams
  gen.generate();

  const auto &graphs = gen.get_unique_graphs();
  EXPECT_EQ(graphs.size(), 1);
}

/*----- Order 4 Diagrams ------*/
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
    int V                          = (int)std::sqrt((double)m.size());
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

/*----- Order 5 Diagrams ------*/
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
    int V                          = (int)std::sqrt((double)m.size());
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
    int V                          = (int)std::sqrt((double)m.size());
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

TEST(GenerateDiagramsTest, Order1IsEmpty) {
  // No order-1 vacuum diagram exists in the pure-hopping branch.
  VacuumDiagramGenerator gen(1);
  gen.generate();
  EXPECT_EQ(gen.get_unique_graphs().size(), 0);
}

TEST(GenerateDiagramsTest, Order8DiagramsBipartite) {
  VacuumDiagramGenerator gen(8);
  gen.generate();

  const auto &graphs = gen.get_unique_graphs();
  EXPECT_EQ(graphs.size(), 32);
}

// TEST(GenerateDiagramsTest, Order10DiagramsBipartite) {
//   VacuumDiagramGenerator gen(10);
//   gen.generate();

//   const auto &graphs = gen.get_unique_graphs();
//   std::cout << "Order 10 bipartite diagrams: " << graphs.size() << std::endl;
// }

/*============================================================================*
 * Rooted diagrams (Track 2)
 *============================================================================*/

// Order 2 with one mark: the only vacuum graph is the digon. There is a
// single rooted form (mark on either vertex, related by the vertex-swap
// automorphism). Stabilizer of one marked vertex is trivial → rooted sym
// factor = 1. Anchor the mark at the origin; the unmarked vertex must be
// NN of the origin → 4 valid embeddings, all bucketed at d²=0 (single-mark
// convention).
TEST(RootedDiagramsTest, Order2OneMark) {
  RootedDiagramGenerator gen(2, /*n_marks=*/1);
  gen.generate();
  auto const &rg = gen.get_unique_rooted_graphs();
  EXPECT_EQ(rg.size(), 1u);
  EXPECT_EQ(rg[0].get_n_marks(), 1);
  EXPECT_EQ(rg[0].get_rooted_symmetry_factor(), 1.0);
  auto const &shell = rg[0].get_shell_multiplicity();
  ASSERT_EQ(shell.size(), 1u);
  EXPECT_EQ(shell.at(0), 4);
}

// Order 2 with two marks: the digon has both vertices marked. They are
// interchangeable (same color), so the swap is a colored automorphism →
// rooted sym factor = 2. Anchor mark[0] at the origin; mark[1] must be NN
// (the digon edges enforce adjacency) → 4 valid embeddings, all at d²=1.
TEST(RootedDiagramsTest, Order2TwoMarks) {
  RootedDiagramGenerator gen(2, /*n_marks=*/2);
  gen.generate();
  auto const &rg = gen.get_unique_rooted_graphs();
  EXPECT_EQ(rg.size(), 1u);
  EXPECT_EQ(rg[0].get_n_marks(), 2);
  EXPECT_EQ(rg[0].get_rooted_symmetry_factor(), 2.0);
  auto const &shell = rg[0].get_shell_multiplicity();
  ASSERT_EQ(shell.size(), 1u);
  EXPECT_EQ(shell.at(1), 4);
}

// Order 4 with two marks, restricted to the "fork" D4b (V=3, vertex 0
// is the degree-4 center, vertices 1 and 2 are the degree-2 leaves
// connected to the center by bidirectional edge pairs: 0↔1, 0↔2).
// Aut(G) is generated by the leaf swap → C_2 of order 2. Mark orbits:
//   - {leaf, leaf}: swap fixes the unordered pair → rooted sym = 2.
//     Both leaves share a sublattice with each other (opposite from the
//     center), so accessible distances are d² ∈ {0, 2, 4}.
//     Anchor one leaf at the origin: center has 4 NN positions;
//     for each, the other leaf has 4 NN-of-center positions. 16 total
//     embeddings, distributed as {0: 4, 2: 8, 4: 4}.
//   - {leaf, center}: swap maps the pair to {other leaf, center} ≠ same
//     pair → rooted sym = 1. The two marked vertices are graph-adjacent,
//     hence lattice-NN, so d² = 1 always. The unmarked leaf must be NN
//     of the center → 4 × 4 = 16 embeddings, all at d² = 1.
TEST(RootedDiagramsTest, Order4ForkTwoMarks) {
  RootedDiagramGenerator gen(4, /*n_marks=*/2, /*bipartite_only=*/true);
  gen.generate();
  auto const &all = gen.get_unique_rooted_graphs();

  std::vector<RootedGraph const *> fork_rooted;
  for (auto const &r : all)
    if (r.get_V() == 3) fork_rooted.push_back(&r);

  ASSERT_EQ(fork_rooted.size(), 2u);

  RootedGraph const *two_leaves   = nullptr;
  RootedGraph const *leaf_center = nullptr;
  for (auto *r : fork_rooted) {
    if (r->get_rooted_symmetry_factor() == 2.0) two_leaves = r;
    else if (r->get_rooted_symmetry_factor() == 1.0) leaf_center = r;
  }
  ASSERT_NE(two_leaves, nullptr);
  ASSERT_NE(leaf_center, nullptr);

  auto const &h_two_leaves  = two_leaves->get_shell_multiplicity();
  auto const &h_leaf_center = leaf_center->get_shell_multiplicity();

  EXPECT_EQ(h_two_leaves.size(), 3u);
  EXPECT_EQ(h_two_leaves.at(0), 4);
  EXPECT_EQ(h_two_leaves.at(2), 8);
  EXPECT_EQ(h_two_leaves.at(4), 4);

  EXPECT_EQ(h_leaf_center.size(), 1u);
  EXPECT_EQ(h_leaf_center.at(1), 16);
}

// Order 4 with two marks, restricted to the 4-cycle (V=4). Hopping lines
// are directed in this codebase, so Aut(C_4) = C_4 (rotations only, order
// 4), not D_4. Unordered mark pairs split into two orbits:
//   - adjacent (e.g. {0,1}): |Stab| = 1 (only identity fixes the pair;
//     rotation by 2 maps {0,1} → {2,3})
//   - diagonal (e.g. {0,2}): |Stab| = 2 (identity + 180° rotation which
//     maps 0↔2, 1↔3)
// So we expect 2 unique rooted graphs with rooted symmetry factors {1, 2}
// (no multi-edge corrections — all mults are 1).
TEST(RootedDiagramsTest, Order4FourCycleTwoMarks) {
  RootedDiagramGenerator gen(4, /*n_marks=*/2, /*bipartite_only=*/true);
  gen.generate();
  auto const &all = gen.get_unique_rooted_graphs();

  std::vector<double> four_cycle_factors;
  for (auto const &r : all)
    if (r.get_V() == 4) four_cycle_factors.push_back(r.get_rooted_symmetry_factor());

  EXPECT_EQ(four_cycle_factors.size(), 2u);
  std::sort(four_cycle_factors.begin(), four_cycle_factors.end());
  EXPECT_EQ(four_cycle_factors[0], 1.0);
  EXPECT_EQ(four_cycle_factors[1], 2.0);
}
