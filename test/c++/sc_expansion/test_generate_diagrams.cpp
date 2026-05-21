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

// Order 2 with two marks on the digon. Two unique rooted forms:
//   - {0,1} distinct vertices: both vertices colored 1, swap is a colored
//     automorphism → rooted sym = 2; mark[1] is graph-NN of mark[0] so all
//     4 NN embeddings live at d²=1.
//   - {0,0}~{1,1} same vertex: colors {2,0} and {0,2} fold under the
//     graph swap → 1 unique form. Vertex-swap stabilizer of a single
//     colored-2 vertex is trivial → rooted sym = 1; mark[1] coincides
//     with mark[0] at the origin, so shell histogram has its mass at d²=0
//     equal to the single-mark embedding count (4, the unmarked vertex's
//     NN positions).
TEST(RootedDiagramsTest, Order2TwoMarks) {
  RootedDiagramGenerator gen(2, /*n_marks=*/2);
  gen.generate();
  auto const &rg = gen.get_unique_rooted_graphs();
  ASSERT_EQ(rg.size(), 2u);

  RootedGraph const *distinct   = nullptr;
  RootedGraph const *coincident = nullptr;
  for (auto const &r : rg) {
    auto const &m = r.get_marks();
    if (m[0] == m[1])
      coincident = &r;
    else
      distinct = &r;
  }
  ASSERT_NE(distinct, nullptr);
  ASSERT_NE(coincident, nullptr);

  EXPECT_EQ(distinct->get_rooted_symmetry_factor(), 2.0);
  auto const &h_distinct = distinct->get_shell_multiplicity();
  ASSERT_EQ(h_distinct.size(), 1u);
  EXPECT_EQ(h_distinct.at(1), 4);

  EXPECT_EQ(coincident->get_rooted_symmetry_factor(), 1.0);
  auto const &h_coincident = coincident->get_shell_multiplicity();
  ASSERT_EQ(h_coincident.size(), 1u);
  EXPECT_EQ(h_coincident.at(0), 4);
}

// Order 4 with two marks, restricted to the "fork" D4b (V=3, vertex 0
// is the degree-4 center, vertices 1 and 2 are the degree-2 leaves
// connected to the center by bidirectional edge pairs: 0↔1, 0↔2).
// Aut(G) is generated by the leaf swap → C_2 of order 2. Four orbits:
//   - {leaf, leaf} (distinct): swap fixes the unordered pair → rooted
//     sym = 2. d² ∈ {0, 2, 4}: anchor one leaf at origin; center has 4
//     NN positions; for each, the other leaf has 4 NN-of-center positions.
//     16 total embeddings, {0: 4, 2: 8, 4: 4}.
//   - {leaf, center} (distinct): swap maps the pair to {other leaf,
//     center} ≠ same pair → rooted sym = 1. Marks are graph-adjacent,
//     hence lattice-NN, so d² = 1 always. Unmarked leaf must be NN of
//     center → 4 × 4 = 16 embeddings, all at d² = 1.
//   - {leaf, leaf} coincident (e.g. {1,1} ~ {2,2}): leaf swap folds the
//     two same-vertex placements together → rooted sym = 1 (stabilizer
//     of a single colored-2 leaf is trivial). Mark anchored at origin;
//     center NN (4 choices); other leaf NN-of-center (4 choices) → 16
//     embeddings, all at d² = 0.
//   - {center, center} coincident: leaf swap fixes the center → rooted
//     sym = 2. Anchor at origin; each leaf NN of center → 4 × 4 = 16
//     embeddings, all at d² = 0.
TEST(RootedDiagramsTest, Order4ForkTwoMarks) {
  RootedDiagramGenerator gen(4, /*n_marks=*/2, /*bipartite_only=*/true);
  gen.generate();
  auto const &all = gen.get_unique_rooted_graphs();

  std::vector<RootedGraph const *> fork_rooted;
  for (auto const &r : all)
    if (r.get_V() == 3) fork_rooted.push_back(&r);

  ASSERT_EQ(fork_rooted.size(), 4u);

  RootedGraph const *two_leaves   = nullptr;
  RootedGraph const *leaf_center  = nullptr;
  RootedGraph const *leaf_coinc   = nullptr;
  RootedGraph const *center_coinc = nullptr;
  for (auto *r : fork_rooted) {
    auto const &m   = r->get_marks();
    bool coincident = (m[0] == m[1]);
    double sym      = r->get_rooted_symmetry_factor();
    if (!coincident && sym == 2.0)
      two_leaves = r;
    else if (!coincident && sym == 1.0)
      leaf_center = r;
    else if (coincident && sym == 1.0)
      leaf_coinc = r;
    else if (coincident && sym == 2.0)
      center_coinc = r;
  }
  ASSERT_NE(two_leaves, nullptr);
  ASSERT_NE(leaf_center, nullptr);
  ASSERT_NE(leaf_coinc, nullptr);
  ASSERT_NE(center_coinc, nullptr);

  auto const &h_two_leaves = two_leaves->get_shell_multiplicity();
  EXPECT_EQ(h_two_leaves.size(), 3u);
  EXPECT_EQ(h_two_leaves.at(0), 4);
  EXPECT_EQ(h_two_leaves.at(2), 8);
  EXPECT_EQ(h_two_leaves.at(4), 4);

  auto const &h_leaf_center = leaf_center->get_shell_multiplicity();
  EXPECT_EQ(h_leaf_center.size(), 1u);
  EXPECT_EQ(h_leaf_center.at(1), 16);

  auto const &h_leaf_coinc = leaf_coinc->get_shell_multiplicity();
  EXPECT_EQ(h_leaf_coinc.size(), 1u);
  EXPECT_EQ(h_leaf_coinc.at(0), 16);

  auto const &h_center_coinc = center_coinc->get_shell_multiplicity();
  EXPECT_EQ(h_center_coinc.size(), 1u);
  EXPECT_EQ(h_center_coinc.at(0), 16);
}

// Order 4 with two marks, restricted to the 4-cycle (V=4). Hopping lines
// are directed in this codebase, so Aut(C_4) = C_4 (rotations only, order
// 4), not D_4. Three unordered mark-multiset orbits:
//   - adjacent (e.g. {0,1}): |Stab| = 1 (only identity fixes the pair;
//     rotation by 2 maps {0,1} → {2,3}) → rooted sym = 1
//   - diagonal (e.g. {0,2}): |Stab| = 2 (identity + 180° rotation which
//     maps 0↔2, 1↔3) → rooted sym = 2
//   - coincident ({v,v}): all four placements equivalent under C_4 →
//     |Stab| of one colored-2 vertex is trivial → rooted sym = 1
// So we expect 3 unique rooted graphs with rooted symmetry factors
// {1, 1, 2} (no multi-edge corrections — all mults are 1).
TEST(RootedDiagramsTest, Order4FourCycleTwoMarks) {
  RootedDiagramGenerator gen(4, /*n_marks=*/2, /*bipartite_only=*/true);
  gen.generate();
  auto const &all = gen.get_unique_rooted_graphs();

  std::vector<double> four_cycle_factors;
  for (auto const &r : all)
    if (r.get_V() == 4) four_cycle_factors.push_back(r.get_rooted_symmetry_factor());

  ASSERT_EQ(four_cycle_factors.size(), 3u);
  std::sort(four_cycle_factors.begin(), four_cycle_factors.end());
  EXPECT_EQ(four_cycle_factors[0], 1.0);
  EXPECT_EQ(four_cycle_factors[1], 1.0);
  EXPECT_EQ(four_cycle_factors[2], 2.0);
}

// Order 4 with two marks, all bipartite topologies combined. Sum across the
// three vacuum diagrams: 4-cycle (V=4) contributes 3 (adjacent, diagonal,
// coincident), fork (V=3) contributes 4 (leaf-leaf, leaf-center, leaf-
// coincident, center-coincident), digon-pair D4c (V=2) contributes 2
// (distinct {0,1} and coincident {v,v}). Total = 9.
TEST(RootedDiagramsTest, Order4TwoMarksTotalCount) {
  RootedDiagramGenerator gen(4, /*n_marks=*/2, /*bipartite_only=*/true);
  gen.generate();
  EXPECT_EQ(gen.get_unique_rooted_graphs().size(), 9u);
}

/*============================================================================*
 * Distance-rooted diagrams
 *============================================================================*/

// r = (2, 0) on the square lattice: |r|_1 = d = 2. At expansion order n = 4
// the contributors live in V in [d+1, n] = [3, 4]:
//   - V = 4: the 4-cycle. Diagonal pairs {0,2} and {1,3} both at graph
//     distance 2; rotation in Aut(C_4) = C_4 collapses them to one rooted
//     topology. Adjacent pairs (dist 1) fall below d and are filtered out;
//     coincident pairs (dist 0) likewise.
//   - V = 3: the fork D4b (path 1—0—2 with bidirectional doublings on
//     both endpoint edges). The leaf pair {1, 2} sits at graph distance 2
//     through the center; the fork embeds on the lattice as leaf at
//     origin, center at (1, 0), other leaf at (2, 0) — Manhattan distance
//     2. The leaf-swap automorphism fixes the pair setwise, so this is
//     one rooted topology.
//   - V = 2 (digon-pair) has diameter 1 < d and is pruned by the diameter
//     step.
// Total: 2 unique rooted topologies.
TEST(DistanceRootedDiagramsTest, R20Order4HasTwoTopologies) {
  DistanceRootedDiagramGenerator gen({2, 0}, /*n=*/4);
  gen.generate();
  auto const &rg = gen.get_rooted_graphs();

  int total = 0;
  for (auto const &[V, rgs] : rg) total += rgs.size();
  EXPECT_EQ(total, 2);
}

// Reporting test (not a strict assertion): print survivor counts for the
// distance-targeted catalog at the chosen (r, n), alongside the full vacuum
// catalog size at the same order. Also report cumulant-degree statistics —
// the max vertex degree in a rooted topology is the largest cumulant order
// the MCMC must recurse into for that topology, so lower is better.
TEST(DistanceRootedDiagramsTest, R30Order10ReportCounts) {
  std::vector<int> r = {4, 0};
  int n              = 8;
  int V_min_report   = 3; // skip V < V_min_report from the breakdown / stats

  DistanceRootedDiagramGenerator dgen(r, n, /*bipartite_only=*/true);
  dgen.generate();
  auto const &rg = dgen.get_rooted_graphs();

  // Cumulant order of vertex i = (sum of row i) + (sum of col i) in the
  // canonical adjacency matrix = in_deg + out_deg = number of hopping-line
  // legs attached to that cumulant vertex.
  auto vertex_cumulant_orders = [](RootedGraph const &g) {
    int V                          = g.get_V();
    std::vector<uint8_t> const adj = g.get_canonical_form();
    std::vector<int> orders(static_cast<size_t>(V), 0);
    for (int i = 0; i < V; ++i) {
      for (int j = 0; j < V; ++j) {
        orders[static_cast<size_t>(i)] += adj[static_cast<size_t>(i) * V + j]; // out-edges of i
        orders[static_cast<size_t>(i)] += adj[static_cast<size_t>(j) * V + i]; // in-edges of i
      }
    }
    return orders;
  };

  int total_rooted           = 0;
  int global_max_cum         = 0;
  long long sum_vertex_cum   = 0;
  long long total_vertices   = 0;
  long long sum_topology_max = 0;

  std::cout << "\n[report] DistanceRootedDiagramGenerator r=(" << r[0] << "," << r[1] << "), n=" << n << ", bipartite_only=true, V>=" << V_min_report
            << "\n";
  for (auto const &[V, rgs] : rg) {
    if (V < V_min_report) continue;
    int V_max_here           = 0;
    long long V_sum_topo_max = 0;
    long long V_sum_vert     = 0;
    for (auto const &g : rgs) {
      auto orders          = vertex_cumulant_orders(g);
      int per_topology_max = *std::max_element(orders.begin(), orders.end());
      V_max_here           = std::max(V_max_here, per_topology_max);
      V_sum_topo_max += per_topology_max;
      for (int o : orders) V_sum_vert += o;

      global_max_cum = std::max(global_max_cum, per_topology_max);
      sum_topology_max += per_topology_max;
      for (int o : orders) sum_vertex_cum += o;
      total_vertices += V;
    }
    std::cout << "  V = " << V << ": " << rgs.size() << " topologies, "
              << "max-cum (this V) = " << V_max_here << ", "
              << "avg per-topo-max (this V) = " << (rgs.empty() ? 0.0 : (double)V_sum_topo_max / rgs.size()) << ", "
              << "avg vertex cum-order (this V) = " << (rgs.empty() ? 0.0 : (double)V_sum_vert / (V * (long long)rgs.size())) << "\n";
    total_rooted += static_cast<int>(rgs.size());
  }
  std::cout << "  TOTAL: " << total_rooted << " rooted topologies\n";
  std::cout << "  MAX cumulant order anywhere in catalog: " << global_max_cum << "\n";
  std::cout << "  AVG per-topology max cumulant order:    " << (total_rooted == 0 ? 0.0 : (double)sum_topology_max / total_rooted) << "\n";
  std::cout << "  AVG vertex cumulant order:              " << (total_vertices == 0 ? 0.0 : (double)sum_vertex_cum / total_vertices) << "\n";

  VacuumDiagramGenerator vgen(n, /*bipartite_only=*/true);
  vgen.generate();
  auto const &vacuum_graphs = vgen.get_unique_graphs();
  int vacuum_total          = static_cast<int>(vacuum_graphs.size());
  int vacuum_filtered       = 0;
  for (auto const &g : vacuum_graphs)
    if (g.get_V() >= V_min_report) ++vacuum_filtered;
  std::cout << "[report] VacuumDiagramGenerator order=" << n << ", bipartite_only=true: " << vacuum_total << " unique vacuum graphs ("
            << vacuum_filtered << " with V>=" << V_min_report << ")\n"
            << std::endl;

  EXPECT_TRUE(true);
}
