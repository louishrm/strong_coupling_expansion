#include <gtest/gtest.h>
#include <limits>
#include <memory>
#include <random>
#include <vector>

#include "sc_expansion/atomic/diagram.hpp"
#include "sc_expansion/atomic/vertex.hpp"
#include "sc_expansion/graph.hpp"
#include "sc_expansion/hubbard_solver.hpp"

using namespace sc_expansion;

namespace {
  constexpr double kTol = 1e-12;

  // 4-cycle as 4 directed hopping lines: 0 -> 1 -> 2 -> 3 -> 0 (order=4, V=4).
  std::vector<uint8_t> four_cycle_adjacency() { return {0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0}; }

  template <typename T> std::vector<std::unique_ptr<atomic::VertexType<T>>> make_vertex_types(int max_co) {
    std::vector<std::unique_ptr<atomic::VertexType<T>>> out;
    for (int k = 1; k <= max_co; ++k) out.push_back(std::make_unique<atomic::VertexType<T>>(2 * k));
    return out;
  }
  template <typename T> std::vector<atomic::VertexType<T> *> to_ptrs(std::vector<std::unique_ptr<atomic::VertexType<T>>> &owned) {
    std::vector<atomic::VertexType<T> *> ptrs;
    for (auto &up : owned) ptrs.push_back(up.get());
    return ptrs;
  }
} // namespace

TEST(StaticDensityFreeMultiplicity, FourCycle) {
  Graph cycle(four_cycle_adjacency(), /*V=*/4, /*bipartite_only=*/true);

  // shell_multiplicity is sparse: d^2 values with no embeddings are absent
  // from the map (.at() would throw). Check absence with .count() == 0.

  // (a) Both marks on the same vertex: only d^2=0 is reachable.
  RootedGraph same_vertex(cycle, /*marks=*/{0, 0});
  auto const &shell_a = same_vertex.get_shell_multiplicity();
  ASSERT_EQ(shell_a.count(0), 1u);
  EXPECT_EQ(shell_a.at(0), 36);
  EXPECT_EQ(shell_a.count(1), 0u);
  EXPECT_EQ(shell_a.count(2), 0u);

  // (b) Marks on opposite vertices (graph-distance 2, bipartite -> even d^2).
  //     d^2=0: 16  (v0=v2=origin; v1,v3 each NN of origin)
  //     d^2=2: 16  (v2 in {(+/-1,+/-1)}: 4 positions, 2*2 placements each)
  //     d^2=4: 4   (v2 in {(+/-2,0),(0,+/-2)}: 4 positions, 1*1 each)
  RootedGraph opposite_pair(cycle, /*marks=*/{0, 2});
  auto const &shell_b = opposite_pair.get_shell_multiplicity();
  ASSERT_EQ(shell_b.count(0), 1u);
  EXPECT_EQ(shell_b.at(0), 16);
  EXPECT_EQ(shell_b.count(1), 0u);
  ASSERT_EQ(shell_b.count(2), 1u);
  EXPECT_EQ(shell_b.at(2), 16);
  ASSERT_EQ(shell_b.count(4), 1u);
  EXPECT_EQ(shell_b.at(4), 4);

  // (c) Marks on adjacent vertices: only d^2=1 is reachable (v1 must be a
  //     NN of v0=origin), no embeddings at d^2=0.
  RootedGraph adjacent_pair(cycle, /*marks=*/{0, 1});
  auto const &shell_c = adjacent_pair.get_shell_multiplicity();
  EXPECT_EQ(shell_c.count(0), 0u);
  ASSERT_EQ(shell_c.count(1), 1u);
  EXPECT_EQ(shell_c.at(1), 36);
  EXPECT_EQ(shell_c.count(2), 0u);
}

// ----------------------------------------------------------------------------
// Order-4 V=3 path graph: undirected 0-1-2 with every edge doubled (one
// directed line each way), so the adjacency matrix is the symmetric
// {{0,1,0},{1,0,1},{0,1,0}} (4 directed lines total). Middle vertex has
// degree 4, endpoints have degree 2. Vacuum free multiplicity anchored at
// v0 = 16 (v1 has 4 NN choices, then v2 has 4).
//
// Mark coverage at squared lattice distances:
//   {0,0}: both insertions on endpoint v0 -> d^2=0 always, mult 16.
//   {1,1}: both insertions on center v1 -> d^2=0 always, mult 16.
//   {0,1}: endpoint + center, graph-distance 1 -> only d^2=1, mult 16.
//   {0,2}: the two endpoints (graph-distance 2) split by d^2 of v2:
//     d^2=0: 4   (v2 returns to origin)
//     d^2=2: 8   (v2 on the diagonal (+/-1, +/-1))
//     d^2=4: 4   (v2 on the axes (+/-2, 0) or (0, +/-2))
//   |r|_2 = sqrt(2) and |r|_2 = 2 sit in distinct d^2 buckets because the
//   shell_multiplicity histogram keys on Euclidean d^2 = x^2 + y^2.
// ----------------------------------------------------------------------------
TEST(StaticDensityFreeMultiplicity, V3Path) {
  std::vector<uint8_t> adj = {0, 1, 0,
                              1, 0, 1,
                              0, 1, 0};
  Graph path(adj, /*V=*/3, /*bipartite_only=*/true);
  ASSERT_EQ(path.get_free_multiplicity(), 16);

  RootedGraph end_double(path, /*marks=*/{0, 0});
  auto const &shell_end = end_double.get_shell_multiplicity();
  ASSERT_EQ(shell_end.count(0), 1u);
  EXPECT_EQ(shell_end.at(0), 16);
  EXPECT_EQ(shell_end.count(1), 0u);
  EXPECT_EQ(shell_end.count(2), 0u);

  RootedGraph center_double(path, /*marks=*/{1, 1});
  auto const &shell_center = center_double.get_shell_multiplicity();
  ASSERT_EQ(shell_center.count(0), 1u);
  EXPECT_EQ(shell_center.at(0), 16);
  EXPECT_EQ(shell_center.count(1), 0u);

  RootedGraph end_center(path, /*marks=*/{0, 1});
  auto const &shell_ec = end_center.get_shell_multiplicity();
  EXPECT_EQ(shell_ec.count(0), 0u);
  ASSERT_EQ(shell_ec.count(1), 1u);
  EXPECT_EQ(shell_ec.at(1), 16);
  EXPECT_EQ(shell_ec.count(2), 0u);

  RootedGraph endpoints(path, /*marks=*/{0, 2});
  auto const &shell_ee = endpoints.get_shell_multiplicity();
  ASSERT_EQ(shell_ee.count(0), 1u);
  EXPECT_EQ(shell_ee.at(0), 4);
  EXPECT_EQ(shell_ee.count(1), 0u);
  ASSERT_EQ(shell_ee.count(2), 1u);
  EXPECT_EQ(shell_ee.at(2), 8);
  ASSERT_EQ(shell_ee.count(4), 1u);
  EXPECT_EQ(shell_ee.at(4), 4);
}

// D4 graph automorphism of the 4-cycle relates the four edge mark pairs
// {0,1}, {1,2}, {2,3}, {0,3} — the four "1st-NN" choices. atomic::Diagram
// indexes taus by hopping-line index (row-major sweep of the adjacency
// matrix: line 0 = v0->v1, line 1 = v1->v2, line 2 = v2->v3, line 3 =
// v3->v0). The rotation marks_k = {k, k+1 mod 4} cycles the line indices
// by k as well, so a genuine D4-invariance check has to apply the same
// cyclic shift to taus. With the matched shift, all four evaluations
// hit the same multiset of (cumulant order, time pair) per vertex and
// must produce exactly the same number.
TEST(StaticDensityDiagram, FourCycleDistanceSymmetric) {
  const double U = 4.0, beta = 1.0, mu = 1.2;
  Parameters<double> params{U, beta, mu, 0.0, true};
  HubbardSolver<1, double> solver(params);

  auto owned = make_vertex_types<double>(/*max_co=*/2);
  auto vts   = to_ptrs(owned);

  // Rooted-graph factor convention: free_multiplicity = 1 (lattice
  // multiplier lives in the measurement layer) and symmetry_factor =
  // rooted automorphism count = 2 (reflection through the marked edge).
  auto make_rooted_cycle = []() {
    return Graph(four_cycle_adjacency(), /*V=*/4, /*automorphism_count=*/2,
                 /*symmetry_factor=*/2, /*free_multiplicity=*/1, /*bipartite_only=*/true);
  };

  // For marks {k, k+1 mod 4}, the rotated line at index j corresponds
  // physically to the original line at index (j - k) mod 4, so the matching
  // tau shift is -k (mod 4), not +k.
  struct Case {
    std::pair<int, int> r;
    std::vector<int> marks;
    int k; // mark rotation index
  };
  const std::vector<Case> cases = {
     {{1, 0},  {0, 1}, 0},
     {{0, 1},  {1, 2}, 1},
     {{-1, 0}, {2, 3}, 2},
     {{0, -1}, {0, 3}, 3},
  };
  const std::vector<int> mark_spins  = {1, 1};
  const std::vector<double> taus_ref = {0.2, 0.114, 0.007, 0.78};
  const int n_lines                  = 4;

  double reference = std::numeric_limits<double>::quiet_NaN();
  for (auto const &c : cases) {
    std::vector<double> taus(n_lines + 1);
    for (int i = 0; i < n_lines; ++i) taus[i] = taus_ref[(i - c.k + n_lines) % n_lines];
    taus[n_lines] = 0.0;

    Graph g = make_rooted_cycle();
    atomic::Diagram<double> diag(g, vts, c.marks, mark_spins);
    diag.mark_all_dirty();
    double val = diag.evaluate(taus, solver, /*infinite_U=*/false);
    if (std::isnan(reference)) reference = val;
    else EXPECT_DOUBLE_EQ(val, reference)
        << "r=(" << c.r.first << "," << c.r.second << ") marks={" << c.marks[0] << "," << c.marks[1] << "}";
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
