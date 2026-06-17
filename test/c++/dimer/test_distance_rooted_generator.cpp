#include <gtest/gtest.h>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>
#include "../c++/sc_expansion/generate_diagrams.hpp"

using namespace sc_expansion;

// =============================================================================
//  DimerDistanceRootedDiagramGenerator — minimal smoke test at order 2.
//
//  At r = (3, 0) and n = 2, the catalog should contain exactly one rooted
//  topology: the directed digon (a₀₁ = a₁₀ = 1) with marks at the two
//  vertices and within-dimer sites (0, 1). Rationale (per
//  project_dimer_rooted_constraints memory):
//
//    - r-parity (rx + ry) = 1 (odd) ⇒ allowed sectors are (0,1) and (1,0).
//    - Sector (0,1): Δû = 0 - 1 + 3 = 2, Δv = 0 ⇒ d_super = 1.
//      V_min = 2, n_min = 2 → admits V=2 graphs.
//    - Sector (1,0): Δû = 1 - 0 + 3 = 4, Δv = 0 ⇒ d_super = 2.
//      V_min = 3, n_min = 4 → no graph at V=2 (and n=2 is below n_min for
//      this sector anyway).
//    - VacuumDiagramGenerator at order=2, V=2 emits the n-cycle (2-cycle =
//      digon) and nothing else; partition (1,1) is filtered out as the
//      all-ones case.
//    - Marks: only (i,j) = (0,1) survives the dist_G(i,j) ≥ d_super = 1
//      filter (i = j has dist 0).
//
//  So the V=2 bucket has one entry; the V=3 bucket is empty (n=2 caps V
//  at 2). Canonical adj after the colored canonicalisation (color 1 at
//  the site-0 mark vertex, color 2 at the site-1 mark vertex) is the same
//  {0,1,1,0}: bliss puts the color-1 vertex at canonical position 0, and
//  the digon's adjacency is unchanged under this labelling.
// =============================================================================

TEST(DimerDistanceRootedGenerator, Order2Topologies) {
  std::vector<uint8_t> expected_adj = {0, 1, 1, 0};
  DimerDistanceRootedDiagramGenerator gen({3, 0}, /*n=*/2);
  gen.generate();
  auto const &by_V = gen.get_rooted_graphs();

  // V=2 is the only populated bucket; nothing else fits at n=2.
  ASSERT_EQ(by_V.size(), 1u);
  ASSERT_TRUE(by_V.count(2) == 1);

  auto const &topologies = by_V.at(2);
  ASSERT_EQ(topologies.size(), 1u);

  auto const &rg = topologies[0];

  EXPECT_EQ(rg.get_V(), 2);
  EXPECT_EQ(rg.get_order(), 2);
  EXPECT_EQ(rg.get_canonical_form(), expected_adj);

  // Marks at both vertices, in canonical labelling.
  ASSERT_EQ(rg.get_marks().size(), 2u);
  EXPECT_EQ(rg.get_marks()[0], 0);
  EXPECT_EQ(rg.get_marks()[1], 1);

  // Sector (0,1): canonical mark 0 has site 0, canonical mark 1 has site 1.
  ASSERT_EQ(rg.get_sites().size(), 2u);
  EXPECT_EQ(rg.get_sites()[0], 0);
  EXPECT_EQ(rg.get_sites()[1], 1);
}

// d_super formula spot-checks for the (3,0) case used above and a couple of
// neighbouring points, mirroring the table in the dimer-rooted-constraints
// memory.
TEST(DimerDistanceRootedGenerator, DSuperFormulaSpotChecks) {
  // r = (3, 0): odd parity, sectors (0,1) and (1,0).
  EXPECT_EQ(DimerDistanceRootedDiagramGenerator::d_super(3, 0, 0, 1), 1);
  EXPECT_EQ(DimerDistanceRootedDiagramGenerator::d_super(3, 0, 1, 0), 2);

  // r = (2, 0): even parity, sectors (0,0) and (1,1); both d_super = 1.
  EXPECT_EQ(DimerDistanceRootedDiagramGenerator::d_super(2, 0, 0, 0), 1);
  EXPECT_EQ(DimerDistanceRootedDiagramGenerator::d_super(2, 0, 1, 1), 1);

  // r = (0, 2): even parity, both sectors d_super = 2 (no advantage along ŷ).
  EXPECT_EQ(DimerDistanceRootedDiagramGenerator::d_super(0, 2, 0, 0), 2);
  EXPECT_EQ(DimerDistanceRootedDiagramGenerator::d_super(0, 2, 1, 1), 2);

  // r = (1, 0): sector (0,1) is the intra-dimer correlator with d_super = 0.
  EXPECT_EQ(DimerDistanceRootedDiagramGenerator::d_super(1, 0, 0, 1), 0);
  EXPECT_EQ(DimerDistanceRootedDiagramGenerator::d_super(1, 0, 1, 0), 1);
}

namespace {

  // The adjacency stores *directed* line counts: order = sum of all entries,
  // and an undirected "double" bond is a 2-cycle of two unit directed lines.
  // The structural predicates below identify a topology independently of the
  // (bliss-internal) canonical vertex labeling.

  // Single directed V-cycle: every row and column sums to exactly 1, and
  // following the permutation from vertex 0 returns to 0 only after V steps.
  // This is the n-cycle topology (the only V == order vacuum graph).
  bool is_single_directed_cycle(std::vector<uint8_t> const &m, int V) {
    if (static_cast<int>(m.size()) != V * V) return false;
    std::vector<int> next(V, -1);
    for (int i = 0; i < V; ++i) {
      int row_sum = 0, col_sum = 0;
      for (int j = 0; j < V; ++j) {
        if (m[i * V + j] != 0 && m[i * V + j] != 1) return false;
        row_sum += m[i * V + j];
        col_sum += m[j * V + i];
        if (m[i * V + j] == 1) next[i] = j;
      }
      if (row_sum != 1 || col_sum != 1) return false;
    }
    // Walk the cycle: must visit all V vertices before closing.
    int cur = 0;
    for (int steps = 0; steps < V; ++steps) {
      cur = next[cur];
      if (cur == 0) return steps == V - 1; // closed exactly at length V
    }
    return false;
  }

  // V == 3 double-digon ("banana"): two undirected double-bonds sharing the
  // center vertex, i.e. the symmetric path P3. Degree sequence (row sums)
  // {2,1,1}, all unit entries, 4 nonzero entries total (= order 4).
  bool is_double_digon(std::vector<uint8_t> const &m, int V) {
    if (V != 3 || m.size() != 9) return false;
    int ones = 0;
    std::vector<int> row_sum(3, 0);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        uint8_t e = m[i * 3 + j];
        if (e != m[j * 3 + i]) return false; // symmetric
        if (e != 0 && e != 1) return false;  // unit entries
        if (e == 1) { ++ones; ++row_sum[i]; }
        if (i == j && e != 0) return false; // no self-line
      }
    }
    if (ones != 4) return false;
    std::sort(row_sum.begin(), row_sum.end());
    return row_sum[0] == 1 && row_sum[1] == 1 && row_sum[2] == 2;
  }

} // namespace

// At r = (5, 0), order 4 (the maximal x-reach at this order: odd-parity
// sector (0,1) has d_super = 2, requiring 2 NN-dimer hops), the catalog
// contains exactly two topologies, each in its own V bucket:
//   - V = 4: the directed 4-cycle (the n-cycle), marks on opposite vertices.
//   - V = 3: the double-digon "banana" {{0,1,0},{1,0,1},{0,1,0}}, marks on
//            the two degree-2 leaves (graph distance 2 via the center).
// The V = 2 banana (diameter 1) cannot host d_super = 2 and is absent.
TEST(DimerDistanceRootedGenerator, Order4_R5_0_NCycleAndBanana) {
  DimerDistanceRootedDiagramGenerator gen({5, 0}, /*n=*/4);
  gen.generate();
  auto const &by_V = gen.get_rooted_graphs();

  // Exactly the V = 3 and V = 4 buckets are populated, one topology each.
  ASSERT_EQ(by_V.size(), 2u);
  ASSERT_EQ(by_V.count(3), 1u);
  ASSERT_EQ(by_V.count(4), 1u);
  ASSERT_EQ(by_V.at(3).size(), 1u);
  ASSERT_EQ(by_V.at(4).size(), 1u);

  // V = 4 entry: the directed 4-cycle (n-cycle).
  auto const &ncycle = by_V.at(4)[0];
  EXPECT_EQ(ncycle.get_V(), 4);
  EXPECT_EQ(ncycle.get_order(), 4);
  EXPECT_TRUE(is_single_directed_cycle(ncycle.get_canonical_form(), 4))
     << "V=4 topology should be the single directed 4-cycle";

  // V = 3 entry: the double-digon "banana".
  auto const &banana = by_V.at(3)[0];
  EXPECT_EQ(banana.get_V(), 3);
  EXPECT_EQ(banana.get_order(), 4);
  EXPECT_TRUE(is_double_digon(banana.get_canonical_form(), 3))
     << "V=3 topology should be the double-digon {{0,1,0},{1,0,1},{0,1,0}}";
}

namespace {

  // V == 2 banana: the only order-4 two-vertex vacuum graph, two directed
  // lines each way ({{0,2},{2,0}}). Marks never change the adjacency, so the
  // canonical form is fixed.
  bool is_banana_v2(std::vector<uint8_t> const &m) {
    return m == std::vector<uint8_t>{0, 2, 2, 0};
  }

  // Smaller of the two parity-allowed sectors' d_super values for r = (x, y).
  // This single number determines the whole order-4 catalog (see test below).
  int dmin_over_sectors(int x, int y) {
    int parity = ((x + y) % 2 + 2) % 2;
    int dmin   = -1;
    for (int a = 0; a < 2; ++a) {
      for (int b = 0; b < 2; ++b) {
        if (((a + b) % 2 + 2) % 2 != parity) continue;
        int d = DimerDistanceRootedDiagramGenerator::d_super(x, y, a, b);
        if (dmin < 0 || d < dmin) dmin = d;
      }
    }
    return dmin;
  }

} // namespace

// Full order-4 sweep over all displacements r = (x, y). The catalog is
// completely determined by Dmin = min over the two parity-allowed sectors of
// d_super:
//
//   Dmin = 0  -> 9 entries  (V2:2, V3:4, V4:3) — only r = (1,0) and images.
//   Dmin = 1  -> 5 entries  (V2:1, V3:2, V4:2).
//   Dmin = 2  -> 2 entries  (V2:0, V3:1, V4:1).
//   Dmin >= 3 -> unreachable at order 4: constructor throws (n < 2*Dmin).
//
// Rationale: the two allowed sectors are always a dimer-inversion pair, so the
// inversion fold collapses them and each (graph, mark-pair class) contributes
// at most one entry, present iff Dmin <= min(dij, V-1). Sweeping that threshold
// over the three order-4 vacuum graphs (banana V2, double-digon V3, n-cycle V4)
// and their mark-pair classes gives the counts above. This test also checks
// that every emitted topology is structurally the expected graph for its V
// bucket and carries order 4.
TEST(DimerDistanceRootedGenerator, Order4FullSweep) {
  constexpr int n = 4;

  for (int x = -7; x <= 7; ++x) {
    for (int y = -4; y <= 4; ++y) {
      if (x == 0 && y == 0) continue; // on-site (V=1) is handled separately

      int dmin = dmin_over_sectors(x, y);
      SCOPED_TRACE("r = (" + std::to_string(x) + ", " + std::to_string(y) + "), Dmin = "
                   + std::to_string(dmin));

      // Out of order-4 reach: the constructor must reject it.
      if (dmin > 2) {
        EXPECT_THROW(DimerDistanceRootedDiagramGenerator({x, y}, n), std::invalid_argument);
        continue;
      }

      DimerDistanceRootedDiagramGenerator gen({x, y}, n);
      gen.generate();
      auto const &by_V = gen.get_rooted_graphs();

      // Expected per-bucket entry counts as a function of Dmin.
      int exp_v2 = (dmin == 0) ? 2 : (dmin == 1) ? 1 : 0;
      int exp_v3 = (dmin == 0) ? 4 : (dmin == 1) ? 2 : 1;
      int exp_v4 = (dmin == 0) ? 3 : (dmin == 1) ? 2 : 1;

      auto bucket_size = [&](int V) -> int {
        auto it = by_V.find(V);
        return it == by_V.end() ? 0 : static_cast<int>(it->second.size());
      };
      EXPECT_EQ(bucket_size(2), exp_v2);
      EXPECT_EQ(bucket_size(3), exp_v3);
      EXPECT_EQ(bucket_size(4), exp_v4);

      // No buckets outside V in {2,3,4}.
      for (auto const &[V, vec] : by_V) {
        EXPECT_GE(V, 2);
        EXPECT_LE(V, 4);
      }

      // Every entry is the expected structural topology for its bucket, order 4.
      if (auto it = by_V.find(2); it != by_V.end())
        for (auto const &rg : it->second) {
          EXPECT_EQ(rg.get_order(), 4);
          EXPECT_TRUE(is_banana_v2(rg.get_canonical_form())) << "V=2 entry should be the banana";
        }
      if (auto it = by_V.find(3); it != by_V.end())
        for (auto const &rg : it->second) {
          EXPECT_EQ(rg.get_order(), 4);
          EXPECT_TRUE(is_double_digon(rg.get_canonical_form(), 3)) << "V=3 entry should be the double-digon";
        }
      if (auto it = by_V.find(4); it != by_V.end())
        for (auto const &rg : it->second) {
          EXPECT_EQ(rg.get_order(), 4);
          EXPECT_TRUE(is_single_directed_cycle(rg.get_canonical_form(), 4)) << "V=4 entry should be the n-cycle";
        }
    }
  }
}
