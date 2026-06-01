#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "../c++/sc_expansion/cumulant.hpp"
#include "../c++/sc_expansion/generate_diagrams.hpp"
#include "../c++/sc_expansion/graph.hpp"
#include "../c++/sc_expansion/hubbard_solver.hpp"
#include "../c++/sc_expansion/dimer/diagram.hpp"
#include "../c++/sc_expansion/dimer/sum_diagrams.hpp"
#include "../c++/sc_expansion/dimer/vertex.hpp"

using namespace sc_expansion;
using namespace sc_expansion::dimer;

// =============================================================================
//  Rooted dimer::Diagram (task 3): the mark-constrained spatial embedding
//  (compute_spatial_configurations_rooted) + the static-density decoration of
//  the marked vertices' cumulants. See dimer_density_density_correlator-03.md.
//
//  Orbital encoding (N_sites=2): orbital = site + spin*2, i.e.
//    0=site0↓, 1=site1↓, 2=site0↑, 3=site1↑.
// =============================================================================

namespace {
  // vt[k] is the cumulant-order-(k+1) VertexType. Sized to max_co so a marked
  // degree-0 vertex (mark_bonus 2 ⇒ index 1) still resolves to a non-null type.
  template <typename T> std::vector<std::unique_ptr<VertexType<T>>> make_vertex_types(int max_co) {
    std::vector<std::unique_ptr<VertexType<T>>> out;
    for (int k = 1; k <= max_co; ++k) out.push_back(std::make_unique<VertexType<T>>(2 * k));
    return out;
  }
  template <typename T> std::vector<VertexType<T> *> to_ptrs(std::vector<std::unique_ptr<VertexType<T>>> &owned) {
    std::vector<VertexType<T> *> ptrs;
    for (auto &up : owned) ptrs.push_back(up.get());
    return ptrs;
  }
} // namespace

// -----------------------------------------------------------------------------
//  On-site correlator ⟨n(0)n(0)⟩ at order 2. Generating the rooted catalog at
//  r=(0,0), n=2 yields a single vacuum topology — the directed digon, vacuum
//  adjmat {{0,1},{1,0}}. The on-site entry is the coincident-mark one (both
//  density insertions on the same dimer at the same within-dimer site).
//
//  Wiring that entry into a rooted Diagram, the mark-constrained embedding pins
//  the marked dimer at (0,0) and sweeps the digon's partner over its 6 NN dimers,
//  producing exactly TWO spatial configurations: the hopping exchange either
//  leaves/enters the marked dimer at the SAME site the density sits on, or at the
//  OTHER within-dimer site. For {0,1,1,0} the line layout is line0 = 0→1,
//  line1 = 1→0, so with both marks on vertex 0 at site 0:
//    - dirs {1,0}: both lines touch vertex 0 at site 0  → density ON the hopping
//      site (same site as the exchange);
//    - dirs {0,1}: both lines touch vertex 0 at site 1  → density on the OTHER
//      site (exchange does not touch the density site).
//  Each pattern folds 3 of the 6 NN-dimer directions ⇒ weight 3, total 6.
// -----------------------------------------------------------------------------
TEST(RootedDimerDiagram, OnSiteR00_Order2_DigonTwoEmbeddings) {
  // (a) Catalog: the only order-2 vacuum topology for the on-site correlator is
  //     the digon, and exactly one catalog entry is the coincident (on-site) one.
  DimerDistanceRootedDiagramGenerator gen({0, 0}, /*n=*/2);
  gen.generate();
  auto const &by_V = gen.get_rooted_graphs();

  // Only V=2 fits at order 2 (the V=1 on-site zeroth-order term is separate).
  ASSERT_EQ(by_V.size(), 1u);
  ASSERT_EQ(by_V.count(2), 1u);

  std::set<std::vector<uint8_t>> topologies;
  int n_coincident = 0;
  for (auto const &rg : by_V.at(2)) {
    topologies.insert(rg.get_canonical_form());
    EXPECT_EQ(rg.get_order(), 2);
    if (rg.get_marks()[0] == rg.get_marks()[1]) { // the coincident on-site entry
      ++n_coincident;
      EXPECT_EQ(rg.get_sites(), (std::vector<int>{0, 0})); // both densities on one site
    }
  }
  ASSERT_EQ(topologies.size(), 1u); // one topology: the digon
  EXPECT_EQ(*topologies.begin(), (std::vector<uint8_t>{0, 1, 1, 0}));
  EXPECT_EQ(n_coincident, 1);

  // (b) Embedding: build the rooted Diagram for the on-site digon, both marks on
  //     vertex 0 at within-dimer site 0. Two spatial embeddings, each weight 3.
  Graph g({0, 1, 1, 0}, 2, /*bipartite_only=*/false);
  auto owned = make_vertex_types<double>(/*max_co=*/2);
  auto vt    = to_ptrs(owned);

  Diagram<double> diag(g, vt, /*marks=*/{0, 0}, /*sites=*/{0, 0}, /*mark_spins=*/{0, 1}, /*r=*/{0, 0});

  auto const &cfgs = diag.get_spatial_configurations();
  ASSERT_EQ(cfgs.size(), 2u);

  std::map<std::vector<uint8_t>, double> weight_of;
  for (auto const &c : cfgs) weight_of[c.directions] = c.weight;

  std::vector<uint8_t> const same_site = {1, 0}; // hopping touches the marked site 0
  std::vector<uint8_t> const off_site  = {0, 1}; // hopping touches the other site 1
  ASSERT_EQ(weight_of.count(same_site), 1u) << "missing the same-site embedding";
  ASSERT_EQ(weight_of.count(off_site), 1u) << "missing the off-site embedding";
  EXPECT_DOUBLE_EQ(weight_of[same_site], 3.0);
  EXPECT_DOUBLE_EQ(weight_of[off_site], 3.0);
  EXPECT_DOUBLE_EQ(diag.get_free_multiplicity(), 6.0);
}

// -----------------------------------------------------------------------------
//  Intra-dimer displacement r=(1,0) at order 2. r is parity-odd, so the two
//  parity-allowed sectors are (0,1) and (1,0); the catalog folds them into TWO
//  entries, BOTH on the same digon vacuum topology {{0,1},{1,0}} and both with
//  within-dimer sites {0,1} (see project_dimer_rooted_constraints):
//
//   (1) COINCIDENT entry — both densities on ONE dimer at sites {0,1} (the
//       order-2 hopping correction to the intra-dimer ⟨n_0 n_1⟩). The marked
//       dimer is pinned and the digon's partner sweeps its 6 NN cells, splitting
//       by which site of the marked dimer the exchange touches:
//         dirs {0,1} (exchange uses site 1) and dirs {1,0} (exchange uses site 0),
//       each folding 3 of the 6 NN directions ⇒ weight 3, total free-mult 6.
//
//       These TWO configs are the UNFOLDED dimer-inversion orbit: site inversion
//       (0↔1) maps {0,1}↔{1,0} while leaving the density product n_0 n_1
//       invariant, so they carry EQUAL value and evaluate() equals the folded
//       "1 config × weight 6" form. The rooted canonicalisation deliberately
//       drops lattice inversion (include_inversion=false) — turning it on would
//       also wrongly fold the r=(0,0) on-site {0,0} orbit (the test above), which
//       is NOT inversion-symmetric. So the orbit stays unfolded here; this is
//       open question #1 in dimer_density_density_correlator-03.md, and a pure
//       representation choice (the correlator value is unaffected). This test
//       pins the current 2×3 form on purpose.
//
//   (2) DIFFERENT-dimer entry — densities on two distinct dimers: exactly ONE
//       spatial config, dirs {1,0}, weight 1. With mark0 pinned at site A(0) of
//       the origin dimer, the +r anchoring forces mark1's dimer onto the origin
//       cell (a vertex collision) and is dropped; only the -r anchoring survives,
//       placing mark1 at site B(1) of the dimer one NN step away. That is the
//       inversion/anchor partner of the "site B of dimer 0 ↔ site A of the dimer
//       directly to the right" picture — the same physical placement read from
//       the opposite endpoint (⟨n(r)n(0)⟩ = ⟨n(-r)n(0)⟩).
// -----------------------------------------------------------------------------
TEST(RootedDimerDiagram, R10_Order2_CoincidentAndDifferentDimer) {
  DimerDistanceRootedDiagramGenerator gen({1, 0}, /*n=*/2);
  gen.generate();
  auto const &by_V = gen.get_rooted_graphs();

  // Two catalog entries, both at V=2 and both on the one digon topology.
  ASSERT_EQ(by_V.size(), 1u);
  ASSERT_EQ(by_V.count(2), 1u);
  ASSERT_EQ(by_V.at(2).size(), 2u);

  std::vector<uint8_t> const digon = {0, 1, 1, 0};
  std::vector<uint8_t> const d01   = {0, 1};
  std::vector<uint8_t> const d10   = {1, 0};

  int n_coincident = 0, n_different = 0;
  for (auto const &rg : by_V.at(2)) {
    EXPECT_EQ(rg.get_canonical_form(), digon);
    EXPECT_EQ(rg.get_order(), 2);
    EXPECT_EQ(rg.get_sites(), (std::vector<int>{0, 1})); // sites {0,1} in both entries

    Graph g(rg.get_canonical_form(), rg.get_V(), /*bipartite_only=*/false);
    auto owned = make_vertex_types<double>(/*max_co=*/2);
    auto vt    = to_ptrs(owned);
    Diagram<double> diag(g, vt, rg.get_marks(), rg.get_sites(), /*mark_spins=*/{0, 1}, /*r=*/{1, 0});

    auto const &cfgs = diag.get_spatial_configurations();
    std::map<std::vector<uint8_t>, double> weight_of;
    for (auto const &c : cfgs) weight_of[c.directions] = c.weight;

    if (rg.get_marks()[0] == rg.get_marks()[1]) {
      // Coincident (intra-dimer): the unfolded inversion orbit — 2 configs, each
      // weight 3, summing to 6 (== the folded 1×6 value).
      ++n_coincident;
      EXPECT_EQ(cfgs.size(), 2u);
      ASSERT_EQ(weight_of.count(d01), 1u) << "missing the site-1-exchange config";
      ASSERT_EQ(weight_of.count(d10), 1u) << "missing the site-0-exchange config";
      EXPECT_DOUBLE_EQ(weight_of[d01], 3.0);
      EXPECT_DOUBLE_EQ(weight_of[d10], 3.0);
      EXPECT_DOUBLE_EQ(diag.get_free_multiplicity(), 6.0);
    } else {
      // Different dimers: exactly one embedding.
      ++n_different;
      EXPECT_EQ(cfgs.size(), 1u);
      ASSERT_EQ(weight_of.count(d10), 1u) << "the single embedding should be dirs {1,0}";
      EXPECT_DOUBLE_EQ(weight_of[d10], 1.0);
      EXPECT_DOUBLE_EQ(diag.get_free_multiplicity(), 1.0);
    }
  }
  EXPECT_EQ(n_coincident, 1);
  EXPECT_EQ(n_different, 1);
}

// -----------------------------------------------------------------------------
//  Embedding weight (definition of done a). The directed digon with marks on
//  the two vertices at within-dimer sites (0, 1) embeds at physical r=(3,0) in
//  exactly ONE way: mark0's dimer at superlattice (0,0) [physical sites 0,1],
//  mark1's forced to (1,0) [physical 2,3], which are NN dimers. The -r=(-3,0)
//  enumeration forces a non-adjacent dimer and contributes nothing at order 2.
//  Pointwise mark-stabiliser is trivial and there are no multi-edges, so the
//  rooted symmetry factor is 1 ⇒ total embedding weight 1.
// -----------------------------------------------------------------------------
TEST(RootedDimerDiagram, DigonR30_EmbeddingWeight) {
  Graph g({0, 1, 1, 0}, 2, /*bipartite_only=*/false);
  auto owned = make_vertex_types<double>(/*max_co=*/2);
  auto vt    = to_ptrs(owned);

  Diagram<double> diag(g, vt, /*marks=*/{0, 1}, /*sites=*/{0, 1}, /*mark_spins=*/{0, 1}, /*r=*/{3, 0});

  EXPECT_EQ(diag.get_spatial_configurations().size(), 1u);
  EXPECT_DOUBLE_EQ(diag.get_free_multiplicity(), 1.0);

  // Spin sum: digon vertex balance forces the two hopping lines to share a
  // spin, so 2 valid (spin-resolved) configs survive, each weight 1 (no spin-
  // flip fold for the rooted path; distinct mark sites ⇒ n_mark_orbit = 1).
  auto const &configs = diag.get_valid_configurations();
  EXPECT_EQ(configs.size(), 2u);
  for (auto const &c : configs) EXPECT_DOUBLE_EQ(c.weight, 1.0);
}

// -----------------------------------------------------------------------------
//  Zero-embedding entry. The digon cannot span r=(5,0) (needs two NN-dimer
//  hops, i.e. ≥ 2 lines), so both the +r and -r pinnings force non-adjacent
//  dimers ⇒ no spatial configs, no valid configs, evaluate() == 0. Mirrors the
//  atomic embedding_counts[i] == 0 skip.
// -----------------------------------------------------------------------------
TEST(RootedDimerDiagram, DigonR50_NoEmbedding) {
  double U = 8.0, beta = 1.0, mu = 2.0, t = 1.0;
  Parameters<double> params{U, beta, mu, t, true};
  HubbardSolver<2, double> solver(params);

  Graph g({0, 1, 1, 0}, 2, /*bipartite_only=*/false);
  auto owned = make_vertex_types<double>(/*max_co=*/2);
  auto vt    = to_ptrs(owned);

  Diagram<double> diag(g, vt, {0, 1}, {0, 1}, {0, 1}, /*r=*/{5, 0});

  EXPECT_EQ(diag.get_spatial_configurations().size(), 0u);
  EXPECT_DOUBLE_EQ(diag.get_free_multiplicity(), 0.0);
  EXPECT_DOUBLE_EQ(diag.evaluate({0.0, 0.0}, solver), 0.0);
}

// -----------------------------------------------------------------------------
//  Order-0 intra-dimer correlator (the simplest density-decorated value). At
//  r=(1,0) two coincident marks sit on a single V=1 vertex at sites (0,1):
//  there are no hopping lines, one embedding (the dimer itself), and the vertex
//  cumulant is the connected intra-dimer density-density correlator
//  ⟨n_{0,s1} n_{1,s2}⟩_c. This pins the density decoration + factored-path
//  wiring against an independent CumulantSolver reference. The rooted prefactor
//  convention is the physical one: ⟨n(r)n(0)⟩_c = -β ∂²Ω/∂J², whose -β cancels
//  the Ω diagram's -1/β, so at order 0 the value is exactly sign·κ₂ (no 1/β).
// -----------------------------------------------------------------------------
TEST(RootedDimerDiagram, IntraDimerR10_MatchesConnectedCorrelator) {
  double U = 4.0, beta = 1.3, mu = 1.1, t = 1.0;
  Parameters<double> params{U, beta, mu, t, true};
  HubbardSolver<2, double> solver(params);

  // V=1 single vertex, no edges. Both marks coincident on vertex 0.
  Graph g(std::vector<uint8_t>{0}, /*V=*/1, /*automorphism_count=*/1, /*symmetry_factor=*/1,
          /*free_multiplicity=*/1, /*bipartite_only=*/false);
  auto owned = make_vertex_types<double>(/*max_co=*/2); // mark_bonus 2 ⇒ needs index 1
  auto vt    = to_ptrs(owned);

  for (int s1 = 0; s1 <= 1; ++s1)
    for (int s2 = 0; s2 <= 1; ++s2) {
      Diagram<double> diag(g, vt, /*marks=*/{0, 0}, /*sites=*/{0, 1}, /*mark_spins=*/{s1, s2}, /*r=*/{1, 0});

      ASSERT_EQ(diag.get_spatial_configurations().size(), 1u);
      EXPECT_DOUBLE_EQ(diag.get_free_multiplicity(), 1.0);

      double val = diag.evaluate(/*taus=*/{}, solver);

      // Reference: connected 2nd cumulant κ₂(n_{site0,s1}, n_{site1,s2}).
      Args<2, double> u({}, {}), p({}, {});
      CumulantSolver<2, double> cs(u, p, solver, /*infinite_U=*/false);
      cs.add_static_density(0 + s1 * 2); // site 0
      cs.add_static_density(1 + s2 * 2); // site 1
      double kappa = cs.compute_cumulant_decomposition();

      double prefactor = diag.get_diagram_sign(); // (-t)^0 = +1; -β cancels the -1/β
      EXPECT_NEAR(val, prefactor * kappa, 1e-12) << "s1=" << s1 << " s2=" << s2;
    }
}

// -----------------------------------------------------------------------------
//  End-to-end finiteness at order 2 (digon, r=(3,0)) with hopping legs present:
//  the density-decorated factored path produces a finite, non-zero value and is
//  cache-stable across repeated / dirty-marked evaluations.
// -----------------------------------------------------------------------------
TEST(RootedDimerDiagram, DigonR30_EvaluatesFiniteAndStable) {
  double U = 8.0, beta = 1.0, mu = 2.0, t = 1.0;
  Parameters<double> params{U, beta, mu, t, true};
  HubbardSolver<2, double> solver(params);

  Graph g({0, 1, 1, 0}, 2, /*bipartite_only=*/false);
  auto owned = make_vertex_types<double>(/*max_co=*/2);
  auto vt    = to_ptrs(owned);

  Diagram<double> diag(g, vt, {0, 1}, {0, 1}, /*mark_spins=*/{0, 1}, /*r=*/{3, 0});

  std::vector<double> taus = {0.4, 0.2};
  double v1                = diag.evaluate(taus, solver);
  EXPECT_TRUE(std::isfinite(v1));
  EXPECT_NE(v1, 0.0);

  EXPECT_DOUBLE_EQ(diag.evaluate(taus, solver), v1); // cache-stable
  diag.mark_all_dirty();
  EXPECT_DOUBLE_EQ(diag.evaluate(taus, solver), v1); // recompute reproduces
}

// -----------------------------------------------------------------------------
// =============================================================================
//  Finite-cluster rooted embedding (compute_spatial_configurations_rooted_cluster).
//  Same 3-dimer triangle used by the free-energy MCMC (test_mcmc_dimer.cpp):
//    A=(0,0), B=(1,0), C=(0,1) on the staggered superlattice. Physical sites:
//    A:(0,0),(1,0)  B:(2,0),(3,0)  C:(1,1),(2,1).
// =============================================================================

// -----------------------------------------------------------------------------
//  Cluster embedding weight (brute-force hand count). Digon, marks {0,1} at
//  sites {0,0}, r=(2,0): mark0 at site0=phys(0,0), mark1 at site0=phys(2,0), so
//  the displacement is (2,0). On the triangle exactly two (anchor, ±r) pinnings
//  keep mark1's forced dimer on the cluster:
//    (+r, anchor A): mark0@A, mark1@B  → dirs {0,1}
//    (-r, anchor B): mark0@B, mark1@A  → dirs {1,0}   (value-equal site1-site1
//                                                       inversion partner)
//  every other anchor×disp forces mark1 onto a non-cluster dimer. Pointwise
//  mark-stabiliser is trivial and there are no multi-edges ⇒ rooted_sym_factor 1;
//  divide the 2 raw embeddings by n_cluster_sites = 3 ⇒ total weight 2/3 across
//  two distinct direction patterns.
// -----------------------------------------------------------------------------
TEST(RootedDimerDiagram, ClusterDigonR20_EmbeddingWeight) {
  std::vector<std::pair<int, int>> cluster = {{0, 0}, {1, 0}, {0, 1}};

  Graph g({0, 1, 1, 0}, 2, /*bipartite_only=*/false);
  auto owned = make_vertex_types<double>(/*max_co=*/2);
  auto vt    = to_ptrs(owned);

  Diagram<double> diag(g, vt, /*marks=*/{0, 1}, /*sites=*/{0, 0}, /*mark_spins=*/{0, 1}, /*r=*/{2, 0}, cluster, /*n_cluster_sites=*/3);

  EXPECT_EQ(diag.get_spatial_configurations().size(), 2u);
  EXPECT_DOUBLE_EQ(diag.get_free_multiplicity(), 2.0 / 3.0);
}

// -----------------------------------------------------------------------------
//  Cluster ⟶ infinite agreement for the intra-dimer (V=1) correlator. The single
//  marked dimer is translation-invariant, so sweeping mark0's home over the three
//  cluster dimers and dividing by n_cluster_sites=3 must reproduce the infinite-
//  lattice value exactly (free_multiplicity 1.0, identical evaluate()). Pins the
//  cluster normalisation against the already-validated infinite path, no ED needed.
// -----------------------------------------------------------------------------
TEST(RootedDimerDiagram, ClusterIntraDimerR10_MatchesInfinite) {
  double U = 4.0, beta = 1.3, mu = 1.1, t = 1.0;
  Parameters<double> params{U, beta, mu, t, true};
  HubbardSolver<2, double> solver(params);

  std::vector<std::pair<int, int>> cluster = {{0, 0}, {1, 0}, {0, 1}};
  Graph g(std::vector<uint8_t>{0}, /*V=*/1, /*automorphism_count=*/1, /*symmetry_factor=*/1,
          /*free_multiplicity=*/1, /*bipartite_only=*/false);

  for (int s1 = 0; s1 <= 1; ++s1)
    for (int s2 = 0; s2 <= 1; ++s2) {
      auto owned_inf = make_vertex_types<double>(2);
      auto vt_inf    = to_ptrs(owned_inf);
      auto owned_cl  = make_vertex_types<double>(2);
      auto vt_cl     = to_ptrs(owned_cl);

      Diagram<double> inf(g, vt_inf, {0, 0}, {0, 1}, {s1, s2}, {1, 0});
      Diagram<double> cl(g, vt_cl, {0, 0}, {0, 1}, {s1, s2}, {1, 0}, cluster, /*n_cluster_sites=*/3);

      EXPECT_DOUBLE_EQ(cl.get_free_multiplicity(), 1.0);
      EXPECT_NEAR(cl.evaluate({}, solver), inf.evaluate({}, solver), 1e-12) << "s1=" << s1 << " s2=" << s2;
    }
}

// -----------------------------------------------------------------------------
//  SumDiagrams cluster rooted constructor smoke test: the owner builds the rooted
//  catalog for r and constructs one cluster-restricted Diagram per entry, each
//  flagged rooted with the ctor's mark spins and target r. density_density() runs
//  end-to-end on the cluster and is finite + cache-stable.
// -----------------------------------------------------------------------------
TEST(RootedDimerDiagram, ClusterSumDiagramsConstructsAndEvaluates) {
  double U = 8.0, beta = 1.0, mu = 2.0, t = 1.0;
  Parameters<double> params{U, beta, mu, t, /*bipartite=*/false};
  std::vector<std::pair<int, int>> cluster = {{0, 0}, {1, 0}, {0, 1}};
  std::vector<int> const r                 = {2, 0};
  int const order                          = 2;

  SumDiagrams<double> sd(params, order, r, /*s1=*/0, /*s2=*/1, cluster, /*n_cluster_sites=*/3);

  EXPECT_TRUE(sd.is_density_density_mode());
  EXPECT_EQ(sd.get_target_r(), r);
  EXPECT_GE(sd.get_n_diagrams(), 1);
  for (auto const &d : sd.get_diagrams()) {
    EXPECT_TRUE(d.is_rooted_diagram());
    EXPECT_EQ(d.get_mark_spins(), (std::vector<int>{0, 1}));
    EXPECT_EQ(d.get_target_r(), r);
  }

  std::vector<double> taus(order, 0.5);
  double v1 = sd.density_density(taus);
  EXPECT_TRUE(std::isfinite(v1));
  EXPECT_EQ(sd.density_density(taus), v1); // cache-stable
  sd.mark_all_dirty();
  EXPECT_EQ(sd.density_density(taus), v1); // recompute reproduces
}

// =============================================================================
//  Independent embedding-weight reference for the INFINITE (thermodynamic-limit)
//  rooted path. Re-derives the staggered-dimer geometry directly in physical Z²
//  domino coordinates, with NO dependence on diagram.cpp's (u,v) superlattice
//  tables (tri_offsets / tri_bond_label / solve_dimer_embedding) — so comparing
//  it against Diagram::get_spatial_configurations() validates that geometry
//  rather than re-asserting it. Extends the order-2 hand checks to the order-4
//  catalog, where bigger finite clusters get expensive and still do not probe
//  the TD limit.
//
//  Physical model (independent re-derivation, cross-checked against the (u,v)
//  tables documented in diagram.cpp):
//    * The underlying lattice is the square lattice Z²; hopping connects
//      nearest-neighbour sites (|Δx| + |Δy| == 1).
//    * A dimer is a horizontal domino identified by its LEFT cell (X,Y) with
//      X ≡ Y (mod 2): it covers site0 = (X,Y) and site1 = (X+1,Y). This is the
//      image of superlattice (u,v) under X = 2u + v%2, Y = v.
//    * Two dimers are NN iff a site of one is a square-NN of a site of the
//      other; each NN pair has exactly ONE such external bond, connecting
//      site (1−L) of the source to site L of the destination, where the bond
//      label L equals the destination within-dimer site. A dimer's 6 NN dominoes
//      are the 8 square-NN of its two sites minus the 2 internal site0–site1
//      bonds.
//    * Lattice sums are UNRESTRICTED (see project_unrestricted_lattice_sums):
//      distinct vertices MAY share a dimer; only EDGE-LINKED vertices cannot
//      (their intra-dimer "bond" carries no label ⇒ that embedding is rejected).
//    * Both ±r anchorings are summed and mark0's dimer is pinned with left cell
//      (0,0), exactly as compute_spatial_configurations_rooted.
//
//  Compared quantities: free_multiplicity (total raw embeddings) and the
//  per-direction weights aggregated by the SORTED label multiset. Sorting is
//  what makes the comparison robust to the code's rooted-automorphism
//  canonicalisation: that step only permutes line positions (no label flip in
//  the rooted path), so every raw embedding's label multiset is preserved within
//  its orbit, and summing the code's configs by sorted(directions) must equal the
//  reference's raw counts by sorted(directions).
// =============================================================================
namespace {

  struct RefEmbedding {
    double free_mult = 0.0;
    std::map<std::vector<int>, double> by_sorted_dirs;
  };

  // The 6 physical NN dominoes (left cells) of the domino with left cell (X,Y).
  std::array<std::pair<int, int>, 6> phys_nn_dimers(int X, int Y) {
    return {{{X - 2, Y}, {X - 1, Y + 1}, {X - 1, Y - 1}, {X + 2, Y}, {X + 1, Y + 1}, {X + 1, Y - 1}}};
  }

  // Bond label of directed line a→b (== destination within-dimer site), or -1 if
  // the two dominoes are not NN (including the same domino: an intra-dimer bond,
  // which has no unique external site-pair).
  int phys_bond_label(int Xa, int Ya, int Xb, int Yb) {
    if (Xa == Xb && Ya == Yb) return -1;
    int label = -1, n = 0;
    for (int sa = 0; sa < 2; ++sa)
      for (int sb = 0; sb < 2; ++sb)
        if (std::abs((Xa + sa) - (Xb + sb)) + std::abs(Ya - Yb) == 1) {
          label = sb;
          ++n;
        }
    return n == 1 ? label : -1;
  }

  // Place every still-unplaced vertex at an NN domino of an already-placed
  // graph-neighbour (the complete candidate set, since edges require NN dimers),
  // filtered by NN-consistency with all placed neighbours. On a full placement,
  // read each directed line's label; reject if any edge maps to a non-NN pair.
  void ref_dfs(int V, std::vector<std::vector<int>> const &adj, std::vector<std::pair<int, int>> const &lines, std::vector<int> &X,
               std::vector<int> &Y, std::vector<bool> &placed, int placed_count, RefEmbedding &out) {
    if (placed_count == V) {
      std::vector<int> dirs;
      dirs.reserve(lines.size());
      for (auto const &ln : lines) {
        int lbl = phys_bond_label(X[ln.first], Y[ln.first], X[ln.second], Y[ln.second]);
        if (lbl < 0) return; // an edge maps to a non-NN (or intra-dimer) pair ⇒ invalid
        dirs.push_back(lbl);
      }
      out.free_mult += 1.0;
      std::sort(dirs.begin(), dirs.end());
      out.by_sorted_dirs[dirs] += 1.0;
      return;
    }

    int target = -1, anchor = -1;
    for (int c = 0; c < V && target < 0; ++c)
      if (!placed[c])
        for (int p = 0; p < V; ++p)
          if (placed[p] && adj[c][p] + adj[p][c] > 0) {
            target = c;
            anchor = p;
            break;
          }
    if (target < 0) return; // remaining vertices are disconnected from the placed set

    for (auto const &d : phys_nn_dimers(X[anchor], Y[anchor])) {
      bool ok = true;
      for (int i = 0; i < V && ok; ++i)
        if (placed[i] && adj[target][i] + adj[i][target] > 0 && phys_bond_label(d.first, d.second, X[i], Y[i]) < 0) ok = false;
      if (!ok) continue;
      X[target]      = d.first;
      Y[target]      = d.second;
      placed[target] = true;
      ref_dfs(V, adj, lines, X, Y, placed, placed_count + 1, out);
      placed[target] = false;
    }
  }

  RefEmbedding ref_embeddings(std::vector<std::vector<int>> const &adj, int V, int m0, int m1, int s0, int s1, int rx, int ry) {
    RefEmbedding out;

    std::vector<std::pair<int, int>> lines; // directed hopping lines (i→j with multiplicity)
    for (int i = 0; i < V; ++i)
      for (int j = 0; j < V; ++j)
        for (int k = 0; k < adj[i][j]; ++k) lines.push_back({i, j});

    std::vector<std::pair<int, int>> disps = {{rx, ry}};
    if (!(rx == 0 && ry == 0)) disps.push_back({-rx, -ry});

    bool coincident = (m0 == m1);
    for (auto const &disp : disps) {
      int dx = disp.first, dy = disp.second;
      std::vector<int> X(V, 0), Y(V, 0); // mark0's domino left cell is (0,0)
      std::vector<bool> placed(V, false);
      int placed_count = 0;

      if (coincident) {
        if (!(dy == 0 && dx == s1 - s0)) continue; // both densities on one domino ⇒ r = (s1−s0, 0)
        placed[m0]   = true;
        placed_count = 1;
      } else {
        int Xm1 = s0 + dx - s1, Ym1 = dy;             // mark1's forced left cell (physical pinning)
        if (((Xm1 - Ym1) % 2 + 2) % 2 != 0) continue; // not a valid domino left cell ⇒ no embedding
        placed[m0]   = true;
        X[m1]        = Xm1;
        Y[m1]        = Ym1;
        placed[m1]   = true;
        placed_count = 2;
      }
      ref_dfs(V, adj, lines, X, Y, placed, placed_count, out);
    }
    return out;
  }

  // Code's spatial configs collapsed by the sorted label multiset (see header).
  std::map<std::vector<int>, double> code_by_sorted(Diagram<double> const &diag) {
    std::map<std::vector<int>, double> m;
    for (auto const &sc : diag.get_spatial_configurations()) {
      std::vector<int> key(sc.directions.begin(), sc.directions.end());
      std::sort(key.begin(), key.end());
      m[key] += sc.weight;
    }
    return m;
  }

  void expect_ref_matches_code(std::vector<uint8_t> const &cf, int V, std::vector<int> const &marks, std::vector<int> const &sites,
                               std::vector<int> const &r) {
    std::vector<std::vector<int>> adj(V, std::vector<int>(V, 0));
    for (int i = 0; i < V; ++i)
      for (int j = 0; j < V; ++j) adj[i][j] = cf[i * V + j];

    RefEmbedding ref = ref_embeddings(adj, V, marks[0], marks[1], sites[0], sites[1], r[0], r[1]);

    Graph g(cf, V, /*bipartite_only=*/false);
    auto owned = make_vertex_types<double>(/*max_co=*/6);
    auto vt    = to_ptrs(owned);
    Diagram<double> diag(g, vt, marks, sites, /*mark_spins=*/{0, 1}, r);
    auto code = code_by_sorted(diag);

    std::ostringstream ctx;
    ctx << "V=" << V << " marks={" << marks[0] << "," << marks[1] << "} sites={" << sites[0] << "," << sites[1] << "} r=(" << r[0]
        << "," << r[1] << ") adj=";
    for (auto b : cf) ctx << (int)b;
    std::string const where = ctx.str();

    EXPECT_DOUBLE_EQ(diag.get_free_multiplicity(), ref.free_mult) << where;
    EXPECT_EQ(code.size(), ref.by_sorted_dirs.size()) << where;
    for (auto const &kv : ref.by_sorted_dirs) {
      ASSERT_EQ(code.count(kv.first), 1u) << where << " (reference direction class absent from code)";
      EXPECT_DOUBLE_EQ(code[kv.first], kv.second) << where;
    }
  }

} // namespace

// Lock the independent reference itself to the order-2 values established by the
// hand-counted tests above, so a later mismatch in the catalog sweep is
// unambiguous (a code regression, not a reference bug).
TEST(RootedDimerReference, ReproducesKnownOrder2ByHand) {
  std::vector<std::vector<int>> const digon = {{0, 1}, {1, 0}};
  std::vector<int> const both               = {0, 1}; // every digon embedding's label multiset

  // on-site coincident ⟨n(0)n(0)⟩: the digon partner sweeps its 6 NN dominoes.
  RefEmbedding onsite = ref_embeddings(digon, 2, /*m*/ 0, 0, /*s*/ 0, 0, /*r*/ 0, 0);
  EXPECT_DOUBLE_EQ(onsite.free_mult, 6.0);
  EXPECT_EQ(onsite.by_sorted_dirs.size(), 1u);
  EXPECT_DOUBLE_EQ(onsite.by_sorted_dirs.at(both), 6.0);

  // intra-dimer coincident r=(1,0): the same 6 embeddings.
  RefEmbedding intra = ref_embeddings(digon, 2, 0, 0, 0, 1, 1, 0);
  EXPECT_DOUBLE_EQ(intra.free_mult, 6.0);
  EXPECT_DOUBLE_EQ(intra.by_sorted_dirs.at(both), 6.0);

  // intra-dimer different-vertex r=(1,0): one embedding (only the −r anchor; the
  // +r anchor collides the two edge-linked marks on one domino and is rejected).
  RefEmbedding diff = ref_embeddings(digon, 2, 0, 1, 0, 1, 1, 0);
  EXPECT_DOUBLE_EQ(diff.free_mult, 1.0);
  EXPECT_DOUBLE_EQ(diff.by_sorted_dirs.at(both), 1.0);

  // distinct-dimer r=(3,0): one embedding (NN dominoes 2 cells apart).
  RefEmbedding r30 = ref_embeddings(digon, 2, 0, 1, 0, 1, 3, 0);
  EXPECT_DOUBLE_EQ(r30.free_mult, 1.0);
  EXPECT_DOUBLE_EQ(r30.by_sorted_dirs.at(both), 1.0);

  // distinct-dimer r=(5,0): the digon cannot span two NN steps ⇒ no embedding.
  RefEmbedding r50 = ref_embeddings(digon, 2, 0, 1, 0, 1, 5, 0);
  EXPECT_DOUBLE_EQ(r50.free_mult, 0.0);
  EXPECT_TRUE(r50.by_sorted_dirs.empty());
}

// The TD-limit validation: for every rooted catalog topology at orders 2 and 4
// and a spread of displacements, the production embedding (Diagram's (u,v)
// solve_dimer_embedding + tri_offsets/tri_bond_label) must reproduce the
// independent physical-coordinate reference — both the total multiplicity and
// the per-direction-class weights.
TEST(RootedDimerReference, CodeMatchesReference_Orders2And4) {
  std::vector<std::vector<int>> const rs = {{0, 0}, {1, 0}, {2, 0}, {1, 1}, {3, 0}, {2, 1}};
  int n_entries                          = 0;
  for (int n : {2, 4}) {
    for (auto const &r : rs) {
      DimerDistanceRootedDiagramGenerator gen({r[0], r[1]}, n);
      gen.generate();
      for (auto const &vg : gen.get_rooted_graphs())
        for (auto const &rg : vg.second) {
          expect_ref_matches_code(rg.get_canonical_form(), rg.get_V(), rg.get_marks(), rg.get_sites(), r);
          ++n_entries;
        }
    }
  }
  EXPECT_GT(n_entries, 0); // guard: the sweep actually exercised catalog entries
}
