#include <gtest/gtest.h>
#include <cmath>
#include <set>
#include <vector>
#include "../c++/sc_expansion/dimer/sum_diagrams.hpp"
#include "../c++/sc_expansion/generate_diagrams.hpp"
#include "../c++/sc_expansion/graph.hpp"

using namespace sc_expansion;
using sc_expansion::dimer::build_rooted_catalog;
using sc_expansion::dimer::SumDiagrams;

// =============================================================================
//  dimer::SumDiagrams (task 2) — structural contract.
//
//  These tests pin the diagram-list OWNER: catalog flattening, the two
//  constructors (catalog-building and prebuilt), accessor wiring, and
//  dirty-marking. They deliberately do NOT check density_density VALUES: the
//  rooted dimer::Diagram embedding + density decoration land in task 3
//  (dimer_density_density_correlator-03.md); until then the rooted Diagram is a
//  well-formed but non-contributing stub (empty spatial configs ⇒ evaluate()
//  returns 0). The value coverage lives with that task. What we lock in here is
//  the contract task 3 must keep satisfying.
// =============================================================================

namespace {
  // (marks, sites) pair as a comparable key, for multiset comparisons that are
  // independent of the V-descending diagram ordering inside SumDiagrams.
  std::vector<int> mark_site_key(std::vector<int> const &marks, std::vector<int> const &sites) {
    std::vector<int> key = marks;
    key.insert(key.end(), sites.begin(), sites.end());
    return key;
  }
} // namespace

// -----------------------------------------------------------------------------
//  build_rooted_catalog: exact flattening at r = (3, 0), order 2.
//
//  The generator (test_distance_rooted_generator.cpp) yields exactly one rooted
//  topology here: the directed digon with marks at vertices {0, 1} and sites
//  {0, 1}. build_rooted_catalog must surface that as parallel single-entry
//  (graphs, marks, sites) vectors, with the graph carrying free_multiplicity = 1.
// -----------------------------------------------------------------------------
TEST(DimerSumDiagrams, BuildRootedCatalogExactAt30Order2) {
  std::vector<Graph> graphs;
  std::vector<std::vector<int>> marks;
  std::vector<std::vector<int>> sites;
  build_rooted_catalog(/*order=*/2, {3, 0}, graphs, marks, sites);

  ASSERT_EQ(graphs.size(), 1u);
  ASSERT_EQ(marks.size(), 1u);
  ASSERT_EQ(sites.size(), 1u);

  EXPECT_EQ(graphs[0].get_V(), 2);
  EXPECT_EQ(graphs[0].get_canonical_form(), (std::vector<uint8_t>{0, 1, 1, 0}));
  // The Graph override ctor pins free_multiplicity = 1 for every rooted entry.
  EXPECT_DOUBLE_EQ(graphs[0].get_free_multiplicity(), 1.0);

  EXPECT_EQ(marks[0], (std::vector<int>{0, 1}));
  EXPECT_EQ(sites[0], (std::vector<int>{0, 1}));
}

// build_rooted_catalog must agree, entry-for-entry, with the underlying
// generator output it flattens (same graphs/marks/sites, just unbucketed).
TEST(DimerSumDiagrams, BuildRootedCatalogMatchesGenerator) {
  std::vector<int> const r = {1, 0};
  int const order          = 4;

  std::vector<Graph> graphs;
  std::vector<std::vector<int>> marks;
  std::vector<std::vector<int>> sites;
  build_rooted_catalog(order, r, graphs, marks, sites);

  ASSERT_EQ(graphs.size(), marks.size());
  ASSERT_EQ(graphs.size(), sites.size());

  DimerDistanceRootedDiagramGenerator gen(r, order);
  gen.generate();
  size_t gen_total = 0;
  std::multiset<std::vector<int>> gen_keys;
  for (auto const &[V, bucket] : gen.get_rooted_graphs()) {
    (void)V;
    gen_total += bucket.size();
    for (auto const &rg : bucket) gen_keys.insert(mark_site_key(rg.get_marks(), rg.get_sites()));
  }

  ASSERT_GT(gen_total, 0u) << "fixture expects a non-trivial catalog";
  EXPECT_EQ(graphs.size(), gen_total);

  std::multiset<std::vector<int>> catalog_keys;
  for (size_t i = 0; i < marks.size(); ++i) catalog_keys.insert(mark_site_key(marks[i], sites[i]));
  EXPECT_EQ(catalog_keys, gen_keys);

  for (auto const &g : graphs) EXPECT_DOUBLE_EQ(g.get_free_multiplicity(), 1.0);
}

// -----------------------------------------------------------------------------
//  Rooted constructor wiring: density-density mode flag, target r, per-diagram
//  rooted state, and that the diagram count matches the flattened catalog.
// -----------------------------------------------------------------------------
TEST(DimerSumDiagrams, RootedConstructorWiring) {
  double const U = 8.0, beta = 1.0, mu = 2.0, t = 1.0;
  Parameters<double> params{U, beta, mu, t, true};
  std::vector<int> const r = {1, 0};
  int const order          = 4;
  int const s1 = 0, s2 = 1; // mark spins: down, up

  std::vector<Graph> cat_graphs;
  std::vector<std::vector<int>> cat_marks;
  std::vector<std::vector<int>> cat_sites;
  build_rooted_catalog(order, r, cat_graphs, cat_marks, cat_sites);
  ASSERT_GE(cat_graphs.size(), 1u) << "fixture expects a non-empty rooted catalog";

  SumDiagrams<double> sd(params, order, r, s1, s2);

  EXPECT_TRUE(sd.is_density_density_mode());
  EXPECT_EQ(sd.get_target_r(), r);
  EXPECT_EQ(sd.get_n_diagrams(), (int)cat_graphs.size());
  EXPECT_EQ((int)sd.get_diagrams().size(), (int)cat_graphs.size());
  EXPECT_EQ((int)sd.get_graphs().size(), (int)cat_graphs.size());

  // Diagrams are stored V-descending; verify the ordering invariant.
  auto const &diagrams = sd.get_diagrams();
  for (size_t i = 1; i < diagrams.size(); ++i) {
    EXPECT_GE(diagrams[i - 1].get_graph().get_V(), diagrams[i].get_graph().get_V());
  }

  // Every diagram is rooted, carries the ctor's mark spins and target r, and has
  // two marks / two sites. The (marks, sites) multiset must equal the catalog's
  // (independent of the V-descending reordering).
  std::multiset<std::vector<int>> diag_keys, cat_keys;
  for (auto const &d : diagrams) {
    EXPECT_TRUE(d.is_rooted_diagram());
    EXPECT_EQ(d.get_mark_spins(), (std::vector<int>{s1, s2}));
    EXPECT_EQ(d.get_target_r(), r);
    ASSERT_EQ(d.get_marks().size(), 2u);
    ASSERT_EQ(d.get_sites().size(), 2u);
    diag_keys.insert(mark_site_key(d.get_marks(), d.get_sites()));
  }
  for (size_t i = 0; i < cat_marks.size(); ++i) cat_keys.insert(mark_site_key(cat_marks[i], cat_sites[i]));
  EXPECT_EQ(diag_keys, cat_keys);
}

// -----------------------------------------------------------------------------
//  Prebuilt-catalog constructor (the MPI path) yields the same diagram list as
//  the catalog-building constructor for the same (r, order, spins).
// -----------------------------------------------------------------------------
TEST(DimerSumDiagrams, PrebuiltConstructorMatchesBuilt) {
  double const U = 8.0, beta = 1.0, mu = 2.0, t = 1.0;
  Parameters<double> params{U, beta, mu, t, true};
  std::vector<int> const r = {1, 0};
  int const order          = 4;
  int const s1 = 1, s2 = 1;

  std::vector<Graph> graphs;
  std::vector<std::vector<int>> marks;
  std::vector<std::vector<int>> sites;
  build_rooted_catalog(order, r, graphs, marks, sites);

  SumDiagrams<double> built(params, order, r, s1, s2);
  SumDiagrams<double> prebuilt(params, order, graphs, marks, sites, r, s1, s2);

  ASSERT_EQ(built.get_n_diagrams(), prebuilt.get_n_diagrams());
  EXPECT_EQ(prebuilt.get_target_r(), r);
  EXPECT_TRUE(prebuilt.is_density_density_mode());

  // Same V-descending sequence and same per-diagram rooted state.
  auto const &b = built.get_diagrams();
  auto const &p = prebuilt.get_diagrams();
  for (size_t i = 0; i < b.size(); ++i) {
    EXPECT_EQ(b[i].get_graph().get_V(), p[i].get_graph().get_V());
    EXPECT_EQ(b[i].get_marks(), p[i].get_marks());
    EXPECT_EQ(b[i].get_sites(), p[i].get_sites());
    EXPECT_EQ(p[i].get_mark_spins(), (std::vector<int>{s1, s2}));
  }
}

// -----------------------------------------------------------------------------
//  density_density evaluator: runs end-to-end over the rooted diagram list and
//  returns a finite scalar, stable across repeated calls and after dirty-
//  marking. (Under the task-2 stub the value is 0 because the rooted Diagram
//  embedding/decoration is not yet wired — task 3 supplies the physics; the
//  point here is that the OWNER drives evaluate() correctly without crashing.)
// -----------------------------------------------------------------------------
TEST(DimerSumDiagrams, DensityDensityEvaluatesAndIsStable) {
  double const U = 8.0, beta = 1.0, mu = 2.0, t = 1.0;
  Parameters<double> params{U, beta, mu, t, true};
  std::vector<int> const r = {1, 0};
  int const order          = 4;

  SumDiagrams<double> sd(params, order, r, /*s1=*/0, /*s2=*/1);
  std::vector<double> taus(order, 0.5);

  double v1 = sd.density_density(taus);
  EXPECT_TRUE(std::isfinite(v1));

  // Repeated evaluation is deterministic (cache-stable).
  double v2 = sd.density_density(taus);
  EXPECT_EQ(v1, v2);

  // Dirty-marking forces recomputation but must reproduce the same value.
  sd.mark_all_dirty();
  EXPECT_EQ(sd.density_density(taus), v1);

  for (int ti = 0; ti < order; ++ti) sd.mark_tau_dirty(ti);
  EXPECT_EQ(sd.density_density(taus), v1);
}

// The Dual instantiation (used for the ∂/∂μ derivative) constructs and
// evaluates through the same owner without issue.
TEST(DimerSumDiagrams, DualInstantiationConstructsAndEvaluates) {
  Dual U{8.0, 0.0}, beta{1.0, 0.0}, mu{2.0, 1.0}, t{1.0, 0.0};
  Parameters<Dual> params{U, beta, mu, t, true};
  std::vector<int> const r = {1, 0};
  int const order          = 4;

  SumDiagrams<Dual> sd(params, order, r, /*s1=*/0, /*s2=*/1);
  EXPECT_EQ(sd.get_n_diagrams(), (int)sd.get_diagrams().size());

  std::vector<double> taus(order, 0.5);
  Dual v = sd.density_density(taus);
  EXPECT_TRUE(std::isfinite(v.value));
  EXPECT_TRUE(std::isfinite(v.derivative));
}
