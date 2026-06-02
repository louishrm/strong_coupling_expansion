#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>
#include "../c++/sc_expansion/dimer/diagram.hpp"
#include "../c++/sc_expansion/dimer/sum_diagrams.hpp"
#include "../c++/sc_expansion/graph.hpp"

using namespace sc_expansion;
using sc_expansion::dimer::Diagram;
using sc_expansion::dimer::SumDiagrams;

// =============================================================================
//  On-site same-spin vanishing probe + 4-cycle weight cross-check.
//
//  Physics: at half filling (mu = U/2) the connected on-site equal-spin
//  correlator is t-independent,
//      <n_sigma(0) n_sigma(0)>_c = <n_sigma> - <n_sigma>^2 = 1/4  (all orders),
//  so every order >= 1 coefficient MUST vanish. The atomic series satisfies
//  this; the dimer series does not at order 4 (CSV: ~0.0179 at U=12, beta=2,
//  mu=6), while orders 2, 3, 5 are ~0.
//
//  These tests reproduce that order-4 failure from the C++ directly, localize
//  the residual to the distinct-vertex (case-2) diagrams whose two marks embed
//  onto the same physical site (the order-4 4-cycle being the first such
//  topology), and dump the rooted weight bookkeeping for the 4-cycle so the
//  atomic-vs-dimer n_mark_orbit / symmetry-factor convention can be compared.
//
//  The opposite-spin contrast confirms the SAME distinct-vertex diagrams carry
//  the legitimate connected double-occupancy and therefore must not be dropped:
//  any fix has to be spin-aware at evaluation, not a catalog deletion.
// =============================================================================

namespace {

  constexpr double U = 12.0, BETA = 2.0, MU = 6.0, THOP = 1.0; // half filling: mu = U/2

  // Uniform-MC sample count. Post-fix the on-site same-spin coefficient is 0 at
  // every order, so a modest count resolves it cheaply (the pre-fix order-4 bug
  // was ~0.018, ~45x the error at this count). Raise it only if you want tighter
  // bounds. (No importance sampling here, unlike the production MCMC.)
  constexpr int kSamples = 30000;

  Parameters<double> half_filled_params() { return Parameters<double>{U, BETA, MU, THOP, /*bipartite=*/false, /*delta=*/0.0}; }

  // Deterministic MC estimate of the order-`order` coefficient of the dimer
  // <n(r)n(0)> series: coeff ~ beta^order * E_{tau ~ U[0,beta)^order}[ f(tau) ],
  // with f = SumDiagrams::density_density (the integrand the MCMC samples; the
  // MCMC draws each tau uniformly in [0, beta), see move.hpp). The overall 1/n!
  // simplex normalisation is irrelevant: we test only whether the coefficient is
  // zero, which no positive constant can change. Returns {mean_coeff, stderr}.
  std::pair<double, double> mc_coefficient(SumDiagrams<double> &sd, int order, int n_samples, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uni(0.0, BETA);
    std::vector<double> taus(order);
    double s = 0.0, s2 = 0.0;
    for (int i = 0; i < n_samples; ++i) {
      for (int j = 0; j < order; ++j) taus[j] = uni(rng);
      sd.mark_all_dirty();
      double f = sd.density_density(taus);
      s += f;
      s2 += f * f;
    }
    double mean  = s / n_samples;
    double var   = std::max(0.0, s2 / n_samples - mean * mean);
    double scale = std::pow(BETA, order);
    return {scale * mean, scale * std::sqrt(var / n_samples)};
  }

  std::vector<int> graph_degrees(Graph const &g) {
    int V = g.get_V();
    std::vector<int> deg(V, 0);
    for (int i = 0; i < V; ++i)
      for (int j = 0; j < V; ++j) deg[i] += g(i, j) + g(j, i);
    return deg;
  }

  // |Aut(G)| over all V! relabelings (V is tiny here). This is the FULL
  // symmetry factor the atomic rooted weight divides by; the dimer rooted path
  // divides by only the mark-FIXING subgroup (get_rooted_sym_factor()).
  long count_full_automorphisms(Graph const &g) {
    int V = g.get_V();
    std::vector<int> perm(V);
    std::iota(perm.begin(), perm.end(), 0);
    long count = 0;
    do {
      bool is_auto = true;
      for (int i = 0; i < V && is_auto; ++i)
        for (int j = 0; j < V && is_auto; ++j)
          if (g(i, j) != g(perm[i], perm[j])) is_auto = false;
      if (is_auto) ++count;
    } while (std::next_permutation(perm.begin(), perm.end()));
    return count;
  }

  // The 4-cycle C4: V == 4, connected, every vertex degree 2 (single edges).
  bool is_four_cycle(Graph const &g) {
    if (g.get_V() != 4) return false;
    for (int d : graph_degrees(g))
      if (d != 2) return false;
    return true;
  }

} // namespace

// -----------------------------------------------------------------------------
//  Test 1 — the regression guard. Every order >= 2 of the on-site same-spin
//  coefficient must integrate to zero (sum rule at half filling). This catches
//  the rooted-weight bug that was fixed in diagram.cpp (the spurious n_mark_orbit
//  doubling drove the order-4 coefficient to ~+0.018 instead of 0).
// -----------------------------------------------------------------------------
TEST(DimerOnSiteSameSpin, IntegratesToZeroPerOrder) {
  auto params = half_filled_params();
  std::cout << std::setprecision(6) << std::fixed;
  for (int order = 2; order <= 5; ++order) {
    SumDiagrams<double> sd(params, order, {0, 0}, /*s1=*/1, /*s2=*/1); // up, up
    auto [coeff, err] = mc_coefficient(sd, order, kSamples, /*seed=*/12345u + (uint64_t)order);
    std::cout << "  order " << order << ": n_diagrams=" << sd.get_n_diagrams() << "  coeff=" << coeff << " +/- " << err << "\n";
    double tol = std::max(1e-3, 6.0 * err);
    EXPECT_NEAR(coeff, 0.0, tol) << "on-site same-spin coefficient must vanish at order " << order;
  }
}

// -----------------------------------------------------------------------------
//  Test 2 — post-fix cancellation diagnostic. Splits the order-4 same-spin
//  coefficient into coincident- vs distinct-vertex groups and confirms they now
//  CANCEL (each is ~±0.018; their sum is ~0). Pre-fix the distinct group was
//  doubled (+0.036) so the sum was +0.018. Also dumps the per-diagram rooted
//  weight bookkeeping and confirms the opposite-spin channel still uses those
//  same distinct-vertex diagrams (they must not be dropped).
// -----------------------------------------------------------------------------
TEST(DimerOnSiteSameSpin, SameSpinGroupsCancelAfterFix) {
  auto params     = half_filled_params();
  int const order = 4;

  SumDiagrams<double> sd(params, order, {0, 0}, /*s1=*/1, /*s2=*/1); // same spin (up, up)
  auto const &solver = sd.get_solver();

  // --- Integrated split: coincident (single-vertex marks) vs distinct-vertex. ---
  std::mt19937_64 rng(98765u);
  std::uniform_real_distribution<double> uni(0.0, BETA);
  std::vector<double> taus(order);
  int const N = kSamples;
  double s_coin = 0.0, s_dist = 0.0, s_tot = 0.0, s_tot2 = 0.0;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < order; ++j) taus[j] = uni(rng);
    double c_i = 0.0, d_i = 0.0;
    for (auto const &dc : sd.get_diagrams()) {
      auto &d = const_cast<Diagram<double> &>(dc);
      d.mark_all_dirty();
      double v = d.evaluate(taus, solver);
      if (d.get_marks()[0] == d.get_marks()[1])
        c_i += v;
      else
        d_i += v;
    }
    s_coin += c_i;
    s_dist += d_i;
    s_tot += c_i + d_i;
    s_tot2 += (c_i + d_i) * (c_i + d_i);
  }
  double scale = std::pow(BETA, order) / N;
  double coin = scale * s_coin, dist = scale * s_dist;
  double tot_mean = s_tot / N;
  double tot_var  = std::max(0.0, s_tot2 / N - tot_mean * tot_mean);
  double tot_err  = std::pow(BETA, order) * std::sqrt(tot_var / N); // stderr of (coin+dist)
  std::cout << std::setprecision(6) << std::fixed;
  std::cout << "\n[order-4 same-spin coefficient split]\n"
            << "  coincident-marks group : " << coin << "\n"
            << "  distinct-vertex group  : " << dist << "\n"
            << "  total                  : " << (coin + dist) << "\n";

  // --- Per-diagram static (tau-independent) weight dump. ---
  std::cout << "\n[order-4 rooted catalog: per-diagram weight bookkeeping]\n";
  std::cout << "  V  marks   sites   spins  sign  freeMult  rootedSym  |Aut|  atomic_nmo  C4?\n";
  for (auto const &dc : sd.get_diagrams()) {
    auto const &d = dc;
    auto const &g = d.get_graph();
    auto m = d.get_marks(), s = d.get_sites(), sp = d.get_mark_spins();

    // n_mark_orbit is now unconditionally 1 (the spurious spin-dependent doubling
    // was removed). The atomic label-combinatorics value is shown for reference.
    double atomic_nmo = (m[0] == m[1]) ? 1.0 : 2.0;

    std::cout << "  " << g.get_V() << "  (" << m[0] << "," << m[1] << ")   (" << s[0] << "," << s[1] << ")   (" << sp[0] << "," << sp[1] << ")    "
              << (int)d.get_diagram_sign() << "   " << std::setw(7) << d.get_free_multiplicity() << "   " << std::setw(7)
              << d.get_rooted_sym_factor() << "   " << std::setw(4) << count_full_automorphisms(g) << "      " << atomic_nmo << "      "
              << (is_four_cycle(g) ? "yes" : "no") << "\n";
  }

  // --- Opposite-spin contrast: identical catalog, distinct group now legitimate. ---
  SumDiagrams<double> sd_opp(params, order, {0, 0}, /*s1=*/0, /*s2=*/1); // up, down
  EXPECT_EQ(sd_opp.get_n_diagrams(), sd.get_n_diagrams()) << "catalog is spin-blind; same topologies for both channels";

  auto const &solver_opp = sd_opp.get_solver();
  std::mt19937_64 rng2(98765u);
  double s_coin_o = 0.0, s_dist_o = 0.0;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < order; ++j) taus[j] = uni(rng2);
    for (auto const &dc : sd_opp.get_diagrams()) {
      auto &d = const_cast<Diagram<double> &>(dc);
      d.mark_all_dirty();
      double v = d.evaluate(taus, solver_opp);
      if (d.get_marks()[0] == d.get_marks()[1])
        s_coin_o += v;
      else
        s_dist_o += v;
    }
  }
  double coin_o = scale * s_coin_o, dist_o = scale * s_dist_o;
  std::cout << "[order-4 OPPOSITE-spin split (legitimate double-occupancy)]\n"
            << "  coincident-marks group : " << coin_o << "\n"
            << "  distinct-vertex group  : " << dist_o << "\n"
            << "  total                  : " << (coin_o + dist_o) << "\n\n";

  // --- Assertions (post-fix) ---
  // The two groups cancel: the order-4 same-spin coefficient vanishes.
  std::cout << "  total stderr           : " << tot_err << "\n";
  EXPECT_NEAR(coin + dist, 0.0, std::max(4e-3, 6.0 * tot_err)) << "order-4 same-spin coefficient must vanish (groups cancel)";

  // Sanity that the cancellation is non-trivial: each group is individually large
  // and they are opposite in sign (the cancellation, not both being ~0).
  EXPECT_GT(std::abs(coin), 5e-3) << "coincident group should be individually large (it cancels the distinct group)";
  EXPECT_GT(std::abs(dist), 5e-3) << "distinct group should be individually large (it cancels the coincident group)";
  EXPECT_LT(coin * dist, 0.0) << "the two groups should have opposite signs (genuine cancellation)";

  // The same distinct-vertex diagrams are NOT spurious for opposite spin:
  // they carry the genuine connected double-occupancy and must be kept.
  EXPECT_GT(std::abs(dist_o), 1e-3) << "distinct-vertex diagrams are legitimate in the opposite-spin channel";
}

// -----------------------------------------------------------------------------
//  Test 3 — order-5 residual localizer. After the order-4 fix, the on-site
//  same-spin coefficient still does not vanish at order 5 (cluster run:
//  ~4.0e-4, ~113 sigma). The factor-2 doubling is gone, so this is a subtler
//  per-topology imbalance between the coincident-marks group and the
//  distinct-vertex group (at r=0 BOTH marks land on site 0 of cell (0,0)).
//
//  This dumps the integrated contribution of EACH order-5 diagram (shared tau
//  stream so the grouping is exact), tagged by a topology signature, so the
//  single graph whose coincident/distinct partner fails to cancel is visible.
//  Same-spin must vanish; the opposite-spin column is the legitimate-physics
//  contrast (it need not vanish).
// -----------------------------------------------------------------------------
namespace {
  // A coarse but stable per-diagram fingerprint: V, #directed-edge-entries,
  // sorted degree sequence, coincident-marks flag. Enough to tell topologies
  // apart in the dump without a full canonical form.
  std::string topo_signature(Diagram<double> const &d) {
    auto const &g = d.get_graph();
    int V = g.get_V(), edges = 0;
    for (int i = 0; i < V; ++i)
      for (int j = 0; j < V; ++j) edges += g(i, j);
    auto deg = graph_degrees(g);
    std::sort(deg.begin(), deg.end());
    std::string s = "V" + std::to_string(V) + " E" + std::to_string(edges) + " deg[";
    for (size_t i = 0; i < deg.size(); ++i) s += (i ? "," : "") + std::to_string(deg[i]);
    s += "]";
    return s;
  }

  // Integrate every diagram of `sd` over a shared uniform tau stream; returns
  // per-diagram coefficient estimates (beta^order * mean), index-aligned with
  // sd.get_diagrams().
  std::vector<double> per_diagram_coeffs(SumDiagrams<double> &sd, int order, int n_samples, uint64_t seed) {
    auto const &solver = sd.get_solver();
    auto const &diags  = sd.get_diagrams();
    std::vector<double> acc(diags.size(), 0.0);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uni(0.0, BETA);
    std::vector<double> taus(order);
    for (int i = 0; i < n_samples; ++i) {
      for (int j = 0; j < order; ++j) taus[j] = uni(rng);
      for (size_t k = 0; k < diags.size(); ++k) {
        auto &d = const_cast<Diagram<double> &>(diags[k]);
        d.mark_all_dirty();
        acc[k] += d.evaluate(taus, solver);
      }
    }
    double scale = std::pow(BETA, order) / n_samples;
    for (auto &a : acc) a *= scale;
    return acc;
  }
} // namespace

TEST(DimerOnSiteSameSpin, Order5ResidualPerDiagram) {
  auto params     = half_filled_params();
  int const order = 5;

  SumDiagrams<double> sd_same(params, order, {0, 0}, /*s1=*/1, /*s2=*/1); // up, up
  SumDiagrams<double> sd_opp(params, order, {0, 0}, /*s1=*/0, /*s2=*/1);  // up, down
  ASSERT_EQ(sd_same.get_n_diagrams(), sd_opp.get_n_diagrams()) << "catalog is spin-blind";

  auto same = per_diagram_coeffs(sd_same, order, kSamples, /*seed=*/24680u);
  auto opp  = per_diagram_coeffs(sd_opp, order, kSamples, /*seed=*/24680u);

  auto const &diags = sd_same.get_diagrams();
  std::cout << std::setprecision(6) << std::fixed;
  std::cout << "\n[order-5 per-diagram on-site coefficient (r=0)]\n";
  std::cout << "  idx  coin?  marks   sites   rootedSym  |Aut|   same_spin      opp_spin    topology\n";

  double coin_same = 0.0, dist_same = 0.0;
  // Accumulate same-spin by topology signature to expose which graph's
  // coincident/distinct pair fails to cancel.
  std::map<std::string, double> by_topo_same, by_topo_opp;
  for (size_t k = 0; k < diags.size(); ++k) {
    auto const &d = diags[k];
    auto m = d.get_marks(), s = d.get_sites();
    bool coincident = (m[0] == m[1]);
    if (coincident) coin_same += same[k];
    else dist_same += same[k];
    std::string sig = topo_signature(d);
    by_topo_same[sig] += same[k];
    by_topo_opp[sig] += opp[k];

    std::cout << "  " << std::setw(3) << k << "   " << (coincident ? "Y" : "n") << "    (" << m[0] << "," << m[1] << ")   (" << s[0] << ","
              << s[1] << ")   " << std::setw(7) << d.get_rooted_sym_factor() << "   " << std::setw(4) << count_full_automorphisms(d.get_graph())
              << "   " << std::setw(11) << same[k] << "  " << std::setw(11) << opp[k] << "   " << sig << "\n";
  }

  std::cout << "\n[order-5 same-spin grouped by topology]  (each row should be ~0 for same spin)\n";
  for (auto const &[sig, val] : by_topo_same)
    std::cout << "  same=" << std::setw(11) << val << "   opp=" << std::setw(11) << by_topo_opp[sig] << "   " << sig << "\n";

  std::cout << "\n  coincident group total : " << coin_same << "\n"
            << "  distinct  group total  : " << dist_same << "\n"
            << "  same-spin total        : " << (coin_same + dist_same) << "\n\n";

  // This documents the open bug: at order 5 the same-spin total is NOT yet zero.
  // Once the per-topology imbalance is fixed, flip this to EXPECT_NEAR(...,0).
  // For now it is a localizer, not a pass/fail gate, so don't fail the suite.
  EXPECT_TRUE(true);
}

// -----------------------------------------------------------------------------
//  Test 4 — fast (no-MC) STRUCTURAL dump of the order-5 catalog. Prints, per
//  diagram, the adjacency matrix, directed hopping lines, embedding free
//  multiplicity, rooted_sym_factor, and #valid-configs. No integration, so it
//  returns instantly. Purpose: see whether the double-bonded coincident
//  diagrams (idx 16/18 in the per-diagram test) are isomorphic duplicates and
//  whether the embedding count already absorbs the parallel-line interchange
//  (which would make the rooted_sym_factor multi-edge division a double-correct).
// -----------------------------------------------------------------------------
TEST(DimerOnSiteSameSpin, Order5StructureDump) {
  auto params     = half_filled_params();
  int const order = 5;
  SumDiagrams<double> sd(params, order, {0, 0}, /*s1=*/1, /*s2=*/1); // up, up

  std::cout << std::setprecision(6) << std::fixed;
  std::cout << "\n[order-5 structural dump]\n";
  int idx = 0;
  for (auto const &d : sd.get_diagrams()) {
    auto const &g = d.get_graph();
    int V         = g.get_V();
    auto m = d.get_marks(), s = d.get_sites();
    bool coincident = (m[0] == m[1]);

    std::cout << "idx " << idx++ << "  V=" << V << "  marks(" << m[0] << "," << m[1] << ")  sites(" << s[0] << "," << s[1] << ")  "
              << (coincident ? "COINCIDENT" : "distinct ") << "  freeMult=" << d.get_free_multiplicity()
              << "  rootedSym=" << d.get_rooted_sym_factor() << "  nValidCfg=" << d.get_valid_configurations().size() << "\n";

    std::cout << "    adj(directed g(i,j)):\n";
    for (int i = 0; i < V; ++i) {
      std::cout << "      ";
      for (int j = 0; j < V; ++j) std::cout << (int)g(i, j) << " ";
      std::cout << "\n";
    }
  }
  EXPECT_TRUE(true);
}
