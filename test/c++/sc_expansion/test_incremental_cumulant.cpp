// Cached vs fresh cumulant evaluation.
//
// evaluate_plan_incremental is the MC hot-path evaluator: it persists per-node
// cumulant values across calls and recomputes a node only when its (transitive)
// line-dependence mask intersects the changed taus. This test drives it over a
// random single-/double-line tau walk -- with occasional in-place reverts that
// mirror an MC reject -- and asserts that at EVERY step the cached value equals
// a fresh full evaluate_plan at the current taus. A clean (reused) node must be
// bit-identical to recomputing it, so any cache bug shows up as an O(1) mismatch.
//
// Covered for both expansions (atom N_sites=1, dimer N_sites=2), both scalar
// types (double, and Dual = the mu-derivative / density observable), finite and
// (atom-only) infinite U, and cumulant orders up to the order-10 "watermelon"
// vertex (5 creators + 5 annihilators) that dominates the dimer free energy.
#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <vector>
#include <memory>
#include <cstdint>
#include "../c++/sc_expansion/hubbard_solver.hpp"
#include "../c++/sc_expansion/cumulant.hpp"

using namespace sc_expansion;

namespace {

  // A cache bug reuses a stale node -> the cumulant changes by an O(1) amount,
  // far above this tolerance; bit-identical reuse passes trivially.
  constexpr double kTol = 1e-12;

  template <int N_sites, typename T> FermionOperator<N_sites, T> make_op(int orbital, bool creation) {
    uint8_t id = (uint8_t)orbital;
    if (creation) id |= FermionOperator<N_sites, T>::ACTION_BIT;
    return FermionOperator<N_sites, T>(id);
  }

  // unprimed = destructions c_{orbs_u[i]} at taus_u[i]; primed = creations
  // c^dag_{orbs_p[i]} at taus_p[i]. Matches CumulantSolver / Args::split_from_raw.
  template <int N_sites, typename T>
  std::pair<Args<N_sites, T>, Args<N_sites, T>> make_args(std::vector<double> const &taus_u, std::vector<int> const &orbs_u,
                                                          std::vector<double> const &taus_p, std::vector<int> const &orbs_p) {
    std::vector<FermionOperator<N_sites, T>> ops_u, ops_p;
    ops_u.reserve(orbs_u.size());
    ops_p.reserve(orbs_p.size());
    for (int o : orbs_u) ops_u.push_back(make_op<N_sites, T>(o, /*creation=*/false));
    for (int o : orbs_p) ops_p.push_back(make_op<N_sites, T>(o, /*creation=*/true));
    return {Args<N_sites, T>(taus_u, std::move(ops_u)), Args<N_sites, T>(taus_p, std::move(ops_p))};
  }

  inline double comp_diff(double a, double b) { return std::abs(a - b); }
  // Dual carries (value, mu-derivative); both components must match.
  inline double comp_diff(Dual a, Dual b) { return std::max(std::abs(a.value - b.value), std::abs(a.derivative - b.derivative)); }

  // Per-node transitive line-dependence mask, re-derived here exactly as
  // Diagram::build_node_cache_metadata does it (independent cross-check of that
  // logic): each line is one operator -- unprimed stable index j -> line j,
  // primed stable index j -> line (k_u + j). A node depends on a line if any of
  // its leaf operators sit on it, or any descendant sub-cumulant does (children
  // precede parents). Mark / <n_sigma> nodes have no leaf operators -> no line dep.
  std::vector<uint64_t> build_node_line_mask(CumulantPlan const &plan, int k_u) {
    std::vector<uint64_t> mask(plan.nodes.size(), 0);
    std::vector<char> is_mark(plan.nodes.size(), 0);
    for (int mid : plan.mark_node_ids)
      if (mid >= 0 && mid < (int)mask.size()) is_mark[mid] = 1;

    for (size_t i = 0; i < plan.nodes.size(); ++i) {
      uint64_t m = 0;
      if (!is_mark[i]) {
        for (int su : plan.nodes[i].leaf.u_global_idx) m |= (uint64_t(1) << su);
        for (int sp : plan.nodes[i].leaf.p_global_idx) m |= (uint64_t(1) << (k_u + sp));
      }
      for (auto const &term : plan.nodes[i].subtraction_terms)
        for (int fid : term.factor_node_ids) m |= mask[fid];
      mask[i] = m;
    }
    return mask;
  }

  template <int N_sites, typename T>
  void incremental_walk(HubbardSolver<N_sites, T> const &solver, std::vector<int> const &orbs_u, std::vector<int> const &orbs_p, bool infinite_U,
                        double beta, int n_steps, uint64_t seed) {
    int k_u = (int)orbs_u.size(), k_p = (int)orbs_p.size(), n_lines = k_u + k_p;

    // tau-independent plan, built once from dummy taus.
    std::vector<double> du(k_u, 0.5), dp(k_p, 0.5);
    auto [u0, p0] = make_args<N_sites, T>(du, orbs_u, dp, orbs_p);
    CumulantSolver<N_sites, T> builder(u0, p0, solver, infinite_U);
    CumulantPlan plan;
    builder.record_plan(plan);
    auto mask = build_node_line_mask(plan, k_u);

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    auto draw = [&] { return beta * U01(rng); };

    std::vector<double> tu(k_u), tp(k_p);
    for (auto &x : tu) x = draw();
    for (auto &x : tp) x = draw();

    auto cur_args = [&] { return make_args<N_sites, T>(tu, orbs_u, tp, orbs_p); };
    auto set_line = [&](int L, double v) {
      if (L < k_u) tu[L] = v;
      else tp[L - k_u] = v;
    };
    auto get_line = [&](int L) { return L < k_u ? tu[L] : tp[L - k_u]; };

    std::vector<T> cache; // persists across steps

    // First evaluation: full recompute (seeds the cache).
    {
      auto [u, p] = cur_args();
      T inc       = evaluate_plan_incremental(plan, u, p, solver, infinite_U, cache, mask, /*recompute_all=*/true, ~uint64_t(0));
      T ref       = evaluate_plan(plan, u, p, solver, infinite_U);
      EXPECT_LE(comp_diff(inc, ref), kTol) << "[init] N_sites=" << N_sites << " order=" << k_u << " infU=" << infinite_U;
    }

    for (int s = 0; s < n_steps; ++s) {
      uint64_t changed = 0;
      int n_ch         = 1 + (int)(rng() % 2); // 1 or 2 lines per step -> exercises accumulated masks (e.g. across a reject)
      for (int c = 0; c < n_ch; ++c) {
        int L      = (int)(rng() % n_lines);
        double old = get_line(L);
        set_line(L, draw());
        changed |= (uint64_t(1) << L);
        if (rng() % 4 == 0) set_line(L, old); // ~1/4: revert in place (line re-dirtied, tau unchanged) == MC reject
      }
      auto [u, p] = cur_args();
      T inc       = evaluate_plan_incremental(plan, u, p, solver, infinite_U, cache, mask, /*recompute_all=*/false, changed);
      T ref       = evaluate_plan(plan, u, p, solver, infinite_U);
      EXPECT_LE(comp_diff(inc, ref), kTol)
         << "[step " << s << "] N_sites=" << N_sites << " order=" << k_u << " infU=" << infinite_U << " changed=" << changed;
    }
  }

} // namespace

// ============================================================================
//  Atom (N_sites = 1).  Orbital encoding: 0 = down, 1 = up.
// ============================================================================

class AtomIncrementalCacheTest : public ::testing::Test {
  protected:
  Parameters<double> params{12.0, 2.0, 3.0, 0.0, true};            // U, beta, mu, t(unused), bipartite
  Parameters<Dual> params_d{Dual(12.0), Dual(2.0), Dual(3.0, 1.0), // mu seeds the density derivative
                            Dual(0.0), true, Dual(0.0)};
  std::unique_ptr<HubbardSolver<1, double>> solver;
  std::unique_ptr<HubbardSolver<1, Dual>> solver_d;
  void SetUp() override {
    solver   = std::make_unique<HubbardSolver<1, double>>(params);
    solver_d = std::make_unique<HubbardSolver<1, Dual>>(params_d);
  }
};

TEST_F(AtomIncrementalCacheTest, Order2_FiniteU) { incremental_walk<1, double>(*solver, {1, 0}, {1, 0}, false, params.beta, 200, 0xA2u); }
TEST_F(AtomIncrementalCacheTest, Order3_FiniteU) { incremental_walk<1, double>(*solver, {1, 0, 1}, {1, 0, 1}, false, params.beta, 200, 0xA3u); }
TEST_F(AtomIncrementalCacheTest, Order4_FiniteU) { incremental_walk<1, double>(*solver, {1, 0, 1, 0}, {1, 0, 1, 0}, false, params.beta, 150, 0xA4u); }
TEST_F(AtomIncrementalCacheTest, Order3_InfiniteU) { incremental_walk<1, double>(*solver, {1, 0, 1}, {1, 0, 1}, true, params.beta, 150, 0xA5u); }
TEST_F(AtomIncrementalCacheTest, Order3_Dual_Density) { incremental_walk<1, Dual>(*solver_d, {1, 0, 1}, {1, 0, 1}, false, 2.0, 150, 0xA6u); }
TEST_F(AtomIncrementalCacheTest, Order4_Dual_Density) {
  incremental_walk<1, Dual>(*solver_d, {1, 0, 1, 0}, {1, 0, 1, 0}, false, 2.0, 120, 0xA7u);
}

// ============================================================================
//  Dimer (N_sites = 2).  Orbital encoding: 0=site0 down, 1=site1 down, 2=site0 up, 3=site1 up.
// ============================================================================

class DimerIncrementalCacheTest : public ::testing::Test {
  protected:
  Parameters<double> params{12.0, 2.0, 3.0, 1.0, true};
  Parameters<Dual> params_d{Dual(12.0), Dual(2.0), Dual(3.0, 1.0), Dual(1.0), true, Dual(0.0)};
  std::unique_ptr<HubbardSolver<2, double>> solver;
  std::unique_ptr<HubbardSolver<2, Dual>> solver_d;
  void SetUp() override {
    solver   = std::make_unique<HubbardSolver<2, double>>(params);
    solver_d = std::make_unique<HubbardSolver<2, Dual>>(params_d);
  }
};

TEST_F(DimerIncrementalCacheTest, Order2_FiniteU) { incremental_walk<2, double>(*solver, {2, 0}, {2, 0}, false, params.beta, 200, 0xD2u); }
TEST_F(DimerIncrementalCacheTest, Order3_FiniteU) { incremental_walk<2, double>(*solver, {2, 0, 3}, {2, 0, 3}, false, params.beta, 200, 0xD3u); }
TEST_F(DimerIncrementalCacheTest, Order4_FiniteU) {
  incremental_walk<2, double>(*solver, {2, 0, 3, 1}, {2, 0, 3, 1}, false, params.beta, 150, 0xD4u);
}
// Order-5 cumulant = the order-10 "watermelon" vertex (5 c^dag + 5 c) that
// dominates the order-10 dimer free energy -- the whole point of the cache.
TEST_F(DimerIncrementalCacheTest, Order5_FiniteU_Watermelon) {
  incremental_walk<2, double>(*solver, {2, 0, 3, 1, 2}, {2, 0, 3, 1, 2}, false, params.beta, 120, 0xD5u);
}
TEST_F(DimerIncrementalCacheTest, Order4_Dual_Density) {
  incremental_walk<2, Dual>(*solver_d, {2, 0, 3, 1}, {2, 0, 3, 1}, false, 2.0, 120, 0xD6u);
}
TEST_F(DimerIncrementalCacheTest, Order5_Dual_Density_Watermelon) {
  incremental_walk<2, Dual>(*solver_d, {2, 0, 3, 1, 2}, {2, 0, 3, 1, 2}, false, 2.0, 100, 0xD7u);
}
// Permuted operator order: different plan layout, same physical cumulant.
TEST_F(DimerIncrementalCacheTest, Order3_PermutedOpOrder) {
  incremental_walk<2, double>(*solver, {3, 2, 0}, {0, 2, 3}, false, params.beta, 150, 0xDEu);
}
