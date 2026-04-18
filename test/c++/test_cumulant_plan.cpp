// Parity test: CumulantSolver::record_plan / evaluate_plan must reproduce
// CumulantSolver::compute_cumulant_decomposition numerically for all τ.
//
// Covers:
//   - atom (N_sites=1) and dimer (N_sites=2)
//   - vertex orders 1, 2, 3, 4
//   - finite-U and infinite-U paths (G0n vs G0n_infinite_U)
//   - random τ sampling (many draws per config)
//   - plan reuse: build plan at τ₁, evaluate at τ₂ ≠ τ₁ (proves τ-independence)

#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <vector>
#include <memory>
#include "../c++/sc_expansion/hubbard_solver.hpp"
#include "../c++/sc_expansion/cumulant.hpp"
#include "../c++/sc_expansion/cumulant_plan.hpp"

using namespace sc_expansion;

namespace {

  constexpr double kTol = 1e-12;

  // Build FermionOperator with N_sites-aware creation/destruction bit.
  template <int N_sites, typename T>
  FermionOperator<N_sites, T> make_op(int orbital, bool creation) {
    uint8_t id = (uint8_t)orbital;
    if (creation) id |= FermionOperator<N_sites, T>::ACTION_BIT;
    return FermionOperator<N_sites, T>(id);
  }

  // Build (unprimed = destructions, primed = creations) Args for a given
  // list of (τ, orbital) pairs. The convention matches CumulantSolver:
  //   unprimed[i] is a destruction c_{orb_u[i]} at τ_u[i]
  //   primed[i]   is a creation  c†_{orb_p[i]} at τ_p[i]
  template <int N_sites, typename T>
  std::pair<Args<N_sites, T>, Args<N_sites, T>>
  make_args(std::vector<double> const &taus_u, std::vector<int> const &orbs_u,
            std::vector<double> const &taus_p, std::vector<int> const &orbs_p) {
    std::vector<FermionOperator<N_sites, T>> ops_u, ops_p;
    ops_u.reserve(orbs_u.size());
    ops_p.reserve(orbs_p.size());
    for (int o : orbs_u) ops_u.push_back(make_op<N_sites, T>(o, /*creation=*/false));
    for (int o : orbs_p) ops_p.push_back(make_op<N_sites, T>(o, /*creation=*/true));
    return {Args<N_sites, T>(taus_u, std::move(ops_u)),
            Args<N_sites, T>(taus_p, std::move(ops_p))};
  }

  // Compare CumulantSolver::compute_cumulant_decomposition against evaluate_plan
  // for a single (unprimed, primed) configuration. Plan is built from this same
  // configuration. Returns (reference, plan_value).
  template <int N_sites>
  std::pair<double, double>
  both_paths(HubbardSolver<N_sites, double> const &solver,
             Args<N_sites, double> const &unprimed,
             Args<N_sites, double> const &primed,
             bool infinite_U) {
    CumulantSolver<N_sites, double> ref_solver(unprimed, primed, solver, infinite_U);
    double ref = ref_solver.compute_cumulant_decomposition();

    CumulantSolver<N_sites, double> rec_solver(unprimed, primed, solver, infinite_U);
    CumulantPlan plan;
    rec_solver.record_plan(plan);
    double from_plan = evaluate_plan(plan, unprimed, primed, solver, infinite_U);

    return {ref, from_plan};
  }

  // Parity across many random τ draws.
  // orbs_u / orbs_p fix the operator structure; τ values are resampled each draw.
  template <int N_sites>
  void parity_sweep(HubbardSolver<N_sites, double> const &solver,
                    std::vector<int> const &orbs_u, std::vector<int> const &orbs_p,
                    bool infinite_U, double beta, int n_draws, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> U01(0.0, 1.0);

    for (int s = 0; s < n_draws; ++s) {
      std::vector<double> tu(orbs_u.size()), tp(orbs_p.size());
      for (auto &x : tu) x = beta * U01(rng);
      for (auto &x : tp) x = beta * U01(rng);

      auto [u, p] = make_args<N_sites, double>(tu, orbs_u, tp, orbs_p);
      auto [ref, via_plan] = both_paths<N_sites>(solver, u, p, infinite_U);

      EXPECT_NEAR(ref, via_plan, kTol)
          << "N_sites=" << N_sites << " order=" << orbs_u.size()
          << " infU=" << infinite_U << " draw=" << s
          << " ref=" << ref << " plan=" << via_plan
          << " diff=" << (ref - via_plan);
    }
  }

  // Plan reuse: build plan at τ₁, evaluate with τ₂, compare to a fresh solver at τ₂.
  // This is the real proof that CumulantPlan is τ-independent.
  template <int N_sites>
  void plan_reuse_check(HubbardSolver<N_sites, double> const &solver,
                        std::vector<int> const &orbs_u, std::vector<int> const &orbs_p,
                        bool infinite_U, double beta, int n_pairs, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> U01(0.0, 1.0);

    auto draw_taus = [&](int n) {
      std::vector<double> t(n);
      for (auto &x : t) x = beta * U01(rng);
      return t;
    };

    for (int s = 0; s < n_pairs; ++s) {
      // Build plan from τ₁
      auto tu1 = draw_taus((int)orbs_u.size());
      auto tp1 = draw_taus((int)orbs_p.size());
      auto [u1, p1] = make_args<N_sites, double>(tu1, orbs_u, tp1, orbs_p);
      CumulantSolver<N_sites, double> builder(u1, p1, solver, infinite_U);
      CumulantPlan plan;
      builder.record_plan(plan);

      // Evaluate plan at τ₂ (different values, same operators)
      auto tu2 = draw_taus((int)orbs_u.size());
      auto tp2 = draw_taus((int)orbs_p.size());
      auto [u2, p2] = make_args<N_sites, double>(tu2, orbs_u, tp2, orbs_p);
      double via_plan = evaluate_plan(plan, u2, p2, solver, infinite_U);

      // Reference: fresh CumulantSolver at τ₂
      CumulantSolver<N_sites, double> ref_solver(u2, p2, solver, infinite_U);
      double ref = ref_solver.compute_cumulant_decomposition();

      EXPECT_NEAR(ref, via_plan, kTol)
          << "[plan-reuse] N_sites=" << N_sites << " order=" << orbs_u.size()
          << " infU=" << infinite_U << " pair=" << s
          << " ref=" << ref << " plan=" << via_plan;
    }
  }

} // namespace

// ============================================================================
//  Atom (N_sites = 1)
// ============================================================================

class AtomCumulantPlanTest : public ::testing::Test {
  protected:
  Parameters<double> params{4.0, 1.5, 2.0, 0.0, true};
  std::unique_ptr<HubbardSolver<1, double>> solver;
  void SetUp() override { solver = std::make_unique<HubbardSolver<1, double>>(params); }
};

TEST_F(AtomCumulantPlanTest, Order1_FiniteU) {
  // orbital 0 = down, 1 = up. One-particle spin-up cumulant.
  parity_sweep<1>(*solver, {1}, {1}, false, params.beta, 64, 0xA1u);
}
TEST_F(AtomCumulantPlanTest, Order2_FiniteU_SameSpin) {
  parity_sweep<1>(*solver, {1, 1}, {1, 1}, false, params.beta, 64, 0xA2u);
}
TEST_F(AtomCumulantPlanTest, Order2_FiniteU_OpposingSpin) {
  parity_sweep<1>(*solver, {0, 1}, {0, 1}, false, params.beta, 64, 0xA3u);
}
TEST_F(AtomCumulantPlanTest, Order3_FiniteU) {
  parity_sweep<1>(*solver, {0, 1, 1}, {0, 1, 1}, false, params.beta, 48, 0xA4u);
}
TEST_F(AtomCumulantPlanTest, Order4_FiniteU) {
  parity_sweep<1>(*solver, {0, 0, 1, 1}, {0, 0, 1, 1}, false, params.beta, 32, 0xA5u);
}
TEST_F(AtomCumulantPlanTest, Order2_InfiniteU) {
  parity_sweep<1>(*solver, {0, 1}, {0, 1}, true, params.beta, 48, 0xA6u);
}
TEST_F(AtomCumulantPlanTest, Order3_InfiniteU) {
  parity_sweep<1>(*solver, {0, 1, 1}, {0, 1, 1}, true, params.beta, 32, 0xA7u);
}
TEST_F(AtomCumulantPlanTest, PlanIsTauIndependent_Order2) {
  plan_reuse_check<1>(*solver, {0, 1}, {0, 1}, false, params.beta, 48, 0xA80u);
}
TEST_F(AtomCumulantPlanTest, PlanIsTauIndependent_Order3) {
  plan_reuse_check<1>(*solver, {0, 1, 1}, {0, 1, 1}, false, params.beta, 32, 0xA8u);
}
TEST_F(AtomCumulantPlanTest, PlanIsTauIndependent_Order4) {
  plan_reuse_check<1>(*solver, {0, 0, 1, 1}, {0, 0, 1, 1}, false, params.beta, 16, 0xA9u);
}

// ============================================================================
//  Dimer (N_sites = 2)
// ============================================================================

class DimerCumulantPlanTest : public ::testing::Test {
  protected:
  // orbital encoding (N_sites=2): 0=site0 down, 1=site1 down, 2=site0 up, 3=site1 up
  Parameters<double> params{4.0, 1.0, 2.0, 1.0, true};
  std::unique_ptr<HubbardSolver<2, double>> solver;
  void SetUp() override { solver = std::make_unique<HubbardSolver<2, double>>(params); }
};

TEST_F(DimerCumulantPlanTest, Order1_FiniteU_Site0Up) {
  parity_sweep<2>(*solver, {2}, {2}, false, params.beta, 64, 0xD1u);
}
TEST_F(DimerCumulantPlanTest, Order2_FiniteU_OnSite) {
  parity_sweep<2>(*solver, {2, 0}, {2, 0}, false, params.beta, 64, 0xD2u);
}
TEST_F(DimerCumulantPlanTest, Order2_FiniteU_InterSite) {
  parity_sweep<2>(*solver, {2, 1}, {2, 1}, false, params.beta, 64, 0xD3u);
}
TEST_F(DimerCumulantPlanTest, Order2_FiniteU_Mixed) {
  parity_sweep<2>(*solver, {2, 1}, {0, 3}, false, params.beta, 64, 0xD4u);
}
TEST_F(DimerCumulantPlanTest, Order3_FiniteU) {
  parity_sweep<2>(*solver, {2, 0, 3}, {2, 0, 3}, false, params.beta, 48, 0xD5u);
}
TEST_F(DimerCumulantPlanTest, Order4_FiniteU) {
  parity_sweep<2>(*solver, {2, 0, 3, 1}, {2, 0, 3, 1}, false, params.beta, 24, 0xD6u);
}
TEST_F(DimerCumulantPlanTest, Order2_InfiniteU) {
  parity_sweep<2>(*solver, {2, 0}, {2, 0}, true, params.beta, 48, 0xD7u);
}
TEST_F(DimerCumulantPlanTest, Order3_InfiniteU) {
  parity_sweep<2>(*solver, {2, 0, 3}, {2, 0, 3}, true, params.beta, 32, 0xD8u);
}
TEST_F(DimerCumulantPlanTest, PlanIsTauIndependent_Order2) {
  plan_reuse_check<2>(*solver, {2, 0}, {2, 0}, false, params.beta, 48, 0xD90u);
}
TEST_F(DimerCumulantPlanTest, PlanIsTauIndependent_Order3) {
  plan_reuse_check<2>(*solver, {2, 0, 3}, {2, 0, 3}, false, params.beta, 32, 0xD9u);
}
TEST_F(DimerCumulantPlanTest, PlanIsTauIndependent_Order4) {
  plan_reuse_check<2>(*solver, {2, 0, 3, 1}, {2, 0, 3, 1}, false, params.beta, 12, 0xDAu);
}

// ============================================================================
//  Permuted operator order: catches leaf-ordering bugs where G0n's internal
//  τ-sign would fall out of step with the combinatorial sign baked into the plan.
// ============================================================================

TEST_F(DimerCumulantPlanTest, Order3_PermutedOpOrder) {
  // Same operator multiset, different index order → different plan layout but same
  // physical cumulant. Both must match compute_cumulant_decomposition.
  parity_sweep<2>(*solver, {3, 2, 0}, {0, 2, 3}, false, params.beta, 32, 0xDEu);
}
