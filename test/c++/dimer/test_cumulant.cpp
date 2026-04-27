#include <gtest/gtest.h>
#include <array>
#include <random>
#include <vector>
#include <memory>
#include "../c++/sc_expansion/hubbard_solver.hpp"
#include "../c++/sc_expansion/cumulant.hpp"

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
  // list of (τ, orbital) pairs.
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

  // Cumulant value via CumulantSolver (the unit under test).
  template <int N_sites>
  double cumulant_value(HubbardSolver<N_sites, double> const &solver,
                        Args<N_sites, double> const &u, Args<N_sites, double> const &p,
                        bool infinite_U = false) {
    CumulantSolver<N_sites, double> cs(u, p, solver, infinite_U);
    return cs.compute_cumulant_decomposition();
  }

  // Build a combined Args for solver.G0n() with interleaved (c†_p[i], c_u[i]) pairs.
  // This matches the operator layout that CumulantSolver constructs internally for
  // its leaf G0n calls, so the cumulant identities below hold without extra signs.
  template <int N_sites>
  Args<N_sites, double>
  combined_for_g0n(std::vector<std::pair<double, int>> const &creations,
                   std::vector<std::pair<double, int>> const &destructions) {
    int n = (int)creations.size();
    std::vector<double> taus;
    std::vector<FermionOperator<N_sites, double>> ops;
    taus.reserve(2 * n);
    ops.reserve(2 * n);
    for (int i = 0; i < n; ++i) {
      taus.push_back(creations[i].first);
      ops.push_back(make_op<N_sites, double>(creations[i].second, /*creation=*/true));
      taus.push_back(destructions[i].first);
      ops.push_back(make_op<N_sites, double>(destructions[i].second, /*creation=*/false));
    }
    return Args<N_sites, double>(std::move(taus), std::move(ops));
  }

  // Single-particle Green's function G(u | p) = G0n on a 2-operator combined Args.
  template <int N_sites>
  double Gpair(HubbardSolver<N_sites, double> const &solver,
               double tu, int orb_u, double tp, int orb_p) {
    return solver.G0n(combined_for_g0n<N_sites>({{tp, orb_p}}, {{tu, orb_u}}));
  }

} // namespace

// ============================================================================
//  Cumulant identity tests: verify the connected-cumulant decomposition matches
//  its analytic combinatorial definition in terms of solver.G0n calls.
// ============================================================================

class DimerCumulantTest : public ::testing::Test {
  protected:
  // orbital encoding (N_sites=2): 0=site0 down, 1=site1 down, 2=site0 up, 3=site1 up
  Parameters<double> params{4.0, 1.0, 2.0, 1.0, true};
  std::unique_ptr<HubbardSolver<2, double>> solver;
  void SetUp() override { solver = std::make_unique<HubbardSolver<2, double>>(params); }
};

TEST_F(DimerCumulantTest, Order1_EqualsG0n) {
  // C₁(u | p) = G(u | p) — order 1 has no lower-order subtractions.
  std::mt19937_64 rng(0x1010u);
  std::uniform_real_distribution<double> U01(0.0, 1.0);
  for (int draw = 0; draw < 64; ++draw) {
    double tu = params.beta * U01(rng);
    double tp = params.beta * U01(rng);
    int ou    = draw % 4;
    int op    = (draw / 4) % 4;

    auto [u, p] = make_args<2, double>({tu}, {ou}, {tp}, {op});
    double c1   = cumulant_value<2>(*solver, u, p);
    double g1   = Gpair<2>(*solver, tu, ou, tp, op);
    EXPECT_NEAR(c1, g1, kTol)
        << "draw=" << draw << " ou=" << ou << " op=" << op
        << " τu=" << tu << " τp=" << tp << " c1=" << c1 << " g1=" << g1;
  }
}

TEST_F(DimerCumulantTest, Order2_EqualsG02MinusDisconnected) {
  // C₂(1,2 | 1',2') = G02(1,2 | 1',2') − G(1|1') G(2|2') + G(1|2') G(2|1').
  std::mt19937_64 rng(0x2020u);
  std::uniform_real_distribution<double> U01(0.0, 1.0);
  std::vector<std::array<int, 4>> orb_combos = {
      {2, 0, 2, 0}, {2, 1, 2, 1}, {2, 1, 0, 3}, {3, 0, 1, 2}, {2, 3, 0, 1},
  };
  for (auto const &orbs : orb_combos) {
    int ou0 = orbs[0], ou1 = orbs[1], op0 = orbs[2], op1 = orbs[3];
    for (int draw = 0; draw < 16; ++draw) {
      double tu0 = params.beta * U01(rng), tu1 = params.beta * U01(rng);
      double tp0 = params.beta * U01(rng), tp1 = params.beta * U01(rng);

      auto [u, p] = make_args<2, double>({tu0, tu1}, {ou0, ou1}, {tp0, tp1}, {op0, op1});
      double c2   = cumulant_value<2>(*solver, u, p);

      auto g02_args = combined_for_g0n<2>({{tp0, op0}, {tp1, op1}}, {{tu0, ou0}, {tu1, ou1}});
      double G02    = solver->G0n(g02_args);

      double g_00 = Gpair<2>(*solver, tu0, ou0, tp0, op0);
      double g_11 = Gpair<2>(*solver, tu1, ou1, tp1, op1);
      double g_01 = Gpair<2>(*solver, tu0, ou0, tp1, op1);
      double g_10 = Gpair<2>(*solver, tu1, ou1, tp0, op0);

      // CumulantSolver works in the stable (sorted-by-τ) indexing internally and
      // applies (u.perm_sign * p.perm_sign) at the root to convert back to the
      // user-input ordering. Our G0n calls return physical values directly, so
      // the same overall factor must be applied to the expected formula.
      double sign        = u.permutation_sign * p.permutation_sign;
      double c2_expected = sign * (G02 - g_00 * g_11 + g_01 * g_10);
      EXPECT_NEAR(c2, c2_expected, kTol)
          << "orbs=(" << ou0 << "," << ou1 << "|" << op0 << "," << op1 << ")"
          << " draw=" << draw << " c2=" << c2 << " expected=" << c2_expected;
    }
  }
}

TEST_F(DimerCumulantTest, Order3_MatchesAnalyticDecomposition) {
  // C₃ = G03 + (C₂·G₁ partitions) + (G₁·G₁·G₁ partitions), with the signs
  // dictated by fermionic Wick combinatorics.
  std::mt19937_64 rng(0x3030u);
  std::uniform_real_distribution<double> U01(0.0, 1.0);
  std::vector<std::array<int, 6>> orb_combos = {
      {2, 0, 3, 2, 0, 3},
      {2, 1, 0, 2, 1, 0},
      {3, 2, 1, 0, 1, 2},
  };
  for (auto const &orbs : orb_combos) {
    std::array<int, 3> ou = {orbs[0], orbs[1], orbs[2]};
    std::array<int, 3> op = {orbs[3], orbs[4], orbs[5]};
    for (int draw = 0; draw < 8; ++draw) {
      std::array<double, 3> tu = {params.beta * U01(rng), params.beta * U01(rng), params.beta * U01(rng)};
      std::array<double, 3> tp = {params.beta * U01(rng), params.beta * U01(rng), params.beta * U01(rng)};

      auto [u, p] = make_args<2, double>({tu[0], tu[1], tu[2]}, {ou[0], ou[1], ou[2]},
                                          {tp[0], tp[1], tp[2]}, {op[0], op[1], op[2]});
      double c3   = cumulant_value<2>(*solver, u, p);

      auto g03_args = combined_for_g0n<2>(
          {{tp[0], op[0]}, {tp[1], op[1]}, {tp[2], op[2]}},
          {{tu[0], ou[0]}, {tu[1], ou[1]}, {tu[2], ou[2]}});
      double G03 = solver->G0n(g03_args);

      auto G1 = [&](int i, int k) { return Gpair<2>(*solver, tu[i], ou[i], tp[k], op[k]); };
      auto C2 = [&](int i, int j, int k, int l) {
        auto [uu, pp] = make_args<2, double>({tu[i], tu[j]}, {ou[i], ou[j]},
                                              {tp[k], tp[l]}, {op[k], op[l]});
        // Strip the (uu × pp) perm-sign factor cumulant_value applies, so we
        // plug the *physical* C₂ value into the order-3 expansion below.
        return cumulant_value<2>(*solver, uu, pp) * uu.permutation_sign * pp.permutation_sign;
      };

      double C12_12_C33 = -C2(0, 1, 0, 1) * G1(2, 2);
      double C12_13_C32 =  C2(0, 1, 0, 2) * G1(2, 1);
      double C12_23_C31 = -C2(0, 1, 1, 2) * G1(2, 0);

      double C13_12_C32 =  C2(0, 2, 0, 1) * G1(1, 2);
      double C13_23_C21 =  C2(0, 2, 1, 2) * G1(1, 0);
      double C13_13_C22 = -C2(0, 2, 0, 2) * G1(1, 1);

      double C23_12_C13 = -C2(1, 2, 0, 1) * G1(0, 2);
      double C23_23_C11 = -C2(1, 2, 1, 2) * G1(0, 0);
      double C23_13_C12 =  C2(1, 2, 0, 2) * G1(0, 1);

      double C02C01_terms = C12_12_C33 + C12_13_C32 + C12_23_C31
                          + C13_12_C32 + C13_23_C21 + C13_13_C22
                          + C23_12_C13 + C23_23_C11 + C23_13_C12;

      double C11_C22_C33 = -G1(0, 0) * G1(1, 1) * G1(2, 2);
      double C11_C23_C32 =  G1(0, 0) * G1(1, 2) * G1(2, 1);
      double C12_C21_C33 =  G1(0, 1) * G1(1, 0) * G1(2, 2);
      double C12_C23_C31 = -G1(0, 1) * G1(1, 2) * G1(2, 0);
      double C13_C22_C31 =  G1(0, 2) * G1(1, 1) * G1(2, 0);
      double C13_C21_C32 = -G1(0, 2) * G1(1, 0) * G1(2, 1);

      double G1G1G1_terms = C11_C22_C33 + C11_C23_C32 + C12_C21_C33
                          + C12_C23_C31 + C13_C22_C31 + C13_C21_C32;

      // Same master-perm-sign correction as in the order-2 test.
      double sign        = u.permutation_sign * p.permutation_sign;
      double c3_expected = sign * (G03 + C02C01_terms + G1G1G1_terms);
      EXPECT_NEAR(c3, c3_expected, kTol)
          << "orbs=(" << ou[0] << "," << ou[1] << "," << ou[2] << "|"
          << op[0] << "," << op[1] << "," << op[2] << ")"
          << " draw=" << draw << " c3=" << c3 << " expected=" << c3_expected;
    }
  }
}

TEST(DimerCumulantNonInteracting, Order2_VanishesAtUEqualsZero) {
  // U = 0 ⇒ Wick's theorem ⇒ G02 factorizes into G₁ × G₁ exactly,
  // so the connected order-2 cumulant must vanish identically.
  Parameters<double> p0{/*U=*/0.0, /*beta=*/2.0, /*mu=*/0.5, /*t=*/1.0, /*bipartite=*/true};
  HubbardSolver<2, double> solver0(p0);

  std::mt19937_64 rng(0x4040u);
  std::uniform_real_distribution<double> U01(0.0, 1.0);
  std::vector<std::array<int, 4>> orb_combos = {
      {2, 0, 2, 0}, {2, 1, 2, 1}, {2, 1, 0, 3}, {3, 0, 1, 2}, {2, 3, 0, 1},
  };
  for (auto const &orbs : orb_combos) {
    int ou0 = orbs[0], ou1 = orbs[1], op0 = orbs[2], op1 = orbs[3];
    for (int draw = 0; draw < 16; ++draw) {
      double tu0 = p0.beta * U01(rng), tu1 = p0.beta * U01(rng);
      double tp0 = p0.beta * U01(rng), tp1 = p0.beta * U01(rng);
      auto [u, p] = make_args<2, double>({tu0, tu1}, {ou0, ou1}, {tp0, tp1}, {op0, op1});
      double c2   = cumulant_value<2>(solver0, u, p);
      EXPECT_NEAR(c2, 0.0, kTol)
          << "orbs=(" << ou0 << "," << ou1 << "|" << op0 << "," << op1 << ")"
          << " draw=" << draw << " c2=" << c2;
    }
  }
}
