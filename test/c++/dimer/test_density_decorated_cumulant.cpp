#include <gtest/gtest.h>
#include <memory>
#include <utility>
#include <vector>
#include "../c++/sc_expansion/hubbard_solver.hpp"
#include "../c++/sc_expansion/cumulant.hpp"

using namespace sc_expansion;

// ============================================================================
// Dimer density-decorated cumulant decomposition tests. Mirrors
// test/c++/atom/test_density_decorated_cumulant.cpp MINUS the infinite-U
// branches (the dimer has no infinite-U path). These are
// decomposition-consistency checks: CumulantSolver::compute_cumulant_
// decomposition() vs a reference assembled from G0n_with_densities / G0n /
// compute_n_sigma. They run through the generic CumulantSolver<2,double> /
// evaluate_plan path unchanged — only G0n_with_densities for N_sites==2
// needed to become real (Step 1a).
//
// Orbital encoding (N_sites=2): 0=site0↓, 1=site1↓, 2=site0↑, 3=site1↑.
// ============================================================================

namespace {

  constexpr double kTol = 1e-12;

  template <typename T> FermionOperator<2, T> mk_op(int orbital, bool creation) {
    uint8_t id = static_cast<uint8_t>(orbital);
    if (creation) id |= FermionOperator<2, T>::ACTION_BIT;
    return FermionOperator<2, T>(id);
  }

} // namespace

class DimerDensityDecoratedCumulant : public ::testing::Test {
  protected:
  double U = 4.0, beta = 1.0, mu = 1.2, t = 1.0;
  Parameters<double> params{U, beta, mu, t, true};
  std::unique_ptr<HubbardSolver<2, double>> solver;
  void SetUp() override { solver = std::make_unique<HubbardSolver<2, double>>(params); }
};

// Single-density connected cumulant: ⟨n_σ⟩₀,c = ⟨n_σ⟩₀.
TEST_F(DimerDensityDecoratedCumulant, SingleDensity_MatchesOccupation) {
  for (int sigma = 0; sigma < 4; ++sigma) {
    Args<2, double> u({}, {});
    Args<2, double> p({}, {});
    CumulantSolver<2, double> cs(u, p, *solver, /*infinite_U=*/false);
    cs.add_static_density(sigma);
    double value    = cs.compute_cumulant_decomposition();
    double expected = solver->compute_n_sigma(sigma);
    EXPECT_NEAR(value, expected, kTol) << "σ=" << sigma;
  }
}

// Same-spin double density: ⟨n_σ n_σ⟩₀,c = ⟨n_σ⟩ − ⟨n_σ⟩² (n_σ idempotent).
TEST_F(DimerDensityDecoratedCumulant, SameSpinDoubleDensity_EqualsVariance) {
  for (int sigma = 0; sigma < 4; ++sigma) {
    Args<2, double> u({}, {});
    Args<2, double> p({}, {});
    CumulantSolver<2, double> cs(u, p, *solver, false);
    cs.add_static_density(sigma);
    cs.add_static_density(sigma);
    double value    = cs.compute_cumulant_decomposition();
    double n        = solver->compute_n_sigma(sigma);
    double expected = n - n * n;
    EXPECT_NEAR(value, expected, kTol) << "σ=" << sigma;
  }
}

// On-site opposite spin: ⟨n_{i↑} n_{i↓}⟩₀,c = ⟨n_{i↑} n_{i↓}⟩₀ − ⟨n_{i↑}⟩⟨n_{i↓}⟩.
TEST_F(DimerDensityDecoratedCumulant, OnSiteOppositeSpin_EqualsDoubleOccMinusProduct) {
  // site0: {↓=0, ↑=2}; site1: {↓=1, ↑=3}.
  std::vector<std::pair<int, int>> onsite = {{0, 2}, {1, 3}};
  for (auto const &[dn, up] : onsite) {
    Args<2, double> u({}, {});
    Args<2, double> p({}, {});
    CumulantSolver<2, double> cs(u, p, *solver, false);
    cs.add_static_density(dn);
    cs.add_static_density(up);
    double value = cs.compute_cumulant_decomposition();

    Args<2, double> empty_args({}, {});
    double D        = solver->G0n_with_densities(empty_args, {dn, up});
    double n_dn     = solver->compute_n_sigma(dn);
    double n_up     = solver->compute_n_sigma(up);
    double expected = D - n_dn * n_up;
    EXPECT_NEAR(value, expected, kTol) << "site dn=" << dn << " up=" << up;
  }
}

// Inter-site double density: ⟨n_{0σ} n_{1σ'}⟩₀,c = ⟨n_{0σ} n_{1σ'}⟩₀ −
// ⟨n_{0σ}⟩⟨n_{1σ'}⟩. The zeroth-order (single-vertex) contribution to the
// intra-dimer correlator.
TEST_F(DimerDensityDecoratedCumulant, InterSiteDoubleDensity_EqualsConnectedCorrelator) {
  std::vector<std::pair<int, int>> inter = {{0, 1}, {0, 3}, {2, 1}, {2, 3}};
  for (auto const &[s0, s1] : inter) {
    Args<2, double> u({}, {});
    Args<2, double> p({}, {});
    CumulantSolver<2, double> cs(u, p, *solver, false);
    cs.add_static_density(s0);
    cs.add_static_density(s1);
    double value = cs.compute_cumulant_decomposition();

    Args<2, double> empty_args({}, {});
    double corr     = solver->G0n_with_densities(empty_args, {s0, s1});
    double n0       = solver->compute_n_sigma(s0);
    double n1       = solver->compute_n_sigma(s1);
    double expected = corr - n0 * n1;
    EXPECT_NEAR(value, expected, kTol) << "orbs=(" << s0 << "," << s1 << ")";
  }
}

// One density + one hybridization pair:
//   κ₂({n_σ, P}) = ⟨n_σ P⟩₀ − ⟨n_σ⟩ ⟨P⟩₀.
TEST_F(DimerDensityDecoratedCumulant, OneDensityPlusPair_EqualsDecoratedMinusDisconnected) {
  int const leg = 2; // site0↑
  double const tau_p = 0.2, tau_u = 0.7;

  Args<2, double> u({tau_u}, {mk_op<double>(leg, false)});
  Args<2, double> p_args({tau_p}, {mk_op<double>(leg, true)});
  Args<2, double> pair_args({tau_p, tau_u}, {mk_op<double>(leg, true), mk_op<double>(leg, false)});

  for (int sigma_dens = 0; sigma_dens < 4; ++sigma_dens) {
    CumulantSolver<2, double> cs(u, p_args, *solver, false);
    cs.add_static_density(sigma_dens);
    double value = cs.compute_cumulant_decomposition();

    double decorated = solver->G0n_with_densities(pair_args, {sigma_dens});
    double G_pair    = solver->G0n(pair_args);
    double n         = solver->compute_n_sigma(sigma_dens);
    double expected  = decorated - n * G_pair;
    EXPECT_NEAR(value, expected, kTol) << "σ_dens=" << sigma_dens;
  }
}

// Same-spin two densities + one hybridization pair:
//   κ₃({n_σ, n_σ, P}) = (1 − 2⟨n_σ⟩) · κ₂({n_σ, P}).
TEST_F(DimerDensityDecoratedCumulant, SameSpinTwoDensitiesPlusPair) {
  int const leg = 2;
  double const tau_p = 0.2, tau_u = 0.7;

  Args<2, double> u({tau_u}, {mk_op<double>(leg, false)});
  Args<2, double> p_args({tau_p}, {mk_op<double>(leg, true)});
  Args<2, double> pair_args({tau_p, tau_u}, {mk_op<double>(leg, true), mk_op<double>(leg, false)});

  for (int sigma_dens = 0; sigma_dens < 4; ++sigma_dens) {
    CumulantSolver<2, double> cs(u, p_args, *solver, false);
    cs.add_static_density(sigma_dens);
    cs.add_static_density(sigma_dens);
    double value = cs.compute_cumulant_decomposition();

    double decorated = solver->G0n_with_densities(pair_args, {sigma_dens});
    double G_pair    = solver->G0n(pair_args);
    double n         = solver->compute_n_sigma(sigma_dens);
    double kappa2    = decorated - n * G_pair;
    double expected  = (1.0 - 2.0 * n) * kappa2;
    EXPECT_NEAR(value, expected, kTol) << "σ_dens=" << sigma_dens;
  }
}

// On-site opposite-spin two densities + one hybridization pair:
//   κ₃({n_↑, n_↓, P}) = ⟨n_↑ n_↓ P⟩₀ − D⟨P⟩₀
//                       − ⟨n_↓⟩⟨n_↑ P⟩₀ − ⟨n_↑⟩⟨n_↓ P⟩₀
//                       + 2 ⟨n_↑⟩⟨n_↓⟩⟨P⟩₀.
TEST_F(DimerDensityDecoratedCumulant, DoubleOccupancyPlusPair) {
  // Densities on site0: ↓=0, ↑=2.
  int const dn = 0, up = 2, leg = 2;
  double const tau_p = 0.2, tau_u = 0.7;

  Args<2, double> u({tau_u}, {mk_op<double>(leg, false)});
  Args<2, double> p_args({tau_p}, {mk_op<double>(leg, true)});
  Args<2, double> pair_args({tau_p, tau_u}, {mk_op<double>(leg, true), mk_op<double>(leg, false)});
  Args<2, double> empty_args({}, {});

  CumulantSolver<2, double> cs(u, p_args, *solver, false);
  cs.add_static_density(dn);
  cs.add_static_density(up);
  double value = cs.compute_cumulant_decomposition();

  double n_up   = solver->compute_n_sigma(up);
  double n_dn   = solver->compute_n_sigma(dn);
  double D      = solver->G0n_with_densities(empty_args, {dn, up});
  double full   = solver->G0n_with_densities(pair_args, {dn, up});
  double dec_up = solver->G0n_with_densities(pair_args, {up});
  double dec_dn = solver->G0n_with_densities(pair_args, {dn});
  double G_pair = solver->G0n(pair_args);

  double expected = full - D * G_pair - n_dn * dec_up - n_up * dec_dn + 2.0 * n_up * n_dn * G_pair;
  EXPECT_NEAR(value, expected, kTol) << "value=" << value << " expected=" << expected;
}

// One density + two hybridization pairs (full sign decomposition). The
// Wick-swap structure is a property of the cumulant partition lattice, not the
// cluster, so the decomposition is identical to the atomic case.
TEST_F(DimerDensityDecoratedCumulant, OneDensityPlusTwoPairs_FullSignDecomposition) {
  int const sigma = 2; // site0↑ legs
  double const tu1 = 0.7, tu2 = 0.4;
  double const tp1 = 0.2, tp2 = 0.55;

  Args<2, double> u({tu1, tu2}, {mk_op<double>(sigma, false), mk_op<double>(sigma, false)});
  Args<2, double> p_args({tp1, tp2}, {mk_op<double>(sigma, true), mk_op<double>(sigma, true)});

  auto pair = [&](double tp, double tu) {
    return Args<2, double>({tp, tu}, {mk_op<double>(sigma, true), mk_op<double>(sigma, false)});
  };
  auto twopair_direct = [&]() {
    return Args<2, double>({tp1, tu1, tp2, tu2}, {mk_op<double>(sigma, true), mk_op<double>(sigma, false),
                                                  mk_op<double>(sigma, true), mk_op<double>(sigma, false)});
  };

  for (int sigma_dens = 0; sigma_dens < 4; ++sigma_dens) {
    CumulantSolver<2, double> cs(u, p_args, *solver, false);
    cs.add_static_density(sigma_dens);
    double value = cs.compute_cumulant_decomposition();

    auto G0n  = [&](Args<2, double> const &a) { return solver->G0n(a); };
    auto G0nD = [&](Args<2, double> const &a, std::vector<int> const &ds) { return solver->G0n_with_densities(a, ds); };

    auto kappa2_pairs = [&](double ta_p, double ta_u, double tb_p, double tb_u) {
      Args<2, double> uu({ta_u, tb_u}, {mk_op<double>(sigma, false), mk_op<double>(sigma, false)});
      Args<2, double> pp({ta_p, tb_p}, {mk_op<double>(sigma, true), mk_op<double>(sigma, true)});
      CumulantSolver<2, double> sub(uu, pp, *solver, false);
      return sub.compute_cumulant_decomposition() * uu.permutation_sign * pp.permutation_sign;
    };
    auto kappa2_n_pair = [&](double tp_, double tu_) {
      Args<2, double> uu({tu_}, {mk_op<double>(sigma, false)});
      Args<2, double> pp({tp_}, {mk_op<double>(sigma, true)});
      CumulantSolver<2, double> sub(uu, pp, *solver, false);
      sub.add_static_density(sigma_dens);
      return sub.compute_cumulant_decomposition() * uu.permutation_sign * pp.permutation_sign;
    };

    double n     = G0nD(Args<2, double>({}, {}), {sigma_dens});
    double P1    = G0n(pair(tp1, tu1));
    double P2    = G0n(pair(tp2, tu2));
    double s12   = G0n(pair(tp1, tu2)); // ⟨c†(1')c(2)⟩
    double s21   = G0n(pair(tp2, tu1)); // ⟨c†(2')c(1)⟩
    double nP1P2 = G0nD(twopair_direct(), {sigma_dens});

    double k2_P1P2  = kappa2_pairs(tp1, tu1, tp2, tu2);
    double k2_nP1   = kappa2_n_pair(tp1, tu1);
    double k2_nP2   = kappa2_n_pair(tp2, tu2);
    double k2_n_s12 = kappa2_n_pair(tp1, tu2);
    double k2_n_s21 = kappa2_n_pair(tp2, tu1);

    double sign     = u.permutation_sign * p_args.permutation_sign;
    double expected = sign * (nP1P2 - n * k2_P1P2 - n * P1 * P2 + n * s12 * s21 - k2_nP1 * P2 - P1 * k2_nP2 +
                              k2_n_s12 * s21 + s12 * k2_n_s21);

    EXPECT_NEAR(value, expected, kTol) << "σ_dens=" << sigma_dens << " value=" << value << " expected=" << expected;
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
