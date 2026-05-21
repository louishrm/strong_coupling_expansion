#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include "../c++/sc_expansion/hubbard_solver.hpp"
#include "../c++/sc_expansion/cumulant.hpp"

using namespace sc_expansion;

namespace {
  constexpr double kTol = 1e-12;
} // namespace

// Single-density connected cumulant: ⟨n_σ⟩₀,c = ⟨n_σ⟩₀.
// For the atomic Hubbard model H = U n↑n↓ − μ(n↑+n↓):
//   finite U  : ⟨n_σ⟩ = (e^{βμ} + e^{β(2μ−U)}) / (1 + 2 e^{βμ} + e^{β(2μ−U)})
//   infinite U: ⟨n_σ⟩ = e^{βμ} / (1 + 2 e^{βμ})
TEST(AtomDensityDecoratedCumulant, SingleDensity_MatchesAtomicOccupation) {
  double const U = 4.0, beta = 1.0, mu = 1.2;
  Parameters<double> params{U, beta, mu, 0.0, true};
  HubbardSolver<1, double> solver(params);

  double const eBmu       = std::exp(beta * mu);
  double const eB2mU      = std::exp(beta * (2.0 * mu - U));
  double const n_finite   = (eBmu + eB2mU) / (1.0 + 2.0 * eBmu + eB2mU);
  double const n_infinite = eBmu / (1.0 + 2.0 * eBmu);

  for (int sigma : {0, 1}) {
    for (bool infinite_U : {false, true}) {
      Args<1, double> u({}, {});
      Args<1, double> p({}, {});
      CumulantSolver<1, double> cs(u, p, solver, infinite_U);
      cs.add_static_density(sigma);
      double value = cs.compute_cumulant_decomposition();

      double expected = infinite_U ? n_infinite : n_finite;
      EXPECT_NEAR(value, expected, kTol) << "σ=" << sigma << " infinite_U=" << infinite_U << " value=" << value << " expected=" << expected;
    }
  }
}

// Same-spin double density: ⟨n_σ n_σ⟩₀,c = ⟨n_σ²⟩ − ⟨n_σ⟩² = ⟨n_σ⟩ − ⟨n_σ⟩²
// (n_σ is idempotent, so the connected cumulant collapses to the variance.)
TEST(AtomDensityDecoratedCumulant, SameSpinDoubleDensity_EqualsVariance) {
  double const U = 4.0, beta = 1.0, mu = 1.2;
  Parameters<double> params{U, beta, mu, 0.0, true};
  HubbardSolver<1, double> solver(params);

  double const eBmu       = std::exp(beta * mu);
  double const eB2mU      = std::exp(beta * (2.0 * mu - U));
  double const n_finite   = (eBmu + eB2mU) / (1.0 + 2.0 * eBmu + eB2mU);
  double const n_infinite = eBmu / (1.0 + 2.0 * eBmu);

  for (int sigma : {0, 1}) {
    for (bool infinite_U : {false, true}) {
      Args<1, double> u({}, {});
      Args<1, double> p({}, {});
      CumulantSolver<1, double> cs(u, p, solver, infinite_U);
      cs.add_static_density(sigma);
      cs.add_static_density(sigma);
      double value = cs.compute_cumulant_decomposition();

      double n        = infinite_U ? n_infinite : n_finite;
      double expected = n - n * n;
      EXPECT_NEAR(value, expected, kTol) << "σ=" << sigma << " infinite_U=" << infinite_U << " value=" << value << " expected=" << expected;
    }
  }
}

// Opposite-spin double density: ⟨n↑ n↓⟩₀,c = ⟨n↑ n↓⟩₀ − ⟨n↑⟩₀⟨n↓⟩₀.
//   finite U  : ⟨n↑ n↓⟩₀ = e^{β(2μ−U)} / (1 + 2 e^{βμ} + e^{β(2μ−U)})
//   infinite U: ⟨n↑ n↓⟩₀ = 0 (double occupancy projected out), so the
//               connected cumulant reduces to −⟨n↑⟩⟨n↓⟩.
TEST(AtomDensityDecoratedCumulant, OppositeSpinDoubleDensity_EqualsDoubleOccupancyMinusProduct) {
  double const U = 4.0, beta = 1.0, mu = 1.2;
  Parameters<double> params{U, beta, mu, 0.0, true};
  HubbardSolver<1, double> solver(params);

  double const eBmu        = std::exp(beta * mu);
  double const eB2mU       = std::exp(beta * (2.0 * mu - U));
  double const Z_finite    = 1.0 + 2.0 * eBmu + eB2mU;
  double const n_finite    = (eBmu + eB2mU) / Z_finite;
  double const docc_finite = eB2mU / Z_finite;
  double const n_infinite  = eBmu / (1.0 + 2.0 * eBmu);

  for (bool infinite_U : {false, true}) {
    Args<1, double> u({}, {});
    Args<1, double> p({}, {});
    CumulantSolver<1, double> cs(u, p, solver, infinite_U);
    cs.add_static_density(0); // ↓
    cs.add_static_density(1); // ↑
    double value = cs.compute_cumulant_decomposition();

    double expected = infinite_U ? -(n_infinite * n_infinite) : (docc_finite - n_finite * n_finite);
    EXPECT_NEAR(value, expected, kTol) << "infinite_U=" << infinite_U << " value=" << value << " expected=" << expected;
  }
}

namespace {
  template <typename T> FermionOperator<1, T> mk_op(int orbital, bool creation) {
    uint8_t id = (uint8_t)orbital;
    if (creation) id |= FermionOperator<1, T>::ACTION_BIT;
    return FermionOperator<1, T>(id);
  }
} // namespace

// One density + one hybridization pair: connected decorated cumulant
//   ⟨n_σ T c†(1') c(1)⟩₀,c = ⟨n_σ T c†(1') c(1)⟩₀ − ⟨n_σ⟩ ⟨T c†(1') c(1)⟩₀
// Verified for both finite U and infinite U.
TEST(AtomDensityDecoratedCumulant, OneDensityPlusPair_EqualsDecoratedMinusDisconnected) {
  double const U = 4.0, beta = 1.0, mu = 1.2;
  Parameters<double> params{U, beta, mu, 0.0, true};
  HubbardSolver<1, double> solver(params);

  // Pair: c†_↑(τ_p) c_↑(τ_u). Choose τ_u > τ_p so the pair has a non-trivial
  // (non-equal-time) value.
  int const sigma_up = 1;
  double const tau_p = 0.2, tau_u = 0.7;

  auto [u, p_args] = std::pair{Args<1, double>({tau_u}, {mk_op<double>(sigma_up, /*creation=*/false)}),
                               Args<1, double>({tau_p}, {mk_op<double>(sigma_up, /*creation=*/true)})};

  Args<1, double> pair_args({tau_p, tau_u}, {mk_op<double>(sigma_up, true), mk_op<double>(sigma_up, false)});

  for (int sigma_dens : {0, 1}) {
    for (bool infinite_U : {false, true}) {
      CumulantSolver<1, double> cs(u, p_args, solver, infinite_U);
      cs.add_static_density(sigma_dens);
      double value = cs.compute_cumulant_decomposition();

      double decorated =
         infinite_U ? solver.G0n_with_densities_infinite_U(pair_args, {sigma_dens}) : solver.G0n_with_densities(pair_args, {sigma_dens});
      double G_pair = infinite_U ? solver.G0n_infinite_U(pair_args) : solver.G0n(pair_args);

      double const eBmu  = std::exp(beta * mu);
      double const eB2mU = std::exp(beta * (2.0 * mu - U));
      double n           = infinite_U ? eBmu / (1.0 + 2.0 * eBmu) : (eBmu + eB2mU) / (1.0 + 2.0 * eBmu + eB2mU);

      double expected = decorated - n * G_pair;
      EXPECT_NEAR(value, expected, kTol) << "σ_dens=" << sigma_dens << " infinite_U=" << infinite_U << " value=" << value << " expected=" << expected;
    }
  }
}

// Same-spin two densities + one hybridization pair:
//   κ₃({n_σ, n_σ, P}) = (1 − 2⟨n_σ⟩) · κ₂({n_σ, P})
//                     = (1 − 2⟨n_σ⟩) · (⟨n_σ P⟩₀ − ⟨n_σ⟩⟨P⟩₀)
// Both finite and infinite U.
TEST(AtomDensityDecoratedCumulant, SameSpinTwoDensitiesPlusPair_EqualsRescaledOneDensityCumulant) {
  double const U = 4.0, beta = 1.0, mu = 1.2;
  Parameters<double> params{U, beta, mu, 0.0, true};
  HubbardSolver<1, double> solver(params);

  int const sigma_up = 1;
  double const tau_p = 0.2, tau_u = 0.7;

  Args<1, double> u({tau_u}, {mk_op<double>(sigma_up, false)});
  Args<1, double> p_args({tau_p}, {mk_op<double>(sigma_up, true)});
  Args<1, double> pair_args({tau_p, tau_u}, {mk_op<double>(sigma_up, true), mk_op<double>(sigma_up, false)});

  double const eBmu       = std::exp(beta * mu);
  double const eB2mU      = std::exp(beta * (2.0 * mu - U));
  double const n_finite   = (eBmu + eB2mU) / (1.0 + 2.0 * eBmu + eB2mU);
  double const n_infinite = eBmu / (1.0 + 2.0 * eBmu);

  for (int sigma_dens : {0, 1}) {
    for (bool infinite_U : {false, true}) {
      CumulantSolver<1, double> cs(u, p_args, solver, infinite_U);
      cs.add_static_density(sigma_dens);
      cs.add_static_density(sigma_dens);
      double value = cs.compute_cumulant_decomposition();

      double decorated =
         infinite_U ? solver.G0n_with_densities_infinite_U(pair_args, {sigma_dens}) : solver.G0n_with_densities(pair_args, {sigma_dens});
      double G_pair = infinite_U ? solver.G0n_infinite_U(pair_args) : solver.G0n(pair_args);
      double n      = infinite_U ? n_infinite : n_finite;

      double kappa2   = decorated - n * G_pair;
      double expected = (1.0 - 2.0 * n) * kappa2;
      EXPECT_NEAR(value, expected, kTol) << "σ_dens=" << sigma_dens << " infinite_U=" << infinite_U << " value=" << value << " expected=" << expected;
    }
  }
}

// Opposite-spin two densities + one hybridization pair (double occupancy
// vertex with a hybridization pair):
//   κ₃({n_↑, n_↓, P}) = ⟨n_↑ n_↓ P⟩₀ − D⟨P⟩₀
//                       − ⟨n_↓⟩⟨n_↑ P⟩₀ − ⟨n_↑⟩⟨n_↓ P⟩₀
//                       + 2 ⟨n_↑⟩⟨n_↓⟩⟨P⟩₀
// where D = ⟨n_↑ n_↓⟩₀ = e^{β(2μ−U)}/Z  (finite U), = 0 (infinite U).
TEST(AtomDensityDecoratedCumulant, DoubleOccupancyPlusPair_FullDecomposition) {
  double const U = 4.0, beta = 1.0, mu = 1.2;
  Parameters<double> params{U, beta, mu, 0.0, true};
  HubbardSolver<1, double> solver(params);

  int const sigma_up = 1;
  double const tau_p = 0.2, tau_u = 0.7;

  Args<1, double> u({tau_u}, {mk_op<double>(sigma_up, false)});
  Args<1, double> p_args({tau_p}, {mk_op<double>(sigma_up, true)});
  Args<1, double> pair_args({tau_p, tau_u}, {mk_op<double>(sigma_up, true), mk_op<double>(sigma_up, false)});

  double const eBmu        = std::exp(beta * mu);
  double const eB2mU       = std::exp(beta * (2.0 * mu - U));
  double const Z_finite    = 1.0 + 2.0 * eBmu + eB2mU;
  double const n_finite    = (eBmu + eB2mU) / Z_finite;
  double const docc_finite = eB2mU / Z_finite;
  double const n_infinite  = eBmu / (1.0 + 2.0 * eBmu);

  for (bool infinite_U : {false, true}) {
    CumulantSolver<1, double> cs(u, p_args, solver, infinite_U);
    cs.add_static_density(0); // ↓
    cs.add_static_density(1); // ↑
    double value = cs.compute_cumulant_decomposition();

    double n_up = infinite_U ? n_infinite : n_finite;
    double n_dn = n_up;
    double D    = infinite_U ? 0.0 : docc_finite;

    double full   = infinite_U ? solver.G0n_with_densities_infinite_U(pair_args, {0, 1}) : solver.G0n_with_densities(pair_args, {0, 1});
    double dec_up = infinite_U ? solver.G0n_with_densities_infinite_U(pair_args, {1}) : solver.G0n_with_densities(pair_args, {1});
    double dec_dn = infinite_U ? solver.G0n_with_densities_infinite_U(pair_args, {0}) : solver.G0n_with_densities(pair_args, {0});
    double G_pair = infinite_U ? solver.G0n_infinite_U(pair_args) : solver.G0n(pair_args);

    double expected = full - D * G_pair - n_dn * dec_up - n_up * dec_dn + 2.0 * n_up * n_dn * G_pair;
    EXPECT_NEAR(value, expected, kTol) << "infinite_U=" << infinite_U << " value=" << value << " expected=" << expected;
  }
}

// One density + two hybridization pairs (sign-logic test).
// Composite-atom Möbius inversion on {n_σ, P1, P2}, with Wick swap
// (1'↔2, 2'↔1) generating additional partitions in which n can decorate
// either swap pair:
//   κ_3 = + ⟨n_σ P1 P2⟩₀
//         − ⟨n_σ⟩₀ · κ_2({P1,P2})
//         − ⟨n_σ⟩₀ ⟨P1⟩₀ ⟨P2⟩₀
//         + ⟨n_σ⟩₀ ⟨c†(1')c(2)⟩₀ ⟨c†(2')c(1)⟩₀
//         − κ_2({n_σ,P1}) · ⟨P2⟩₀
//         − ⟨P1⟩₀ · κ_2({n_σ,P2})
//         + κ_2({n_σ, c†(1')c(2)}) · ⟨c†(2')c(1)⟩₀
//         + ⟨c†(1')c(2)⟩₀ · κ_2({n_σ, c†(2')c(1)})
// with:
//   κ_2({P1,P2})    = ⟨P1 P2⟩₀ − ⟨P1⟩₀⟨P2⟩₀ + ⟨c†(1')c(2)⟩₀⟨c†(2')c(1)⟩₀
//   κ_2({n_σ, X})   = ⟨n_σ X⟩₀ − ⟨n_σ⟩₀⟨X⟩₀          (X = direct or swap pair)
// Both finite and infinite U.
TEST(AtomDensityDecoratedCumulant, OneDensityPlusTwoPairs_FullSignDecomposition) {
  double const U = 4.0, beta = 1.0, mu = 1.2;
  Parameters<double> params{U, beta, mu, 0.0, true};
  HubbardSolver<1, double> solver(params);

  int const sigma  = 1; // ↑ for the hybridization legs
  double const tu1 = 0.7, tu2 = 0.4;
  double const tp1 = 0.2, tp2 = 0.55;

  // Master args: 2 destructions (c(1), c(2)), 2 creations (c†(1'), c†(2')).
  Args<1, double> u({tu1, tu2}, {mk_op<double>(sigma, false), mk_op<double>(sigma, false)});
  Args<1, double> p_args({tp1, tp2}, {mk_op<double>(sigma, true), mk_op<double>(sigma, true)});

  // Reference building blocks. All G0n calls use (c†, c) interleaved layout.
  auto pair           = [&](double tp, double tu) { return Args<1, double>({tp, tu}, {mk_op<double>(sigma, true), mk_op<double>(sigma, false)}); };
  auto twopair_direct = [&]() {
    // c†(1') c(1) c†(2') c(2)
    return Args<1, double>({tp1, tu1, tp2, tu2},
                           {mk_op<double>(sigma, true), mk_op<double>(sigma, false), mk_op<double>(sigma, true), mk_op<double>(sigma, false)});
  };

  for (int sigma_dens : {0, 1}) {
    for (bool infinite_U : {false, true}) {
      CumulantSolver<1, double> cs(u, p_args, solver, infinite_U);
      cs.add_static_density(sigma_dens);
      double value = cs.compute_cumulant_decomposition();

      auto G0n  = [&](Args<1, double> const &a) { return infinite_U ? solver.G0n_infinite_U(a) : solver.G0n(a); };
      auto G0nD = [&](Args<1, double> const &a, std::vector<int> const &ds) {
        return infinite_U ? solver.G0n_with_densities_infinite_U(a, ds) : solver.G0n_with_densities(a, ds);
      };

      // Hybridization-only sub-cumulant κ_2({P_a, P_b}) via CumulantSolver.
      // Returned in user-input ordering — strip the perm-sign factor to get
      // the physical sub-cumulant we can plug into the expansion.
      auto kappa2_pairs = [&](double ta_p, double ta_u, double tb_p, double tb_u) {
        Args<1, double> uu({ta_u, tb_u}, {mk_op<double>(sigma, false), mk_op<double>(sigma, false)});
        Args<1, double> pp({ta_p, tb_p}, {mk_op<double>(sigma, true), mk_op<double>(sigma, true)});
        CumulantSolver<1, double> sub(uu, pp, solver, infinite_U);
        return sub.compute_cumulant_decomposition() * uu.permutation_sign * pp.permutation_sign;
      };
      // κ_2({n_σ, c†(tp)c(tu)}) via CumulantSolver — validated in earlier tests.
      auto kappa2_n_pair = [&](double tp_, double tu_) {
        Args<1, double> uu({tu_}, {mk_op<double>(sigma, false)});
        Args<1, double> pp({tp_}, {mk_op<double>(sigma, true)});
        CumulantSolver<1, double> sub(uu, pp, solver, infinite_U);
        sub.add_static_density(sigma_dens);
        return sub.compute_cumulant_decomposition() * uu.permutation_sign * pp.permutation_sign;
      };

      double n     = G0nD(Args<1, double>({}, {}), {sigma_dens});
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

      // CumulantSolver returns the result in the user-input ordering; its
      // internal stable (sorted-by-τ) machinery applies (u.perm × p.perm) at
      // the root. Our reference is built from physical G0n calls, so we apply
      // the same factor.
      double sign     = u.permutation_sign * p_args.permutation_sign;
      double expected = sign * (nP1P2 - n * k2_P1P2 - n * P1 * P2 + n * s12 * s21 - k2_nP1 * P2 - P1 * k2_nP2 + k2_n_s12 * s21 + s12 * k2_n_s21);

      EXPECT_NEAR(value, expected, kTol) << "σ_dens=" << sigma_dens << " infinite_U=" << infinite_U << " value=" << value << " expected=" << expected;
    }
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
