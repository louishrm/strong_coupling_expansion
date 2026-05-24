
#include <gtest/gtest.h>
#include "sc_expansion/atomic/sum_diagrams.hpp"
#include "sc_expansion/hilbert_traits.hpp"
#include <vector>

using namespace sc_expansion;

// Thin wrapper around SumDiagrams::density_density_infinite_U_coefficient.
// Reference coefficients below are computed on a 2-site dimer cluster, where
// every (graph, marks, r) admits exactly one embedding. We bypass the default
// Z² embedding-count path by passing override_lm = 1 so each diagram enters
// with multiplier 1 instead of its actual Z² count — i.e. this is a
// per-diagram, no-lattice-sum check of the SJT evaluator and the rooted
// catalog construction.
//
// spin convention: 0 = ↓, 1 = ↑.
static double compute_infinite_U_coefficient(int order, double beta, double mu, std::vector<int> const &r, int s1, int s2) {
  Parameters<double> params{/*U=*/0.0, beta, mu, 0.0, /*bipartite=*/true};
  atomic::SumDiagrams<double> calculator(params, order, r, s1, s2, /*override_lm=*/1);
  auto coeff_map = calculator.density_density_infinite_U_coefficient();
  return coeff_map.at(calculator.get_target_d_sq()).second;
}

// -----------------------------------------------------------------------
// Tests: 2×2 grid of (r ∈ {0, (1,0)}) × ((s1,s2) ∈ {(↑,↑), (↑,↓)}).
// Reference values are exact infinite-U TDL coefficients supplied by the
// user.
// -----------------------------------------------------------------------

namespace {
  [[maybe_unused]] constexpr int SPIN_DOWN = 0;
  constexpr int SPIN_UP   = 1;
} // namespace

TEST(CorrelatorAtomInfU, OnSiteSameSpinOrder2) {
  // ⟨n_{0,↑} n_{0,↑}⟩ at U=∞, β=2, μ=1, t² coefficient on the bipartite
  // hypercubic lattice. Reference: -3.285401994562779e-03.
  double c           = compute_infinite_U_coefficient(/*order=*/2, /*beta=*/2.0, /*mu=*/1.0,
                                                      /*r=*/{0, 0}, SPIN_UP, SPIN_UP);
  double exact_coeff = -3.285401994562779e-03;
  EXPECT_NEAR(c, exact_coeff, 1e-12) << "Computed " << c << " vs exact " << exact_coeff << " (Δ = " << (c - exact_coeff) << ")";
}

TEST(CorrelatorAtomInfU, OnSiteSameSpinOrder4) {
  // ⟨n_{0,↑} n_{0,↑}⟩ at U=∞, β=2, μ=1, t⁴ coefficient.
  // Reference: -3.002141546904583e-03.
  double c           = compute_infinite_U_coefficient(/*order=*/4, /*beta=*/2.0, /*mu=*/1.0,
                                                      /*r=*/{0, 0}, SPIN_UP, SPIN_UP);
  double exact_coeff = -3.002141546904583e-03;
  EXPECT_NEAR(c, exact_coeff, 1e-12) << "Computed " << c << " vs exact " << exact_coeff << " (Δ = " << (c - exact_coeff) << ")";
}

TEST(CorrelatorAtomInfU, OnSiteOppositeSpinOrder2) {
  // ⟨n_{0,↑} n_{0,↓}⟩ at U=∞, β=2, μ=1, t² coefficient on the bipartite
  // hypercubic lattice. The full value of n_↑n_↓ vanishes at U=∞, but the
  // connected diagrammatic coefficient is nonzero through the cumulant
  // subtraction. Reference: 0.04855203929072596.
  double c           = compute_infinite_U_coefficient(/*order=*/2, /*beta=*/2.0, /*mu=*/1.0,
                                                      /*r=*/{0, 0}, SPIN_UP, SPIN_DOWN);
  double exact_coeff = 0.04855203929072596;
  EXPECT_NEAR(c, exact_coeff, 1e-12) << "Computed " << c << " vs exact " << exact_coeff << " (Δ = " << (c - exact_coeff) << ")";
}

TEST(CorrelatorAtomInfU, OnSiteOppositeSpinOrder4) {
  // ⟨n_{0,↑} n_{0,↓}⟩ at U=∞, β=2, μ=1, t⁴ coefficient.
  // Reference: 0.001968298731591109.
  double c           = compute_infinite_U_coefficient(/*order=*/4, /*beta=*/2.0, /*mu=*/1.0,
                                                      /*r=*/{0, 0}, SPIN_UP, SPIN_DOWN);
  double exact_coeff = 0.001968298731591109;
  EXPECT_NEAR(c, exact_coeff, 1e-12) << "Computed " << c << " vs exact " << exact_coeff << " (Δ = " << (c - exact_coeff) << ")";
}

TEST(CorrelatorAtomInfU, NearestNeighborSameSpinOrder2) {
  // ⟨n_{(1,0),↑} n_{0,↑}⟩ at U=∞, β=2, μ=1, t² coefficient.
  // Reference: -0.0035238528031628622.
  double c           = compute_infinite_U_coefficient(/*order=*/2, /*beta=*/2.0, /*mu=*/1.0,
                                                      /*r=*/{1, 0}, SPIN_UP, SPIN_UP);
  double exact_coeff = -0.0035238528031628622;
  EXPECT_NEAR(c, exact_coeff, 1e-12) << "Computed " << c << " vs exact " << exact_coeff << " (Δ = " << (c - exact_coeff) << ")";
}

TEST(CorrelatorAtomInfU, NearestNeighborSameSpinOrder4) {
  // ⟨n_{(1,0),↑} n_{0,↑}⟩ at U=∞, β=2, μ=1, t⁴ coefficient.
  // Reference: -0.003025005435902937.
  double c           = compute_infinite_U_coefficient(/*order=*/4, /*beta=*/2.0, /*mu=*/1.0,
                                                      /*r=*/{1, 0}, SPIN_UP, SPIN_UP);
  double exact_coeff = -0.003025005435902937;
  EXPECT_NEAR(c, exact_coeff, 1e-12) << "Computed " << c << " vs exact " << exact_coeff << " (Δ = " << (c - exact_coeff) << ")";
}

TEST(CorrelatorAtomInfU, NearestNeighborOppositeSpinOrder2) {
  // ⟨n_{(1,0),↑} n_{0,↓}⟩ at U=∞, β=2, μ=1, t² coefficient.
  // Reference: -0.0035238528031628622.
  double c           = compute_infinite_U_coefficient(/*order=*/2, /*beta=*/2.0, /*mu=*/1.0,
                                                      /*r=*/{1, 0}, SPIN_UP, SPIN_DOWN);
  double exact_coeff = -0.0035238528031628622;
  EXPECT_NEAR(c, exact_coeff, 1e-12) << "Computed " << c << " vs exact " << exact_coeff << " (Δ = " << (c - exact_coeff) << ")";
}

TEST(CorrelatorAtomInfU, NearestNeighborOppositeSpinOrder4) {
  // ⟨n_{(1,0),↑} n_{0,↓}⟩ at U=∞, β=2, μ=1, t⁴ coefficient.
  // Reference: -0.003025005435902937.
  double c           = compute_infinite_U_coefficient(/*order=*/4, /*beta=*/2.0, /*mu=*/1.0,
                                                      /*r=*/{1, 0}, SPIN_UP, SPIN_DOWN);
  double exact_coeff = -0.003025005435902937;
  EXPECT_NEAR(c, exact_coeff, 1e-12) << "Computed " << c << " vs exact " << exact_coeff << " (Δ = " << (c - exact_coeff) << ")";
}

TEST(CorrelatorAtomInfU, NearestNeighborOppositeSpinOrder6) {
  // ⟨n_{(1,0),↑} n_{0,↓}⟩ at U=∞, β=2, μ=1, t⁶ coefficient.
  // Reference: -0.0003127844534866856.
  double c           = compute_infinite_U_coefficient(/*order=*/6, /*beta=*/2.0, /*mu=*/1.0,
                                                      /*r=*/{1, 0}, SPIN_UP, SPIN_DOWN);
  double exact_coeff = -0.0003127844534866856;
  EXPECT_NEAR(c, exact_coeff, 1e-12) << "Computed " << c << " vs exact " << exact_coeff << " (Δ = " << (c - exact_coeff) << ")";
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}