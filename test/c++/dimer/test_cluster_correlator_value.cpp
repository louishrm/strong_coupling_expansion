#include <gtest/gtest.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include "../c++/sc_expansion/dimer/diagram.hpp"
#include "../c++/sc_expansion/dimer/sum_diagrams.hpp"

using namespace sc_expansion;
using sc_expansion::dimer::Diagram;
using sc_expansion::dimer::SumDiagrams;

// =============================================================================
//  Deterministic (serial, no-MPI) cluster value checks for the dimer
//  <n(r)n(0)> correlator, to confirm the rooted-weight fix (removal of the
//  spurious spin-dependent n_mark_orbit doubling in diagram.cpp):
//
//   * Opposite-spin on-site coefficients still match the exact finite-cluster
//     ED references used by test_mcmc_dimer_correlator.cpp. The fix never
//     touched the opp-spin channel (its n_mark_orbit was already 1), so this is
//     the regression guard that the edit left opp-spin intact — and it also
//     re-validates that the deterministic Simpson integral reproduces ED.
//
//   * Same-spin on-site obeys SU(2) (up,up == down,down) and the off-site value
//     is printed for an external ED cross-check (no hardcoded reference yet).
//
//  The integrand and normalisation mirror dump_dimer_correlator_diagnostic in
//  test_mcmc_dimer_correlator.cpp: coeff = ∫_{[0,β]^order} density_density(τ)
//  d^order τ, evaluated by composite Simpson — exactly what the MC ratio
//  estimator targets. The 3-dimer triangular cluster {(0,0),(1,0),(0,1)} with
//  pin_origin=true anchors n(0) at the single pendant reference site (0,0).
// =============================================================================

namespace {

  constexpr int SPIN_DOWN = 0, SPIN_UP = 1;

  double simpson_weight(int i, int n) { return (i == 0 || i == n) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0); }

  // Order-2 connected coefficient on the 3-dimer cluster (pin_origin), by Simpson
  // quadrature over [0,β]^2. n_grid must be even.
  double cluster_coeff_order2(double U, double mu, double beta, std::vector<int> const &r, int s1, int s2, int n_grid = 80) {
    double t                                           = 1.0; // intra-dimer hopping
    std::vector<std::pair<int, int>> cluster_positions = {{0, 0}, {1, 0}, {0, 1}};
    int n_cluster_sites                                = (int)cluster_positions.size();

    Parameters<double> params{U, beta, mu, t, /*bipartite=*/false};
    SumDiagrams<double> calc(params, /*order=*/2, r, s1, s2, cluster_positions, n_cluster_sites, /*pin_origin=*/true);

    double h        = beta / n_grid;
    double integral = 0.0;
    std::vector<double> taus(2);
    for (int i = 0; i <= n_grid; ++i) {
      for (int j = 0; j <= n_grid; ++j) {
        taus[0]   = i * h;
        taus[1]   = j * h;
        double sw = simpson_weight(i, n_grid) * simpson_weight(j, n_grid);
        calc.mark_all_dirty();
        integral += sw * calc.density_density(taus);
      }
    }
    return integral * (h / 3.0) * (h / 3.0);
  }

} // namespace

// -----------------------------------------------------------------------------
//  Opposite-spin on-site, order 2 — must match the exact ED references from
//  test_mcmc_dimer_correlator.cpp (OnSiteUpDownOrder2 / DiagnosticOnSiteOrder2).
//  Regression guard: the fix must not perturb the opposite-spin channel.
// -----------------------------------------------------------------------------
TEST(DimerClusterCorrelator, OppSpinOnSiteMatchesED_U12) {
  double coeff = cluster_coeff_order2(/*U=*/12.0, /*mu=*/3.0, /*beta=*/1.0, {0, 0}, SPIN_UP, SPIN_DOWN);
  double ed    = 0.000675371799; // exact finite-cluster ED (single-site, pinned)
  std::cout << std::setprecision(10) << "U=12 opp-spin on-site order-2: Simpson=" << coeff << "  ED=" << ed << "  rel=" << std::abs(coeff - ed) / ed
            << "\n";
  EXPECT_NEAR(coeff, ed, 0.02 * std::abs(ed)) << "opposite-spin on-site coefficient must still match ED after the fix";
}

TEST(DimerClusterCorrelator, OppSpinOnSiteMatchesED_U6) {
  double coeff = cluster_coeff_order2(/*U=*/6.0, /*mu=*/3.0, /*beta=*/1.0, {0, 0}, SPIN_UP, SPIN_DOWN);
  double ed    = 0.001334786422;
  std::cout << std::setprecision(10) << "U=6 opp-spin on-site order-2: Simpson=" << coeff << "  ED=" << ed << "  rel=" << std::abs(coeff - ed) / ed
            << "\n";
  EXPECT_NEAR(coeff, ed, 0.02 * std::abs(ed)) << "opposite-spin on-site coefficient must still match ED after the fix";
}

// -----------------------------------------------------------------------------
//  Same-spin channel (the one the fix changed). SU(2) requires up,up == down,down.
//  The off-site value is printed for an external ED cross-check.
// -----------------------------------------------------------------------------
TEST(DimerClusterCorrelator, SameSpinObeysSU2AndPrintsValues) {
  // On-site same-spin must be spin-flip symmetric.
  double uu = cluster_coeff_order2(12.0, 3.0, 1.0, {0, 0}, SPIN_UP, SPIN_UP);
  double dd = cluster_coeff_order2(12.0, 3.0, 1.0, {0, 0}, SPIN_DOWN, SPIN_DOWN);
  std::cout << std::setprecision(10) << "same-spin on-site order-2:  up,up=" << uu << "  down,down=" << dd << "\n";
  EXPECT_NEAR(uu, dd, 1e-9) << "SU(2): up,up and down,down on-site coefficients must be equal";

  // Off-site (intra-dimer neighbor) same-spin coefficient — print for an external
  // ED cross-check (compute the order-2 t-coefficient of <n_{(1,0),up} n_{(0,0),up}>_c
  // on the 3-dimer cluster from your Python HubbardED and compare).
  double off_uu = cluster_coeff_order2(12.0, 3.0, 1.0, {1, 0}, SPIN_UP, SPIN_UP);
  std::cout << std::setprecision(10) << "same-spin r=(1,0) order-2 (compare to ED): " << off_uu << "\n";
  // TODO(ED): replace with EXPECT_NEAR(off_uu, <ED value>, tol) once the reference
  // is computed. The opp-spin matches above already validate the Simpson machinery.
  EXPECT_TRUE(std::isfinite(off_uu));
}
