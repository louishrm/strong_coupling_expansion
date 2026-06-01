#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>
#include "../c++/sc_expansion/hubbard_solver.hpp"
#include "../c++/sc_expansion/dual.hpp"
#include "../c++/sc_expansion/args.hpp"

using namespace sc_expansion;

// ============================================================================
// Dimer G0n_with_densities (Step 1 of the staggered-expansion density–density
// correlator). Mirrors test/c++/atom/test_hubbard_solver.cpp MINUS the
// infinite-U tests and MINUS the closed-form G01 tests (the dimer has no
// closed-form G01). The 16-state dimer spectrum has no compact closed form, so
// we check the trace-loop result against a direct spectral evaluation.
//
// Orbital encoding (N_sites=2): 0=site0↓, 1=site1↓, 2=site0↑, 3=site1↑.
// Fock state bit b is occupied ⇔ (state >> b) & 1.
// ============================================================================

namespace {

  constexpr int NS = HubbardSolver<2, double>::N_STATES; // 16

  template <typename T> FermionOperator<2, T> mk_op(int orbital, bool creation) {
    uint8_t id = static_cast<uint8_t>(orbital);
    if (creation) id |= FermionOperator<2, T>::ACTION_BIT;
    return FermionOperator<2, T>(id);
  }

  // Independent ED reference for a static (τ=0) product of number operators
  // ⟨Π_k n_{σ_k}⟩ = Σ_α (e^{-βE_α}/Z) Σ_{f∈α} |c_{α,f}|² Π_k occ_{σ_k}(f).
  // Built purely from eigenstate coefficients + Fock occupancy — it never
  // touches the operator-matrix machinery used by the implementation, so it is
  // a genuine cross-check.
  double ed_static_expectation(HubbardSolver<2, double> const &solver, std::vector<int> const &orbitals) {
    double result = 0.0;
    for (int a = 0; a < NS; ++a) {
      Eigenstate<double> const &es = solver.get_eigenstate(a);
      double diag                  = 0.0;
      for (auto const &[fock, coeff] : es.coefficients) {
        bool all_occ = true;
        for (int orb : orbitals) {
          if (!((fock >> orb) & 1)) {
            all_occ = false;
            break;
          }
        }
        if (all_occ) diag += coeff * coeff;
      }
      result += solver.get_exp_beta_E(a) * diag;
    }
    return result / solver.get_Z();
  }

  using Mat = std::array<std::array<double, NS>, NS>;

  Mat mat_zero() {
    Mat m{};
    return m;
  }

  Mat mat_mul(Mat const &A, Mat const &B) {
    Mat C = mat_zero();
    for (int i = 0; i < NS; ++i)
      for (int k = 0; k < NS; ++k) {
        if (A[i][k] == 0.0) continue;
        for (int j = 0; j < NS; ++j) C[i][j] += A[i][k] * B[k][j];
      }
    return C;
  }

  // ⟨T_τ Π_k n_{σ_k}(0) · (time-ordered hybridization ops)⟩ / Z via explicit
  // dense eigenbasis matrix products. `timed_ops` is a list of (op_id, τ) in
  // INPUT order (not yet τ-sorted); the fermionic time-ordering sign is the
  // sign of the permutation that sorts them into descending τ. The density
  // sits at τ=0 (right-most slot). This is a separate, hand-written evaluator
  // (dense products + explicit ordering) used to pin G0n_with_densities.
  double ed_time_ordered(HubbardSolver<2, double> const &solver, std::vector<std::pair<int, double>> const &timed_ops,
                         std::vector<int> const &density_orbitals) {
    std::array<double, NS> E{};
    for (int a = 0; a < NS; ++a) E[a] = solver.get_eigenstate(a).energy;

    auto dense_op = [&](int op_id, double tau) {
      Mat M                       = mat_zero();
      SparseMatrix<double> const &sm = solver.get_operator_matrix(static_cast<uint8_t>(op_id));
      for (auto const &e : sm.entries) M[e.row][e.col] = e.value * std::exp(tau * (E[e.row] - E[e.col]));
      return M;
    };

    // Density matrix N = Π_k n_{σ_k} in the eigenbasis (τ=0, no evolution).
    constexpr uint8_t ACTION = FermionOperator<2, double>::ACTION_BIT;
    Mat Nmat                 = mat_zero();
    for (int i = 0; i < NS; ++i) Nmat[i][i] = 1.0;
    for (int orb : density_orbitals) {
      Mat cdag = dense_op(ACTION | orb, 0.0);
      Mat c    = dense_op(orb, 0.0);
      Nmat     = mat_mul(mat_mul(cdag, c), Nmat);
    }

    // Sort ops into descending τ, counting inversions for the fermion sign.
    std::vector<int> order(timed_ops.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = static_cast<int>(i);
    int inversions = 0;
    for (size_t i = 0; i < timed_ops.size(); ++i)
      for (size_t j = i + 1; j < timed_ops.size(); ++j)
        if (timed_ops[i].second < timed_ops[j].second) inversions++;
    std::stable_sort(order.begin(), order.end(), [&](int x, int y) { return timed_ops[x].second > timed_ops[y].second; });
    double sign = (inversions % 2 == 0) ? 1.0 : -1.0;

    // Product (descending τ, left→right) · N.
    Mat prod = Nmat;
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
      auto const &[op_id, tau] = timed_ops[*it];
      prod                     = mat_mul(dense_op(op_id, tau), prod);
    }

    double trace = 0.0;
    for (int a = 0; a < NS; ++a) trace += solver.get_exp_beta_E(a) * prod[a][a];
    return sign * trace / solver.get_Z();
  }

} // namespace

class DimerHubbardSolverTest : public ::testing::Test {
  protected:
  double U    = 4.0;
  double beta = 1.0;
  double mu   = 2.0;
  double t    = 1.0;
  Parameters<double> params{U, beta, mu, t, true};
  std::unique_ptr<HubbardSolver<2, double>> solver;
  void SetUp() override { solver = std::make_unique<HubbardSolver<2, double>>(params); }
};

TEST_F(DimerHubbardSolverTest, DensityDecoratedEmptyOrbitalsMatchesG0n) {
  // {} decoration must recover plain G0n exactly, with and without legs.
  {
    Args<2, double> empty_args({}, {});
    EXPECT_NEAR(solver->G0n_with_densities(empty_args, {}), solver->G0n(empty_args), 1e-12);
  }
  {
    // One hybridization pair c†_{site0↑}(0.2) c_{site0↑}(0.7).
    std::vector<double> taus = {0.2, 0.7};
    std::vector<FermionOperator<2, double>> ops = {mk_op<double>(2, true), mk_op<double>(2, false)};
    Args<2, double> args(taus, ops);
    EXPECT_NEAR(solver->G0n_with_densities(args, {}), solver->G0n(args), 1e-12);
  }
}

TEST_F(DimerHubbardSolverTest, DensityDecoratedSingleMatchesOccupation) {
  // ⟨n_σ⟩ via G0n_with_densities({}, {σ}) must match compute_n_sigma(σ).
  Args<2, double> empty_args({}, {});
  for (int sigma = 0; sigma < 4; ++sigma) {
    double via_dec = solver->G0n_with_densities(empty_args, {sigma});
    double via_ref = solver->compute_n_sigma(sigma);
    EXPECT_NEAR(via_dec, via_ref, 1e-12) << "sigma=" << sigma;
  }
}

TEST_F(DimerHubbardSolverTest, DensityDecoratedSameSpinIdempotent) {
  // n_σ² = n_σ ⇒ {σ,σ} == {σ}.
  Args<2, double> empty_args({}, {});
  for (int sigma = 0; sigma < 4; ++sigma) {
    double squared = solver->G0n_with_densities(empty_args, {sigma, sigma});
    double single  = solver->compute_n_sigma(sigma);
    EXPECT_NEAR(squared, single, 1e-12) << "sigma=" << sigma;
  }
}

TEST_F(DimerHubbardSolverTest, DensityDecoratedOnSiteDoubleOcc) {
  // On-site double occupancy ⟨n_{i↑} n_{i↓}⟩ vs the independent ED reference.
  // site0: {↓=0, ↑=2}; site1: {↓=1, ↑=3}.
  std::vector<std::pair<int, int>> onsite = {{0, 2}, {1, 3}};
  Args<2, double> empty_args({}, {});
  for (auto const &[dn, up] : onsite) {
    double via_dec = solver->G0n_with_densities(empty_args, {dn, up});
    double via_ref = ed_static_expectation(*solver, {dn, up});
    EXPECT_NEAR(via_dec, via_ref, 1e-12) << "site dn=" << dn << " up=" << up;
  }
}

TEST_F(DimerHubbardSolverTest, DensityDecoratedInterSite) {
  // One density per site ⟨n_{0σ} n_{1σ'}⟩ — the intra-dimer NN correlator the
  // expansion ultimately needs. Pins the τ=0 matrix insertion across sites.
  std::vector<std::pair<int, int>> inter = {{0, 1}, {0, 3}, {2, 1}, {2, 3}};
  Args<2, double> empty_args({}, {});
  for (auto const &[s0, s1] : inter) {
    double via_dec = solver->G0n_with_densities(empty_args, {s0, s1});
    double via_ref = ed_static_expectation(*solver, {s0, s1});
    EXPECT_NEAR(via_dec, via_ref, 1e-12) << "orbs=(" << s0 << "," << s1 << ")";
  }
}

TEST_F(DimerHubbardSolverTest, DensityDecoratedWithHybridizationPair) {
  // With one (c†_σ(τ_p), c_σ(τ_u)) pair, G0n_with_densities(args, {σ_dens})
  // equals ⟨T_τ n_{σ_dens}(0) c†(τ_p) c(τ_u)⟩ (validates τ-ordering + matrix
  // insert together). Legs on site0↑ (orbital 2).
  int const leg_sigma = 2;
  double const tau_p = 0.2, tau_u = 0.7;
  std::vector<double> taus                    = {tau_p, tau_u};
  std::vector<FermionOperator<2, double>> ops = {mk_op<double>(leg_sigma, true), mk_op<double>(leg_sigma, false)};
  Args<2, double> pair_args(taus, ops);

  std::vector<std::pair<int, double>> timed = {{static_cast<int>(mk_op<double>(leg_sigma, true).op), tau_p},
                                               {static_cast<int>(mk_op<double>(leg_sigma, false).op), tau_u}};

  for (int sigma_dens = 0; sigma_dens < 4; ++sigma_dens) {
    double via_dec = solver->G0n_with_densities(pair_args, {sigma_dens});
    double via_ref = ed_time_ordered(*solver, timed, {sigma_dens});
    EXPECT_NEAR(via_dec, via_ref, 1e-12) << "sigma_dens=" << sigma_dens;
  }

  // Also check the reversed time ordering τ_u < τ_p path.
  double const tau_p2 = 0.7, tau_u2 = 0.2;
  std::vector<double> taus2                    = {tau_p2, tau_u2};
  std::vector<FermionOperator<2, double>> ops2 = {mk_op<double>(leg_sigma, true), mk_op<double>(leg_sigma, false)};
  Args<2, double> pair_args2(taus2, ops2);
  std::vector<std::pair<int, double>> timed2 = {{static_cast<int>(mk_op<double>(leg_sigma, true).op), tau_p2},
                                                {static_cast<int>(mk_op<double>(leg_sigma, false).op), tau_u2}};
  for (int sigma_dens = 0; sigma_dens < 4; ++sigma_dens) {
    double via_dec = solver->G0n_with_densities(pair_args2, {sigma_dens});
    double via_ref = ed_time_ordered(*solver, timed2, {sigma_dens});
    EXPECT_NEAR(via_dec, via_ref, 1e-12) << "reversed sigma_dens=" << sigma_dens;
  }
}

TEST_F(DimerHubbardSolverTest, DensityDecoratedDualDerivativeMu) {
  // Dual value + ∂/∂μ of G0n_with_densities against a numerically
  // differentiated ED reference, confirming the matrix path is Dual-clean.
  int const leg_sigma = 2;
  double const tau_p = 0.2, tau_u = 0.7;

  Dual U_d(U, 0.0), beta_d(beta, 0.0), mu_d(mu, 1.0), t_d(t, 0.0);
  Parameters<Dual> params_d{U_d, beta_d, mu_d, t_d, true};
  HubbardSolver<2, Dual> solver_d(params_d);

  std::vector<double> taus                  = {tau_p, tau_u};
  std::vector<FermionOperator<2, Dual>> ops = {mk_op<Dual>(leg_sigma, true), mk_op<Dual>(leg_sigma, false)};
  Args<2, Dual> pair_args(taus, ops);

  std::vector<std::pair<int, double>> timed = {{static_cast<int>(mk_op<double>(leg_sigma, true).op), tau_p},
                                               {static_cast<int>(mk_op<double>(leg_sigma, false).op), tau_u}};

  double const h = 1e-6;
  for (int sigma_dens = 0; sigma_dens < 4; ++sigma_dens) {
    Dual via_dec = solver_d.G0n_with_densities(pair_args, {sigma_dens});

    // Numerical ∂/∂μ via central difference of the double ED reference.
    Parameters<double> p_plus{U, beta, mu + h, t, true};
    Parameters<double> p_minus{U, beta, mu - h, t, true};
    HubbardSolver<2, double> s_plus(p_plus), s_minus(p_minus);
    double d_num = (ed_time_ordered(s_plus, timed, {sigma_dens}) - ed_time_ordered(s_minus, timed, {sigma_dens})) / (2.0 * h);

    double value = ed_time_ordered(*solver, timed, {sigma_dens});
    EXPECT_NEAR(via_dec.value, value, 1e-12) << "value sigma_dens=" << sigma_dens;
    EXPECT_NEAR(via_dec.derivative, d_num, 1e-6) << "deriv sigma_dens=" << sigma_dens;
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
