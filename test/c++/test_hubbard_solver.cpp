#include <gtest/gtest.h>
#include <cmath>
#include "../c++/sc_expansion/hubbard_solver.hpp"
#include "../c++/sc_expansion/dual.hpp"
#include "../c++/sc_expansion/args.hpp"
#include <iostream>
#include <vector>
#include <memory>

using namespace sc_expansion;

class HubbardSolverTest : public ::testing::Test {
  protected:
  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;
  Parameters<double> params{U, beta, mu, true};

  std::unique_ptr<HubbardSolver<double>> solver;

  void SetUp() override { solver = std::make_unique<HubbardSolver<double>>(params); }
};

TEST_F(HubbardSolverTest, AtomG01FiniteU) {
  double tau = 0.5;

  // cdag_up(0.0), c_up(0.5)
  // Orbital index: up is 1, down is 0
  // Action: create is 1, destroy is 0
  // FermionOperator(uint8_t op) where Bit N_sites (bit 1 for N_sites=1) is action, Bit 0 is orbital index
  // cdag_up: (1 << 1) | 1 = 3
  // c_up:    (0 << 1) | 1 = 1

  std::vector<double> taus                    = {0.0, 0.5};
  std::vector<FermionOperator<double>> ops = {FermionOperator<double>(3), FermionOperator<double>(1)};

  Args<double> args(taus, ops);
  double g0 = solver->G0n(args);

  double Z_exact = 1 + 2 * std::exp(beta * mu) + std::exp(beta * (2 * mu - U));
  // Exact G0(tau|0) = < T cdag_up(0) c_up(tau) >
  // For 0 < tau < beta: G0(tau|0) = - < c_up(tau) cdag_up(0) > = - 1/Z Tr( e^{-(beta-tau)H} c_up e^{-tau H} cdag_up )
  double G0_exact = -1.0 / Z_exact * (std::exp(tau * mu) + std::exp(beta * mu) * std::exp(-tau * (U - mu)));

  EXPECT_NEAR(g0, G0_exact, 1e-10);
}

TEST_F(HubbardSolverTest, AtomG01MatchesG0n) {
  // Regression test: closed-form G01 must agree with the trace-loop G0n across
  // a range of (tau_1, tau_1') values and both input operator orderings.
  struct Case {
    double t1, t1p;
  };
  std::vector<Case> cases = {{0.0, 0.5}, {0.5, 0.0}, {0.2, 0.8}, {0.9, 0.1}, {0.3, 0.3001}, {0.4, 0.4}};

  for (auto const &c : cases) {
    // Input order [c^dag(t1'), c(t1)]
    {
      std::vector<double> taus                    = {c.t1p, c.t1};
      std::vector<FermionOperator<double>> ops = {FermionOperator<double>(3), FermionOperator<double>(1)};
      Args<double> args(taus, ops);
      double g0  = solver->G0n(args);
      double g01 = solver->G01(args);
      EXPECT_NEAR(g01, g0, 1e-12) << "mismatch at (t1,t1')=(" << c.t1 << "," << c.t1p << ") order [cdag,c]";
    }
    // Input order [c(t1), c^dag(t1')]
    {
      std::vector<double> taus                    = {c.t1, c.t1p};
      std::vector<FermionOperator<double>> ops = {FermionOperator<double>(1), FermionOperator<double>(3)};
      Args<double> args(taus, ops);
      double g0  = solver->G0n(args);
      double g01 = solver->G01(args);
      EXPECT_NEAR(g01, g0, 1e-12) << "mismatch at (t1,t1')=(" << c.t1 << "," << c.t1p << ") order [c,cdag]";
    }
  }
}

TEST_F(HubbardSolverTest, AtomG01DualMatchesG0n) {
  // Dual version: check both value and derivative match between G01 and G0n.
  Dual U_d(U, 0.0);
  Dual beta_d(beta, 0.0);
  Dual mu_d(mu, 1.0); // derivative wrt mu
  Parameters<Dual> params_d{U_d, beta_d, mu_d, true};
  HubbardSolver<Dual> solver_d(params_d);

  std::vector<std::pair<double, double>> cases = {{0.0, 0.5}, {0.5, 0.0}, {0.2, 0.8}, {0.9, 0.1}};
  for (auto const &[t1, t1p] : cases) {
    std::vector<double> taus                  = {t1p, t1};
    std::vector<FermionOperator<Dual>> ops = {FermionOperator<Dual>(3), FermionOperator<Dual>(1)};
    Args<Dual> args(taus, ops);
    Dual g0  = solver_d.G0n(args);
    Dual g01 = solver_d.G01(args);
    EXPECT_NEAR(g01.value, g0.value, 1e-12) << "value mismatch at (" << t1 << "," << t1p << ")";
    EXPECT_NEAR(g01.derivative, g0.derivative, 1e-12) << "derivative mismatch at (" << t1 << "," << t1p << ")";
  }
}

TEST_F(HubbardSolverTest, AtomG0DerivativeMu) {
  double tau = 0.4;

  Dual U_d(U, 0.0);
  Dual beta_d(beta, 0.0);
  Dual mu_d(mu, 1.0); // derivative wrt mu

  Parameters<Dual> params_d{U_d, beta_d, mu_d, true};
  HubbardSolver<Dual> solver_d(params_d);

  std::vector<double> taus                  = {0.0, tau};
  std::vector<FermionOperator<Dual>> ops = {FermionOperator<Dual>(3), FermionOperator<Dual>(1)};

  Dual g0 = solver_d.G0n(Args<Dual>(taus, ops));

  // Analytical derivative
  double exp_tau_mu    = std::exp(tau * mu);
  double exp_beta_mu   = std::exp(beta * mu);
  double exp_beta_2muU = std::exp(beta * (2 * mu - U));
  double exp_tau_Umu   = std::exp(-tau * (U - mu));

  double Z    = 1.0 + 2.0 * exp_beta_mu + exp_beta_2muU;
  double dZ_dmu = 2.0 * beta * exp_beta_mu + 2.0 * beta * exp_beta_2muU;

  double G_num = -(exp_tau_mu + exp_beta_mu * exp_tau_Umu);
  double dG_num_dmu = -(tau * exp_tau_mu + (beta + tau) * exp_beta_mu * exp_tau_Umu);

  double dG_dmu_analytical = dG_num_dmu / Z - (G_num / (Z * Z)) * dZ_dmu;

  EXPECT_NEAR(g0.value, G_num / Z, 1e-10);
  EXPECT_NEAR(g0.derivative, dG_dmu_analytical, 1e-10);
}

TEST_F(HubbardSolverTest, AtomPHSymmetryDouble) {
  double U_ph    = 8.0;
  double beta_ph = 1.0;
  double mu_ph   = U_ph / 2.0;
  Parameters<double> ph_params{U_ph, beta_ph, mu_ph, true};
  HubbardSolver<double> ph_solver(ph_params);

  double tau = 0.3;
  // G(tau) = - < T c_up(tau) cdag_up(0) >
  std::vector<double> taus1                    = {0.0, tau};
  std::vector<FermionOperator<double>> ops1 = {FermionOperator<double>(3), FermionOperator<double>(1)};
  double g_tau                                 = ph_solver.G0n(Args<double>(taus1, ops1));

  // G(beta-tau)
  std::vector<double> taus2 = {0.0, beta_ph - tau};
  double g_beta_minus_tau   = ph_solver.G0n(Args<double>(taus2, ops1));

  // For mu = U/2, G(tau) = G(beta-tau) (Hubbard atom is PH symmetric)
  EXPECT_NEAR(g_tau, g_beta_minus_tau, 1e-10);
}

TEST_F(HubbardSolverTest, AtomPHSymmetryDual) {
  double U_val    = 8.0;
  double beta_val = 1.0;
  double mu_val   = U_val / 2.0;

  Dual U(U_val, 0.0);
  Dual beta(beta_val, 0.0);
  Dual mu(mu_val, 1.0); // derivative wrt mu

  Parameters<Dual> ph_params{U, beta, mu, true};
  HubbardSolver<Dual> ph_solver(ph_params);

  double tau = 0.3;
  // G(tau) = - < T c_up(tau) cdag_up(0) >
  std::vector<double> taus1                  = {0.0, tau};
  std::vector<FermionOperator<Dual>> ops1 = {FermionOperator<Dual>(3), FermionOperator<Dual>(1)};
  Dual g_tau                                 = ph_solver.G0n(Args<Dual>(taus1, ops1));

  // G(beta-tau)
  std::vector<double> taus2 = {0.0, beta_val - tau};
  Dual g_beta_minus_tau     = ph_solver.G0n(Args<Dual>(taus2, ops1));

  EXPECT_NEAR(g_tau.value, g_beta_minus_tau.value, 1e-10);
  EXPECT_NEAR(g_tau.derivative, -g_beta_minus_tau.derivative, 1e-10);
}

TEST_F(HubbardSolverTest, AtomDensity) {
  // <cdag_up(0) c_up(0)>
  // In G0n, with taus = {0.0, 0.0}, sorting preserves order {cdag_up, c_up}
  // The loop runs i=1 (c_up), then i=0 (cdag_up), so it computes Tr(e^-bH cdag c)
  std::vector<double> taus                    = {0.0, 0.0};
  std::vector<FermionOperator<double>> ops = {FermionOperator<double>(3), FermionOperator<double>(1)};

  double n_up = solver->G0n(Args<double>(taus, ops));

  double Z_exact    = 1 + 2 * std::exp(beta * mu) + std::exp(beta * (2 * mu - U));
  double n_up_exact = (std::exp(beta * mu) + std::exp(beta * (2 * mu - U))) / Z_exact;

  EXPECT_NEAR(n_up, n_up_exact, 1e-10);
}

TEST_F(HubbardSolverTest, AtomG04InfiniteUForbidden) {
  // Case 1: Incorrect operator sequence for infinite U (hits state 3)
  // cdag_up(0.3), cdag_dn(0.2), c_dn(0.1), c_up(0.0)
  // Sorting: cdag_up, cdag_dn, c_dn, c_up
  // Bits for N_sites=1: cdag_up=3, cdag_dn=2, c_dn=0, c_up=1
  std::vector<double> taus                  = {0.3, 0.2, 0.1, 0.0};
  std::vector<FermionOperator<double>> ops = {FermionOperator<double>(3), FermionOperator<double>(2), FermionOperator<double>(0), FermionOperator<double>(1)};

  Args<double> args1(taus, ops);
  EXPECT_FALSE(args1.operator_sequence_is_valid_infinite_U());
  EXPECT_DOUBLE_EQ(solver->G0n_infinite_U(args1), 0.0);

  // Case 2: Valid for finite U but hits state 3 in infinite U
  // Sequence: c_up(0.3), cdag_up(0.2), cdag_dn(0.1), c_dn(0.0)
  // Start state |up>: |up> -> |0> -> |up> -> |up down> (FORBIDDEN) -> |up>
  std::vector<double> taus2                  = {0.3, 0.2, 0.1, 0.0};
  std::vector<FermionOperator<double>> ops2 = {FermionOperator<double>(1), FermionOperator<double>(3), FermionOperator<double>(2), FermionOperator<double>(0)};

  Args<double> args2(taus2, ops2);
  EXPECT_TRUE(args2.operator_sequence_is_valid());
  EXPECT_FALSE(args2.operator_sequence_is_valid_infinite_U());
  EXPECT_DOUBLE_EQ(solver->G0n_infinite_U(args2), 0.0);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
