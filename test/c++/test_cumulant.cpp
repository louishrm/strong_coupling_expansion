#include <gtest/gtest.h>
#include <cmath>
#include "../c++/sc_expansion/hubbard_solver.hpp"
#include "../c++/sc_expansion/cumulant.hpp"
#include "../c++/sc_expansion/dual.hpp"

using namespace sc_expansion;

// Helper to call G0 with ArgList
double call_G0_helper(HubbardSolver<1, double> const &solver, ArgList const &unprimed, ArgList const &primed, bool infinite_U = false) {
  int n = unprimed.size();
  std::vector<double> taus;
  std::vector<FermionOperator<1, double>> ops;
  taus.reserve(2 * n);
  ops.reserve(2 * n);

  for (int i = 0; i < n; ++i) {
    // Primed (cdag)
    taus.push_back(primed[i].first);
    ops.push_back(FermionOperator<1, double>((1 << 1) | primed[i].second));
    // Unprimed (c)
    taus.push_back(unprimed[i].first);
    ops.push_back(FermionOperator<1, double>(unprimed[i].second));
  }

  Args<1, double> args(std::move(taus), std::move(ops));
  return infinite_U ? solver.G0n_infinite_U(args) : solver.G0n(args);
}

// The Test Fixture: Sets up the data common to all tests
class AtomCumulantTest : public ::testing::Test {
  protected:
  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;
  Parameters<double> params{U, beta, mu, 0.0, true};

  std::unique_ptr<HubbardSolver<1, double>> solver;

  void SetUp() override { solver = std::make_unique<HubbardSolver<1, double>>(params); }
};

TEST_F(AtomCumulantTest, CumulantOrderOneMatchesExactResult) {

  double tau       = 0.5;
  ArgList unprimed = {{tau, 0}};
  ArgList primed   = {{0, 0}};

  double G01 = call_G0_helper(*solver, unprimed, primed);
  double C01 = compute_cumulant_decomposition<1, double>(unprimed, primed, *solver);

  EXPECT_NEAR(G01, C01, 1e-12);
}

TEST_F(AtomCumulantTest, CumulantOrderTwoMatchesExactResult) {

  ArgList unprimed_args = {{0.5, 1}, {0.8, 0}};
  ArgList primed_args   = {{0.0, 0}, {0.3, 1}};
  double G02            = call_G0_helper(*solver, unprimed_args, primed_args); //G0(1,2|1',2')

  double G0_11 = call_G0_helper(*solver, {unprimed_args[0]}, {primed_args[0]}); //G(1|3)
  double G0_22 = call_G0_helper(*solver, {unprimed_args[1]}, {primed_args[1]}); //G(2|4)

  double G0_12 = call_G0_helper(*solver, {unprimed_args[0]}, {primed_args[1]}); //G(1|4)
  double G0_21 = call_G0_helper(*solver, {unprimed_args[1]}, {primed_args[0]}); //G(2|3)

  double C02_exact = G02 - G0_11 * G0_22 + G0_12 * G0_21;

  double C02 = compute_cumulant_decomposition<1, double>(unprimed_args, primed_args, *solver);

  EXPECT_NEAR(C02, C02_exact, 1e-12);
}

TEST_F(AtomCumulantTest, CumulantOrderThreeMatchesExactResult) {

  ArgList unprimed_args = {{0.2, 0}, {0.56, 0}, {0.32, 0}};
  ArgList primed_args   = {{0.07, 0}, {0.262, 0}, {0.651, 0}};
  double G03            = call_G0_helper(*solver, unprimed_args, primed_args); //G03

  double C12_12_C33 = -compute_cumulant_decomposition<1, double>({unprimed_args[0], unprimed_args[1]}, {primed_args[0], primed_args[1]}, *solver)
     * call_G0_helper(*solver, {unprimed_args[2]}, {primed_args[2]}); //-C(1,2|1',2')*G(3|3')

  double C12_13_C32 = compute_cumulant_decomposition<1, double>({unprimed_args[0], unprimed_args[1]}, {primed_args[0], primed_args[2]}, *solver)
     * call_G0_helper(*solver, {unprimed_args[2]}, {primed_args[1]}); //C(1,2|1',3')*G(3|2')

  double C12_23_C31 = -compute_cumulant_decomposition<1, double>({unprimed_args[0], unprimed_args[1]}, {primed_args[1], primed_args[2]}, *solver)
     * call_G0_helper(*solver, {unprimed_args[2]}, {primed_args[0]}); //-C(1,2|2',3')*G(3|1')

  double C13_12_C32 = compute_cumulant_decomposition<1, double>({unprimed_args[0], unprimed_args[2]}, {primed_args[0], primed_args[1]}, *solver)
     * call_G0_helper(*solver, {unprimed_args[1]}, {primed_args[2]}); //C(1,3|1',2')*G(2|3')

  double C13_23_C21 = compute_cumulant_decomposition<1, double>({unprimed_args[0], unprimed_args[2]}, {primed_args[1], primed_args[2]}, *solver)
     * call_G0_helper(*solver, {unprimed_args[1]}, {primed_args[0]}); //C(1,3|2',3')*G(2|1')

  double C13_13_C22 = -compute_cumulant_decomposition<1, double>({unprimed_args[0], unprimed_args[2]}, {primed_args[0], primed_args[2]}, *solver)
     * call_G0_helper(*solver, {unprimed_args[1]}, {primed_args[1]}); //-C(1,3|1',3')*G(2|2')

  double C23_12_C13 = -compute_cumulant_decomposition<1, double>({unprimed_args[1], unprimed_args[2]}, {primed_args[0], primed_args[1]}, *solver)
     * call_G0_helper(*solver, {unprimed_args[0]}, {primed_args[2]}); //-C(2,3|1',2')*G(1|3')

  double C23_23_C11 = -compute_cumulant_decomposition<1, double>({unprimed_args[1], unprimed_args[2]}, {primed_args[1], primed_args[2]}, *solver)
     * call_G0_helper(*solver, {unprimed_args[0]}, {primed_args[0]}); //-C(2,3|2',3')*G(1|1')

  double C23_13_C12 = compute_cumulant_decomposition<1, double>({unprimed_args[1], unprimed_args[2]}, {primed_args[0], primed_args[2]}, *solver)
     * call_G0_helper(*solver, {unprimed_args[0]}, {primed_args[1]}); //C(2,3|1',3')*G(1|2')

  double C02C01_terms = C12_12_C33 + C12_13_C32 + C12_23_C31 + C13_12_C32 + C13_23_C21 + C13_13_C22 + C23_12_C13 + C23_23_C11 + C23_13_C12;

  double C11_C22_C33 = -call_G0_helper(*solver, {unprimed_args[0]}, {primed_args[0]}) * call_G0_helper(*solver, {unprimed_args[1]}, {primed_args[1]})
     * call_G0_helper(*solver, {unprimed_args[2]}, {primed_args[2]}); //-G(1|1')*G(2|2')*G(3|3')

  double C11_C23_C32 = call_G0_helper(*solver, {unprimed_args[0]}, {primed_args[0]}) * call_G0_helper(*solver, {unprimed_args[1]}, {primed_args[2]})
     * call_G0_helper(*solver, {unprimed_args[2]}, {primed_args[1]}); //G(1|1')*G(2|3')*G(3|2')

  double C12_C21_C33 = call_G0_helper(*solver, {unprimed_args[0]}, {primed_args[1]}) * call_G0_helper(*solver, {unprimed_args[1]}, {primed_args[0]})
     * call_G0_helper(*solver, {unprimed_args[2]}, {primed_args[2]}); //G(1|2')*G(2|1')*G(3|3')

  double C12_C23_C31 = -call_G0_helper(*solver, {unprimed_args[0]}, {primed_args[1]}) * call_G0_helper(*solver, {unprimed_args[1]}, {primed_args[2]})
     * call_G0_helper(*solver, {unprimed_args[2]}, {primed_args[0]}); //  -G(1|2')*G(2|3')*G(3|1')

  double C13_C22_C31 = call_G0_helper(*solver, {unprimed_args[0]}, {primed_args[2]}) * call_G0_helper(*solver, {unprimed_args[1]}, {primed_args[1]})
     * call_G0_helper(*solver, {unprimed_args[2]}, {primed_args[0]}); //G(1|3')*G(2|2')*G(3|1')

  double C13_C21_C32 = -call_G0_helper(*solver, {unprimed_args[0]}, {primed_args[2]}) * call_G0_helper(*solver, {unprimed_args[1]}, {primed_args[0]})
     * call_G0_helper(*solver, {unprimed_args[2]}, {primed_args[1]}); //-G(1|3')*G(2|1')*G(3|2')

  double G1G1G1_terms = C11_C22_C33 + C11_C23_C32 + C12_C21_C33 + C12_C23_C31 + C13_C22_C31 + C13_C21_C32;

  double C03_exact = G03 + C02C01_terms + G1G1G1_terms;

  double C03 = compute_cumulant_decomposition<1, double>(unprimed_args, primed_args, *solver, false, true);

  EXPECT_NEAR(C03, C03_exact, 1e-12);
}

TEST_F(AtomCumulantTest, SpinConservationOfCumulant) {

  ArgList unprimed_args1 = {{0.5, 1}, {0.8, 1}};
  ArgList primed_args1   = {{0.0, 0}, {0.3, 0}};

  double C02_1 = compute_cumulant_decomposition<1, double>(unprimed_args1, primed_args1, *solver);

  EXPECT_NEAR(C02_1, 0.0, 1e-12);

  ArgList unprimed_args2 = {{0.5, 1}, {0.8, 0}};
  ArgList primed_args2   = {{0.0, 0}, {0.3, 0}};

  double C02_2 = compute_cumulant_decomposition<1, double>(unprimed_args2, primed_args2, *solver);
  EXPECT_NEAR(C02_2, 0.0, 1e-12);
}

TEST_F(AtomCumulantTest, InfiniteUCumulantOrder2MatchesExactResult) {

  ArgList unprimed_args = {{0.8, 0}, {0.4, 0}};
  ArgList primed_args   = {{0.6, 0}, {0.2, 0}};
  double G02            = call_G0_helper(*solver, unprimed_args, primed_args, true); //G0(1,2|1',2')

  double G0_11 = call_G0_helper(*solver, {unprimed_args[0]}, {primed_args[0]}, true); //G(1|3)
  double G0_22 = call_G0_helper(*solver, {unprimed_args[1]}, {primed_args[1]}, true); //G(2|4)

  double G0_12 = call_G0_helper(*solver, {unprimed_args[0]}, {primed_args[1]}, true); //G(1|4)
  double G0_21 = call_G0_helper(*solver, {unprimed_args[1]}, {primed_args[0]}, true); //G(2|3)

  double C02_exact = G02 - G0_11 * G0_22 + G0_12 * G0_21;

  double C02 = compute_cumulant_decomposition<1, double>(unprimed_args, primed_args, *solver, true);

  EXPECT_NEAR(C02, C02_exact, 1e-12);
}

TEST(AtomCumulantDualTest, ParticleHoleSymmetryCumulant) {
  double U_val    = 8.0;
  double beta_val = 1.0;
  double delta    = 0.1;

  Dual mu1_val(U_val / 2.0 + delta, 1.0);
  Dual mu2_val(U_val / 2.0 - delta, 1.0);

  Parameters<Dual> params1{Dual(U_val, 0.0), Dual(beta_val, 0.0), mu1_val, Dual(0.0, 0.0), true};
  Parameters<Dual> params2{Dual(U_val, 0.0), Dual(beta_val, 0.0), mu2_val, Dual(0.0, 0.0), true};

  HubbardSolver<1, Dual> solver1(params1);
  HubbardSolver<1, Dual> solver2(params2);

  ArgList u1 = {{0.8, 0}, {0.4, 1}};
  ArgList p1 = {{0.6, 0}, {0.2, 1}};

  Dual c1 = compute_cumulant_decomposition<1, Dual>(u1, p1, solver1);
  Dual c2 = compute_cumulant_decomposition<1, Dual>(p1, u1, solver2);

  EXPECT_NEAR(c1.value, c2.value, 1e-12);
  EXPECT_NEAR(c1.derivative, -c2.derivative, 1e-12);
}

TEST(AtomCumulantDualTest, DualMu_VanishInNonInteractingLimit) {
  double U_val    = 0.0;
  double beta_val = 1.0;
  Dual mu_val(0.5, 1.0); // mu = 0.5, dmu/dmu = 1

  Parameters<Dual> params{Dual(U_val, 0.0), Dual(beta_val, 0.0), mu_val, Dual(0.0, 0.0), true};
  HubbardSolver<1, Dual> solver(params);

  // 2-particle cumulant (4 operators)
  ArgList unprimed = {{0.8, 0}, {0.4, 1}};
  ArgList primed   = {{0.6, 0}, {0.2, 1}};

  Dual c = compute_cumulant_decomposition<1, Dual>(unprimed, primed, solver);

  // At U=0, the cumulant should be zero (Gaussian),
  // and its derivative w.r.t mu should also be zero.
  EXPECT_NEAR(c.value, 0.0, 1e-12);
  EXPECT_NEAR(c.derivative, 0.0, 1e-12);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
