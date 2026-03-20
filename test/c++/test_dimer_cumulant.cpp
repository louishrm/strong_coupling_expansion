#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include "../c++/sc_expansion/hubbard_solver.hpp"
#include "../c++/sc_expansion/cumulant.hpp"

using namespace sc_expansion;

// Helper to call G0 with ArgList for N_sites=2
double call_G0_helper_dimer(HubbardSolver<2, double> const &solver, ArgList const &unprimed, ArgList const &primed) {
  int n = unprimed.size();
  std::vector<double> taus;
  std::vector<FermionOperator<2, double>> ops;
  taus.reserve(2 * n);
  ops.reserve(2 * n);

  for (int i = 0; i < n; ++i) {
    // Primed (cdag)
    taus.push_back(primed[i].first);
    ops.push_back(FermionOperator<2, double>((1 << 2) | primed[i].second));
    // Unprimed (c)
    taus.push_back(unprimed[i].first);
    ops.push_back(FermionOperator<2, double>(unprimed[i].second));
  }

  Args<2, double> args(std::move(taus), std::move(ops));
  return solver.G0n(args);
}

class DimerCumulantTest : public ::testing::Test {
  protected:
  double U    = 4.0;
  double beta = 1.0;
  double mu   = 2.0;
  double t    = 1.0;
  Parameters<double> params{U, beta, mu, t, true};

  std::unique_ptr<HubbardSolver<2, double>> solver;

  void SetUp() override { solver = std::make_unique<HubbardSolver<2, double>>(params); }
};

TEST_F(DimerCumulantTest, CumulantOrderOneMatchesExactResult) {
  // 2-point cumulant on site 0, spin up (index 2)
  // orbital indices for N_sites=2:
  // 0: site 0 down, 1: site 1 down, 2: site 0 up, 3: site 1 up
  int orbital = 2;

  double tau       = 0.5;
  ArgList unprimed = {{tau, orbital}};
  ArgList primed   = {{0.0, orbital}};

  double G01 = call_G0_helper_dimer(*solver, unprimed, primed);
  double C01 = compute_cumulant_decomposition<2, double>(unprimed, primed, *solver);

  EXPECT_NEAR(G01, C01, 1e-12);
}

TEST_F(DimerCumulantTest, CumulantOrderTwoMatchesExactResult) {
  // 4-point cumulant on site 0
  // unprimed: (tau1, up), (tau2, down)
  // primed:   (tau3, up), (tau4, down)
  int up   = 2;
  int down = 0;

  ArgList unprimed_args = {{0.5, up}, {0.8, down}};
  ArgList primed_args   = {{0.3, up}, {0.0, down}};

  double G02 = call_G0_helper_dimer(*solver, unprimed_args, primed_args);

  // Exact 2nd order cumulant: C(1,2|1',2') = G(1,2|1',2') - G(1|1')G(2|2') + G(1|2')G(2|1')
  double G_11 = call_G0_helper_dimer(*solver, {unprimed_args[0]}, {primed_args[0]});
  double G_22 = call_G0_helper_dimer(*solver, {unprimed_args[1]}, {primed_args[1]});
  double G_12 = call_G0_helper_dimer(*solver, {unprimed_args[0]}, {primed_args[1]});
  double G_21 = call_G0_helper_dimer(*solver, {unprimed_args[1]}, {primed_args[0]});

  double C02_exact = G02 - G_11 * G_22 + G_12 * G_21;
  double C02       = compute_cumulant_decomposition<2, double>(unprimed_args, primed_args, *solver);

  EXPECT_NEAR(C02, C02_exact, 1e-12);
}

TEST_F(DimerCumulantTest, CumulantOrderThreeMatchesExactResult) {
  // 6-point cumulant on site 0, all spin up
  int up = 2;

  ArgList unprimed_args = {{0.2, up}, {0.56, up}, {0.32, up}};
  ArgList primed_args   = {{0.07, up}, {0.262, up}, {0.651, up}};

  double G03 = call_G0_helper_dimer(*solver, unprimed_args, primed_args);

  // Manual decomposition for C03
  auto C1 = [&](int i, int j) { return compute_cumulant_decomposition<2, double>({unprimed_args[i]}, {primed_args[j]}, *solver); };
  auto C2 = [&](int i1, int i2, int j1, int j2) {
    return compute_cumulant_decomposition<2, double>({unprimed_args[i1], unprimed_args[i2]}, {primed_args[j1], primed_args[j2]}, *solver);
  };

  // C3 = G3 - sum(C2*C1) - sum(C1*C1*C1)
  // Term sum(C2*C1) has 3*3 = 9 terms
  double C2C1_terms = 0;
  // Partitions of {1,2,3}|{1',2',3'} into {a,b}|{a'} and {c}|{c'}
  // 12|1'2' * 3|3', 12|1'3' * 3|2', 12|2'3' * 3|1' ...
  // Signs must be consistent with compute_cumulant_decomposition

  // Actually, we can just use the definition G3 = C3 + sum(C2*C1) + sum(C1*C1*C1)
  // to verify compute_cumulant_decomposition returns the correct value.
  // The logic in test_cumulant.cpp for C03_exact is:
  // C03_exact = G03 + C02C01_terms + G1G1G1_terms
  // wait, the signs in test_cumulant.cpp are:
  // C03 = G03 - sum(C2*C1) - sum(C1*C1*C1) ? No, the recursive relation is G = sum(products of C).

  // Let's use the same exact calculation as in test_cumulant.cpp but adapted for dimer
  auto G1 = [&](int i, int j) { return call_G0_helper_dimer(*solver, {unprimed_args[i]}, {primed_args[j]}); };

  double C12_12_C33 = -C2(0, 1, 0, 1) * G1(2, 2);
  double C12_13_C32 = C2(0, 1, 0, 2) * G1(2, 1);
  double C12_23_C31 = -C2(0, 1, 1, 2) * G1(2, 0);

  double C13_12_C32 = C2(0, 2, 0, 1) * G1(1, 2);
  double C13_23_C21 = C2(0, 2, 1, 2) * G1(1, 0);
  double C13_13_C22 = -C2(0, 2, 0, 2) * G1(1, 1);

  double C23_12_C13 = -C2(1, 2, 0, 1) * G1(0, 2);
  double C23_23_C11 = -C2(1, 2, 1, 2) * G1(0, 0);
  double C23_13_C12 = C2(1, 2, 0, 2) * G1(0, 1);

  double C02C01_terms = C12_12_C33 + C12_13_C32 + C12_23_C31 + C13_12_C32 + C13_23_C21 + C13_13_C22 + C23_12_C13 + C23_23_C11 + C23_13_C12;

  double C11_C22_C33 = -G1(0, 0) * G1(1, 1) * G1(2, 2);
  double C11_C23_C32 = G1(0, 0) * G1(1, 2) * G1(2, 1);
  double C12_C21_C33 = G1(0, 1) * G1(1, 0) * G1(2, 2);
  double C12_C23_C31 = -G1(0, 1) * G1(1, 2) * G1(2, 0);
  double C13_C22_C31 = G1(0, 2) * G1(1, 1) * G1(2, 0);
  double C13_C21_C32 = -G1(0, 2) * G1(1, 0) * G1(2, 1);

  double G1G1G1_terms = C11_C22_C33 + C11_C23_C32 + C12_C21_C33 + C12_C23_C31 + C13_C22_C31 + C13_C21_C32;

  double C03_exact = G03 + C02C01_terms + G1G1G1_terms;
  double C03       = compute_cumulant_decomposition<2, double>(unprimed_args, primed_args, *solver);

  EXPECT_NEAR(C03, C03_exact, 1e-12);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
