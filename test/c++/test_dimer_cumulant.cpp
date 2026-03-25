#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include "../c++/sc_expansion/hubbard_solver.hpp"
#include "../c++/sc_expansion/cumulant.hpp"

using namespace sc_expansion;

// Helper to call G0 with Args for N_sites=2
template <int N_sites, typename T>
T call_G0_helper_dimer(HubbardSolver<N_sites, T> const &solver, Args<N_sites, T> const &unprimed, Args<N_sites, T> const &primed) {
  int n = unprimed.order;
  std::vector<double> taus;
  std::vector<FermionOperator<N_sites, T>> ops;
  taus.reserve(2 * n);
  ops.reserve(2 * n);

  for (int i = 0; i < n; ++i) {
    // Primed (cdag)
    taus.push_back(primed.taus[i]);
    ops.push_back(primed.ops[i]);
    // Unprimed (c)
    taus.push_back(unprimed.taus[i]);
    ops.push_back(unprimed.ops[i]);
  }

  Args<N_sites, T> args(std::move(taus), std::move(ops));
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
  Args<2, double> unprimed({tau}, {FermionOperator<2, double>(orbital)});
  Args<2, double> primed({0.0}, {FermionOperator<2, double>((1 << 2) | orbital)});

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

  Args<2, double> unprimed_args({0.5, 0.8}, {FermionOperator<2, double>(up), FermionOperator<2, double>(down)});
  Args<2, double> primed_args({0.3, 0.0}, {FermionOperator<2, double>((1 << 2) | up), FermionOperator<2, double>((1 << 2) | down)});

  double G02 = call_G0_helper_dimer(*solver, unprimed_args, primed_args);

  auto get_sub = [&](const Args<2, double>& a, std::vector<int> idxs) {
      std::vector<double> ts;
      std::vector<FermionOperator<2, double>> os;
      for (int i : idxs) {
          ts.push_back(a.taus[i]);
          os.push_back(a.ops[i]);
      }
      return Args<2, double>(ts, os);
  };

  // Exact 2nd order cumulant: C(1,2|1',2') = G(1,2|1',2') - G(1|1')G(2|2') + G(1|2')G(2|1')
  double G_11 = call_G0_helper_dimer(*solver, get_sub(unprimed_args, {0}), get_sub(primed_args, {0}));
  double G_22 = call_G0_helper_dimer(*solver, get_sub(unprimed_args, {1}), get_sub(primed_args, {1}));
  double G_12 = call_G0_helper_dimer(*solver, get_sub(unprimed_args, {0}), get_sub(primed_args, {1}));
  double G_21 = call_G0_helper_dimer(*solver, get_sub(unprimed_args, {1}), get_sub(primed_args, {0}));

  double C02_exact = G02 - G_11 * G_22 + G_12 * G_21;
  double C02       = compute_cumulant_decomposition<2, double>(unprimed_args, primed_args, *solver);

  EXPECT_NEAR(C02, C02_exact, 1e-12);
}

TEST_F(DimerCumulantTest, CumulantOrderThreeMatchesExactResult) {
  // 6-point cumulant on site 0, all spin up
  int up = 2;

  Args<2, double> unprimed_args({0.2, 0.56, 0.32}, {FermionOperator<2, double>(up), FermionOperator<2, double>(up), FermionOperator<2, double>(up)});
  Args<2, double> primed_args({0.07, 0.262, 0.651}, {FermionOperator<2, double>((1 << 2) | up), FermionOperator<2, double>((1 << 2) | up), FermionOperator<2, double>((1 << 2) | up)});

  double G03 = call_G0_helper_dimer(*solver, unprimed_args, primed_args);

  auto get_sub = [&](const Args<2, double>& a, std::vector<int> idxs) {
      std::vector<double> ts;
      std::vector<FermionOperator<2, double>> os;
      for (int i : idxs) {
          ts.push_back(a.taus[i]);
          os.push_back(a.ops[i]);
      }
      return Args<2, double>(ts, os);
  };

  // Manual decomposition for C03
  auto C1 = [&](int i, int j) { return compute_cumulant_decomposition<2, double>(get_sub(unprimed_args, {i}), get_sub(primed_args, {j}), *solver); };
  auto C2 = [&](int i1, int i2, int j1, int j2) {
    return compute_cumulant_decomposition<2, double>(get_sub(unprimed_args, {i1, i2}), get_sub(primed_args, {j1, j2}), *solver);
  };

  // Let's use the same exact calculation as in test_cumulant.cpp but adapted for dimer
  auto G1 = [&](int i, int j) { return call_G0_helper_dimer(*solver, get_sub(unprimed_args, {i}), get_sub(primed_args, {j})); };

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
