#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include "../c++/sc_expansion/hubbard_solver.hpp"
#include "../c++/sc_expansion/cumulant.hpp"
#include "../c++/sc_expansion/dual.hpp"
#include "../c++/sc_expansion/args.hpp"

using namespace sc_expansion;

// Helper to call G0 with local lists
double call_G0_helper(HubbardSolver<1, double> const &solver, 
                      const std::vector<std::pair<double, int>>& unprimed, 
                      const std::vector<std::pair<double, int>>& primed, 
                      bool infinite_U = false) {
  std::vector<double> taus;
  std::vector<FermionOperator<1, double>> ops;
  
  // Unprimed are destruction (action 0)
  for (auto const& p : unprimed) {
    taus.push_back(p.first);
    ops.push_back(FermionOperator<1, double>(p.second)); // Action 0 | orbital
  }
  // Primed are creation (action 1)
  for (auto const& p : primed) {
    taus.push_back(p.first);
    ops.push_back(FermionOperator<1, double>((1 << 1) | p.second)); // Action 1 | orbital
  }

  Args<1, double> args(std::move(taus), std::move(ops));
  return infinite_U ? solver.G0n_infinite_U(args) : solver.G0n(args);
}

double compute_cumul_helper(HubbardSolver<1, double> const &solver, 
                            const std::vector<std::pair<double, int>>& unprimed, 
                            const std::vector<std::pair<double, int>>& primed, 
                            bool infinite_U = false) {
  std::vector<double> taus;
  std::vector<FermionOperator<1, double>> ops;
  for (auto const& p : unprimed) {
    taus.push_back(p.first);
    ops.push_back(FermionOperator<1, double>(p.second));
  }
  for (auto const& p : primed) {
    taus.push_back(p.first);
    ops.push_back(FermionOperator<1, double>((1 << 1) | p.second));
  }
  Args<1, double> args(std::move(taus), std::move(ops));
  return compute_cumulant_decomposition<1, double>(args, solver, infinite_U);
}

class HubbardAtomTest : public ::testing::Test {
  protected:
  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;
  Parameters<double> params{U, beta, mu};
  std::unique_ptr<HubbardSolver<1, double>> solver;

  void SetUp() override { solver = std::make_unique<HubbardSolver<1, double>>(params); }
};

TEST_F(HubbardAtomTest, CumulantOrderOneMatchesExactResult) {
  double tau = 0.5;
  std::vector<std::pair<double, int>> unprimed = {{tau, 0}};
  std::vector<std::pair<double, int>> primed   = {{0, 0}};

  double G01 = call_G0_helper(*solver, unprimed, primed);
  double C01 = compute_cumul_helper(*solver, unprimed, primed);

  EXPECT_NEAR(G01, C01, 1e-12);
}

TEST_F(HubbardAtomTest, CumulantOrderTwoMatchesExactResult) {
  std::vector<std::pair<double, int>> unprimed_args = {{0.5, 1}, {0.8, 0}};
  std::vector<std::pair<double, int>> primed_args   = {{0.0, 0}, {0.3, 1}};
  
  double G02 = call_G0_helper(*solver, unprimed_args, primed_args);

  double G0_11 = call_G0_helper(*solver, {unprimed_args[0]}, {primed_args[0]});
  double G0_22 = call_G0_helper(*solver, {unprimed_args[1]}, {primed_args[1]});
  double G0_12 = call_G0_helper(*solver, {unprimed_args[0]}, {primed_args[1]});
  double G0_21 = call_G0_helper(*solver, {unprimed_args[1]}, {primed_args[0]});

  double C02_exact = G02 - G0_11 * G0_22 + G0_12 * G0_21;
  double C02 = compute_cumul_helper(*solver, unprimed_args, primed_args);

  EXPECT_NEAR(C02, C02_exact, 1e-12);
}

TEST_F(HubbardAtomTest, CumulantOrderThreeMatchesExactResult) {
  std::vector<std::pair<double, int>> unprimed_args = {{0.2, 0}, {0.56, 0}, {0.32, 0}};
  std::vector<std::pair<double, int>> primed_args   = {{0.07, 0}, {0.262, 0}, {0.651, 0}};
  
  double G03 = call_G0_helper(*solver, unprimed_args, primed_args);

  auto C2 = [&](int i1, int i2, int j1, int j2) {
    return compute_cumul_helper(*solver, {unprimed_args[i1], unprimed_args[i2]}, {primed_args[j1], primed_args[j2]});
  };
  auto G1 = [&](int i, int j) {
    return call_G0_helper(*solver, {unprimed_args[i]}, {primed_args[j]});
  };

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
  double C03 = compute_cumul_helper(*solver, unprimed_args, primed_args, false);

  EXPECT_NEAR(C03, C03_exact, 1e-12);
}

TEST_F(HubbardAtomTest, SpinConservationOfCumulant) {
  std::vector<std::pair<double, int>> unprimed_args1 = {{0.5, 1}, {0.8, 1}};
  std::vector<std::pair<double, int>> primed_args1   = {{0.0, 0}, {0.3, 0}};
  double C02_1 = compute_cumul_helper(*solver, unprimed_args1, primed_args1);
  EXPECT_NEAR(C02_1, 0.0, 1e-12);

  std::vector<std::pair<double, int>> unprimed_args2 = {{0.5, 1}, {0.8, 0}};
  std::vector<std::pair<double, int>> primed_args2   = {{0.0, 0}, {0.3, 0}};
  double C02_2 = compute_cumul_helper(*solver, unprimed_args2, primed_args2);
  EXPECT_NEAR(C02_2, 0.0, 1e-12);
}

TEST_F(HubbardAtomTest, InfiniteUCumulantOrder2MatchesExactResult) {
  std::vector<std::pair<double, int>> unprimed_args = {{0.8, 0}, {0.4, 0}};
  std::vector<std::pair<double, int>> primed_args   = {{0.6, 0}, {0.2, 0}};
  
  double G02 = call_G0_helper(*solver, unprimed_args, primed_args, true);

  double G0_11 = call_G0_helper(*solver, {unprimed_args[0]}, {primed_args[0]}, true);
  double G0_22 = call_G0_helper(*solver, {unprimed_args[1]}, {primed_args[1]}, true);
  double G0_12 = call_G0_helper(*solver, {unprimed_args[0]}, {primed_args[1]}, true);
  double G0_21 = call_G0_helper(*solver, {unprimed_args[1]}, {primed_args[0]}, true);

  double C02_exact = G02 - G0_11 * G0_22 + G0_12 * G0_21;
  double C02 = compute_cumul_helper(*solver, unprimed_args, primed_args, true);

  EXPECT_NEAR(C02, C02_exact, 1e-12);
}

TEST(HubbardAtomDualTest, ParticleHoleSymmetryCumulant) {
  double U_val    = 8.0;
  double beta_val = 1.0;
  double delta    = 0.1;

  Dual mu1_val(U_val / 2.0 + delta, 1.0);
  Dual mu2_val(U_val / 2.0 - delta, 1.0);

  Parameters<Dual> params1{Dual(U_val, 0.0), Dual(beta_val, 0.0), mu1_val};
  Parameters<Dual> params2{Dual(U_val, 0.0), Dual(beta_val, 0.0), mu2_val};

  HubbardSolver<1, Dual> solver1(params1);
  HubbardSolver<1, Dual> solver2(params2);

  std::vector<std::pair<double, int>> u1 = {{0.8, 0}, {0.4, 1}};
  std::vector<std::pair<double, int>> p1 = {{0.6, 0}, {0.2, 1}};

  auto get_dual_cumul = [](auto const& s, auto const& u, auto const& p) {
    std::vector<double> taus;
    std::vector<FermionOperator<1, Dual>> ops;
    for (auto const& pair : u) { taus.push_back(pair.first); ops.push_back(FermionOperator<1, Dual>(pair.second)); }
    for (auto const& pair : p) { taus.push_back(pair.first); ops.push_back(FermionOperator<1, Dual>((1 << 1) | pair.second)); }
    Args<1, Dual> args(std::move(taus), std::move(ops));
    return compute_cumulant_decomposition<1, Dual>(args, s);
  };

  Dual c1 = get_dual_cumul(solver1, u1, p1);
  Dual c2 = get_dual_cumul(solver2, p1, u1);

  EXPECT_NEAR(c1.value, c2.value, 1e-12);
  EXPECT_NEAR(c1.derivative, -c2.derivative, 1e-12);
}

TEST(HubbardAtomDualTest, DualMu_VanishInNonInteractingLimit) {
  double U_val    = 0.0;
  double beta_val = 1.0;
  Dual mu_val(0.5, 1.0);

  Parameters<Dual> params{Dual(U_val, 0.0), Dual(beta_val, 0.0), mu_val};
  HubbardSolver<1, Dual> solver(params);

  std::vector<std::pair<double, int>> unprimed = {{0.8, 0}, {0.4, 1}};
  std::vector<std::pair<double, int>> primed   = {{0.6, 0}, {0.2, 1}};

  std::vector<double> taus;
  std::vector<FermionOperator<1, Dual>> ops;
  for (auto const& p : unprimed) { taus.push_back(p.first); ops.push_back(FermionOperator<1, Dual>(p.second)); }
  for (auto const& p : primed) { taus.push_back(p.first); ops.push_back(FermionOperator<1, Dual>((1 << 1) | p.second)); }
  Args<1, Dual> args(std::move(taus), std::move(ops));

  Dual c = compute_cumulant_decomposition<1, Dual>(args, solver);

  EXPECT_NEAR(c.value, 0.0, 1e-12);
  EXPECT_NEAR(c.derivative, 0.0, 1e-12);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
