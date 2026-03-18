#include <gtest/gtest.h>
#include <cmath>
#include "../c++/sc_expansion/hubbard_solver.hpp"
#include "../c++/sc_expansion/cumulant.hpp"
#include "../c++/sc_expansion/dual.hpp"

using namespace sc_expansion;

// Helper to call G0 with separated unprimed/primed args and build the interleaved format
double call_G0_helper(HubbardAtom<double> const &atom, std::vector<std::pair<double, int>> const &unprimed, std::vector<std::pair<double, int>> const &primed, bool infinite_U = false) {
  int n = unprimed.size();
  std::vector<double> taus(2 * n);
  std::vector<int> ops(2 * n);
  for (int i = 0; i < n; ++i) {
    taus[2 * i]      = primed[i].first;
    ops[2 * i]       = (primed[i].second << 1) | 1; // action 1
    taus[2 * i + 1]  = unprimed[i].first;
    ops[2 * i + 1]   = (unprimed[i].second << 1) | 0; // action 0
  }
  return infinite_U ? atom.G0_infinite_U(taus, ops) : atom.G0(taus, ops);
}

// Wrapper to match old API for cumulant
template <typename T>
T compute_cumulant_decomposition_wrapper(std::vector<std::pair<double, int>> const &unprimed, std::vector<std::pair<double, int>> const &primed, HubbardAtom<T> const &atom, bool infinite_U = false, bool verbose = false) {
  int n = unprimed.size();
  std::vector<double> taus(2 * n);
  std::vector<int> ops(2 * n);
  for (int i = 0; i < n; ++i) {
    taus[2 * i]      = primed[i].first;
    ops[2 * i]       = (primed[i].second << 1) | 1; // action 1
    taus[2 * i + 1]  = unprimed[i].first;
    ops[2 * i + 1]   = (unprimed[i].second << 1) | 0; // action 0
  }
  return compute_cumulant_decomposition(taus, ops, atom, infinite_U, verbose);
}


// The Test Fixture: Sets up the data common to all tests
class HubbardAtomTest : public ::testing::Test {
  protected:
  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;
  Parameters<double> params{U, beta, mu, 0.0, true};

  std::unique_ptr<HubbardAtom<double>> atom;

  void SetUp() override { atom = std::make_unique<HubbardAtom<double>>(params); }
};

TEST_F(HubbardAtomTest, CumulantOrderOneMatchesExactResult) {

  double tau       = 0.5;
  std::vector<std::pair<double, int>> unprimed = {{tau, 0}};
  std::vector<std::pair<double, int>> primed   = {{0, 0}};

  double G01 = call_G0_helper(*atom, unprimed, primed);
  double C01 = compute_cumulant_decomposition_wrapper<double>(unprimed, primed, *atom);

  EXPECT_NEAR(G01, C01, 1e-12);
}

TEST_F(HubbardAtomTest, CumulantOrderTwoMatchesExactResult) {

  std::vector<std::pair<double, int>> unprimed_args = {{0.5, 1}, {0.8, 0}};
  std::vector<std::pair<double, int>> primed_args   = {{0.0, 0}, {0.3, 1}};
  double G02            = call_G0_helper(*atom, unprimed_args, primed_args); //G0(1,2|1',2')

  double G0_11 = call_G0_helper(*atom, {unprimed_args[0]}, {primed_args[0]}); //G(1|3)
  double G0_22 = call_G0_helper(*atom, {unprimed_args[1]}, {primed_args[1]}); //G(2|4)

  double G0_12 = call_G0_helper(*atom, {unprimed_args[0]}, {primed_args[1]}); //G(1|4)
  double G0_21 = call_G0_helper(*atom, {unprimed_args[1]}, {primed_args[0]}); //G(2|3)

  double C02_exact = G02 - G0_11 * G0_22 + G0_12 * G0_21;

  double C02 = compute_cumulant_decomposition_wrapper<double>(unprimed_args, primed_args, *atom);

  EXPECT_NEAR(C02, C02_exact, 1e-12);
}

TEST_F(HubbardAtomTest, CumulantOrderThreeMatchesExactResult) {

  std::vector<std::pair<double, int>> unprimed_args = {{0.2, 0}, {0.56, 0}, {0.32, 0}};
  std::vector<std::pair<double, int>> primed_args   = {{0.07, 0}, {0.262, 0}, {0.651, 0}};
  double G03            = call_G0_helper(*atom, unprimed_args, primed_args); //G03

  double C12_12_C33 = -compute_cumulant_decomposition_wrapper<double>({unprimed_args[0], unprimed_args[1]}, {primed_args[0], primed_args[1]}, *atom)
     * call_G0_helper(*atom, {unprimed_args[2]}, {primed_args[2]}); //-C(1,2|1',2')*G(3|3')

  double C12_13_C32 = compute_cumulant_decomposition_wrapper<double>({unprimed_args[0], unprimed_args[1]}, {primed_args[0], primed_args[2]}, *atom)
     * call_G0_helper(*atom, {unprimed_args[2]}, {primed_args[1]}); //C(1,2|1',3')*G(3|2')

  double C12_23_C31 = -compute_cumulant_decomposition_wrapper<double>({unprimed_args[0], unprimed_args[1]}, {primed_args[1], primed_args[2]}, *atom)
     * call_G0_helper(*atom, {unprimed_args[2]}, {primed_args[0]}); //-C(1,2|2',3')*G(3|1')

  double C13_12_C32 = compute_cumulant_decomposition_wrapper<double>({unprimed_args[0], unprimed_args[2]}, {primed_args[0], primed_args[1]}, *atom)
     * call_G0_helper(*atom, {unprimed_args[1]}, {primed_args[2]}); //C(1,3|1',2')*G(2|3')

  double C13_23_C21 = compute_cumulant_decomposition_wrapper<double>({unprimed_args[0], unprimed_args[2]}, {primed_args[1], primed_args[2]}, *atom)
     * call_G0_helper(*atom, {unprimed_args[1]}, {primed_args[0]}); //C(1,3|2',3')*G(2|1')

  double C13_13_C22 = -compute_cumulant_decomposition_wrapper<double>({unprimed_args[0], unprimed_args[2]}, {primed_args[0], primed_args[2]}, *atom)
     * call_G0_helper(*atom, {unprimed_args[1]}, {primed_args[1]}); //-C(1,3|1',3')*G(2|2')

  double C23_12_C13 = -compute_cumulant_decomposition_wrapper<double>({unprimed_args[1], unprimed_args[2]}, {primed_args[0], primed_args[1]}, *atom)
     * call_G0_helper(*atom, {unprimed_args[0]}, {primed_args[2]}); //-C(2,3|1',2')*G(1|3')

  double C23_23_C11 = -compute_cumulant_decomposition_wrapper<double>({unprimed_args[1], unprimed_args[2]}, {primed_args[1], primed_args[2]}, *atom)
     * call_G0_helper(*atom, {unprimed_args[0]}, {primed_args[0]}); //-C(2,3|2',3')*G(1|1')

  double C23_13_C12 = compute_cumulant_decomposition_wrapper<double>({unprimed_args[1], unprimed_args[2]}, {primed_args[0], primed_args[2]}, *atom)
     * call_G0_helper(*atom, {unprimed_args[0]}, {primed_args[1]}); //C(2,3|1',3')*G(1|2')

  double C02C01_terms = C12_12_C33 + C12_13_C32 + C12_23_C31 + C13_12_C32 + C13_23_C21 + C13_13_C22 + C23_12_C13 + C23_23_C11 + C23_13_C12;

  double C11_C22_C33 = -call_G0_helper(*atom, {unprimed_args[0]}, {primed_args[0]}) * call_G0_helper(*atom, {unprimed_args[1]}, {primed_args[1]})
     * call_G0_helper(*atom, {unprimed_args[2]}, {primed_args[2]}); //-G(1|1')*G(2|2')*G(3|3')

  double C11_C23_C32 = call_G0_helper(*atom, {unprimed_args[0]}, {primed_args[0]}) * call_G0_helper(*atom, {unprimed_args[1]}, {primed_args[2]})
     * call_G0_helper(*atom, {unprimed_args[2]}, {primed_args[1]}); //G(1|1')*G(2|3')*G(3|2')

  double C12_C21_C33 = call_G0_helper(*atom, {unprimed_args[0]}, {primed_args[1]}) * call_G0_helper(*atom, {unprimed_args[1]}, {primed_args[0]})
     * call_G0_helper(*atom, {unprimed_args[2]}, {primed_args[2]}); //G(1|2')*G(2|1')*G(3|3')

  double C12_C23_C31 = -call_G0_helper(*atom, {unprimed_args[0]}, {primed_args[1]}) * call_G0_helper(*atom, {unprimed_args[1]}, {primed_args[2]})
     * call_G0_helper(*atom, {unprimed_args[2]}, {primed_args[0]}); //  -G(1|2')*G(2|3')*G(3|1')

  double C13_C22_C31 = call_G0_helper(*atom, {unprimed_args[0]}, {primed_args[2]}) * call_G0_helper(*atom, {unprimed_args[1]}, {primed_args[1]})
     * call_G0_helper(*atom, {unprimed_args[2]}, {primed_args[0]}); //G(1|3')*G(2|2')*G(3|1')

  double C13_C21_C32 = -call_G0_helper(*atom, {unprimed_args[0]}, {primed_args[2]}) * call_G0_helper(*atom, {unprimed_args[1]}, {primed_args[0]})
     * call_G0_helper(*atom, {unprimed_args[2]}, {primed_args[1]}); //-G(1|3')*G(2|1')*G(3|2')

  double G1G1G1_terms = C11_C22_C33 + C11_C23_C32 + C12_C21_C33 + C12_C23_C31 + C13_C22_C31 + C13_C21_C32;

  double C03_exact = G03 + C02C01_terms + G1G1G1_terms;

  double C03 = compute_cumulant_decomposition_wrapper<double>(unprimed_args, primed_args, *atom, false, true);

  EXPECT_NEAR(C03, C03_exact, 1e-12);
}

TEST_F(HubbardAtomTest, SpinConservationOfCumulant) {

  std::vector<std::pair<double, int>> unprimed_args1 = {{0.5, 1}, {0.8, 1}};
  std::vector<std::pair<double, int>> primed_args1   = {{0.0, 0}, {0.3, 0}};

  double C02_1 = compute_cumulant_decomposition_wrapper<double>(unprimed_args1, primed_args1, *atom);

  EXPECT_NEAR(C02_1, 0.0, 1e-12);

  std::vector<std::pair<double, int>> unprimed_args2 = {{0.5, 1}, {0.8, 0}};
  std::vector<std::pair<double, int>> primed_args2   = {{0.0, 0}, {0.3, 0}};

  double C02_2 = compute_cumulant_decomposition_wrapper<double>(unprimed_args2, primed_args2, *atom);
  EXPECT_NEAR(C02_2, 0.0, 1e-12);
}

TEST_F(HubbardAtomTest, InfiniteUCumulantOrder2MatchesExactResult) {

  std::vector<std::pair<double, int>> unprimed_args = {{0.8, 0}, {0.4, 0}};
  std::vector<std::pair<double, int>> primed_args   = {{0.6, 0}, {0.2, 0}};
  double G02            = call_G0_helper(*atom, unprimed_args, primed_args, true); //G0(1,2|1',2')

  double G0_11 = call_G0_helper(*atom, {unprimed_args[0]}, {primed_args[0]}, true); //G(1|3)
  double G0_22 = call_G0_helper(*atom, {unprimed_args[1]}, {primed_args[1]}, true); //G(2|4)

  double G0_12 = call_G0_helper(*atom, {unprimed_args[0]}, {primed_args[1]}, true); //G(1|4)
  double G0_21 = call_G0_helper(*atom, {unprimed_args[1]}, {primed_args[0]}, true); //G(2|3)

  double C02_exact = G02 - G0_11 * G0_22 + G0_12 * G0_21;

  double C02 = compute_cumulant_decomposition_wrapper<double>(unprimed_args, primed_args, *atom, true);

  EXPECT_NEAR(C02, C02_exact, 1e-12);
}

TEST(HubbardAtomDualTest, ParticleHoleSymmetryCumulant) {
  double U_val    = 8.0;
  double beta_val = 1.0;
  double delta    = 0.1;

  Dual mu1_val(U_val / 2.0 + delta, 1.0);
  Dual mu2_val(U_val / 2.0 - delta, 1.0);

  Parameters<Dual> params1{Dual(U_val, 0.0), Dual(beta_val, 0.0), mu1_val, Dual(0.0, 0.0), true};
  Parameters<Dual> params2{Dual(U_val, 0.0), Dual(beta_val, 0.0), mu2_val, Dual(0.0, 0.0), true};

  HubbardAtom<Dual> atom1(params1);
  HubbardAtom<Dual> atom2(params2);

  std::vector<std::pair<double, int>> u1 = {{0.8, 0}, {0.4, 1}};
  std::vector<std::pair<double, int>> p1 = {{0.6, 0}, {0.2, 1}};

  Dual c1 = compute_cumulant_decomposition_wrapper<Dual>(u1, p1, atom1);
  Dual c2 = compute_cumulant_decomposition_wrapper<Dual>(p1, u1, atom2);

  EXPECT_NEAR(c1.value, c2.value, 1e-12);
  EXPECT_NEAR(c1.derivative, -c2.derivative, 1e-12);
}

TEST(HubbardAtomDualTest, DualMu_VanishInNonInteractingLimit) {
  double U_val    = 0.0;
  double beta_val = 1.0;
  Dual mu_val(0.5, 1.0); // mu = 0.5, dmu/dmu = 1

  Parameters<Dual> params{Dual(U_val, 0.0), Dual(beta_val, 0.0), mu_val, Dual(0.0, 0.0), true};
  HubbardAtom<Dual> atom(params);

  // 2-particle cumulant (4 operators)
  std::vector<std::pair<double, int>> unprimed = {{0.8, 0}, {0.4, 1}};
  std::vector<std::pair<double, int>> primed   = {{0.6, 0}, {0.2, 1}};

  Dual c = compute_cumulant_decomposition_wrapper<Dual>(unprimed, primed, atom);

  // At U=0, the cumulant should be zero (Gaussian),
  // and its derivative w.r.t mu should also be zero.
  EXPECT_NEAR(c.value, 0.0, 1e-12);
  EXPECT_NEAR(c.derivative, 0.0, 1e-12);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
