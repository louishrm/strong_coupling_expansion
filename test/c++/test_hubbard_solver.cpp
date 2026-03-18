#include <gtest/gtest.h>
#include <cmath>
#include "../c++/sc_expansion/hubbard_solver.hpp"
#include "../c++/sc_expansion/dual.hpp"
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace sc_expansion;

// --- HubbardAtom Tests ---

class HubbardAtomTest : public ::testing::Test {
  protected:
  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;
  Parameters<double> params{U, beta, mu, 0.0, true};

  std::unique_ptr<HubbardAtom<double>> atom;

  void SetUp() override { atom = std::make_unique<HubbardAtom<double>>(params); }
};

TEST_F(HubbardAtomTest, ConstructorStateEnergies) {
  // E = {0, -mu, -mu, 2*mu - U}
  EXPECT_DOUBLE_EQ(atom->E[0], 0.0);
  EXPECT_DOUBLE_EQ(atom->E[1], -2.0);
  EXPECT_DOUBLE_EQ(atom->E[2], -2.0);
  EXPECT_DOUBLE_EQ(atom->E[3], 8.0 - 2 * 2.0); // -4.0
}

TEST_F(HubbardAtomTest, PartitionFunctions) {
  double Z_expected          = 1 + 2 * std::exp(beta * mu) + std::exp(beta * (2 * mu - U));
  double Z_infinite_expected = 1 + 2 * std::exp(beta * mu);

  EXPECT_NEAR(atom->Z, Z_expected, 1e-10);
  EXPECT_NEAR(atom->Z_infinite_U, Z_infinite_expected, 1e-10);
}

TEST_F(HubbardAtomTest, LookupTableVerify) {

  // op = 2 (c_up: Action=0, Spin=1). state 0 -> state 0 + op 2 = index 2
  EXPECT_EQ(HubbardAtom<double>::lookup_table[(0 << 2) | 2].connected_state, 1);
  EXPECT_DOUBLE_EQ(HubbardAtom<double>::lookup_table[(0 << 2) | 2].matrix_element, 1.0);

  // op = 3 (cdag_up: Action=1, Spin=1). state 1 (|dn>) + op 3 = index 7
  EXPECT_EQ(HubbardAtom<double>::lookup_table[(1 << 2) | 3].connected_state, 3);
  EXPECT_DOUBLE_EQ(HubbardAtom<double>::lookup_table[(1 << 2) | 3].matrix_element, -1.0);

  // op = 0 (c_dn: Action=0, Spin=0). state 1 (|dn>) + op 0 = index 4
  EXPECT_EQ(HubbardAtom<double>::lookup_table[(1 << 2) | 0].connected_state, 0);
  EXPECT_DOUBLE_EQ(HubbardAtom<double>::lookup_table[(1 << 2) | 0].matrix_element, 1.0);

  // op = 1 (cdag_dn: Action=1, Spin=0). state 3 (|up dn>) + op 1 = index 13
  EXPECT_EQ(HubbardAtom<double>::lookup_table[(3 << 2) | 1].connected_state, 1);
  EXPECT_DOUBLE_EQ(HubbardAtom<double>::lookup_table[(3 << 2) | 1].matrix_element, -1.0);

  // op = 2 (c_up: Action=0, Spin=1). state 1 (|dn>) + op 2 = index 6
  EXPECT_DOUBLE_EQ(HubbardAtom<double>::lookup_table[(1 << 2) | 2].matrix_element, 0.0);
}

TEST_F(HubbardAtomTest, GreenFunctionG0_FiniteU_Order1) {
  double tau = 0.5;

  std::vector<double> taus = {0.0, 0.5};
  std::vector<int> ops     = {3, 2}; // cdag_up, c_up
  double g0                = atom->G0(taus, ops);

  double Z_exact  = 1 + 2 * std::exp(beta * mu) + std::exp(beta * (2 * mu - U));
  double G0_exact = -1.0 / Z_exact * (std::exp(tau * mu) + std::exp(beta * mu) * std::exp(-tau * (U - mu)));
  EXPECT_NEAR(g0, G0_exact, 1e-10);
}

TEST_F(HubbardAtomTest, GreenFunctionG0_InfiniteU_Order1) {

  std::vector<double> taus = {0.0, 0.5};
  std::vector<int> ops     = {3, 2}; // cdag_up, c_up
  double g0_inf            = atom->G0_infinite_U(taus, ops); //<T_tau c^dag_up(0.0) cup(0.5)> = <cup(0.5) c^dag_up(0.0)>

  double Z_inf_exact = 1 + 2 * std::exp(beta * mu);
  double expected    = -1.0 / Z_inf_exact;

  EXPECT_NEAR(g0_inf, expected, 1e-10);
}

TEST(HubbardAtomDualTest, G0IsParticleHoleSymmetric) {
  double U     = 8.0;
  double beta  = 1.0;
  double delta = 0.1;

  Dual mu1_val(U / 2.0 + delta, 1.0);
  Dual mu2_val(U / 2.0 - delta, 1.0);

  Parameters<Dual> params1{Dual(U, 0.0), Dual(beta, 0.0), mu1_val, Dual(0.0, 0.0), true};
  Parameters<Dual> params2{Dual(U, 0.0), Dual(beta, 0.0), mu2_val, Dual(0.0, 0.0), true};

  HubbardAtom<Dual> atom1(params1);
  HubbardAtom<Dual> atom2(params2);

  std::vector<double> taus = {0.8, 0.6, 0.4, 0.2};
  std::vector<int> ops     = {3, 2, 1, 0}; // cdag_up, c_up, cdag_dn, c_dn

  std::vector<double> taus_ph = {taus[1], taus[0], taus[3], taus[2]};
  std::vector<int> ops_ph     = {2, 3, 0, 1}; // Swap creation and annihilation

  Dual g1 = atom1.G0(taus, ops);
  Dual g2 = atom2.G0(taus_ph, ops_ph); // Pass the swapped times to mu2

  // Now, the physics is perfectly mirrored, and these will pass!
  EXPECT_NEAR(g1.value, g2.value, 1e-12);

  // The chain rule dictates that d/dmu G(mu) = d/dmu G(U - mu) * (-1)
  EXPECT_NEAR(g1.derivative, -g2.derivative, 1e-12);
}

TEST(HubbardDimerTest, LookupTable) {

  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;
  double t    = 1.0;

  Parameters<double> params{U, beta, mu, t, true};
  HubbardDimer<double> dimer(params);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
