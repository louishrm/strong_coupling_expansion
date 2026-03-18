#include <gtest/gtest.h>
#include <cmath>
#include "../c++/sc_expansion/hubbard_solver.hpp"
#include "../c++/sc_expansion/dual.hpp"
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace sc_expansion;

// --- Args Tests ---

TEST(ArgsTest, ConstructorAndSorting) {

  std::vector<double> taus = {0.2, 0.8};
  std::vector<int> spins   = {0, 1};
  Args args(taus, spins);

  EXPECT_EQ(args.order, 2);

  EXPECT_DOUBLE_EQ(args.taus[0], 0.8);
  EXPECT_DOUBLE_EQ(args.taus[1], 0.2);

  EXPECT_EQ(args.ops[0], 1);
  EXPECT_EQ(args.ops[1], 2);

  EXPECT_DOUBLE_EQ(args.permutation_sign, -1.0);
}

TEST(ArgsTest, PermutationSignComplex) {
  // 4 operators, multiple swaps
  std::vector<double> taus = {0.1, 0.4, 0.2, 0.3}; // argsort: {1, 3, 2, 0}
  std::vector<int> spins   = {0, 0, 1, 1};
  Args args(taus, spins);

  EXPECT_DOUBLE_EQ(args.permutation_sign, 1.0);
}

TEST(ArgsTest, VerifyConsecutiveTermsInfiniteU) {

  {
    Args args({0.8, 0.2}, {0, 0});
    EXPECT_TRUE(args.verify_consecutive_terms_infinite_U());
  }

  {
    Args args({0.8, 0.2}, {0, 1});
    EXPECT_FALSE(args.verify_consecutive_terms_infinite_U());
  }

  {
    Args args({0.8, 0.2}, {1, 0});
    EXPECT_FALSE(args.verify_consecutive_terms_infinite_U());
  }

  {
    Args args({0.7, 0.9, 0.3, 0.5}, {0, 0, 1, 1});
    EXPECT_TRUE(args.verify_consecutive_terms_infinite_U());
  }
}

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

  EXPECT_EQ(HubbardAtom<double>::lookup_table[(0 << 2) | 2].connected_state, 1);
  EXPECT_DOUBLE_EQ(HubbardAtom<double>::lookup_table[(0 << 2) | 2].matrix_element, 1.0);

  EXPECT_EQ(HubbardAtom<double>::lookup_table[(1 << 2) | 3].connected_state, 3);
  EXPECT_DOUBLE_EQ(HubbardAtom<double>::lookup_table[(1 << 2) | 3].matrix_element, -1.0);

  EXPECT_EQ(HubbardAtom<double>::lookup_table[(1 << 2) | 0].connected_state, 0);
  EXPECT_DOUBLE_EQ(HubbardAtom<double>::lookup_table[(1 << 2) | 0].matrix_element, 1.0);

  EXPECT_EQ(HubbardAtom<double>::lookup_table[(3 << 2) | 1].connected_state, 1);
  EXPECT_DOUBLE_EQ(HubbardAtom<double>::lookup_table[(3 << 2) | 1].matrix_element, -1.0);

  EXPECT_DOUBLE_EQ(HubbardAtom<double>::lookup_table[(1 << 2) | 2].matrix_element, 0.0);
}

TEST_F(HubbardAtomTest, GreenFunctionG0_FiniteU_Order1) {
  double tau = 0.5;

  std::vector<double> taus = {0.0, 0.5};
  std::vector<int> spins   = {0, 0};
  double g0                = atom->G0(taus, spins);

  double Z_exact  = 1 + 2 * std::exp(beta * mu) + std::exp(beta * (2 * mu - U));
  double G0_exact = -1.0 / Z_exact * (std::exp(tau * mu) + std::exp(beta * mu) * std::exp(-tau * (U - mu)));
  EXPECT_NEAR(g0, G0_exact, 1e-10);
}

TEST_F(HubbardAtomTest, GreenFunctionG0_InfiniteU_Order1) {

  std::vector<double> taus = {0.0, 0.5};
  std::vector<int> spins   = {0, 0};
  double g0_inf            = atom->G0_infinite_U(taus, spins); //<T_tau c^dag_up(0.0) cup(0.5)> = <cup(0.5) c^dag_up(0.0)>

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
  std::vector<int> spins   = {0, 0, 1, 1}; // up, up, down, down

  std::vector<double> taus_ph = {taus[1], taus[0], taus[3], taus[2]};

  Dual g1 = atom1.G0(taus, spins);
  Dual g2 = atom2.G0(taus_ph, spins); // Pass the swapped times to mu2

  // Now, the physics is perfectly mirrored, and these will pass!
  EXPECT_NEAR(g1.value, g2.value, 1e-12);

  // The chain rule dictates that d/dmu G(mu) = d/dmu G(U - mu) * (-1)
  EXPECT_NEAR(g1.derivative, -g2.derivative, 1e-12);
}

// --- HubbardDimer Tests ---

class HubbardDimerTest : public ::testing::Test {
  protected:
  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;
  double t    = 1.0;
  Parameters<double> params{U, beta, mu, t, true};

  std::unique_ptr<HubbardDimer<double>> dimer;

  void SetUp() override { dimer = std::make_unique<HubbardDimer<double>>(params, t); }
};

TEST_F(HubbardDimerTest, TransitionTableVacuumToN1) {
  // op_index: 0: 1down destroy, 1: 2down destroy, 2: 1up destroy, 3: 2up destroy
  //           4: 1down create, 5: 2down create, 6: 1up create, 7: 2up create

  // Eigenstate 0 is vacuum |0>
  // cdag_1down (op 4) acting on |0>:
  // |0> is basis state 0.
  // cdag_1down (op 4) |0> -> basis state 1 (|down, 0>) with mel 1.0.
  // Eigenstates with basis state 1:
  // state 1: 1/sqrt(2) (|1> + |2>) -> mel 1/sqrt(2)
  // state 2: 1/sqrt(2) (|1> - |2>) -> mel 1/sqrt(2)

  const auto &transitions = dimer->transition_table[4][0].transitions;
  EXPECT_EQ(transitions.size(), 2);

  // Check state 1
  auto it1 = std::find_if(transitions.begin(), transitions.end(), [](const auto &t) { return t.connected_state == 1; });
  ASSERT_NE(it1, transitions.end());
  EXPECT_NEAR(it1->matrix_element, SQRT2_INV, 1e-12);

  // Check state 2
  auto it2 = std::find_if(transitions.begin(), transitions.end(), [](const auto &t) { return t.connected_state == 2; });
  ASSERT_NE(it2, transitions.end());
  EXPECT_NEAR(it2->matrix_element, SQRT2_INV, 1e-12);
}

TEST_F(HubbardDimerTest, TransitionTableSpinConservation) {
  // c_1up (op 2) acting on |0> should have no transitions
  const auto &transitions = dimer->transition_table[2][0].transitions;
  EXPECT_EQ(transitions.size(), 0);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
