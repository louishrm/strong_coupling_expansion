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
  // Scenario: t1 < t2, but input as {t1, t2}. Should be sorted to {t2, t1}.
  // Index 0: spin 0 (up) -> creation (op=2)
  // Index 1: spin 1 (down) -> annihilation (op=1)
  std::vector<double> taus = {0.2, 0.8};
  std::vector<int> spins   = {0, 1};
  Args args(taus, spins);

  EXPECT_EQ(args.order, 2);
  // Sorted descending:
  EXPECT_DOUBLE_EQ(args.taus[0], 0.8);
  EXPECT_DOUBLE_EQ(args.taus[1], 0.2);

  // Initial ops were: ops[0]=2 (cdag_up), ops[1]=1 (cdn)
  // After sorting (argsort={1, 0}):
  EXPECT_EQ(args.ops[0], 1); // cdn
  EXPECT_EQ(args.ops[1], 2); // cdag_up

  // One swap: permutation sign should be -1.0
  EXPECT_DOUBLE_EQ(args.permutation_sign, -1.0);
}

TEST(ArgsTest, PermutationSignComplex) {
  // 4 operators, multiple swaps
  std::vector<double> taus = {0.1, 0.4, 0.2, 0.3}; // argsort: {1, 3, 2, 0}
  std::vector<int> spins   = {0, 0, 1, 1};
  Args args(taus, spins);

  // argsort = {1, 3, 2, 0}
  // 1 3 2 0
  // Initial: 0 1 2 3
  // Swap(0, 1) -> 1 0 2 3 (sign -1)
  // Swap(1, 3) -> 1 3 2 0 (sign +1) -> Wait, 1 3 2 0 is not 1 swap from 1 0 2 3.
  // 0 1 2 3 -> 1 0 2 3 (1) -> 1 2 0 3 (2) -> 1 3 0 2 (3) -> 1 3 2 0 (4)
  // Or: inversions in {1, 3, 2, 0}:
  // (1,0), (3,2), (3,0), (2,0) -> 4 inversions. Sign = +1.0.
  EXPECT_DOUBLE_EQ(args.permutation_sign, 1.0);
}

TEST(ArgsTest, VerifyConsecutiveTermsInfiniteU) {
  // Test valid sequence: cdag_up(0.8), cup(0.2)
  // taus = {0.8, 0.2}, spins = {0, 0} -> ops = {2, 0}.
  // i=0: ops[0]=2 (late), i+1=1: ops[1]=0 (early)
  // early_is_create=false, late_is_create=true. early_is_create != late_is_create (True).
  // No Rule B check because early is not create.
  // Beta check: op_beta = 2, op_zero = 0. beta_is_create = true. op_zero == 2-2=0. (True).
  {
    Args args({0.8, 0.2}, {0, 0});
    EXPECT_TRUE(args.verify_consecutive_terms_infinite_U());
  }

  // Invalid: Consecutive creations
  // taus = {0.8, 0.2}, spins = {0, 1} -> ops = {2, 3}.
  // early_is_create=true, late_is_create=true -> False.
  {
    Args args({0.8, 0.2}, {0, 1});
    EXPECT_FALSE(args.verify_consecutive_terms_infinite_U());
  }

  // Invalid: Spin mismatch (creation of up, destruction of down)
  // taus = {0.8, 0.2}, spins = {1, 0} -> ops = {3, 0}.
  // i=0: op_late=3 (cdag_dn), i+1=1: op_early=0 (cup).
  // early_is_create=false, late_is_create=true.
  // Rule B check: None.
  // Beta check: op_beta=3 (cdag_dn), op_zero=0 (cup). beta_is_create=true.
  // op_zero != 3-2 (0 != 1) -> False.
  {
    Args args({0.8, 0.2}, {1, 0});
    EXPECT_FALSE(args.verify_consecutive_terms_infinite_U());
  }

  // Valid: D(0.9), C(0.7), D(0.5), C(0.3)
  // taus={0.9, 0.7, 0.5, 0.3}, spins={0, 0, 1, 1}
  // ops={0, 2, 1, 3} (Wait, i=0: cdag, i=1: c ... so this is D, C, D, C?)
  // No, i=0: cdag, i=1: c.
  // input: taus={0.7, 0.9, 0.3, 0.5}, spins={0, 0, 1, 1}
  // indices: 0: cdag_up, 1: cup, 2: cdag_dn, 3: cdn
  // sorted: 1(0.9), 0(0.7), 3(0.5), 2(0.3)
  // sorted_ops: ops[1]=0, ops[0]=2, ops[3]=1, ops[2]=3
  // D(0.9), C(0.7), D(0.5), C(0.3)
  // Pairs: (0, 2) -> early=2, late=0. early_is_create=true. op_late = 0 == 2-2. OK.
  // (2, 1) -> early=1, late=2. early_is_create=false. OK.
  // (1, 3) -> early=3, late=1. early_is_create=true. op_late = 1 == 3-2. OK.
  // Boundary: op_beta=0, op_zero=3. beta_is_create=false. OK.
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
  // index = (state << 2) | op
  // op 0: cup, 1: cdn, 2: cdag_up, 3: cdag_dn
  // States: 0:|0>, 1:|up>, 2:|dn>, 3:|up dn>

  // cdag_up on |0> (state 0, op 2) -> |up> (state 1), mel 1.0
  EXPECT_EQ(HubbardAtom<double>::lookup_table[(0 << 2) | 2].connected_state, 1);
  EXPECT_DOUBLE_EQ(HubbardAtom<double>::lookup_table[(0 << 2) | 2].matrix_element, 1.0);

  // cdag_dn on |up> (state 1, op 3) -> |up dn> (state 3), mel -1.0 (jump over up)
  EXPECT_EQ(HubbardAtom<double>::lookup_table[(1 << 2) | 3].connected_state, 3);
  EXPECT_DOUBLE_EQ(HubbardAtom<double>::lookup_table[(1 << 2) | 3].matrix_element, -1.0);

  // cup on |up> (state 1, op 0) -> |0> (state 0), mel 1.0
  EXPECT_EQ(HubbardAtom<double>::lookup_table[(1 << 2) | 0].connected_state, 0);
  EXPECT_DOUBLE_EQ(HubbardAtom<double>::lookup_table[(1 << 2) | 0].matrix_element, 1.0);

  // cdn on |up dn> (state 3, op 1) -> |up> (state 1), mel -1.0 (jump over up)
  EXPECT_EQ(HubbardAtom<double>::lookup_table[(3 << 2) | 1].connected_state, 1);
  EXPECT_DOUBLE_EQ(HubbardAtom<double>::lookup_table[(3 << 2) | 1].matrix_element, -1.0);

  // Pauli exclusion: cdag_up on |up> (state 1, op 2) -> mel 0.0
  EXPECT_DOUBLE_EQ(HubbardAtom<double>::lookup_table[(1 << 2) | 2].matrix_element, 0.0);
}

TEST_F(HubbardAtomTest, GreenFunctionG0_FiniteU_Order1) {
  double tau = 0.5;
  // Compute G0(tau, 0) = - < T c_up(tau) cdag_up(0) >
  // For tau > 0, G0 = - < c_up(tau) cdag_up(0) >
  // Input to G0 method: {tau, 0} for c, cdag.
  // Wait, the code expects {tau_cdag, tau_c}.
  // Let's check G0 implementation again.
  // It computes Tr(rho O_n(tau_n) ... O_1(tau_1)) where tau_n > ... > tau_1.
  // To get < c_up(tau) cdag_up(0) >, we need tau_c = tau, tau_cdag = 0.
  // Since tau > 0, tau_c is the LATEST operator.
  // Args({tau_c, tau_cdag}, {0, 0}) -> taus={0.5, 0.0}, spins={0, 0}.
  // ops[0]=2 (cdag_up), ops[1]=0 (cup).
  // Sorted: 0.5, 0.0. sorted_ops: ops[0], ops[1] -> 2, 0.
  // Wait, if taus={0.5, 0.0}, argsort={0, 1}. sorted_ops = {ops[0], ops[1]} = {2, 0}.
  // Loop i=1 (tau=0.0): op=0 (cup). Loop i=0 (tau=0.5): op=2 (cdag_up).
  // This computes Tr(rho cdag_up(0.5) cup(0)).

  // To get Tr(rho cup(0.5) cdag_up(0)):
  // We need op_late = cup, op_early = cdag_up.
  // Since tau_late > tau_early, we need tau_cup=0.5, tau_cdag=0.0.
  // In Args(taus, spins):
  // ops[i] = 2+spin if i even, spin if i odd.
  // So we need i=1 for cup, i=0 for cdag_up.
  // taus = {0.0, 0.5}, spins = {0, 0}.
  // ops[0] = 2 (cdag_up), ops[1] = 0 (cup).
  // Sorted: 0.5(i=1), 0.0(i=0). sorted_ops: ops[1], ops[0] -> 0, 2.
  // Permutation sign: swap (0, 1) -> -1.0.
  // Result = -1.0 * Tr(rho cup(0.5) cdag_up(0)).
  // This is exactly G0(0.5, 0)!

  std::vector<double> taus = {0.0, 0.5};
  std::vector<int> spins   = {0, 0};
  double g0                = atom->G0(taus, spins);

  // Exact G0(tau) = -1/Z * (exp(-beta*E0) * exp(tau*(E0-E1)) * 1 * 1 + exp(-beta*E2) * exp(tau*(E2-E3)) * 1 * 1)
  // G0(tau) = -1/Z * (exp(tau*mu) + exp(-beta*(-mu)) * exp(tau*(-mu - (2*mu-U))))
  // G0(tau) = -1/Z * (exp(tau*mu) + exp(beta*mu) * exp(tau*(U-mu)))
  // Wait, my signs...
  // Tr(rho cup(tau) cdag_up(0)) = 1/Z * Sum_a <a| e^{-beta H} e^{tau H} cup e^{-tau H} cdag |a>
  // a=0: <0| e^{-beta E0} e^{tau E0} cup e^{-tau E1} cdag |0> = exp(-beta*0) * exp(tau*0) * 1 * exp(-tau*(-mu)) * 1 = exp(tau*mu)
  // a=2: <2| e^{-beta E2} e^{tau E2} cup e^{-tau E3} cdag |2> = exp(-beta*-mu) * exp(tau*-mu) * (-1) * exp(-tau*(2*mu-U)) * (-1) = exp(beta*mu) * exp(tau*(U-mu))
  // G0 = -1/Z * (exp(tau*mu) + exp(beta*mu) * exp(tau*(U-mu)))

  double Z_exact  = 1 + 2 * std::exp(beta * mu) + std::exp(beta * (2 * mu - U));
  double G0_exact = -1.0 / Z_exact * (std::exp(tau * mu) + std::exp(beta * mu) * std::exp(-tau * (U - mu)));
  EXPECT_NEAR(g0, G0_exact, 1e-10);
}

TEST_F(HubbardAtomTest, GreenFunctionG0_InfiniteU_Order1) {
  // taus = {0.0, 0.5}, spins = {0, 0}
  // G0_infinite_U implementation:
  // args.ops = {0, 2} sorted -> {0, 2} (wait, 0.5 is at index 1).
  // sorted_ops = {0, 2}. first_op = 0.
  // Result = args.permutation_sign * (-exp(beta*mu) / Z_infinite_U)
  // permutation_sign = -1.0. Result = exp(beta*mu) / Z_infinite_U.

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

  // FIX: Swap the creation and annihilation times to test the "hole" counterpart!
  // Old creation (0, 2) become new annihilation. Old annihilation (1, 3) become new creation.
  std::vector<double> taus_ph = {taus[1], taus[0], taus[3], taus[2]};

  Dual g1 = atom1.G0(taus, spins);
  Dual g2 = atom2.G0(taus_ph, spins); // Pass the swapped times to mu2

  // Now, the physics is perfectly mirrored, and these will pass!
  EXPECT_NEAR(g1.value, g2.value, 1e-12);

  // The chain rule dictates that d/dmu G(mu) = d/dmu G(U - mu) * (-1)
  EXPECT_NEAR(g1.derivative, -g2.derivative, 1e-12);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
