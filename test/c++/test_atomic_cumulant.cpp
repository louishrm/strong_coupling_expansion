#include <gtest/gtest.h>
#include <cmath>
#include "../c++/sc_expansion/hubbard_atom.hpp"
#include "../c++/sc_expansion/cumulant.hpp"

using namespace sc_expansion;

// The Test Fixture: Sets up the data common to all tests
class HubbardAtomTest : public ::testing::Test {
  protected:
  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;

  // We declare the atom_diag object here, but initialize it in SetUp or constructor
  // Depending on TRIQS version, smart pointers (std::shared_ptr) might be safer/easier here
  // to avoid default constructor issues, but assuming copy/move works:
  std::unique_ptr<triqs::atom_diag::atom_diag<false>> ad;

  void SetUp() override {
    auto fops = hubbard_atom::make_fops();
    auto H0   = hubbard_atom::make_H0(U, mu);

    // This allows you to keep using {} because it bypasses template deduction
    ad.reset(new triqs::atom_diag::atom_diag<false>(H0, fops, {}));
  }
};

TEST_F(HubbardAtomTest, CumulantOrderOneMatchesExactResult) {

  double tau                        = 0.5;
  hubbard_atom::cumul_args unprimed = {{tau, 0}};
  hubbard_atom::cumul_args primed   = {{0, 0}};

  double G01 = hubbard_atom::G0(*ad, beta, unprimed, primed);
  double C01 = compute_cumulant_decomposition(unprimed, primed, *ad, beta);

  // EXPECT_NEAR is better than 'abs > 1e-10' because it handles output formatting
  EXPECT_NEAR(G01, C01, 1e-12);
}

TEST_F(HubbardAtomTest, CumulantOrderTwoMatchesExactResult) {

  hubbard_atom::cumul_args unprimed_args = {{0.5, 1}, {0.8, 0}};
  hubbard_atom::cumul_args primed_args   = {{0.0, 0}, {0.3, 1}};
  double G02                             = hubbard_atom::G0(*ad, beta, unprimed_args, primed_args); //G0(1,2|1',2')

  double G0_11 = hubbard_atom::G0(*ad, beta, {unprimed_args[0]}, {primed_args[0]}); //G(1|3)
  double G0_22 = hubbard_atom::G0(*ad, beta, {unprimed_args[1]}, {primed_args[1]}); //G(2|4)

  double G0_12 = hubbard_atom::G0(*ad, beta, {unprimed_args[0]}, {primed_args[1]}); //G(1|4)
  double G0_21 = hubbard_atom::G0(*ad, beta, {unprimed_args[1]}, {primed_args[0]}); //G(2|3)

  double C02_exact = G02 - G0_11 * G0_22 + G0_12 * G0_21;

  double C02 = compute_cumulant_decomposition(unprimed_args, primed_args, *ad, beta);

  EXPECT_NEAR(C02, C02_exact, 1e-12);
}

TEST_F(HubbardAtomTest, CumulantOrderThreeMatchesExactResult) {

  hubbard_atom::cumul_args unprimed_args = {{0.5, 1}, {0.8, 0}, {0.11, 0}};
  hubbard_atom::cumul_args primed_args   = {{0.0, 0}, {0.3, 1}, {0.65, 1}};
  double G03                             = hubbard_atom::G0(*ad, beta, unprimed_args, primed_args); //G03

  double C12_12_C33 = -compute_cumulant_decomposition({unprimed_args[0], unprimed_args[1]}, {primed_args[0], primed_args[1]}, *ad, beta)
     * hubbard_atom::G0(*ad, beta, {unprimed_args[2]}, {primed_args[2]}); //-C(1,2|1',2')*G(3|3')

  double C12_13_C32 = compute_cumulant_decomposition({unprimed_args[0], unprimed_args[1]}, {primed_args[0], primed_args[2]}, *ad, beta)
     * hubbard_atom::G0(*ad, beta, {unprimed_args[2]}, {primed_args[1]}); //C(1,2|1',3')*G(3|2')

  double C12_23_C31 = -compute_cumulant_decomposition({unprimed_args[0], unprimed_args[1]}, {primed_args[1], primed_args[2]}, *ad, beta)
     * hubbard_atom::G0(*ad, beta, {unprimed_args[2]}, {primed_args[0]}); //-C(1,2|2',3')*G(3|1')

  double C13_12_C32 = compute_cumulant_decomposition({unprimed_args[0], unprimed_args[2]}, {primed_args[0], primed_args[1]}, *ad, beta)
     * hubbard_atom::G0(*ad, beta, {unprimed_args[1]}, {primed_args[2]}); //C(1,3|1',2')*G(2|3')

  double C13_23_C21 = compute_cumulant_decomposition({unprimed_args[0], unprimed_args[2]}, {primed_args[1], primed_args[2]}, *ad, beta)
     * hubbard_atom::G0(*ad, beta, {unprimed_args[1]}, {primed_args[0]}); //C(1,3|2',3')*G(2|1')

  double C13_13_C22 = -compute_cumulant_decomposition({unprimed_args[0], unprimed_args[2]}, {primed_args[0], primed_args[2]}, *ad, beta)
     * hubbard_atom::G0(*ad, beta, {unprimed_args[1]}, {primed_args[1]}); //-C(1,3|1',3')*G(2|2')

  double C23_12_C13 = -compute_cumulant_decomposition({unprimed_args[1], unprimed_args[2]}, {primed_args[0], primed_args[1]}, *ad, beta)
     * hubbard_atom::G0(*ad, beta, {unprimed_args[0]}, {primed_args[2]}); //-C(2,3|1',2')*G(1|3')

  double C23_23_C11 = -compute_cumulant_decomposition({unprimed_args[1], unprimed_args[2]}, {primed_args[1], primed_args[2]}, *ad, beta)
     * hubbard_atom::G0(*ad, beta, {unprimed_args[0]}, {primed_args[0]}); //-C(2,3|2',3')*G(1|1')

  double C23_13_C12 = compute_cumulant_decomposition({unprimed_args[1], unprimed_args[2]}, {primed_args[0], primed_args[2]}, *ad, beta)
     * hubbard_atom::G0(*ad, beta, {unprimed_args[0]}, {primed_args[1]}); //C(2,3|1',3')*G(1|2')

  double C02C01_terms = C12_12_C33 + C12_13_C32 + C12_23_C31 + C13_12_C32 + C13_23_C21 + C13_13_C22 + C23_12_C13 + C23_23_C11 + C23_13_C12;

  double C11_C22_C33 = -hubbard_atom::G0(*ad, beta, {unprimed_args[0]}, {primed_args[0]})
     * hubbard_atom::G0(*ad, beta, {unprimed_args[1]}, {primed_args[1]})
     * hubbard_atom::G0(*ad, beta, {unprimed_args[2]}, {primed_args[2]}); //-G(1|1')*G(2|2')*G(3|3')

  double C11_C23_C32 = hubbard_atom::G0(*ad, beta, {unprimed_args[0]}, {primed_args[0]})
     * hubbard_atom::G0(*ad, beta, {unprimed_args[1]}, {primed_args[2]})
     * hubbard_atom::G0(*ad, beta, {unprimed_args[2]}, {primed_args[1]}); //G(1|1')*G(2|3')*G(3|2')

  double C12_C21_C33 = hubbard_atom::G0(*ad, beta, {unprimed_args[0]}, {primed_args[1]})
     * hubbard_atom::G0(*ad, beta, {unprimed_args[1]}, {primed_args[0]})
     * hubbard_atom::G0(*ad, beta, {unprimed_args[2]}, {primed_args[2]}); //G(1|2')*G(2|1')*G(3|3')

  double C12_C23_C31 = -hubbard_atom::G0(*ad, beta, {unprimed_args[0]}, {primed_args[1]})
     * hubbard_atom::G0(*ad, beta, {unprimed_args[1]}, {primed_args[2]})
     * hubbard_atom::G0(*ad, beta, {unprimed_args[2]}, {primed_args[0]}); //  -G(1|2')*G(2|3')*G(3|1')

  double C13_C22_C31 = hubbard_atom::G0(*ad, beta, {unprimed_args[0]}, {primed_args[2]})
     * hubbard_atom::G0(*ad, beta, {unprimed_args[1]}, {primed_args[1]})
     * hubbard_atom::G0(*ad, beta, {unprimed_args[2]}, {primed_args[0]}); //G(1|3')*G(2|2')*G(3|1')

  double C13_C21_C32 = -hubbard_atom::G0(*ad, beta, {unprimed_args[0]}, {primed_args[2]})
     * hubbard_atom::G0(*ad, beta, {unprimed_args[1]}, {primed_args[0]})
     * hubbard_atom::G0(*ad, beta, {unprimed_args[2]}, {primed_args[1]}); //-G(1|3')*G(2|1')*G(3|2')

  double G1G1G1_terms = C11_C22_C33 + C11_C23_C32 + C12_C21_C33 + C12_C23_C31 + C13_C22_C31 + C13_C21_C32;

  double C03_exact = G03 + C02C01_terms + G1G1G1_terms;

  double C03 = compute_cumulant_decomposition(unprimed_args, primed_args, *ad, beta);

  EXPECT_NEAR(C03, C03_exact, 1e-12);
}

TEST_F(HubbardAtomTest, SpinConservationOfCumulant) {

  hubbard_atom::cumul_args unprimed_args1 = {{0.5, 1}, {0.8, 1}};
  hubbard_atom::cumul_args primed_args1   = {{0.0, 0}, {0.3, 0}};

  double C02_1 = compute_cumulant_decomposition(unprimed_args1, primed_args1, *ad, beta);

  EXPECT_NEAR(C02_1, 0.0, 1e-12);

  hubbard_atom::cumul_args unprimed_args2 = {{0.5, 1}, {0.8, 0}};
  hubbard_atom::cumul_args primed_args2   = {{0.0, 0}, {0.3, 0}};

  double C02_2 = compute_cumulant_decomposition(unprimed_args2, primed_args2, *ad, beta);
  EXPECT_NEAR(C02_2, 0.0, 1e-12);
}
