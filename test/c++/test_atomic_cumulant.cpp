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

  std::unique_ptr<HubbardAtom> atom;

  void SetUp() override { atom = std::make_unique<HubbardAtom>(U, beta, mu); }
};

TEST_F(HubbardAtomTest, CumulantOrderOneMatchesExactResult) {

  double tau                       = 0.5;
  HubbardAtom::cumul_args unprimed = {{tau, 0}};
  HubbardAtom::cumul_args primed   = {{0, 0}};

  double G01 = atom->G0(unprimed, primed);
  double C01 = compute_cumulant_decomposition(unprimed, primed, *atom, false);

  // EXPECT_NEAR is better than 'abs > 1e-10' because it handles output formatting
  EXPECT_NEAR(G01, C01, 1e-12);
}

TEST_F(HubbardAtomTest, CumulantOrderTwoMatchesExactResult) {

  HubbardAtom::cumul_args unprimed_args = {{0.5, 1}, {0.8, 0}};
  HubbardAtom::cumul_args primed_args   = {{0.0, 0}, {0.3, 1}};
  double G02                            = atom->G0(unprimed_args, primed_args); //G0(1,2|1',2')

  double G0_11 = atom->G0({unprimed_args[0]}, {primed_args[0]}); //G(1|3)
  double G0_22 = atom->G0({unprimed_args[1]}, {primed_args[1]}); //G(2|4)

  double G0_12 = atom->G0({unprimed_args[0]}, {primed_args[1]}); //G(1|4)
  double G0_21 = atom->G0({unprimed_args[1]}, {primed_args[0]}); //G(2|3)

  double C02_exact = G02 - G0_11 * G0_22 + G0_12 * G0_21;

  double C02 = compute_cumulant_decomposition(unprimed_args, primed_args, *atom, false);

  EXPECT_NEAR(C02, C02_exact, 1e-12);
}

TEST_F(HubbardAtomTest, CumulantOrderThreeMatchesExactResult) {

  HubbardAtom::cumul_args unprimed_args = {{0.5, 1}, {0.8, 0}, {0.11, 0}};
  HubbardAtom::cumul_args primed_args   = {{0.0, 0}, {0.3, 1}, {0.65, 0}};
  double G03                            = atom->G0(unprimed_args, primed_args); //G03

  double C12_12_C33 = -compute_cumulant_decomposition({unprimed_args[0], unprimed_args[1]}, {primed_args[0], primed_args[1]}, *atom, false)
     * atom->G0({unprimed_args[2]}, {primed_args[2]}); //-C(1,2|1',2')*G(3|3')

  double C12_13_C32 = compute_cumulant_decomposition({unprimed_args[0], unprimed_args[1]}, {primed_args[0], primed_args[2]}, *atom, false)
     * atom->G0({unprimed_args[2]}, {primed_args[1]}); //C(1,2|1',3')*G(3|2')

  double C12_23_C31 = -compute_cumulant_decomposition({unprimed_args[0], unprimed_args[1]}, {primed_args[1], primed_args[2]}, *atom, false)
     * atom->G0({unprimed_args[2]}, {primed_args[0]}); //-C(1,2|2',3')*G(3|1')

  double C13_12_C32 = compute_cumulant_decomposition({unprimed_args[0], unprimed_args[2]}, {primed_args[0], primed_args[1]}, *atom, false)
     * atom->G0({unprimed_args[1]}, {primed_args[2]}); //C(1,3|1',2')*G(2|3')

  double C13_23_C21 = compute_cumulant_decomposition({unprimed_args[0], unprimed_args[2]}, {primed_args[1], primed_args[2]}, *atom, false)
     * atom->G0({unprimed_args[1]}, {primed_args[0]}); //C(1,3|2',3')*G(2|1')

  double C13_13_C22 = -compute_cumulant_decomposition({unprimed_args[0], unprimed_args[2]}, {primed_args[0], primed_args[2]}, *atom, false)
     * atom->G0({unprimed_args[1]}, {primed_args[1]}); //-C(1,3|1',3')*G(2|2')

  double C23_12_C13 = -compute_cumulant_decomposition({unprimed_args[1], unprimed_args[2]}, {primed_args[0], primed_args[1]}, *atom, false)
     * atom->G0({unprimed_args[0]}, {primed_args[2]}); //-C(2,3|1',2')*G(1|3')

  double C23_23_C11 = -compute_cumulant_decomposition({unprimed_args[1], unprimed_args[2]}, {primed_args[1], primed_args[2]}, *atom, false)
     * atom->G0({unprimed_args[0]}, {primed_args[0]}); //-C(2,3|2',3')*G(1|1')

  double C23_13_C12 = compute_cumulant_decomposition({unprimed_args[1], unprimed_args[2]}, {primed_args[0], primed_args[2]}, *atom, false)
     * atom->G0({unprimed_args[0]}, {primed_args[1]}); //C(2,3|1',3')*G(1|2')

  double C02C01_terms = C12_12_C33 + C12_13_C32 + C12_23_C31 + C13_12_C32 + C13_23_C21 + C13_13_C22 + C23_12_C13 + C23_23_C11 + C23_13_C12;

  double C11_C22_C33 = -atom->G0({unprimed_args[0]}, {primed_args[0]}) * atom->G0({unprimed_args[1]}, {primed_args[1]})
     * atom->G0({unprimed_args[2]}, {primed_args[2]}); //-G(1|1')*G(2|2')*G(3|3')

  double C11_C23_C32 = atom->G0({unprimed_args[0]}, {primed_args[0]}) * atom->G0({unprimed_args[1]}, {primed_args[2]})
     * atom->G0({unprimed_args[2]}, {primed_args[1]}); //G(1|1')*G(2|3')*G(3|2')

  double C12_C21_C33 = atom->G0({unprimed_args[0]}, {primed_args[1]}) * atom->G0({unprimed_args[1]}, {primed_args[0]})
     * atom->G0({unprimed_args[2]}, {primed_args[2]}); //G(1|2')*G(2|1')*G(3|3')

  double C12_C23_C31 = -atom->G0({unprimed_args[0]}, {primed_args[1]}) * atom->G0({unprimed_args[1]}, {primed_args[2]})
     * atom->G0({unprimed_args[2]}, {primed_args[0]}); //  -G(1|2')*G(2|3')*G(3|1')

  double C13_C22_C31 = atom->G0({unprimed_args[0]}, {primed_args[2]}) * atom->G0({unprimed_args[1]}, {primed_args[1]})
     * atom->G0({unprimed_args[2]}, {primed_args[0]}); //G(1|3')*G(2|2')*G(3|1')

  double C13_C21_C32 = -atom->G0({unprimed_args[0]}, {primed_args[2]}) * atom->G0({unprimed_args[1]}, {primed_args[0]})
     * atom->G0({unprimed_args[2]}, {primed_args[1]}); //-G(1|3')*G(2|1')*G(3|2')

  double G1G1G1_terms = C11_C22_C33 + C11_C23_C32 + C12_C21_C33 + C12_C23_C31 + C13_C22_C31 + C13_C21_C32;

  double C03_exact = G03 + C02C01_terms + G1G1G1_terms;

  double C03 = compute_cumulant_decomposition(unprimed_args, primed_args, *atom, true);

  EXPECT_NEAR(C03, C03_exact, 1e-12);
}

TEST_F(HubbardAtomTest, SpinConservationOfCumulant) {

  HubbardAtom::cumul_args unprimed_args1 = {{0.5, 1}, {0.8, 1}};
  HubbardAtom::cumul_args primed_args1   = {{0.0, 0}, {0.3, 0}};

  double C02_1 = compute_cumulant_decomposition(unprimed_args1, primed_args1, *atom, false);

  EXPECT_NEAR(C02_1, 0.0, 1e-12);

  HubbardAtom::cumul_args unprimed_args2 = {{0.5, 1}, {0.8, 0}};
  HubbardAtom::cumul_args primed_args2   = {{0.0, 0}, {0.3, 0}};

  double C02_2 = compute_cumulant_decomposition(unprimed_args2, primed_args2, *atom, false);
  EXPECT_NEAR(C02_2, 0.0, 1e-12);
}