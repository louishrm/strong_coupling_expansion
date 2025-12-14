#include <gtest/gtest.h>
#include <cmath>
#include "../c++/sc_expansion/hubbard_atom.hpp"

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

TEST_F(HubbardAtomTest, PartitionFunctionMatchesExactResult) {
  double Z_exact = 1 + 2 * std::exp(beta * mu) + std::exp(-beta * (U - 2 * mu));
  double Z       = hubbard_atom::partition_function(*ad, beta);

  // EXPECT_NEAR is better than 'abs > 1e-10' because it handles output formatting
  EXPECT_NEAR(Z, Z_exact, 1e-10);
}

TEST_F(HubbardAtomTest, GreenFunctionG01MatchesExactResult) {
  double tau = 0.5;

  double Z_exact   = 1 + 2 * std::exp(beta * mu) + std::exp(-beta * (U - 2 * mu));
  double G01_exact = 1.0 / Z_exact * (std::exp(tau * mu) + std::exp(beta * mu) * std::exp(-tau * (U - mu)));

  hubbard_atom::cumul_args unprimed = {{tau, 0}};
  hubbard_atom::cumul_args primed   = {{0, 0}};

  double G01 = hubbard_atom::G0(*ad, beta, unprimed, primed);

  EXPECT_NEAR(G01, G01_exact, 1e-10);
}

// Main is usually provided by gtest_main, but if you need a custom one:
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}