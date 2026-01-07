#include <gtest/gtest.h>
#include <cmath>
#include "../c++/sc_expansion/hubbard_atom.hpp"
#include <iostream>

// The Test Fixture: Sets up the data common to all tests
class HubbardAtomTest : public ::testing::Test {
  protected:
  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;

  triqs::atom_diag::atom_diag<false> ad;

  void SetUp() override {
    auto fops = hubbard_atom::make_fops();
    auto H0   = hubbard_atom::make_H0(U, mu);

    // This allows you to keep using {} because it bypasses template deduction
    ad = triqs::atom_diag::atom_diag<false>(H0, fops, {});
  }
};

TEST_F(HubbardAtomTest, PartitionFunctionMatchesExactResult) {
  double Z_exact = 1 + 2 * std::exp(beta * mu) + std::exp(-beta * (U - 2 * mu));
  double gs      = ad.get_gs_energy();
  Z_exact *= std::exp(beta * gs); //H -> H-E0 shift for numerical stability
  double Z = hubbard_atom::_partition_function(ad, beta);

  // EXPECT_NEAR is better than 'abs > 1e-10' because it handles output formatting
  EXPECT_NEAR(Z, Z_exact, 1e-10);
}

TEST_F(HubbardAtomTest, AtomicDensityMatrixMatchesExactResult) {
  double Z_exact = 1 + 2 * std::exp(beta * mu) + std::exp(-beta * (U - 2 * mu));
  double gs      = ad.get_gs_energy();
  Z_exact *= std::exp(beta * gs); //H -> H-E0 shift for numerical stability

  auto rho = triqs::atom_diag::atomic_density_matrix(ad, beta)[0];

  std::vector<double> rho_entries_exact = {
     std::exp(-beta * (0 - gs)) / Z_exact,
     std::exp(-beta * (-mu - gs)) / Z_exact,
     std::exp(-beta * (-mu - gs)) / Z_exact,
     std::exp(-beta * (U - 2 * mu - gs)) / Z_exact,
  };

  std::sort(rho_entries_exact.begin(), rho_entries_exact.end(), std::greater<double>());

  EXPECT_NEAR(rho(0, 0), rho_entries_exact[0], 1e-10);
  EXPECT_NEAR(rho(1, 1), rho_entries_exact[1], 1e-10);
  EXPECT_NEAR(rho(2, 2), rho_entries_exact[2], 1e-10);
  EXPECT_NEAR(rho(3, 3), rho_entries_exact[3], 1e-10);
}

TEST_F(HubbardAtomTest, ImaginaryTimeEvolutionMatchesExactResult) {

  double tau = 0.5;
  double gs  = ad.get_gs_energy();

  auto left = triqs::atom_diag::atomic_density_matrix(ad, -tau)[0];
  left *= hubbard_atom::_partition_function(ad, -tau); //Multiply by Z0 to get unnormalized left operator

  std::vector<double> left_exact = {std::exp(tau * (0 - gs)), std::exp(tau * (-mu - gs)), std::exp(tau * (-mu - gs)),
                                    std::exp(tau * (U - 2 * mu - gs))};

  std::sort(left_exact.begin(), left_exact.end());
  EXPECT_NEAR(left(0, 0), left_exact[0], 1e-10);
  EXPECT_NEAR(left(1, 1), left_exact[1], 1e-10);
  EXPECT_NEAR(left(2, 2), left_exact[2], 1e-10);
  EXPECT_NEAR(left(3, 3), left_exact[3], 1e-10);
}

TEST_F(HubbardAtomTest, GreenFunctionG01MatchesExactResult) {
  //also checks shift invariance of G0 with respect to H -> H-E0
  double tau = 0.5;

  double Z_exact   = 1 + 2 * std::exp(beta * mu) + std::exp(-beta * (U - 2 * mu));
  double G01_exact = 1.0 / Z_exact * (std::exp(tau * mu) + std::exp(beta * mu) * std::exp(-tau * (U - mu)));

  hubbard_atom::cumul_args unprimed = {{tau, 0}};
  hubbard_atom::cumul_args primed   = {{0, 0}};

  double G01 = hubbard_atom::G0(ad, beta, unprimed, primed);

  EXPECT_NEAR(G01, G01_exact, 1e-10);
}

// Main is usually provided by gtest_main, but if you need a custom one:
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
