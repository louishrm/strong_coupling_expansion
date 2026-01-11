#include <gtest/gtest.h>
#include <cmath>
#include "../c++/sc_expansion/hubbard_atom.hpp"
#include <iostream>

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

TEST_F(HubbardAtomTest, PartitionFunctionMatchesExactResult) {
  double Z_exact = 1 + 2 * std::exp(beta * mu) + std::exp(-beta * (U - 2 * mu));
  double gs      = atom->ad.get_gs_energy();
  Z_exact *= std::exp(beta * gs); //H -> H-E0 shift for numerical stability
  double Z = triqs::atom_diag::partition_function(atom->ad, beta);

  // EXPECT_NEAR is better than 'abs > 1e-10' because it handles output formatting
  EXPECT_NEAR(Z, Z_exact, 1e-10);
}

TEST_F(HubbardAtomTest, AtomicDensityMatrixMatchesExactResult) {
  double Z_exact = 1 + 2 * std::exp(beta * mu) + std::exp(-beta * (U - 2 * mu));
  double gs      = atom->ad.get_gs_energy();
  Z_exact *= std::exp(beta * gs); //H -> H-E0 shift for numerical stability

  auto rho = triqs::atom_diag::atomic_density_matrix(atom->ad, beta)[0];

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
  double gs  = atom->ad.get_gs_energy();

  auto left = triqs::atom_diag::atomic_density_matrix(atom->ad, -tau)[0];
  left *= triqs::atom_diag::partition_function(atom->ad, -tau); //Multiply by Z0 to get unnormalized left operator

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

  HubbardAtom::cumul_args unprimed = {{tau, 0}};
  HubbardAtom::cumul_args primed   = {{0, 0}};

  double G01 = atom->G0(unprimed, primed);

  EXPECT_NEAR(G01, G01_exact, 1e-10);
}

TEST_F(HubbardAtomTest, GreenFunctionG02VanishConsecutiveTimes) {

  HubbardAtom::cumul_args unprimed = {{0.1, 0}, {0.2, 0}};
  HubbardAtom::cumul_args primed   = {{0.3, 0}, {0.4, 0}};

  double G02 = atom->G0(unprimed, primed);

  EXPECT_NEAR(G02, 0.0, 1e-10);
}

// Main is usually provided by gtest_main, but if you need a custom one:
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}