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
  double G01_exact = -1.0 / Z_exact * (std::exp(tau * mu) + std::exp(beta * mu) * std::exp(-tau * (U - mu)));

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

TEST_F(HubbardAtomTest, GreenFunctionG01InfiniteUMatchesExactResult) {
  double tau = 0.5;

  double Z_infinite_U = 1 + 2 * std::exp(beta * mu);

  //tau = destroy
  double G01_exact_1 = 1.0 / Z_infinite_U; //* (std::exp(tau * mu));

  double G01_exact_2 = -G01_exact_1 * std::exp(beta * mu);

  HubbardAtom::cumul_args unprimed = {{tau, 0}};
  HubbardAtom::cumul_args primed   = {{0, 0}};

  double G01_1 = atom->G0_infinite_U(unprimed, primed);
  double G01_2 = atom->G0_infinite_U(primed, unprimed);

  EXPECT_NEAR(G01_1, G01_exact_1, 1e-12);
  EXPECT_NEAR(G01_2, G01_exact_2, 1e-12);
}

TEST_F(HubbardAtomTest, GreenFunctionG02InfiniteUMatchesExactResult) {

  double Z_infinite_U = 1 + 2 * std::exp(beta * mu);

  //tau = destroy
  double G02_exact_1 = -1.0 / Z_infinite_U; //* (std::exp(tau * mu));

  double G02_exact_2 = G02_exact_1 * std::exp(beta * mu);

  HubbardAtom::cumul_args unprimed = {{0.7, 0}, {0.3, 0}};
  HubbardAtom::cumul_args primed   = {{0.5, 0}, {0.1, 0}};

  double G02_1 = atom->G0_infinite_U(unprimed, primed);
  double G02_2 = atom->G0_infinite_U(primed, unprimed);

  EXPECT_NEAR(G02_1, G02_exact_1, 1e-12);
  EXPECT_NEAR(G02_2, G02_exact_2, 1e-12);
}

TEST_F(HubbardAtomTest, VerifyConsecutiveTerms) {
  // 1. Three consecutive creates/destroys -> False
  EXPECT_FALSE(HubbardAtom::verify_consecutive_terms({0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 1, 1})); //3 consecutive creates
  EXPECT_FALSE(HubbardAtom::verify_consecutive_terms({0, 0, 0, 0, 0, 0}, {1, 1, 1, 0, 0, 0})); //3 consecutive destroys

  // 2. Two consecutive creates/destroys with same spin -> False
  EXPECT_FALSE(HubbardAtom::verify_consecutive_terms({0, 0, 1, 1}, {0, 0, 1, 1})); //creates with same spin
  EXPECT_FALSE(HubbardAtom::verify_consecutive_terms({1, 1, 0, 0}, {1, 1, 0, 0})); //destroys with same spin

  // Two consecutive create with a destroy of opposite spin in between -> false
  EXPECT_FALSE(HubbardAtom::verify_consecutive_terms({0, 1, 0, 1}, {0, 1, 0, 1})); //up+ dn - up+ dn -

  //Valid alternating sequence
  EXPECT_TRUE(HubbardAtom::verify_consecutive_terms({0, 1, 0, 1}, {1, 0, 0, 1}));

  //complex sequence
  std::vector<int> complex_spins = {0, 1, 0, 1, 0, 0};
  std::vector<int> complex_flags = {1, 0, 0, 1, 1, 0};
  EXPECT_TRUE(HubbardAtom::verify_consecutive_terms(complex_spins, complex_flags));
}
TEST_F(HubbardAtomTest, VerifyConsecutiveTermsInfiniteU) {
  // 1. Consecutive Identical Flags -> False
  EXPECT_FALSE(HubbardAtom::verify_consecutive_terms_infinite_U({0, 0}, {0, 0}, 1)); // 0, 0
  EXPECT_FALSE(HubbardAtom::verify_consecutive_terms_infinite_U({0, 0}, {1, 1}, 1)); // 1, 1

  // 2. Destroy(t0), Create(t1) [t0 > t1] -> Particle in (t1, t0). Spins MUST match.
  // flags = {0, 1}
  EXPECT_TRUE(HubbardAtom::verify_consecutive_terms_infinite_U({0, 0}, {0, 1}, 1));  // Match
  EXPECT_FALSE(HubbardAtom::verify_consecutive_terms_infinite_U({0, 1}, {0, 1}, 1)); // Mismatch

  // 3. Create(t0), Destroy(t1) [t0 > t1] -> Particle in (t0, beta) U (0, t1). Spins MUST match (Boundary).
  // flags = {1, 0}
  EXPECT_TRUE(HubbardAtom::verify_consecutive_terms_infinite_U({0, 0}, {1, 0}, 1));  // Match
  EXPECT_FALSE(HubbardAtom::verify_consecutive_terms_infinite_U({0, 1}, {1, 0}, 1)); // Mismatch

  // 4. Complex Chain: D(t0), C(t1), D(t2), C(t3).
  // (t1, t0) -> Particle. (t3, t2) -> Particle.
  // flags = {0, 1, 0, 1}
  // spins = {0, 0, 1, 1} -> OK.
  EXPECT_TRUE(HubbardAtom::verify_consecutive_terms_infinite_U({0, 0, 1, 1}, {0, 1, 0, 1}, 2));

  // spins = {0, 1, 1, 1} -> Fail at first pair (0 vs 1).
  EXPECT_FALSE(HubbardAtom::verify_consecutive_terms_infinite_U({0, 1, 1, 1}, {0, 1, 0, 1}, 2));
}

// Main is usually provided by gtest_main, but if you need a custom one:
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}