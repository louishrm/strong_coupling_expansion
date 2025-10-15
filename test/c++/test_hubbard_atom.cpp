#include "../c++/sc_expansion/hubbard_atom.hpp"

void test_partition_function(double U, double mu, double beta, triqs::atom_diag::atom_diag<false> ad) {

  double Z_exact = 1 + 2 * exp(beta * mu) + exp(-beta * (U - 2 * mu));
  double Z       = hubbard_atom::partition_function(ad, beta);

  if (std::abs(Z - Z_exact) > 1e-10) {
    std::cout << "partition function test failed" << std::endl;
  } else {
    std::cout << "partition function test passed" << std::endl;
  }
}

void test_G01(double U, double mu, double beta, double tau, triqs::atom_diag::atom_diag<false> ad) {

  double Z_exact   = 1 + 2 * exp(beta * mu) + exp(-beta * (U - 2 * mu));
  double G01_exact = 1 / Z_exact * (exp(tau * mu) + exp(beta * mu) * exp(-tau * (U - mu)));

  hubbard_atom::cumul_args unprimed = {{tau, 0}};
  hubbard_atom::cumul_args primed   = {{0, 0}};
  double G01                        = hubbard_atom::G0(ad, beta, unprimed, primed);

  if (std::abs(G01 - G01_exact) > 1e-10) {
    std::cout << "G01 test failed" << std::endl;
  } else {
    std::cout << "G01 test passed" << std::endl;
  }
}

int main() {

  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;

  triqs::hilbert_space::fundamental_operator_set fops = hubbard_atom::make_fops();

  triqs::operators::many_body_operator_generic<double> H0 = hubbard_atom::make_H0(U, mu);

  triqs::atom_diag::atom_diag<false> ad(H0, fops, {}); // atom_diag object

  test_partition_function(U, mu, beta, ad);
  test_G01(U, mu, beta, 0.5, ad);

  return 0;
}