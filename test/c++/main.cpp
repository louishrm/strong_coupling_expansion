//#include "sc_expansion/sc_expansion.hpp"
#include "../c++/sc_expansion/hubbard_atom.hpp"

#include <cmath>

double o2_exact(double U, double mu, double beta) {
  double term1 = std::exp(beta * mu);
  double term2 = beta * (1 + std::exp(-beta * (U - 2 * mu)));
  double term3 = (4 / U) * std::exp(-beta * (U / 2 - mu)) * std::sinh(beta * U / 2);

  double expr = 2 * term1 * beta * (term2 + term3);
  return expr;
}

double exact_one_body_GF(double U, double mu, double beta, double Z0, double tau_1, double tau_2) {
  double expr = 1 / Z0 * std::exp(-mu * (tau_1 - tau_2)) * (std::exp(beta * mu) + std::exp(U * (tau_1 - tau_2)) * std::exp(-beta * (U - 2 * mu)));

  return expr;
}

double free_energy_o2_trimer_contrib(auto ad, double beta, double tau1, double tau2) {

  std::vector<int> spins   = {0, 0};
  std::vector<int> flags_1 = {0, 1};
  std::vector<int> flags_2 = {1, 0};

  double G01 = hubbard_atom::G0(ad, beta, {tau1, tau2}, spins, flags_1);
  double G02 = hubbard_atom::G0(ad, beta, {tau1, tau2}, spins, flags_2);

  double symmetry_factor = 2.0 * 6.0; // 0.5 from series, 2 for spin, 6 for permutations of sites

  double contrib = symmetry_factor * (G01 * G02);

  return contrib;
}

double free_energy_o2_trimer(auto ad, double beta, double delta) {

  double integral = 0.0;
  for (double tau_1 = 0.0; tau_1 < beta; tau_1 += delta) {

    for (double tau_2 = 0.0; tau_2 < tau_1; tau_2 += delta) {

      double value = free_energy_o2_trimer_contrib(ad, beta, tau_1, tau_2);
      integral += value * delta * delta;
    }
  }

  return integral;
}

int main(int argc, char *argv[]) {

  double U    = 4.0;
  double mu   = 5.0;
  double beta = 0.5;

  triqs::hilbert_space::fundamental_operator_set fops = hubbard_atom::make_fops();

  triqs::operators::many_body_operator_generic<double> H0 = hubbard_atom::make_H0(U, mu);

  triqs::atom_diag::atom_diag<false> ad(H0, fops, {}); // atom_diag object

  //double Z0 = triqs::atom_diag::partition_function(ad, beta); // Z0
  double delta       = beta / 500.0;
  double free_energy = free_energy_o2_trimer(ad, beta, delta);

  std::cout << "Trimer second order for U = " << U << ", mu = " << mu << ", beta = " << beta << ": " << free_energy << std::endl;
  return 0;
}