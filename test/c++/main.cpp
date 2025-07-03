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

int main(int argc, char *argv[]) {

  double U    = 4.0;
  double mu   = 1.0;
  double beta = 1.0;

  triqs::hilbert_space::fundamental_operator_set fops = hubbard_atom::make_fops();

  triqs::operators::many_body_operator_generic<double> H0 = hubbard_atom::make_H0(U, mu);

  triqs::atom_diag::atom_diag<false> ad(H0, fops, {}); // atom_diag object

  //double Z0 = triqs::atom_diag::partition_function(ad, beta); // Z0

  double Z0 = 1 + 2 * std::exp(beta * mu) + std::exp(-beta * (U - 2 * mu));
  //Z0 *= std::exp(-beta * 2 * mu);
  double Z0test = hubbard_atom::partition_function(ad, beta); // Z0

  //double Z0test = triqs::atom_diag::partition_function(ad, beta); // Z0 from atom_diag

  std::vector<double> tau_test = {1.0, 0.3};
  std::vector<int> spins_test  = {0, 0};
  std::vector<int> flags_test  = {0, 1};

  std::vector<int> spins   = {0, 0};
  std::vector<int> flags_1 = {0, 1};
  std::vector<int> flags_2 = {1, 0};

  //std::vector<std::vector<double>> grid_vals;
  double delta    = 5e-4; // Time step for integration
  double integral = 0.0;

  for (double tau_1 = 0.0; tau_1 < beta; tau_1 += delta) {
    //std::vector<double> tau_2_row;

    for (double tau_2 = 0.0; tau_2 < tau_1; tau_2 += delta) {

      double G011 = hubbard_atom::G0(ad, beta, {tau_1, tau_2}, spins, flags_1);
      double G012 = hubbard_atom::G0(ad, beta, {tau_1, tau_2}, spins, flags_2);

      double value = G011 * G012;
      integral += value * delta * delta;

      //tau_2_row.push_back(value);
    }

    //grid_vals.push_back(tau_2_row);
  }

  integral *= 4 * Z0 * Z0; // Normalize by Z0^2

  double exact = o2_exact(U, mu, beta);
  std::cout << "Integral = " << integral << std::endl;
  std::cout << "Exact = " << exact << std::endl;

  return 0;
}