//#include "sc_expansion/sc_expansion.hpp"
#include "../c++/sc_expansion/hubbard_atom.hpp"

#include <cmath>
#include <iostream>
#include <xfac/tensor/tensor_ci_2.h>

double dimer_Omega4a(auto ad, double beta, std::vector<double> tau) {

  std::vector<int> spins = {0, 0};
  std::vector<int> flags = {1, 0};

  double G01_14 = hubbard_atom::G0(ad, beta, {tau[0], tau[3]}, spins, flags); //G(1|4)
  double G01_21 = hubbard_atom::G0(ad, beta, {tau[1], tau[0]}, spins, flags); //G(2|1)
  double G01_32 = hubbard_atom::G0(ad, beta, {tau[2], tau[1]}, spins, flags); //G(3|2)
  double G01_43 = hubbard_atom::G0(ad, beta, {tau[3], tau[2]}, spins, flags); //G(4|3)

  double sign              = -1.0;
  double symmetry_factor   = 1 / 4.0;
  double free_multiplicity = 2.0;
  double spin_factor       = 2.0; // for spin degeneracy
  double prefactor         = sign * symmetry_factor * free_multiplicity * spin_factor;

  return prefactor * G01_14 * G01_21 * G01_32 * G01_43;
}

double dimer_Omega4b(auto ad, double beta, std::vector<double> tau) {
  double result          = 0.0;
  std::vector<int> flags = {1, 1, 0, 0};
  std::vector<int> spins = {0, 0, 0, 0};

  double C02_up    = hubbard_atom::C02(ad, beta, {tau[1], tau[3], tau[0], tau[2]}, spins, flags); //C(2 4 | 1 3)
  double G01_14_up = hubbard_atom::G0(ad, beta, {tau[0], tau[1]}, {0, 0}, {1, 0});                //G(1|4)
  double G01_32_up = hubbard_atom::G0(ad, beta, {tau[2], tau[3]}, {0, 0}, {1, 0});                //G(4|1)
  result += 2 * C02_up * G01_14_up * G01_32_up;

  double C02_updown  = hubbard_atom::C02(ad, beta, {tau[1], tau[3], tau[0], tau[2]}, {0, 1, 0, 1}, flags); //C(2 4 | 1 3)
  double G01_14_up_2 = hubbard_atom::G0(ad, beta, {tau[0], tau[1]}, {0, 0}, {1, 0});                       //G(1|4)
  double G01_32_down = hubbard_atom::G0(ad, beta, {tau[2], tau[3]}, {1, 1}, {1, 0});                       //G(3|2)
  result += 2 * C02_updown * G01_14_up_2 * G01_32_down;

  double sign              = 1.0;
  double symmetry_factor   = 1.0 / 2.0;
  double free_multiplicity = 2.0;
  double prefactor         = sign * symmetry_factor * free_multiplicity;

  return prefactor * result;
}

double dimer_Omega4c(auto ad, double beta, std::vector<double> tau) {

  double result          = 0.0;
  std::vector<int> flags = {1, 1, 0, 0};

  result += 2 * hubbard_atom::C02(ad, beta, {tau[1], tau[3], tau[0], tau[2]}, {0, 0, 0, 0}, flags)
     * hubbard_atom::C02(ad, beta, {tau[0], tau[2], tau[1], tau[3]}, {0, 0, 0, 0}, flags);

  result += 2 * hubbard_atom::C02(ad, beta, {tau[1], tau[3], tau[0], tau[2]}, {0, 1, 0, 1}, flags)
     * hubbard_atom::C02(ad, beta, {tau[0], tau[2], tau[1], tau[3]}, {0, 1, 0, 1}, flags);

  result += 2 * hubbard_atom::C02(ad, beta, {tau[1], tau[3], tau[0], tau[2]}, {1, 0, 0, 1}, flags)
     * hubbard_atom::C02(ad, beta, {tau[0], tau[2], tau[1], tau[3]}, {0, 1, 1, 0}, flags);

  double sign              = 1.0;
  double symmetry_factor   = 1 / 8.0;
  double free_multiplicity = 2.0;
  double prefactor         = sign * symmetry_factor * free_multiplicity;

  return prefactor * result;
}

double dimer_Omega2a(auto ad, double beta, std::vector<double> tau) {

  std::vector<int> spins = {0, 0};
  std::vector<int> flags = {1, 0};

  double G01_12 = hubbard_atom::G0(ad, beta, {tau[0], tau[1]}, spins, flags); //G(1|2)
  double G01_21 = hubbard_atom::G0(ad, beta, {tau[1], tau[0]}, spins, flags); //G(2|1)

  double sign              = 1.0;
  double symmetry_factor   = 1 / 2.0;
  double free_multiplicity = 2.0;
  double spin_factor       = 2.0; // for spin degeneracy
  double prefactor         = sign * symmetry_factor * free_multiplicity * spin_factor;

  return prefactor * G01_12 * G01_21;
}

double dimer_Omega2(auto ad, double beta, double delta) {

  double result = 0.0;

  for (double tau2 = 0.0; tau2 <= beta; tau2 += delta) {
    std::vector<double> tau = {0, tau2};
    double result_a         = dimer_Omega2a(ad, beta, tau);
    result += result_a * delta;
  }
  return result;
}

double dimer_Omega4(auto ad, double beta, double delta) {

  double result = 0.0;

  for (double tau2 = 0.0; tau2 <= beta; tau2 += delta) {
    for (double tau3 = 0.0; tau3 <= beta; tau3 += delta) {
      for (double tau4 = 0.0; tau4 <= beta; tau4 += delta) {

        std::vector<double> tau = {0, tau2, tau3, tau4};

        double result_a = dimer_Omega4a(ad, beta, tau);
        double result_b = dimer_Omega4b(ad, beta, tau);
        double result_c = dimer_Omega4c(ad, beta, tau);
        result += (result_a + result_b + result_c) * delta * delta * delta;
      }
    }
  }
  return result;
}

double Omega8a(auto ad, double beta, std::vector<double> tau) {

  //8th order polygon diagram with all vertices having same spin

  std::vector<int> spins = {0, 0};
  std::vector<int> flags = {1, 0};
  double result          = 1.0;

  //tau pairs = {tau[0], tau[-1]}, {tau[2], tau[0]}
  int n = tau.size();
  std::vector<std::vector<double>> tau_pairs;
  for (int i = 0; i < n; i++) {
    std::vector<double> pair = {tau[i], tau[(i + n - 1) % n]};
    tau_pairs.push_back(pair);

    result *= hubbard_atom::G0(ad, beta, pair, spins, flags); //G(i|i-1)
  }

  return result;
}

int main(int argc, char *argv[]) {

  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;

  triqs::hilbert_space::fundamental_operator_set fops = hubbard_atom::make_fops();

  triqs::operators::many_body_operator_generic<double> H0 = hubbard_atom::make_H0(U, mu);

  triqs::atom_diag::atom_diag<false> ad(H0, fops, {}); // atom_diag object

  xfac::TensorCI2Param params;
  std::cout << params.bondDim << std::endl;

  // double delta = beta / 100.0;

  // double free_energy = dimer_Omega4(ad, beta, delta);
  // // double free_energy = dimer_Omega2(ad, beta, delta);
  // free_energy *= beta;
  // std::cout << "Free energy order 4 dimer for U = " << U << ", mu = " << mu << ", beta = " << beta << ": " << free_energy << std::endl;

  std::vector<double> tau = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
  double result           = Omega8a(ad, beta, tau);
  std::cout << "Omega8a: " << result << std::endl;

  return 0;
}