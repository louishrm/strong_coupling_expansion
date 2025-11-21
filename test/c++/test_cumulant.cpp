#include <vector>
#include <algorithm>
#include <iostream>
#include <set>
#include <sstream>
#include "../c++/sc_expansion/hubbard_atom.hpp"
#include "../c++/sc_expansion/cumulant.hpp"
#include <numeric>
#include <chrono>

double dimer_Omega2a(auto ad, double beta, std::vector<double> tau) {

  double G0_12 = compute_cumulant_decomposition({{tau[0], 0}}, {{tau[1], 0}}, ad, beta); //G(1|2)
  double G0_21 = compute_cumulant_decomposition({{tau[1], 0}}, {{tau[0], 0}}, ad, beta); //G(2|1)

  double sign              = 1.0;
  double symmetry_factor   = 1 / 2.0;
  double free_multiplicity = 2.0;
  double spin_factor       = 2.0; // for spin degeneracy
  double prefactor         = sign * symmetry_factor * free_multiplicity * spin_factor;

  return prefactor * G0_12 * G0_21;
}

double dimer_Omega2(auto ad, double beta, double delta) {

  double result = 0.0;

  for (double tau2 = 0.0; tau2 <= beta; tau2 += delta) {

    std::vector<double> tau = {0, tau2};

    double result_a = dimer_Omega2a(ad, beta, tau);
    result += result_a * delta;
  }
  return result;
}

double dimer_Omega4a(auto ad, double beta, std::vector<double> tau) {

  double G01_14 = compute_cumulant_decomposition({{tau[0], 0}}, {{tau[3], 0}}, ad, beta); //G(1|4)
  double G01_21 = compute_cumulant_decomposition({{tau[1], 0}}, {{tau[0], 0}}, ad, beta); //G(2|1)
  double G01_32 = compute_cumulant_decomposition({{tau[2], 0}}, {{tau[1], 0}}, ad, beta); //G(3|2)
  double G01_43 = compute_cumulant_decomposition({{tau[3], 0}}, {{tau[2], 0}}, ad, beta); //G(4|3)

  double sign              = -1.0;
  double symmetry_factor   = 1 / 4.0;
  double free_multiplicity = 2.0;
  double spin_factor       = 2.0; // for spin degeneracy
  double prefactor         = sign * symmetry_factor * free_multiplicity * spin_factor;

  return prefactor * G01_14 * G01_21 * G01_32 * G01_43;
}

double dimer_Omega4b(auto ad, double beta, std::vector<double> tau) {
  double result = 0.0;

  double C02_up = compute_cumulant_decomposition({{tau[1], 0}, {tau[3], 0}}, {{tau[0], 0}, {tau[2], 0}}, ad, beta); //C(2 4 | 1 3)

  double G01_14_up = compute_cumulant_decomposition({{tau[0], 0}}, {{tau[1], 0}}, ad, beta);
  double G01_32_up = compute_cumulant_decomposition({{tau[2], 0}}, {{tau[3], 0}}, ad, beta);
  result += 2 * C02_up * G01_14_up * G01_32_up;

  double C02_updown = compute_cumulant_decomposition({{tau[1], 0}, {tau[3], 1}}, {{tau[0], 0}, {tau[2], 1}}, ad, beta); //C(2 4 | 1 3)

  double G01_14_up_2 = compute_cumulant_decomposition({{tau[0], 0}}, {{tau[1], 0}}, ad, beta); //G(1|4)
  double G01_32_down = compute_cumulant_decomposition({{tau[2], 1}}, {{tau[3], 1}}, ad, beta); //G(3|2)

  result += 2 * C02_updown * G01_14_up_2 * G01_32_down;

  double sign              = 1.0;
  double symmetry_factor   = 1.0 / 2.0;
  double free_multiplicity = 2.0;
  double prefactor         = sign * symmetry_factor * free_multiplicity;

  return prefactor * result;
}

double dimer_Omega4c(auto ad, double beta, std::vector<double> tau) {

  double result = 0.0;

  result += 2 * compute_cumulant_decomposition({{tau[1], 0}, {tau[3], 0}}, {{tau[0], 0}, {tau[2], 0}}, ad, beta)
     * compute_cumulant_decomposition({{tau[0], 0}, {tau[2], 0}}, {{tau[1], 0}, {tau[3], 0}}, ad, beta);

  result += 2 * compute_cumulant_decomposition({{tau[1], 0}, {tau[3], 1}}, {{tau[0], 0}, {tau[2], 1}}, ad, beta)
     * compute_cumulant_decomposition({{tau[0], 0}, {tau[2], 1}}, {{tau[1], 0}, {tau[3], 1}}, ad, beta);

  result += 2 * compute_cumulant_decomposition({{tau[1], 1}, {tau[3], 0}}, {{tau[0], 0}, {tau[2], 1}}, ad, beta)
     * compute_cumulant_decomposition({{tau[0], 0}, {tau[2], 1}}, {{tau[1], 1}, {tau[3], 0}}, ad, beta);

  double sign              = 1.0;
  double symmetry_factor   = 1 / 8.0;
  double free_multiplicity = 2.0;
  double prefactor         = sign * symmetry_factor * free_multiplicity;

  return prefactor * result;
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

int main() {

  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;

  triqs::hilbert_space::fundamental_operator_set fops = hubbard_atom::make_fops();

  triqs::operators::many_body_operator_generic<double> H0 = hubbard_atom::make_H0(U, mu);

  triqs::atom_diag::atom_diag<false> ad(H0, fops, {}); // atom_diag object

  auto start              = std::chrono::high_resolution_clock::now();
  std::vector<double> tau = {0.0, beta / 3.0, 2.0 * beta / 3.0, beta};

  double G01_14 = 0.0;
  for (int i = 0; i < 100; i++) {
    G01_14 += compute_cumulant_decomposition({{tau[0], 0}, {tau[2], 0}}, {{tau[1], 0}, {tau[3], 0}}, ad, beta); //G(1|4)
  }
  G01_14 /= 100.0;
  std::cout << "Res 4th order cumulant " << G01_14 << std::endl;
  auto end                              = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Elapsed time: " << elapsed.count() << std::endl;
  std::cout << "4th order cumulant time: " << elapsed.count() / 100.0 << " seconds" << std::endl;

  // double delta = beta / 1000.0;
  // double F2    = dimer_Omega2(ad, beta, delta);

  // std::cout << "Free energy 2nd order for U: " << U << ", mu: " << mu << ", beta: " << beta << " is " << F2 * beta << std::endl;

  //   double delta2 = beta / 100.0;
  //   double F4     = beta * dimer_Omega4(ad, beta, delta2);
  //   std::cout << "Free energy 4th order for U: " << U << ", mu: " << mu << ", beta: " << beta << " is " << F4 * beta << std::endl;
  return 0;
}