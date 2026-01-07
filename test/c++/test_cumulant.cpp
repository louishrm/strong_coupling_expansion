// #include <vector>
// #include <algorithm>
// #include <iostream>
// #include <set>
// #include <sstream>
// #include "../c++/sc_expansion/hubbard_atom.hpp"
// #include "../c++/sc_expansion/cumulant.hpp"
// #include "../c++/sc_expansion/diagram.hpp"
// #include <numeric>
// #include <chrono>

// using namespace sc_expansion;

// double dimer_Omega2a(auto ad, double beta, std::vector<double> tau) {

//   double G0_12 = compute_cumulant_decomposition({{tau[0], 0}}, {{tau[1], 0}}, ad, beta); //G(1|2)
//   double G0_21 = compute_cumulant_decomposition({{tau[1], 0}}, {{tau[0], 0}}, ad, beta); //G(2|1)

//   double sign              = 1.0;
//   double symmetry_factor   = 1 / 2.0;
//   double free_multiplicity = 2.0;
//   double spin_factor       = 2.0; // for spin degeneracy
//   double prefactor         = sign * symmetry_factor * free_multiplicity * spin_factor;

//   return prefactor * G0_12 * G0_21;
// }

// double dimer_Omega2(auto ad, double beta, double delta) {

//   double result = 0.0;

//   for (double tau2 = 0.0; tau2 <= beta; tau2 += delta) {

//     std::vector<double> tau = {0, tau2};

//     double result_a = dimer_Omega2a(ad, beta, tau);
//     result += result_a * delta;
//   }
//   return result;
// }

// double dimer_Omega4a(auto ad, double beta, std::vector<double> tau) {

//   double G01_14 = compute_cumulant_decomposition({{tau[0], 0}}, {{tau[3], 0}}, ad, beta); //G(1|4)
//   double G01_21 = compute_cumulant_decomposition({{tau[1], 0}}, {{tau[0], 0}}, ad, beta); //G(2|1)
//   double G01_32 = compute_cumulant_decomposition({{tau[2], 0}}, {{tau[1], 0}}, ad, beta); //G(3|2)
//   double G01_43 = compute_cumulant_decomposition({{tau[3], 0}}, {{tau[2], 0}}, ad, beta); //G(4|3)

//   double sign              = -1.0;
//   double symmetry_factor   = 1 / 4.0;
//   double free_multiplicity = 1.0;
//   double spin_factor       = 2.0; // for spin degeneracy
//   double prefactor         = sign * symmetry_factor * free_multiplicity * spin_factor;

//   return prefactor * G01_14 * G01_21 * G01_32 * G01_43;
// }

// double dimer_Omega4b(auto ad, double beta, std::vector<double> tau) {
//   double result = 0.0;

//   double C02_up = compute_cumulant_decomposition({{tau[0], 0}, {tau[1], 0}}, {{tau[2], 0}, {tau[3], 0}}, ad, beta); //C(2 4 | 1 3)

//   double G01_14_up = compute_cumulant_decomposition({{tau[2], 0}}, {{tau[0], 0}}, ad, beta);
//   double G01_32_up = compute_cumulant_decomposition({{tau[3], 0}}, {{tau[1], 0}}, ad, beta);
//   result += 2 * C02_up * G01_14_up * G01_32_up;

//   double C02_updown = compute_cumulant_decomposition({{tau[0], 0}, {tau[1], 1}}, {{tau[2], 0}, {tau[3], 1}}, ad, beta); //C(2 4 | 1 3)

//   double G01_14_up_2 = compute_cumulant_decomposition({{tau[2], 0}}, {{tau[0], 0}}, ad, beta); //G(1|4)
//   double G01_32_down = compute_cumulant_decomposition({{tau[3], 1}}, {{tau[1], 1}}, ad, beta); //G(3|2)

//   result += 2 * C02_updown * G01_14_up_2 * G01_32_down;

//   double sign              = 1.0;
//   double symmetry_factor   = 1.0 / 2.0;
//   double free_multiplicity = 1.0;
//   double prefactor         = sign * symmetry_factor * free_multiplicity;

//   return prefactor * result;
// }

// double dimer_Omega4c(auto ad, double beta, std::vector<double> tau) {

//   double result = 0.0;

//   result += 2 * compute_cumulant_decomposition({{tau[0], 0}, {tau[1], 0}}, {{tau[2], 0}, {tau[3], 0}}, ad, beta)
//      * compute_cumulant_decomposition({{tau[2], 0}, {tau[3], 0}}, {{tau[0], 0}, {tau[1], 0}}, ad, beta);

//   result += 2 * compute_cumulant_decomposition({{tau[0], 0}, {tau[1], 1}}, {{tau[2], 0}, {tau[3], 1}}, ad, beta)
//      * compute_cumulant_decomposition({{tau[2], 0}, {tau[3], 1}}, {{tau[0], 0}, {tau[1], 1}}, ad, beta);

//   result += 2 * compute_cumulant_decomposition({{tau[0], 1}, {tau[1], 0}}, {{tau[2], 0}, {tau[3], 1}}, ad, beta)
//      * compute_cumulant_decomposition({{tau[2], 0}, {tau[3], 1}}, {{tau[0], 1}, {tau[1], 0}}, ad, beta);

//   double sign              = 1.0;
//   double symmetry_factor   = 1 / 8.0;
//   double free_multiplicity = 1.0;
//   double prefactor         = sign * symmetry_factor * free_multiplicity;

//   return prefactor * result;
// }

// double dimer_Omega4(auto ad, double beta, double delta) {

//   double result = 0.0;

//   for (double tau2 = 0; tau2 < beta; tau2 += delta) {
//     for (double tau3 = 0; tau3 < beta; tau3 += delta) {
//       for (double tau4 = 0; tau4 < beta; tau4 += delta) {

//         std::vector<double> tau = {0, tau2, tau3, tau4};

//         double result_a = dimer_Omega4a(ad, beta, tau);
//         double result_b = dimer_Omega4b(ad, beta, tau);
//         double result_c = dimer_Omega4c(ad, beta, tau);
//         result += (result_a + result_b + result_c) * delta * delta * delta;
//       }
//     }
//   }
//   return result;
// }

// int main() {

//   double U    = 8.0;
//   double beta = 1.0;
//   double mu   = 2.0;

//   triqs::hilbert_space::fundamental_operator_set fops = hubbard_atom::make_fops();

//   triqs::operators::many_body_operator_generic<double> H0 = hubbard_atom::make_H0(U, mu);

//   triqs::atom_diag::atom_diag<false> ad(H0, fops, {}); // atom_diag object

//   double delta2 = beta / 50.0;
//   double F4     = beta * dimer_Omega4(ad, beta, delta2);
//   std::cout << "Free energy 4th order for U: " << U << ", mu: " << mu << ", beta: " << beta << " is " << F4 * beta << std::endl;

//   adjmat D4a   = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}; //4-cycle
//   auto omega4a = Diagram(D4a, U, beta, mu);
//   adjmat D4b   = {{0, 1, 1}, {1, 0, 0}, {1, 0, 0}}; //3-cycle with double lines
//   auto omega4b = Diagram(D4b, U, beta, mu);
//   adjmat D4c   = {{0, 2}, {2, 0}}; //2-cycle with double lines
//   auto omega4c = Diagram(D4c, U, beta, mu);

//   std::vector<double> tau = {0.2, 0.4, 0.1, 0.3};
//   double testa            = dimer_Omega4a(ad, beta, tau);
//   double testa2           = omega4a.evaluate_at_taus(tau);
//   double testb            = dimer_Omega4b(ad, beta, tau);
//   double testb2           = omega4b.evaluate_at_taus(tau);
//   double testc            = dimer_Omega4c(ad, beta, tau);
//   double testc2           = omega4c.evaluate_at_taus(tau);

//   std::cout << "Test Omega4a: " << testa << " " << testa2 << std::endl;
//   std::cout << "Test Omega4b: " << testb << " " << testb2 << std::endl;
//   std::cout << "Test Omega4c: " << testc << " " << testc2 << std::endl;

//   double C02_up = compute_cumulant_decomposition({{0.2, 0}, {0.4, 0}}, {{0.1, 0}, {0.3, 0}}, ad, beta); //C(2 4 | 1 3)

//   double G01_14_up = compute_cumulant_decomposition({{0.1, 0}}, {{0.2, 0}}, ad, beta);
//   double G01_32_up = compute_cumulant_decomposition({{0.3, 0}}, {{0.4, 0}}, ad, beta);
//   double test4b    = C02_up * G01_14_up * G01_32_up;

//   hubbard_atom::cumul_args args4b = {{tau[0], 0}, {tau[1], 0}, {tau[2], 0}, {tau[3], 0}};
//   double test4b2                  = omega4b.evaluate_at_points(args4b);

//   std::cout << "Test Omega4b part: " << test4b << " " << test4b2 << std::endl;

//   return 0;
// }
