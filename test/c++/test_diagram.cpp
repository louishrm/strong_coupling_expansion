#include "../../c++/sc_expansion/diagram.hpp"
#include "../../c++/sc_expansion/hubbard_atom.hpp"
#include "../../c++/sc_expansion/cumulant.hpp"

#include <iostream>

double dimer_Omega2a(auto ad, double beta, std::vector<double> tau) {

  double G0_12 = compute_cumulant_decomposition({{tau[0], 0}}, {{tau[1], 0}}, ad, beta); //G(1|2)
  double G0_21 = compute_cumulant_decomposition({{tau[1], 0}}, {{tau[0], 0}}, ad, beta); //G(2|1)

  return G0_12 * G0_21;
}

int main() {

  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;

  triqs::hilbert_space::fundamental_operator_set fops     = hubbard_atom::make_fops();
  triqs::operators::many_body_operator_generic<double> H0 = hubbard_atom::make_H0(U, mu);
  triqs::atom_diag::atom_diag<false> ad(H0, fops, {}); // atom_diag object

  //adjmat adjmat_test = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}; //4-cycle
  //adjmat adjmat_test = {{0, 1, 1}, {1, 0, 0}, {1, 0, 0}}; //3-cycle with double lines
  adjmat adjmat_test = {{0, 1}, {1, 0}};
  Diagram D(adjmat_test);
  bool test1 = D.is_connected();
  bool test2 = D.is_particle_number_conserving();
  int sf     = D.get_symmetry_factor();
  std::cout << "Connected: " << test1 << ", PNC: " << test2 << ", Symmetry factor: " << sf << std::endl;

  auto hopping_lines = D.get_hopping_lines();

  for (auto line : hopping_lines) { std::cout << "Line from " << line.from_vertex << " to " << line.to_vertex << std::endl; }

  double exact = dimer_Omega2a(ad, beta, {0.0, 0.5});
  std::cout << "Exact Omega2a at (0,0.5): " << exact << std::endl;

  hubbard_atom::cumul_args args = {{0.0, 0}, {0.5, 0}};
  double eval                   = D.evaluate_at_points(ad, beta, args);
  std::cout << "Diagram evaluated at same points: " << eval << std::endl;
  return 0;
}