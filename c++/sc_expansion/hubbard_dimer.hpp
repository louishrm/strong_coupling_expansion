#pragma once
#include <triqs/atom_diag.hpp>
#include <triqs/atom_diag/functions.hpp>

namespace hubbard_dimer {

  //unperturbed Hamiltonian
  auto make_H0(double t, double mu);

  auto boltz_weights_and_imag_time(auto energies, double tau1, double tau2);

  double unperturbed_green_function(double U, double mu, double beta, std::vector<std::pair<int, int>> path, std::vector<double> tau);
  //expectation value of the occupation operator at site `site_idx`
  double occupation(double U, double mu, double beta);

} // namespace hubbard_dimer