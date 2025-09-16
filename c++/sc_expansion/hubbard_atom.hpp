#pragma once
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/atom_diag/functions.hpp>
#include <vector>
#include <nda/nda.hpp>
#include <utility>
#include <iostream>
#include <algorithm>
#include <numeric>
using namespace triqs::operators;

namespace hubbard_atom {

  using cumul_args = std::vector<std::pair<double, int>>; //vector of pairs of imaginary time and spin

  //unperturbed Hamiltonian
  triqs::hilbert_space::fundamental_operator_set make_fops();

  triqs::operators::many_body_operator_generic<double> make_H0(double U, double mu);

  double partition_function(triqs::atom_diag::atom_diag<false> ad, double beta);

  int calculate_permutation_sign(const std::vector<int> &p);

  std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, int>
  sort_operators(const std::vector<double> &times, const std::vector<int> &spins, const std::vector<int> &flags);

  nda::matrix<double> make_interaction_picture_destroy_op(triqs::atom_diag::atom_diag<false> ad, double tau, int state_index);

  nda::matrix<double> make_interaction_picture_create_op(triqs::atom_diag::atom_diag<false> ad, double tau, int state_index);

  double G0(triqs::atom_diag::atom_diag<false> ad, double beta, cumul_args unprimed_args, cumul_args primed_args);

  // double C02(triqs::atom_diag::atom_diag<false> ad, double beta, std::vector<double> times, std::vector<int> spins, std::vector<int> flags);

} // namespace hubbard_atom