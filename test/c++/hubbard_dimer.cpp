#include <triqs/atom_diag.hpp>
#include <iostream>
#include <triqs/atom_diag/functions.hpp>
#include "./hubbard_dimer.hpp"

// Prepare funcdamental operator set
triqs::hilbert_space::fundamental_operator_set make_fops() {
  triqs::hilbert_space::fundamental_operator_set fops;
  for (int o : range(1)) fops.insert("dn", o);
  for (int o : range(1)) fops.insert("up", o);
  return fops;
}

int dimer_groundstate() {

  //triqs::hilbert_space::fundamental_operator_set fops = make_fops();
  double t = 1.0;
  double U = 1.0;
  //std::cout << "t = " << t << ", U = " << U << endl;
  return 0;
}