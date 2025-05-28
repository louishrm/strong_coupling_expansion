#pragma once
#include <triqs/atom_diag.hpp>
#include <triqs/atom_diag/functions.hpp>

namespace hubbard_dimer {

  // Prepare fundamental operator set
  triqs::hilbert_space::fundamental_operator_set make_fops();

  // Main function to compute dimer ground state
  int dimer_groundstate();

} // namespace hubbard_dimer