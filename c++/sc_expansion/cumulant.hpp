#ifndef CUMULANT_HPP
#define CUMULANT_HPP

#include "hubbard_atom.hpp"
#include <vector>

double compute_cumulant_decomposition(const hubbard_atom::cumul_args &unprimed, const hubbard_atom::cumul_args &primed,
                                      const triqs::atom_diag::atom_diag<false> &ad, double beta);

#endif
