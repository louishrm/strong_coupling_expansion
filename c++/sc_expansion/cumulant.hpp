#ifndef CUMULANT_HPP
#define CUMULANT_HPP

#include "hubbard_atom.hpp"
#include <vector>

namespace sc_expansion {

  double compute_cumulant_decomposition(HubbardAtom::cumul_args const &unprimed, HubbardAtom::cumul_args const &primed, HubbardAtom const &atom);

}

#endif