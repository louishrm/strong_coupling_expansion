#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/atom_diag/functions.hpp>
#include <vector>
#include <nda/nda.hpp>
#include <utility>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <tuple>

namespace sc_expansion {

  template <typename T> class HubbardDimer {

    public:
    HubbardDimer(T t, T U, T beta, T mu);

    T Z;
    T G0(std::vector<double> const &taus, std::vector<int> const &spins, std::vector<int> const &sites) const;
  };
} // namespace sc_expansion