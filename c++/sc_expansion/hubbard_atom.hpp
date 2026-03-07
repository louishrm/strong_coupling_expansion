#pragma once

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

  using Arg     = std::pair<double, int>;
  using ArgList = std::vector<Arg>;

  template <typename T>
  struct Parameters {
    T U;
    T beta;
    T mu;
  };

  struct Args {

    std::vector<double> taus;
    std::vector<int> spins;
    std::vector<int> ops; // 0: cup, 1: cdn, 2: cdag_up, 3: cdag_dn
    double permutation_sign;
    int order;

    Args(std::vector<double> taus, std::vector<int> spins);

    void sort_args();
    bool verify_consecutive_terms_infinite_U() const;
  };

  struct Transition {

    int connected_state;
    double matrix_element;
  };

  template <typename T>
  class HubbardAtom {

    public:
    HubbardAtom(T U, T beta, T mu);

    T G0(std::vector<double> const &taus, std::vector<int> const &spins) const;
    T G0_infinite_U(std::vector<double> const &taus, std::vector<int> const &spins) const;

    // Internal members made public for testing and cumulant solver efficiency
    T U;
    T beta;
    T mu;

    T Z;
    T Z_infinite_U;
    std::array<T, 4> E;

    static const std::array<Transition, 16> lookup_table; 

    static constexpr int valid_start_states[4][2] = {
       {1, 3}, // op 0 (c_up):       Needs state 1 (|up>) or 3 (|up down>)
       {2, 3}, // op 1 (c_down):     Needs state 2 (|down>) or 3 (|up down>)
       {0, 2}, // op 2 (cdag_up):    Needs state 0 (|0>) or 2 (|down>)
       {0, 1}  // op 3 (cdag_down):  Needs state 0 (|0>) or 1 (|up>)
    };
  };
} // namespace sc_expansion
