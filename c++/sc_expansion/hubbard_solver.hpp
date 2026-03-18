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
#include <cmath>
#include "dual.hpp"

namespace sc_expansion {

  using Arg     = std::pair<double, int>;
  using ArgList = std::vector<Arg>;

  template <typename T> struct Parameters {
    T U;
    T beta;
    T mu;
    T t            = T(0.0); // Only used for the dimer, but included here for convenience
    bool bipartite = true;
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

  template <typename T> struct TransitionT {
    int connected_state;
    T matrix_element;
  };

  template <typename T> struct TransitionList {
    // transition list the dimer, operators connect to linear combination of eigenstates (max 6)
    std::vector<TransitionT<T>> transitions;
  };

  template <typename T> class HubbardAtom {

    public:
    static constexpr int N_STATES = 4;
    static constexpr int N_OPS    = 4;

    HubbardAtom(Parameters<T> const &params);

    T G0(std::vector<double> const &taus, std::vector<int> const &spins) const;
    T G0_infinite_U(std::vector<double> const &taus, std::vector<int> const &spins) const;

    // Internal members made public for testing and cumulant solver efficiency
    Parameters<T> const &params;

    T Z;
    T Z_infinite_U;
    std::array<T, N_STATES> E;

    static const std::array<Transition, N_STATES * N_OPS> lookup_table;

    static constexpr int valid_start_states[N_OPS][2] = {
       {1, 3}, // op 0 (c_up):       Needs state 1 (|up>) or 3 (|up down>)
       {2, 3}, // op 1 (c_down):     Needs state 2 (|down>) or 3 (|up down>)
       {0, 2}, // op 2 (cdag_up):    Needs state 0 (|0>) or 2 (|down>)
       {0, 1}  // op 3 (cdag_down):  Needs state 0 (|0>) or 1 (|up>)
    };
  };

  template <typename T> std::pair<T, int> fermion_operator_act(int op_index, int state);

  static constexpr double SQRT2_INV = 0.70710678118654752440;

  template <typename T> T Eplus(T t, T U, T mu) { return T(0.5) * (U + sqrt(U * U + 16.0 * t * t)) - 2.0 * mu; }
  template <typename T> T Eminus(T t, T U, T mu) { return T(0.5) * (U - sqrt(U * U + 16.0 * t * t)) - 2.0 * mu; }

  template <typename T> struct Eigenstate {
    std::vector<std::pair<int, T>> coefficients; // List of (basis state index, coefficient) pairs
    T energy;
    int particle_number;
  };

  template <typename T> class HubbardDimer {
    public:
    static constexpr int N_STATES = 16;
    static constexpr int N_OPS    = 8;
    HubbardDimer(Parameters<T> const &params, T t);

    // transition_table[op_index][eigenv_ket_index] -> list of {eigenv_bra_index, matrix_element}
    std::array<std::array<TransitionList<T>, N_STATES>, N_OPS> transition_table;

    private:
    Parameters<T> const &params;
    T t;
    T Z;
    std::array<Eigenstate<T>, N_STATES> all_eigenstates;

    void compute_eigenstates();
    void compute_transition_table();
  };

} // namespace sc_expansion
