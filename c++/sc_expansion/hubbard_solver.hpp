#pragma once
#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <cmath>
#include "dual.hpp"
#include "args.hpp"
#include "fock_space.hpp"

namespace sc_expansion {

  template <typename T> struct Parameters {
    T U;
    T beta;
    T mu;
    T t            = T(0.0); // Only used for the dimer, but included here for convenience
    bool bipartite = true;
  };

  template <typename T> inline T Eplus(T t, T U, T mu) {
    using std::sqrt;
    return T(0.5) * (U + sqrt(U * U + 16.0 * t * t)) - 2.0 * mu;
  }

  template <typename T> inline T Eminus(T t, T U, T mu) {
    using std::sqrt;
    return T(0.5) * (U - sqrt(U * U + 16.0 * t * t)) - 2.0 * mu;
  }

  template <int N_sites, typename T> class HubbardSolver {

    public:
    Parameters<T> const &params;
    HubbardSolver(Parameters<T> const &params);

    T Z;
    T Z_infinite_U;

    T G0n(Args<N_sites, T> const &args) const;
    T G0n_infinite_U(Args<N_sites, T> const &args) const;

    // Accessors for testing
    T get_Z() const { return this->Z; }
    T get_exp_beta_E(int i) const { return this->exp_beta_E[i]; }
    const Eigenstate<T>& get_eigenstate(int i) const { return this->all_eigenstates[i]; }
    const SparseMatrix<T>& get_operator_matrix(int op_idx) const { return this->operator_matrices[op_idx]; }

    private:
    constexpr static int N_STATES     = 1 << (2 * N_sites); // 4 states for atom, 16 for dimer
    constexpr static int N_OPS        = 4 * N_sites;        // 4 ops for atom, 8 for dimer
    static constexpr double SQRT2_INV = 0.70710678118654752440;
    std::array<FermionOperator<N_sites, T>, N_OPS> operators;
    std::array<Eigenstate<T>, N_STATES> all_eigenstates;
    std::array<std::array<TransitionList<T>, N_STATES>, N_OPS> transition_table;
    std::array<SparseMatrix<T>, N_OPS> operator_matrices;
    std::array<T, N_STATES> exp_beta_E;

    void compute_eigenstates();
    void compute_transition_table();
  };

} // namespace sc_expansion
