#pragma once
#include <array>
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
    bool bipartite = true;
    // Shift absorbed into the perturbation: H_pert = -t sum c^dag c + delta sum n.
    // HubbardSolver still diagonalises H_0 with the `mu` stored here; the caller
    // is responsible for passing the already-shifted chemical potential (mu - delta)
    // if they want the atomic reference to match the variational choice.
    T delta = T(0.0);
  };

  template <typename T> class HubbardSolver {

    public:
    Parameters<T> const &params;
    HubbardSolver(Parameters<T> const &params);

    T Z;
    T Z_infinite_U;

    T G0n(Args<T> const &args) const;
    T G0n_infinite_U(Args<T> const &args) const;
    T G01(Args<T> const &args) const; // closed-form one-body propagator

    // Accessors for testing
    T get_Z() const { return this->Z; }
    T get_exp_beta_E(int i) const { return this->exp_beta_E[i]; }
    const Eigenstate<T> &get_eigenstate(int i) const { return this->all_eigenstates[i]; }
    const SparseMatrix<T> &get_operator_matrix(int op_idx) const { return this->operator_matrices[op_idx]; }

    private:
    constexpr static int N_STATES      = 4;
    constexpr static int N_OPS         = 4;
    constexpr static int MAX_G0N_ORDER = 16; // max 2*cumulant_order supported in G0n

    using ExpTable = std::array<std::array<T, N_STATES>, MAX_G0N_ORDER>;

    std::array<FermionOperator<T>, N_OPS> operators;
    std::array<Eigenstate<T>, N_STATES> all_eigenstates;
    std::array<std::array<TransitionList<T>, N_STATES>, N_OPS> transition_table;
    std::array<SparseMatrix<T>, N_OPS> operator_matrices;
    std::array<T, N_STATES> exp_beta_E;
    std::array<T, N_STATES> eigenstate_energies; // E[s] cached contiguously for fast access in G0n

    void compute_eigenstates();
    void compute_transition_table();
    void build_tau_exp_tables(Args<T> const &args, ExpTable &fwd, ExpTable &inv) const;
  };

} // namespace sc_expansion
