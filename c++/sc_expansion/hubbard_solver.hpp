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
#include "hilbert_traits.hpp"

namespace sc_expansion {

  template <int N_sites, typename T> class HubbardSolver {

    public:
    using Traits                       = HilbertTraits<N_sites>;
    static constexpr int N_STATES      = Traits::N_STATES;
    static constexpr int N_OPS         = Traits::N_OPS;
    static constexpr int MAX_G0N_ORDER = 16; // max 2*cumulant_order supported in G0n

    Parameters<T> const &params;
    HubbardSolver(Parameters<T> const &params);

    T Z;
    T Z_infinite_U;

    T G0n(Args<N_sites, T> const &args) const;
    T G0n_infinite_U(Args<N_sites, T> const &args) const;
    // Closed-form fast path exposed for tests. Returns 0 when no closed form
    // applies (e.g. dimer, or higher orders). Production code should call G0n,
    // which dispatches to this when applicable.
    T G01(Args<N_sites, T> const &args) const;

    // Accessors for testing
    T get_Z() const { return this->Z; }
    T get_exp_beta_E(int i) const { return this->exp_beta_E[i]; }
    const Eigenstate<T> &get_eigenstate(int i) const { return this->all_eigenstates[i]; }
    const SparseMatrix<T> &get_operator_matrix(int op_idx) const { return this->operator_matrices[op_idx]; }

    private:
    using ExpTable = std::array<std::array<T, N_STATES>, MAX_G0N_ORDER>;

    std::array<FermionOperator<N_sites, T>, N_OPS> operators;
    std::array<Eigenstate<T>, N_STATES> all_eigenstates;
    std::array<std::array<TransitionList<T>, N_STATES>, N_OPS> transition_table;
    std::array<SparseMatrix<T>, N_OPS> operator_matrices;
    std::array<T, N_STATES> exp_beta_E;
    std::array<T, N_STATES> eigenstate_energies;

    void compute_transition_table();
    void build_tau_exp_tables(Args<N_sites, T> const &args, ExpTable &fwd, ExpTable &inv) const;
  };

} // namespace sc_expansion
