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

    // ⟨T_τ Π_k n_{σ_k}(0) Π_i O_i(τ_i)⟩ / Z, atomic only (N_sites == 1).
    // Each entry of `density_orbitals` is one static n_{σ}(0) insertion at τ=0;
    // 0, 1, or 2 entries supported (callers of the static density-density
    // correlator need at most 2). Since H_0 is diagonal in occupation for the
    // atomic case, n_σ commutes with everything and can sit OUTSIDE the
    // time-ordering — `args` is reused as-is (no re-sort, no extra sign).
    // For N_sites != 1 the implementation is a placeholder that asserts and
    // returns 0; cluster generalization is deferred.
    T G0n_with_densities(Args<N_sites, T> const &args, std::vector<int> const &density_orbitals) const;
    T G0n_with_densities_infinite_U(Args<N_sites, T> const &args, std::vector<int> const &density_orbitals) const;

    // ⟨n_σ⟩ = Tr(e^{-βH} c†_σ c_σ) / Z, computed directly from the
    // eigenbasis. Used by the rooted (correlator) cumulant evaluator to
    // bypass the τ=0 equal-time ambiguity that ⟨T_τ c_σ(0) c†_σ(0)⟩
    // would otherwise introduce when the mark block is folded into
    // the partition lattice as a sub-cumulant.
    // `orbital` is the FermionOperator orbital index (σ ∈ {0..N_ORBITALS}).
    T compute_n_sigma(int orbital) const;
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
    // Dense N_STATES×N_STATES matrix in the energy eigenbasis. Small (4 or 16
    // per side), so density operators are assembled and applied densely.
    using DenseMatrix = std::array<std::array<T, N_STATES>, N_STATES>;

    // Build N = Π_k n_{σ_k} (n_σ = c†_σ c_σ) as a dense matrix in the energy
    // eigenbasis. The factors commute, so build order is irrelevant. Used by
    // the cluster (N_sites != 1) density-insertion path, where n_σ is not
    // diagonal in the hybridized eigenbasis.
    DenseMatrix build_density_matrix(std::vector<int> const &density_orbitals) const;

    // Trace evaluator shared by the density paths: seeds each start state,
    // applies the τ=0 density matrix N (no evolution factor) at the right-most
    // trace slot, then the hybridization operators O_i exactly as G0n does.
    // Correct for any N_sites; the atomic case keeps a diagonal fast path only
    // for performance.
    T G0n_with_density_matrix(Args<N_sites, T> const &args, DenseMatrix const &N_matrix) const;

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
