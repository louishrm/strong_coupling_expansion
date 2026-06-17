#pragma once
#include <array>
#include <cstdint>
#include <vector>
#include <tuple>
#include <utility>

namespace sc_expansion {

  // FockState<N_sites>: bit layout is [low N_sites bits = spin down, next N_sites bits = spin up].
  // For N_sites=1 this matches the legacy single-site convention (bit 0 = down, bit 1 = up).
  template <int N_sites> struct FockState {
    int state;
    FockState(int state_);
    bool is_occupied(uint8_t orbital_index) const;
  };

  template <typename T> struct Transition {
    int connected_state;
    T matrix_element;
    Transition(int connected_state, T mel);
  };

  template <typename T> struct Eigenstate {
    std::vector<std::pair<int, T>> coefficients;
    T energy;
  };

  template <typename T> struct TransitionList {
    std::vector<Transition<T>> transitions;
  };

  template <typename T> struct SparseMatrix {
    struct Entry {
      int row, col;
      T value;
      T deltaE;
    };
    std::vector<Entry> entries;
  };

  // FermionOperator<N_sites, T>: encodes (action, orbital_index) in `op` byte.
  //   ACTION_BIT   = 1 << N_sites (bit set => create)
  //   ORBITAL_MASK = ACTION_BIT - 1  (low bits hold the orbital index 0..2*N_sites-1)
  //   N_STATES     = 4^N_sites      (Hilbert dimension)
  // Note: ACTION_BIT == 1<<N_sites gives N_sites bits for the orbital, which is enough
  // for 2*N_sites orbitals iff 2*N_sites <= 2^N_sites — true for N_sites >= 1.
  template <int N_sites, typename T> struct FermionOperator {

    public:
    uint8_t op;
    FermionOperator() = default;
    FermionOperator(uint8_t op_);
    uint8_t get_action() const;
    uint8_t get_orbital_index() const;
    Transition<T> act_on_state(FockState<N_sites> const &fock_state) const;
    T compute_matrix_element(Eigenstate<T> const &bra, Eigenstate<T> const &ket) const;

    FermionOperator<N_sites, T> apply_spin_flip() const;
    // Site-swap reflection. No-op for N_sites=1; swaps site 0 <-> site 1 for N_sites=2.
    FermionOperator<N_sites, T> apply_reflection() const;

    static constexpr uint8_t N_STATES     = 1 << (2 * N_sites);
    static constexpr uint8_t ACTION_BIT   = 1 << N_sites;
    static constexpr uint8_t ORBITAL_MASK = ACTION_BIT - 1;

    SparseMatrix<T> compute_sparse_matrix(std::array<Eigenstate<T>, N_STATES> const &all_eigenstates) const;
  };
} // namespace sc_expansion
