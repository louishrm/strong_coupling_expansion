#pragma once
#include <array>
#include <cstdint>
#include <vector>
#include <tuple>
#include <utility>

namespace sc_expansion {

  struct FockState {
    int state; // bit 0: spin down, bit 1: spin up
    FockState(int state_);
    bool is_occupied(uint8_t orbital_index) const;
  };

  template <typename T> struct Transition {
    int connected_state;
    T matrix_element;
    Transition(int connected_state, T mel);
  };

  template <typename T> struct Eigenstate {
    std::vector<std::pair<int, T>> coefficients; // List of (basis state index, coefficient) pairs
    T energy;
  };

  template <typename T> struct TransitionList {
    // transition list
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

  template <typename T> struct FermionOperator {

    public:
    uint8_t op; // Bit 1 is the 'action' bit: 0 = destroy, 1 = create
    FermionOperator() = default;
    FermionOperator(uint8_t op_);
    uint8_t get_action() const;
    uint8_t get_orbital_index() const;
    Transition<T> act_on_state(FockState const &fock_state) const;
    T compute_matrix_element(Eigenstate<T> const &bra, Eigenstate<T> const &ket) const;

    FermionOperator<T> apply_spin_flip() const;

    static constexpr uint8_t N_STATES     = 4;
    static constexpr uint8_t ACTION_BIT   = 2;
    static constexpr uint8_t ORBITAL_MASK = ACTION_BIT - 1;

    SparseMatrix<T> compute_sparse_matrix(std::array<Eigenstate<T>, N_STATES> const &all_eigenstates) const;
  };
} // namespace sc_expansion
