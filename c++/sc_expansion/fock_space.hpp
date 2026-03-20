#pragma once
#include <cstdint>
#include <vector>
#include <tuple>
#include <utility>

namespace sc_expansion {

  template <int N_sites> struct FockState {
    int state; // lowest N_sites bits: spin down, next N_sites bits: spin up
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

  template <int N_sites, typename T> struct FermionOperator {

    public:
    uint8_t op; // Bit N_sites is the 'action' bit: 0 = destroy, 1 = create
    FermionOperator() = default;
    FermionOperator(uint8_t op_);
    uint8_t get_action() const;
    uint8_t get_orbital_index() const;
    Transition<T> act_on_state(FockState<N_sites> const &fock_state) const;
    T compute_matrix_element(Eigenstate<T> const &bra, Eigenstate<T> const &ket) const;

    template <typename Container> SparseMatrix<T> compute_sparse_matrix(Container const &eigenstates) const {
      SparseMatrix<T> matrix;
      for (std::size_t i = 0; i < eigenstates.size(); ++i) {
        for (std::size_t j = 0; j < eigenstates.size(); ++j) {
          T mel = this->compute_matrix_element(eigenstates[i], eigenstates[j]);
          if (mel != T(0.0)) { matrix.entries.push_back({static_cast<int>(i), static_cast<int>(j), mel, eigenstates[i].energy - eigenstates[j].energy}); }
        }
      }
      return matrix;
    }

    static constexpr uint8_t ACTION_BIT   = (1 << N_sites);
    static constexpr uint8_t ORBITAL_MASK = ACTION_BIT - 1;
    static constexpr uint8_t N_STATES     = 1 << (2 * N_sites); //4^N_sites
  };

} // namespace sc_expansion
