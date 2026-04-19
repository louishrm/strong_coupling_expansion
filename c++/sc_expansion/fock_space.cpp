#include "fock_space.hpp"
#include "dual.hpp"

namespace sc_expansion {

  // --- FockState ---
  FockState::FockState(int state_) { this->state = state_; }

  bool FockState::is_occupied(uint8_t orbital_index) const { return (this->state >> orbital_index) & 1; }

  // --- Transition ---
  template <typename T> Transition<T>::Transition(int state, T mel) : connected_state(state) { this->matrix_element = mel; }

  // --- FermionOperator ---
  template <typename T> FermionOperator<T>::FermionOperator(uint8_t op_) { this->op = op_; }

  template <typename T> uint8_t FermionOperator<T>::get_action() const { return (this->op >> 1) & 1; }

  template <typename T> uint8_t FermionOperator<T>::get_orbital_index() const { return this->op & ORBITAL_MASK; }

  template <typename T> Transition<T> FermionOperator<T>::act_on_state(FockState const &fock_state) const {
    uint8_t idx    = this->get_orbital_index();
    bool create    = (this->get_action() == 1);
    int occupation = (fock_state.state >> idx) & 1;

    // Pauli exclusion: cannot create if occupied, cannot destroy if empty
    if (create == (bool)occupation) return Transition<T>(-1, 0.0);

    int next_state = fock_state.state ^ (1 << idx);

    // Fermionic sign: count electrons to the right (lower indices)
    int count = __builtin_popcount(fock_state.state & ((1 << idx) - 1));
    T sign    = (count % 2 == 0) ? T(1.0) : T(-1.0);

    return Transition<T>(next_state, sign);
  }

  template <typename T> T FermionOperator<T>::compute_matrix_element(Eigenstate<T> const &bra, Eigenstate<T> const &ket) const {
    T overlap = T(0.0);

    for (const auto &[basis_idx, coeff] : ket.coefficients) {
      Transition<T> transition = this->act_on_state(FockState(basis_idx));
      if (transition.matrix_element == T(0.0)) continue;

      for (const auto &[bra_basis_idx, bra_coeff] : bra.coefficients) {
        if (bra_basis_idx == transition.connected_state) {
          overlap = overlap + coeff * transition.matrix_element * bra_coeff;
          break;
        }
      }
    }
    return overlap;
  }

  template <typename T>
  SparseMatrix<T> FermionOperator<T>::compute_sparse_matrix(std::array<Eigenstate<T>, N_STATES> const &all_eigenstates) const {
    SparseMatrix<T> res;
    auto is_zero = [](auto const &val) {
      if constexpr (std::is_same_v<std::decay_t<decltype(val)>, Dual>) {
        return std::abs(val.value) < 1e-15;
      } else {
        return std::abs(val) < 1e-15;
      }
    };

    for (int ket_idx = 0; ket_idx < (int)N_STATES; ++ket_idx) {
      for (int bra_idx = 0; bra_idx < (int)N_STATES; ++bra_idx) {
        T overlap = this->compute_matrix_element(all_eigenstates[bra_idx], all_eigenstates[ket_idx]);
        if (!is_zero(overlap)) {
          res.entries.push_back({bra_idx, ket_idx, overlap, all_eigenstates[bra_idx].energy - all_eigenstates[ket_idx].energy});
        }
      }
    }
    return res;
  }

  template <typename T> FermionOperator<T> FermionOperator<T>::apply_spin_flip() const {
    uint8_t action          = this->op & ACTION_BIT;
    uint8_t orbital         = this->get_orbital_index();
    uint8_t flipped_orbital = orbital ^ 1;
    return FermionOperator<T>(action | flipped_orbital);
  }

  // --- Explicit Instantiations ---

  template struct Transition<double>;
  template struct Transition<Dual>;

  template struct FermionOperator<double>;
  template struct FermionOperator<Dual>;

} // namespace sc_expansion
