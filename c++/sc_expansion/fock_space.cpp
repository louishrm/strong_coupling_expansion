#include "fock_space.hpp"
#include "dual.hpp"

namespace sc_expansion {

  // --- FockState ---
  template <int N_sites> FockState<N_sites>::FockState(int state_) { this->state = state_; }

  template <int N_sites> bool FockState<N_sites>::is_occupied(uint8_t orbital_index) const { return (this->state >> orbital_index) & 1; }

  // --- Transition ---
  template <typename T> Transition<T>::Transition(int state, T mel) : connected_state(state) { this->matrix_element = mel; }

  // --- FermionOperator ---
  template <int N_sites, typename T> FermionOperator<N_sites, T>::FermionOperator(uint8_t op_) { this->op = op_; }

  template <int N_sites, typename T> uint8_t FermionOperator<N_sites, T>::get_action() const { return (this->op >> N_sites) & 1; }

  template <int N_sites, typename T> uint8_t FermionOperator<N_sites, T>::get_orbital_index() const { return this->op & ORBITAL_MASK; }

  template <int N_sites, typename T> Transition<T> FermionOperator<N_sites, T>::act_on_state(FockState<N_sites> const &fock_state) const {
    uint8_t idx    = this->get_orbital_index();
    bool create    = (this->get_action() == 1);
    int occupation = (fock_state.state >> idx) & 1;

    // Pauli exclusion: cannot create if occupied, cannot destroy if empty
    if (create == (bool)occupation) return Transition<T>(-1, 0.0);

    int next_state = fock_state.state ^ (1 << idx);

    // Fermionic sign: count electrons to the right (lower indices)
    // Basis order: [up_N-1 ... up_0 | down_N-1 ... down_0]
    int count = __builtin_popcount(fock_state.state & ((1 << idx) - 1));
    T sign    = (count % 2 == 0) ? T(1.0) : T(-1.0);

    return Transition<T>(next_state, sign);
  }

  template <int N_sites, typename T> T FermionOperator<N_sites, T>::compute_matrix_element(Eigenstate<T> const &bra, Eigenstate<T> const &ket) const {
    T overlap = T(0.0);

    for (const auto &[basis_idx, coeff] : ket.coefficients) {
      Transition<T> transition = this->act_on_state(FockState<N_sites>(basis_idx));
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

  template <int N_sites, typename T>
  SparseMatrix<T> FermionOperator<N_sites, T>::compute_sparse_matrix(std::array<Eigenstate<T>, N_STATES> const &all_eigenstates) const {
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

  template <int N_sites, typename T> FermionOperator<N_sites, T> FermionOperator<N_sites, T>::apply_spin_flip() const {
    uint8_t action          = this->op & ACTION_BIT;
    uint8_t orbital         = this->get_orbital_index();
    uint8_t flipped_orbital = orbital ^ N_sites;
    return FermionOperator<N_sites, T>(action | flipped_orbital);
  }

  template <int N_sites, typename T> FermionOperator<N_sites, T> FermionOperator<N_sites, T>::apply_reflection() const {
    if constexpr (N_sites == 1) {
      return *this;
    } else {
      //flip the rightmost bit 0<--->1
      uint8_t new_op_index = this->op ^ 1;
      return FermionOperator<N_sites, T>(new_op_index);
    }
  }
  // --- Explicit Instantiations ---

  template struct FockState<1>;
  template struct FockState<2>;

  template struct Transition<double>;
  template struct Transition<Dual>;

  template struct FermionOperator<1, double>;
  template struct FermionOperator<2, double>;

  template struct FermionOperator<1, Dual>;
  template struct FermionOperator<2, Dual>;

} // namespace sc_expansion
