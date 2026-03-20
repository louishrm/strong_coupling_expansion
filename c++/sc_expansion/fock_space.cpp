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
