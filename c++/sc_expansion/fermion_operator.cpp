#include "fermion_operator.hpp"
#include "dual.hpp"

namespace sc_expansion {

  // --- FockState ---
  template <int N_sites> FockState<N_sites>::FockState(int state_) { this->state = state_; }

  template <int N_sites> bool FockState<N_sites>::is_occupied(uint8_t orbital_index) const {
    return (this->state >> orbital_index) & 1;
  }

  // --- Transition ---
  template <int N_sites, typename T>
  Transition<N_sites, T>::Transition(FockState<N_sites> state, T mel) : connected_state(state) {
    this->matrix_element = mel;
  }

  // --- FermionOperator ---
  template <int N_sites> FermionOperator<N_sites>::FermionOperator(uint8_t op_) { this->op = op_; }

  template <int N_sites> uint8_t FermionOperator<N_sites>::get_action() const { return (this->op >> N_sites) & 1; }

  template <int N_sites> uint8_t FermionOperator<N_sites>::get_orbital_index() const { return this->op & ORBITAL_MASK; }

  template <int N_sites> Transition<N_sites, double> FermionOperator<N_sites>::act_on_state(FockState<N_sites> const &fock_state) const {
    uint8_t idx    = this->get_orbital_index();
    bool create    = (this->get_action() == 1);
    int occupation = (fock_state.state >> idx) & 1;

    // Pauli exclusion: cannot create if occupied, cannot destroy if empty
    if (create == (bool)occupation) return Transition<N_sites, double>(FockState<N_sites>(-1), 0.0);

    int next_state = fock_state.state ^ (1 << idx);

    // Fermionic sign: count electrons to the right (lower indices)
    // Basis order: [up_N-1 ... up_0 | down_N-1 ... down_0]
    int count   = __builtin_popcount(fock_state.state & ((1 << idx) - 1));
    double sign = (count % 2 == 0) ? 1.0 : -1.0;

    return Transition<N_sites, double>(FockState<N_sites>(next_state), sign);
  }

  // --- Explicit Instantiations ---
  template struct FockState<1>;
  template struct FockState<2>;

  template struct Transition<1, double>;
  template struct Transition<2, double>;
  template struct Transition<1, Dual>;
  template struct Transition<2, Dual>;

  template struct FermionOperator<1>;
  template struct FermionOperator<2>;

} // namespace sc_expansion
