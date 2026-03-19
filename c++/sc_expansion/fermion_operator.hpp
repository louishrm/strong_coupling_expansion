#pragma once
#include <cstdint>

namespace sc_expansion {

  template <int N_sites> struct FockState {
    int state; // lowest N_sites bits: spin down, next N_sites bits: spin up
    FockState(int state_);
    bool is_occupied(uint8_t orbital_index) const;
  };

  template <int N_sites, typename T> struct Transition {
    FockState<N_sites> connected_state;
    T matrix_element;
    Transition(FockState<N_sites> state, T mel);
  };

  template <int N_sites> struct FermionOperator {
    uint8_t op; // Bit N_sites is the 'action' bit: 0 = destroy, 1 = create
    static constexpr uint8_t ACTION_BIT   = (1 << N_sites);
    static constexpr uint8_t ORBITAL_MASK = ACTION_BIT - 1;

    FermionOperator(uint8_t op_);
    uint8_t get_action() const;
    uint8_t get_orbital_index() const;
    Transition<N_sites, double> act_on_state(FockState<N_sites> const &fock_state) const;
  };

} // namespace sc_expansion
