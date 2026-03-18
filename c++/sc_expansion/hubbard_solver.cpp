#include "./hubbard_solver.hpp"
#include "./combinatorics.hpp"
#include "./dual.hpp"

namespace sc_expansion {

  template <typename T>
  const std::array<Transition, HubbardAtom<T>::N_STATES * HubbardAtom<T>::N_OPS> HubbardAtom<T>::lookup_table = []() {
    std::array<Transition, HubbardAtom<T>::N_STATES * HubbardAtom<T>::N_OPS> table;

    // Iterate through all 4 states (00, 01, 10, 11)
    for (int state = 0; state < HubbardAtom<T>::N_STATES; ++state) {

      // Iterate through all 4 operators (Bit 0: Action, Bit 1: Spin)
      // 0=c_dn, 1=cdag_dn, 2=c_up, 3=cdag_up
      for (int op = 0; op < HubbardAtom<T>::N_OPS; ++op) {

        // The 1D array index: (State * 4) + Operator
        int index = (state * HubbardAtom<T>::N_OPS) | op;

        int next_state = 0;
        double mel     = 0.0;

        // Decode the operator type
        bool is_create = (op & 1) == 1;     // Bit 0 is action (1 = create)
        bool is_down   = ((op >> 1) & 1) == 0; // Bit 1 is spin (0 = down)

        // Decode the current state's orbital occupancies
        // Bit 0 (right) is UP, Bit 1 (left) is DOWN
        int n_up = state & 1;
        int n_dn = (state >> 1) & 1;

        int target_orbital_occ = is_down ? n_dn : n_up;

        // 1. Pauli Exclusion / Annihilation Check
        if (is_create && target_orbital_occ == 1) {
          mel = 0.0; // Cannot create if orbital is already full
        } else if (!is_create && target_orbital_occ == 0) {
          mel = 0.0; // Cannot destroy if orbital is already empty
        } else {
          // The state survives the operator!

          // 2. Calculate New State (Toggle the target bit using XOR)
          int bit_to_flip = is_down ? 2 : 1; // 2 in binary is 10, 1 is 01
          next_state      = state ^ bit_to_flip;

          // 3. Calculate Fermionic Sign
          // Basis rule: |up down> = cdag_up cdag_down |0>
          // Any 'down' operator must jump over the 'up' orbital.
          if (is_down && n_up == 1) {
            mel = -1.0; // Jumped over an occupied up-electron
          } else {
            mel = 1.0; // No jump, or jumped over an empty orbital
          }
        }

        // Save to the table
        table[index] = {next_state, mel};
      }
    }

    return table;
  }();

  template <typename T> HubbardAtom<T>::HubbardAtom(Parameters<T> const &params_) : params(params_) {
    this->Z = T(1.0) + 2.0 * exp(params.beta * params.mu) + exp(params.beta * (2.0 * params.mu - params.U)); // Partition function at finite U
    this->Z_infinite_U = T(1.0) + 2.0 * exp(params.beta * params.mu);                                        // Partition function at infinite U
    this->E            = {T(0.0), -params.mu, -params.mu, params.U - 2.0 * params.mu};
  }

  template <typename T> T HubbardAtom<T>::G0(std::vector<double> const &taus, std::vector<int> const &ops) const {
    OperatorSequence seq(taus, ops);
    T G0_value = T(0.0);

    if (!seq.verify_flavor_alternation(0)) return T(0.0);

    // get the only two states that the list op in the sorted list can act on
    int last_op = seq.ops.back();

    for (int initial_state : this->valid_start_states[last_op]) {

      int current_state   = initial_state;
      T current_trace_val = T(1.0);

      for (int i = seq.order - 1; i >= 0; --i) {

        int table_index = current_state * HubbardAtom<T>::N_OPS + seq.ops[i];
        Transition t    = this->lookup_table[table_index];

        if (t.matrix_element == 0.0) {
          current_trace_val = T(0.0);
          break;
        }
        int next_state    = t.connected_state;
        T energy_diff     = this->E[next_state] - this->E[current_state];
        current_trace_val = current_trace_val * t.matrix_element * exp(seq.taus[i] * energy_diff);
        current_state     = next_state;
      }

      // Add to value for Dual properly
      if (current_state == initial_state) {
        current_trace_val = current_trace_val * exp(-this->params.beta * this->E[current_state]);
        G0_value          = G0_value + current_trace_val;
      }
    }

    return T(1.0) / this->Z * seq.permutation_sign * G0_value;
  }

  template <typename T> T HubbardAtom<T>::G0_infinite_U(std::vector<double> const &taus, std::vector<int> const &ops) const {
    OperatorSequence seq(taus, ops);
    T G0_value = T(0.0);

    if (!seq.verify_infinite_U_constraints()) { return T(0.0); }

    // If the sequence is valid, we can directly compute the contribution from the first operator.
    // The contribution from the rest of the sequence is guaranteed to be 1 due to the strict alternation and spin conservation rules.

    int first_op = seq.ops[0];
    if ((first_op & 1) == 1) { // Action is Bit 0 (1 = create)
      // First operator is a creation operator
      G0_value = exp(this->params.beta * this->params.mu) / this->Z_infinite_U;
    } else {
      // First operator is a destruction operator
      G0_value = T(1.0) / this->Z_infinite_U;
    }

    return seq.permutation_sign * G0_value;
  }

  template <typename T> std::pair<T, int> fermion_operator_act(int op, int state) {
    int action  = (op >> 2) & 1; // 0 destroy, 1 create
    int orbital = op & 3;        // 0: 1down 1: 2down 2: 1up 3: 2up

    // check occupation of the orbital
    int occupation = (state >> orbital) & 1;
    if (occupation == action) { // the bits have to be different for the operator to act non-trivially
      return {T(0.0), -1};      // zero mel
    }

    int new_state = state ^ (1 << orbital);

    unsigned int mask    = (1 << orbital) - 1;
    int electrons_jumped = std::popcount(static_cast<unsigned int>(state) & mask); // Number of occupied orbitals before the target orbital

    return {(electrons_jumped % 2 == 0 ? T(1.0) : T(-1.0)), new_state};
  }

  template <typename T> HubbardDimer<T>::HubbardDimer(Parameters<T> const &params_, T t_) : params(params_), t(t_) {
    this->compute_eigenstates();
    this->compute_transition_table();

    // Compute Z = sum(exp(-beta * E_i))
    this->Z = T(0.0);
    for (const auto &state : this->all_eigenstates) { this->Z = this->Z + exp(-this->params.beta * state.energy); }
  }

  template <typename T> void HubbardDimer<T>::compute_eigenstates() {

    // |n> = |n2u n1u n2d n1d> = bit(n2u, n1u, n2d, n1d)

    this->all_eigenstates[0] = Eigenstate<T>{{{0, T(1.0)}}, T(0.0), 0};

    // N=1, Sz=-1/2 |down,0> ± |0,down>
    this->all_eigenstates[1] = Eigenstate<T>{{{1, T(SQRT2_INV)}, {2, T(SQRT2_INV)}}, -this->params.t - this->params.mu, 1};
    this->all_eigenstates[2] = Eigenstate<T>{{{1, T(SQRT2_INV)}, {2, -T(SQRT2_INV)}}, this->params.t - this->params.mu, 1};

    // N=1, Sz=+1/2 |up,0> ± |0,up>
    this->all_eigenstates[3] = Eigenstate<T>{{{4, T(SQRT2_INV)}, {8, T(SQRT2_INV)}}, -this->params.t - this->params.mu, 1};
    this->all_eigenstates[4] = Eigenstate<T>{{{4, T(SQRT2_INV)}, {8, -T(SQRT2_INV)}}, this->params.t - this->params.mu, 1};

    // N=2, Sz=-1 |down,down>
    this->all_eigenstates[5] = Eigenstate<T>{{{3, T(1.0)}}, -2.0 * this->params.mu, 2};

    // N=2, Sz=0, parity = even
    // a(|down up, 0> + |0, down up>) + b(|down, up> + |up, down>)
    T Ep              = Eplus(this->t, this->params.U, this->params.mu);
    T Em              = Eminus(this->t, this->params.U, this->params.mu);
    T norm_plus       = Ep * SQRT2_INV / (sqrt(Ep * Ep + 16.0 * this->t * this->t));
    T norm_minus      = Em * SQRT2_INV / (sqrt(Em * Em + 16.0 * this->t * this->t));
    T component_plus  = -2.0 * this->t / Ep * norm_plus;
    T component_minus = -2.0 * this->t / Em * norm_minus;

    this->all_eigenstates[6] = Eigenstate<T>{{{5, T(norm_plus)}, {10, T(norm_plus)}, {9, T(component_plus)}, {6, T(component_plus)}}, Ep, 2};
    this->all_eigenstates[7] = Eigenstate<T>{{{5, T(norm_minus)}, {10, T(norm_minus)}, {9, T(component_minus)}, {6, T(component_minus)}}, Em, 2};

    // N=2, Sz=0, parity = odd
    //|down up, 0> - |0, down up> and |up, down> - |down,up>
    this->all_eigenstates[8] = Eigenstate<T>{{{5, T(SQRT2_INV)}, {10, -T(SQRT2_INV)}}, this->params.U - 2.0 * this->params.mu, 2};
    this->all_eigenstates[9] = Eigenstate<T>{{{9, T(SQRT2_INV)}, {6, -T(SQRT2_INV)}}, -2.0 * this->params.mu, 2};

    // N=2, Sz=+1 |up,up>
    this->all_eigenstates[10] = Eigenstate<T>{{{12, T(1.0)}}, -2.0 * this->params.mu, 2};

    // N=3, Sz = -1/2, |down up, down> ± |down, down up>
    this->all_eigenstates[11] = Eigenstate<T>{{{7, T(SQRT2_INV)}, {11, T(SQRT2_INV)}}, this->params.U + this->params.t - 3.0 * this->params.mu, 3};
    this->all_eigenstates[12] =
       Eigenstate<T>{{{7, T(SQRT2_INV)}, {11, -T(SQRT2_INV)}}, this->params.U - this->params.t - 3.0 * this->params.mu, 3};

    // N=3, Sz = +1/2, |down up, up> ± |up, down up>
    this->all_eigenstates[13] =
       Eigenstate<T>{{{13, T(SQRT2_INV)}, {14, T(SQRT2_INV)}}, this->params.U + this->params.t - 3.0 * this->params.mu, 3};
    this->all_eigenstates[14] =
       Eigenstate<T>{{{13, T(SQRT2_INV)}, {14, -T(SQRT2_INV)}}, this->params.U - this->params.t - 3.0 * this->params.mu, 3};

    // N=4, Sz=0 |down up, down up>
    this->all_eigenstates[15] = Eigenstate<T>{{{15, T(1.0)}}, 2.0 * this->params.U - 4.0 * this->params.mu, 4};
  }

  template <typename T> void HubbardDimer<T>::compute_transition_table() {

    for (int op_index = 0; op_index < this->N_OPS; ++op_index) {
      for (int eigenv_ket_index = 0; eigenv_ket_index < this->N_STATES; ++eigenv_ket_index) {

        for (int eigenv_bra_index = 0; eigenv_bra_index < this->N_STATES; ++eigenv_bra_index) {

          T matrix_element = T(0.0);
          for (auto const &[ket_basis_index, ket_coeff] : this->all_eigenstates[eigenv_ket_index].coefficients) {
            auto [mel, new_state] = fermion_operator_act<T>(op_index, ket_basis_index);

            if (new_state != -1) {
              for (auto const &[bra_basis_index, bra_coeff] : this->all_eigenstates[eigenv_bra_index].coefficients) {
                if (bra_basis_index == new_state) { matrix_element = matrix_element + bra_coeff * mel * ket_coeff; }
              }
            }
          }

          auto is_zero = [](auto const &val) {
            if constexpr (std::is_same_v<std::decay_t<decltype(val)>, Dual>) { return std::abs(val.value) < 1e-15; }
            else {
              return std::abs(val) < 1e-15;
            }
          };

          if (!is_zero(matrix_element)) {
            this->transition_table[op_index][eigenv_ket_index].transitions.push_back({eigenv_bra_index, matrix_element});
          }
        }
      }
    }
  }

  template class HubbardAtom<double>;
  template class HubbardAtom<Dual>;

  template class HubbardDimer<double>;
  template class HubbardDimer<Dual>;

} // namespace sc_expansion
