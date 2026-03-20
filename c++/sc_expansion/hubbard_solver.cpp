#include "./hubbard_solver.hpp"
#include "./combinatorics.hpp"
#include "./dual.hpp"

namespace sc_expansion {

  template <int N_sites, typename T> HubbardSolver<N_sites, T>::HubbardSolver(Parameters<T> const &params_) : params(params_) {
    for (int i = 0; i < N_OPS; ++i) { this->operators[i] = FermionOperator<N_sites, T>(static_cast<uint8_t>(i)); }

    this->compute_eigenstates();
    this->compute_transition_table();

    for (int i = 0; i < N_OPS; ++i) { this->operator_matrices[i] = this->operators[i].compute_sparse_matrix(this->all_eigenstates); }

    // Precompute exp(-beta * E_i)
    using std::exp;
    for (int i = 0; i < N_STATES; ++i) { this->exp_beta_E[i] = exp(-this->params.beta * this->all_eigenstates[i].energy); }

    // Compute Z = sum(exp(-beta * E_i))
    this->Z = T(0.0);
    for (int i = 0; i < N_STATES; ++i) { this->Z = this->Z + this->exp_beta_E[i]; }

    // Compute Z_infinite_U
    this->Z_infinite_U = T(0.0);
    if constexpr (N_sites == 1) {
      for (int i = 0; i < 3; ++i) { this->Z_infinite_U = this->Z_infinite_U + this->exp_beta_E[i]; }
    } else {
      this->Z_infinite_U = this->Z; // Fallback
    }
  }

  template <int N_sites, typename T> void HubbardSolver<N_sites, T>::compute_eigenstates() {
    if constexpr (N_sites == 1) {
      // Basis: |0>, |down>, |up>, |up down>
      // Bit order: Bit 0 is DOWN, Bit 1 is UP
      this->all_eigenstates[0] = Eigenstate<T>{{{0, T(1.0)}}, T(0.0)};
      this->all_eigenstates[1] = Eigenstate<T>{{{1, T(1.0)}}, -this->params.mu};
      this->all_eigenstates[2] = Eigenstate<T>{{{2, T(1.0)}}, -this->params.mu};
      this->all_eigenstates[3] = Eigenstate<T>{{{3, T(1.0)}}, this->params.U - 2.0 * this->params.mu};
    } else if constexpr (N_sites == 2) {
      T t  = this->params.t;
      T U  = this->params.U;
      T mu = this->params.mu;

      this->all_eigenstates[0] = Eigenstate<T>{{{0, T(1.0)}}, T(0.0)};

      // N=1, Sz=-1/2 |down,0> ± |0,down>
      this->all_eigenstates[1] = Eigenstate<T>{{{1, T(SQRT2_INV)}, {2, T(SQRT2_INV)}}, -t - mu};
      this->all_eigenstates[2] = Eigenstate<T>{{{1, T(SQRT2_INV)}, {2, -T(SQRT2_INV)}}, t - mu};

      // N=1, Sz=+1/2 |up,0> ± |0,up>
      this->all_eigenstates[3] = Eigenstate<T>{{{4, T(SQRT2_INV)}, {8, T(SQRT2_INV)}}, -t - mu};
      this->all_eigenstates[4] = Eigenstate<T>{{{4, T(SQRT2_INV)}, {8, -T(SQRT2_INV)}}, t - mu};

      // N=2, Sz=-1 |down,down>
      this->all_eigenstates[5] = Eigenstate<T>{{{3, T(1.0)}}, -2.0 * mu};

      // N=2, Sz=0, parity = even
      T Ep = Eplus(t, U, mu);
      T Em = Eminus(t, U, mu);
      using std::sqrt;
      T norm_plus       = Ep * SQRT2_INV / (sqrt(Ep * Ep + 16.0 * t * t));
      T norm_minus      = Em * SQRT2_INV / (sqrt(Em * Em + 16.0 * t * t));
      T component_plus  = -2.0 * t / Ep * norm_plus;
      T component_minus = -2.0 * t / Em * norm_minus;

      this->all_eigenstates[6] = Eigenstate<T>{{{5, T(norm_plus)}, {10, T(norm_plus)}, {9, T(component_plus)}, {6, T(component_plus)}}, Ep};
      this->all_eigenstates[7] = Eigenstate<T>{{{5, T(norm_minus)}, {10, T(norm_minus)}, {9, T(component_minus)}, {6, T(component_minus)}}, Em};

      // N=2, Sz=0, parity = odd
      this->all_eigenstates[8] = Eigenstate<T>{{{5, T(SQRT2_INV)}, {10, -T(SQRT2_INV)}}, U - 2.0 * mu};
      this->all_eigenstates[9] = Eigenstate<T>{{{9, T(SQRT2_INV)}, {6, -T(SQRT2_INV)}}, -2.0 * mu};

      // N=2, Sz=+1 |up,up>
      this->all_eigenstates[10] = Eigenstate<T>{{{12, T(1.0)}}, -2.0 * mu};

      // N=3, Sz = -1/2, |down up, down> ± |down, down up>
      this->all_eigenstates[11] = Eigenstate<T>{{{7, T(SQRT2_INV)}, {11, T(SQRT2_INV)}}, U + t - 3.0 * mu};
      this->all_eigenstates[12] = Eigenstate<T>{{{7, T(SQRT2_INV)}, {11, -T(SQRT2_INV)}}, U - t - 3.0 * mu};

      // N=3, Sz = +1/2, |down up, up> ± |up, down up>
      this->all_eigenstates[13] = Eigenstate<T>{{{13, T(SQRT2_INV)}, {14, T(SQRT2_INV)}}, U + t - 3.0 * mu};
      this->all_eigenstates[14] = Eigenstate<T>{{{13, T(SQRT2_INV)}, {14, -T(SQRT2_INV)}}, U - t - 3.0 * mu};

      // N=4, Sz=0 |down up, down up>
      this->all_eigenstates[15] = Eigenstate<T>{{{15, T(1.0)}}, 2.0 * U - 4.0 * mu};
    }
  }

  template <int N_sites, typename T> void HubbardSolver<N_sites, T>::compute_transition_table() {
    auto is_zero = [](auto const &val) {
      if constexpr (std::is_same_v<std::decay_t<decltype(val)>, Dual>) {
        return std::abs(val.value) < 1e-15;
      } else {
        return std::abs(val) < 1e-15;
      }
    };

    for (int op_idx = 0; op_idx < N_OPS; ++op_idx) {
      auto const &op = this->operators[op_idx];
      for (int ket_idx = 0; ket_idx < N_STATES; ++ket_idx) {
        auto const &ket = this->all_eigenstates[ket_idx];
        for (int bra_idx = 0; bra_idx < N_STATES; ++bra_idx) {
          auto const &bra = this->all_eigenstates[bra_idx];
          T overlap       = T(0.0);

          for (const auto &[basis_idx, coeff] : ket.coefficients) {
            auto transition = op.act_on_state(FockState<N_sites>(basis_idx));
            if (transition.matrix_element == 0.0) continue;

            for (const auto &[bra_basis_idx, bra_coeff] : bra.coefficients) {
              if (bra_basis_idx == transition.connected_state) {
                overlap = overlap + coeff * T(transition.matrix_element) * bra_coeff;
                break;
              }
            }
          }

          if (!is_zero(overlap)) { this->transition_table[op_idx][ket_idx].transitions.push_back({bra_idx, overlap}); }
        }
      }
    }
  }

  template <int N_sites, typename T> T HubbardSolver<N_sites, T>::G0n(Args<N_sites, T> const &args) const {
    T result     = T(0.0);
    auto is_zero = [](auto const &val) {
      if constexpr (std::is_same_v<std::decay_t<decltype(val)>, Dual>) {
        return std::abs(val.value) < 1e-15;
      } else {
        return std::abs(val) < 1e-15;
      }
    };

    using std::exp;

    for (int start_state = 0; start_state < N_STATES; ++start_state) {
      std::array<T, N_STATES> amplitudes;
      amplitudes.fill(T(0.0));
      amplitudes[start_state] = this->exp_beta_E[start_state];

      // Act with operators from earliest to latest time
      for (int i = args.order - 1; i >= 0; --i) {
        std::array<T, N_STATES> next_amplitudes;
        next_amplitudes.fill(T(0.0));
        uint8_t op_idx = args.ops[i].op;

        for (auto const &entry : this->operator_matrices[op_idx].entries) {
          if (is_zero(amplitudes[entry.col])) continue;
          next_amplitudes[entry.row] = next_amplitudes[entry.row] + amplitudes[entry.col] * entry.value * exp(args.taus[i] * entry.deltaE);
        }
        amplitudes = next_amplitudes;
      }
      result = result + amplitudes[start_state];
    }
    return (T(1.0) / this->Z) * args.permutation_sign * result;
  }

  template <int N_sites, typename T> T HubbardSolver<N_sites, T>::G0n_infinite_U(Args<N_sites, T> const &args) const {

    T result = T(0.0);

    using std::exp;

    if (!(args.operator_sequence_is_valid_infinite_U())) { return result; }

    FermionOperator<N_sites, T> first_op = args.ops[0];
    if (first_op.get_action() == 0) {

      return (T(1.0) / this->Z_infinite_U) * args.permutation_sign;
    } else {
      return exp(this->params.beta * this->params.mu) * (T(1.0) / this->Z_infinite_U) * args.permutation_sign;
    }
  }

  template class HubbardSolver<1, double>;
  template class HubbardSolver<1, Dual>;
  template class HubbardSolver<2, double>;
  template class HubbardSolver<2, Dual>;

} // namespace sc_expansion
