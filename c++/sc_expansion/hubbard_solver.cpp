#include "./hubbard_solver.hpp"
#include "./combinatorics.hpp"
#include "./dual.hpp"

namespace sc_expansion {

  namespace detail {
    template <typename V> inline bool is_zero(V const &val) {
      if constexpr (std::is_same_v<std::decay_t<V>, Dual>) {
        return std::abs(val.value) < 1e-15;
      } else {
        return std::abs(val) < 1e-15;
      }
    }
  } // namespace detail

  template <typename T> HubbardSolver<T>::HubbardSolver(Parameters<T> const &params_) : params(params_) {
    for (int i = 0; i < N_OPS; ++i) { this->operators[i] = FermionOperator<T>(static_cast<uint8_t>(i)); }

    this->compute_eigenstates();
    for (int i = 0; i < N_OPS; ++i) { this->operator_matrices[i] = this->operators[i].compute_sparse_matrix(this->all_eigenstates); }
    this->compute_transition_table();

    // Cache eigenstate energies contiguously for fast access in G0n
    for (int i = 0; i < N_STATES; ++i) { this->eigenstate_energies[i] = this->all_eigenstates[i].energy; }

    // Precompute exp(-beta * E_i)
    using std::exp;
    for (int i = 0; i < N_STATES; ++i) { this->exp_beta_E[i] = exp(-this->params.beta * this->eigenstate_energies[i]); }

    // Compute Z = sum(exp(-beta * E_i))
    this->Z = T(0.0);
    for (int i = 0; i < N_STATES; ++i) { this->Z = this->Z + this->exp_beta_E[i]; }

    // Compute Z_infinite_U (excludes doubly-occupied state)
    this->Z_infinite_U = T(0.0);
    for (int i = 0; i < 3; ++i) { this->Z_infinite_U = this->Z_infinite_U + this->exp_beta_E[i]; }
  }

  template <typename T> void HubbardSolver<T>::compute_eigenstates() {
    // Basis: |0>, |down>, |up>, |up down>
    // Bit order: Bit 0 is DOWN, Bit 1 is UP
    this->all_eigenstates[0] = Eigenstate<T>{{{0, T(1.0)}}, T(0.0)};
    this->all_eigenstates[1] = Eigenstate<T>{{{1, T(1.0)}}, -this->params.mu};
    this->all_eigenstates[2] = Eigenstate<T>{{{2, T(1.0)}}, -this->params.mu};
    this->all_eigenstates[3] = Eigenstate<T>{{{3, T(1.0)}}, this->params.U - 2.0 * this->params.mu};
  }

  template <typename T> void HubbardSolver<T>::compute_transition_table() {
    for (int op_idx = 0; op_idx < N_OPS; ++op_idx) {
      for (const auto &entry : this->operator_matrices[op_idx].entries) {
        this->transition_table[op_idx][entry.col].transitions.push_back({entry.row, entry.value});
      }
    }
  }

  template <typename T>
  void HubbardSolver<T>::build_tau_exp_tables(Args<T> const &args, ExpTable &fwd, ExpTable &inv) const {
    using std::exp;
    int n = args.order;
    for (int i = 0; i < n; ++i) {
      double tau_i = args.taus[i];
      for (int s = 0; s < N_STATES; ++s) {
        fwd[i][s] = exp(tau_i * this->eigenstate_energies[s]);
        inv[i][s] = T(1.0) / fwd[i][s];
      }
    }
  }

  template <typename T> T HubbardSolver<T>::G0n(Args<T> const &args) const {
    T result = T(0.0);

    if (!args.operator_sequence_is_valid()) { return result; }

    int n = args.order;

    // Fast path: closed-form one-body propagator.
    if (n == 2) return this->G01(args);

    // exp(tau * deltaE) = exp(tau * E_row) * inv_exp(tau * E_col), avoiding exp() in the inner loop.
    ExpTable exp_tau_E, inv_exp_tau_E;
    this->build_tau_exp_tables(args, exp_tau_E, inv_exp_tau_E);

    // Double-buffer for amplitudes (pointer swap instead of array copy)
    std::array<T, N_STATES> buf_a, buf_b;

    for (int start_state = 0; start_state < N_STATES; ++start_state) {
      T *amplitudes = buf_a.data();
      T *next       = buf_b.data();

      std::fill(amplitudes, amplitudes + N_STATES, T(0.0));
      amplitudes[start_state] = this->exp_beta_E[start_state];

      for (int i = n - 1; i >= 0; --i) {
        std::fill(next, next + N_STATES, T(0.0));
        uint8_t op_idx = args.ops[i].op;

        for (auto const &entry : this->operator_matrices[op_idx].entries) {
          if (detail::is_zero(amplitudes[entry.col])) continue;
          next[entry.row] = next[entry.row] + amplitudes[entry.col] * entry.value * exp_tau_E[i][entry.row] * inv_exp_tau_E[i][entry.col];
        }
        std::swap(amplitudes, next);
      }
      result = result + amplitudes[start_state];
    }
    return (T(1.0) / this->Z) * args.permutation_sign * result;
  }

  template <typename T> T HubbardSolver<T>::G0n_infinite_U(Args<T> const &args) const {

    T result = T(0.0);

    using std::exp;

    if (!(args.operator_sequence_is_valid_infinite_U())) { return result; }

    FermionOperator<T> first_op = args.ops[0];
    if (first_op.get_action() == 0) {

      return (T(1.0) / this->Z_infinite_U) * args.permutation_sign;
    } else {
      return exp(this->params.beta * this->params.mu) * (T(1.0) / this->Z_infinite_U) * args.permutation_sign;
    }
  }

  template <typename T> T HubbardSolver<T>::G01(Args<T> const &args) const {
    /* Closed-form one-body propagator.
       Operates on the (already time-sorted) args: taus[0] >= taus[1], delta = taus[0] - taus[1] >= 0.
       Works out, by inspecting which op sits at the larger tau, which raw trace applies:
         - annihilator first (sorted): Tr[e^{-bH} c(tau_1) c^dag(tau_1')]  (physical tau_1 > tau_1')
         - creator     first (sorted): Tr[e^{-bH} c^dag(tau_1') c(tau_1)]  (physical tau_1' > tau_1)
       The fermionic-input sign is applied via args.permutation_sign, matching G0n's convention. */

    using std::exp;

    if (!args.operator_sequence_is_valid()) return T(0.0);

    T const &mu   = this->params.mu;
    T const &U    = this->params.U;
    T const &beta = this->params.beta;
    double delta  = args.taus[0] - args.taus[1];
    bool c_first  = (args.ops[0].get_action() == 0); // action 0 = annihilator

    T trace;
    if (c_first) {
      // Tr[e^{-bH} c(tau_1) c^dag(tau_1')], delta = tau_1 - tau_1' > 0
      trace = exp(mu * delta) + exp(beta * mu + (mu - U) * delta);
    } else {
      // Tr[e^{-bH} c^dag(tau_1') c(tau_1)], delta = tau_1' - tau_1 > 0
      trace = exp(beta * mu - mu * delta) + exp(-beta * (U - T(2.0) * mu) + (U - mu) * delta);
    }
    return args.permutation_sign * trace / this->Z;
  }

  template class HubbardSolver<double>;
  template class HubbardSolver<Dual>;

} // namespace sc_expansion
