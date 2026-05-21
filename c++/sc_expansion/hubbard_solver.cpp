#include "./hubbard_solver.hpp"
#include "./combinatorics.hpp"
#include "./dual.hpp"

namespace sc_expansion {

  template <int N_sites, typename T> HubbardSolver<N_sites, T>::HubbardSolver(Parameters<T> const &params_) : params(params_) {
    for (int i = 0; i < N_OPS; ++i) { this->operators[i] = FermionOperator<N_sites, T>(static_cast<uint8_t>(i)); }

    Traits::diagonalize(this->params, this->all_eigenstates);

    for (int i = 0; i < N_OPS; ++i) { this->operator_matrices[i] = this->operators[i].compute_sparse_matrix(this->all_eigenstates); }
    this->compute_transition_table();

    for (int i = 0; i < N_STATES; ++i) { this->eigenstate_energies[i] = this->all_eigenstates[i].energy; }

    using std::exp;
    for (int i = 0; i < N_STATES; ++i) { this->exp_beta_E[i] = exp(-this->params.beta * this->eigenstate_energies[i]); }

    this->Z = T(0.0);
    for (int i = 0; i < N_STATES; ++i) { this->Z = this->Z + this->exp_beta_E[i]; }

    this->Z_infinite_U = Traits::compute_Z_infinite_U(this->exp_beta_E);
  }

  template <int N_sites, typename T> void HubbardSolver<N_sites, T>::compute_transition_table() {
    for (int op_idx = 0; op_idx < N_OPS; ++op_idx) {
      for (const auto &entry : this->operator_matrices[op_idx].entries) {
        this->transition_table[op_idx][entry.col].transitions.push_back({entry.row, entry.value});
      }
    }
  }

  template <int N_sites, typename T>
  void HubbardSolver<N_sites, T>::build_tau_exp_tables(Args<N_sites, T> const &args, ExpTable &fwd, ExpTable &inv) const {
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

  template <int N_sites, typename T> T HubbardSolver<N_sites, T>::G0n(Args<N_sites, T> const &args) const {
    if (!args.operator_sequence_is_valid()) return T(0.0);

    // Per-cluster fast paths (e.g. atomic n=2 closed form). Returns false on miss.
    {
      T fast;
      if (Traits::try_closed_form_G0n(args, this->params, this->Z, fast)) return fast;
    }

    int n = args.order;

    // exp(tau * deltaE) = exp(tau * E_row) * inv_exp(tau * E_col), avoiding exp() in the inner loop.
    ExpTable exp_tau_E, inv_exp_tau_E;
    this->build_tau_exp_tables(args, exp_tau_E, inv_exp_tau_E);

    std::array<T, N_STATES> buf_a, buf_b;

    T result = T(0.0);
    for (int start_state = 0; start_state < N_STATES; ++start_state) {
      T *amplitudes = buf_a.data();
      T *next       = buf_b.data();

      std::fill(amplitudes, amplitudes + N_STATES, T(0.0));
      amplitudes[start_state] = this->exp_beta_E[start_state];

      for (int i = n - 1; i >= 0; --i) {
        std::fill(next, next + N_STATES, T(0.0));
        uint8_t op_idx = args.ops[i].op;

        for (auto const &entry : this->operator_matrices[op_idx].entries) {
          if (is_zero(amplitudes[entry.col])) continue;
          next[entry.row] = next[entry.row] + amplitudes[entry.col] * entry.value * exp_tau_E[i][entry.row] * inv_exp_tau_E[i][entry.col];
        }
        std::swap(amplitudes, next);
      }
      result = result + amplitudes[start_state];
    }
    return (T(1.0) / this->Z) * args.permutation_sign * result;
  }

  template <int N_sites, typename T> T HubbardSolver<N_sites, T>::G0n_infinite_U(Args<N_sites, T> const &args) const {
    using std::exp;

    if (!(args.operator_sequence_is_valid_infinite_U())) return T(0.0);

    // Single-site infinite-U short-circuit. For N_sites>=2 the validity check above
    // already returns false, so this branch is dead but compiles uniformly.
    FermionOperator<N_sites, T> first_op = args.ops[0];
    if (first_op.get_action() == 0) {
      return (T(1.0) / this->Z_infinite_U) * args.permutation_sign;
    } else {
      return exp(this->params.beta * this->params.mu) * (T(1.0) / this->Z_infinite_U) * args.permutation_sign;
    }
  }

  template <int N_sites, typename T>
  T HubbardSolver<N_sites, T>::G0n_with_densities(Args<N_sites, T> const &args, std::vector<int> const &density_orbitals) const {
    if constexpr (N_sites != 1) {
      // Cluster (N_sites >= 2) generalization is a separate workstream.
      static_cast<void>(args);
      static_cast<void>(density_orbitals);
      return T(0.0);
    } else {
      if (density_orbitals.empty()) return this->G0n(args);
      if (!args.operator_sequence_is_valid()) return T(0.0);

      int n = args.order;

      // For atomic, n_σ is diagonal in the eigenbasis (states 0..3 are pure
      // occupation bitstrings). Precompute the per-start-state density factor.
      std::array<T, N_STATES> density_factor;
      for (int s = 0; s < N_STATES; ++s) {
        T factor = T(1.0);
        for (int orb : density_orbitals) { factor = factor * T((s & (1 << orb)) ? 1.0 : 0.0); }
        density_factor[s] = factor;
      }

      if (n == 0) {
        // Pure density expectation, no time-evolution kernel: just weight by Boltzmann.
        T result = T(0.0);
        for (int s = 0; s < N_STATES; ++s) { result = result + density_factor[s] * this->exp_beta_E[s]; }
        return (T(1.0) / this->Z) * args.permutation_sign * result;
      }

      ExpTable exp_tau_E, inv_exp_tau_E;
      this->build_tau_exp_tables(args, exp_tau_E, inv_exp_tau_E);

      std::array<T, N_STATES> buf_a, buf_b;
      T result = T(0.0);
      for (int start_state = 0; start_state < N_STATES; ++start_state) {
        if (is_zero(density_factor[start_state])) continue;

        T *amplitudes = buf_a.data();
        T *next       = buf_b.data();

        std::fill(amplitudes, amplitudes + N_STATES, T(0.0));
        amplitudes[start_state] = this->exp_beta_E[start_state];

        for (int i = n - 1; i >= 0; --i) {
          std::fill(next, next + N_STATES, T(0.0));
          uint8_t op_idx = args.ops[i].op;

          for (auto const &entry : this->operator_matrices[op_idx].entries) {
            if (is_zero(amplitudes[entry.col])) continue;
            next[entry.row] = next[entry.row] + amplitudes[entry.col] * entry.value * exp_tau_E[i][entry.row] * inv_exp_tau_E[i][entry.col];
          }
          std::swap(amplitudes, next);
        }
        result = result + density_factor[start_state] * amplitudes[start_state];
      }
      return (T(1.0) / this->Z) * args.permutation_sign * result;
    }
  }

  template <int N_sites, typename T>
  T HubbardSolver<N_sites, T>::G0n_with_densities_infinite_U(Args<N_sites, T> const &args, std::vector<int> const &density_orbitals) const {
    if constexpr (N_sites != 1) {
      static_cast<void>(args);
      static_cast<void>(density_orbitals);
      return T(0.0);
    } else {
      if (density_orbitals.empty()) return this->G0n_infinite_U(args);
      if (!args.operator_sequence_is_valid_infinite_U()) return T(0.0);

      int n = args.order;

      // Same diagonal-occupation density factor as the finite-U path; the
      // doubly-occupied state (s == 3) is then automatically excluded below.
      std::array<T, N_STATES> density_factor;
      for (int s = 0; s < N_STATES; ++s) {
        T factor = T(1.0);
        for (int orb : density_orbitals) { factor = factor * T((s & (1 << orb)) ? 1.0 : 0.0); }
        density_factor[s] = factor;
      }

      // Restrict to the projected sector {|0>, |down>, |up>} = states 0, 1, 2.
      auto in_projected = [](int s) { return s != 3; };

      if (n == 0) {
        T result = T(0.0);
        for (int s = 0; s < N_STATES; ++s) {
          if (!in_projected(s)) continue;
          result = result + density_factor[s] * this->exp_beta_E[s];
        }
        return (T(1.0) / this->Z_infinite_U) * args.permutation_sign * result;
      }

      ExpTable exp_tau_E, inv_exp_tau_E;
      this->build_tau_exp_tables(args, exp_tau_E, inv_exp_tau_E);

      std::array<T, N_STATES> buf_a, buf_b;
      T result = T(0.0);
      for (int start_state = 0; start_state < N_STATES; ++start_state) {
        if (!in_projected(start_state)) continue;
        if (is_zero(density_factor[start_state])) continue;

        T *amplitudes = buf_a.data();
        T *next       = buf_b.data();

        std::fill(amplitudes, amplitudes + N_STATES, T(0.0));
        amplitudes[start_state] = this->exp_beta_E[start_state];

        for (int i = n - 1; i >= 0; --i) {
          std::fill(next, next + N_STATES, T(0.0));
          uint8_t op_idx = args.ops[i].op;

          for (auto const &entry : this->operator_matrices[op_idx].entries) {
            if (!in_projected(entry.row)) continue;
            if (is_zero(amplitudes[entry.col])) continue;
            next[entry.row] = next[entry.row] + amplitudes[entry.col] * entry.value * exp_tau_E[i][entry.row] * inv_exp_tau_E[i][entry.col];
          }
          std::swap(amplitudes, next);
        }
        result = result + density_factor[start_state] * amplitudes[start_state];
      }
      return (T(1.0) / this->Z_infinite_U) * args.permutation_sign * result;
    }
  }

  template <int N_sites, typename T> T HubbardSolver<N_sites, T>::G01(Args<N_sites, T> const &args) const {
    T out = T(0.0);
    Traits::try_closed_form_G0n(args, this->params, this->Z, out);
    return out;
  }

  template <int N_sites, typename T> T HubbardSolver<N_sites, T>::compute_n_sigma(int orbital) const {
    // ⟨n_σ⟩ = Σ_α (e^{-βE_α}/Z) ⟨α| c†_σ c_σ |α⟩
    //       = Σ_α (e^{-βE_α}/Z) Σ_β c_σ_{βα} · c†_σ_{αβ}
    constexpr uint8_t ACTION_BIT_LOCAL = FermionOperator<N_sites, T>::ACTION_BIT;
    uint8_t c_op_id      = static_cast<uint8_t>(orbital);
    uint8_t c_dag_op_id  = static_cast<uint8_t>(ACTION_BIT_LOCAL | orbital);

    T result = T(0.0);
    for (int alpha = 0; alpha < N_STATES; ++alpha) {
      // Apply c_σ to |α⟩; record the resulting amplitudes by state.
      std::array<T, N_STATES> after_c{};
      for (auto const &e : this->operator_matrices[c_op_id].entries) {
        if (e.col == alpha) after_c[e.row] = after_c[e.row] + e.value;
      }
      // ⟨α| c†_σ |β⟩ for the resulting β states.
      T diag = T(0.0);
      for (auto const &e : this->operator_matrices[c_dag_op_id].entries) {
        if (e.row == alpha) diag = diag + e.value * after_c[e.col];
      }
      result = result + this->exp_beta_E[alpha] * diag;
    }
    return result / this->Z;
  }

  template class HubbardSolver<1, double>;
  template class HubbardSolver<1, Dual>;
  template class HubbardSolver<2, double>;
  template class HubbardSolver<2, Dual>;

} // namespace sc_expansion
