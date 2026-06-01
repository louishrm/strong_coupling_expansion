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
  typename HubbardSolver<N_sites, T>::DenseMatrix
  HubbardSolver<N_sites, T>::build_density_matrix(std::vector<int> const &density_orbitals) const {
    constexpr uint8_t ACTION_BIT_LOCAL = FermionOperator<N_sites, T>::ACTION_BIT;

    // Start from the identity; multiply in each n_σ. The n_σ commute, so the
    // accumulation order does not matter.
    DenseMatrix acc{};
    for (int i = 0; i < N_STATES; ++i) { acc[i][i] = T(1.0); }

    for (int orb : density_orbitals) {
      uint8_t c_idx     = static_cast<uint8_t>(orb);
      uint8_t c_dag_idx = static_cast<uint8_t>(ACTION_BIT_LOCAL | orb);

      // n_σ = M(c†_σ) · M(c_σ) in the eigenbasis (both are sparse eigenbasis
      // matrices held by the solver).
      DenseMatrix n_sigma{};
      for (auto const &cd : this->operator_matrices[c_dag_idx].entries) {
        for (auto const &ca : this->operator_matrices[c_idx].entries) {
          if (ca.row == cd.col) { n_sigma[cd.row][ca.col] = n_sigma[cd.row][ca.col] + cd.value * ca.value; }
        }
      }

      // acc ← n_σ · acc.
      DenseMatrix prod{};
      for (int row = 0; row < N_STATES; ++row) {
        for (int k = 0; k < N_STATES; ++k) {
          if (is_zero(n_sigma[row][k])) continue;
          for (int col = 0; col < N_STATES; ++col) {
            if (is_zero(acc[k][col])) continue;
            prod[row][col] = prod[row][col] + n_sigma[row][k] * acc[k][col];
          }
        }
      }
      acc = prod;
    }
    return acc;
  }

  template <int N_sites, typename T>
  T HubbardSolver<N_sites, T>::G0n_with_density_matrix(Args<N_sites, T> const &args, DenseMatrix const &N_matrix) const {
    int n = args.order;

    // exp(tau * deltaE) factors for the hybridization operators (see G0n).
    ExpTable exp_tau_E, inv_exp_tau_E;
    this->build_tau_exp_tables(args, exp_tau_E, inv_exp_tau_E);

    std::array<T, N_STATES> buf_a, buf_b;

    T result = T(0.0);
    for (int start_state = 0; start_state < N_STATES; ++start_state) {
      T *amplitudes = buf_a.data();
      T *next       = buf_b.data();

      std::fill(amplitudes, amplitudes + N_STATES, T(0.0));
      amplitudes[start_state] = this->exp_beta_E[start_state];

      // Apply N at τ=0 (right-most slot of the time-ordered trace, no
      // evolution factor): amplitudes ← N · amplitudes.
      std::fill(next, next + N_STATES, T(0.0));
      for (int col = 0; col < N_STATES; ++col) {
        if (is_zero(amplitudes[col])) continue;
        for (int row = 0; row < N_STATES; ++row) {
          if (is_zero(N_matrix[row][col])) continue;
          next[row] = next[row] + N_matrix[row][col] * amplitudes[col];
        }
      }
      std::swap(amplitudes, next);

      // Then the hybridization operators, exactly as G0n does.
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

  template <int N_sites, typename T>
  T HubbardSolver<N_sites, T>::G0n_with_densities(Args<N_sites, T> const &args, std::vector<int> const &density_orbitals) const {
    // Empty decoration must recover plain G0n exactly (EmptyOrbitalsMatchesG0n).
    if (density_orbitals.empty()) return this->G0n(args);
    if (!args.operator_sequence_is_valid()) return T(0.0);

    if constexpr (N_sites != 1) {
      // Cluster path: H_dimer mixes orbitals, so n_σ is NOT diagonal in the
      // energy eigenbasis and the atomic start-state weighting is wrong.
      // Insert N = Π_k n_{σ_k} as a matrix at the τ=0 trace slot.
      DenseMatrix N_matrix = this->build_density_matrix(density_orbitals);
      return this->G0n_with_density_matrix(args, N_matrix);
    } else {
      // ---- Atomic diagonal fast path ----
      // For atomic, n_σ is diagonal in the eigenbasis (states 0..3 are pure
      // occupation bitstrings). Precompute the per-start-state density factor.
      int n = args.order;

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
    using std::exp;
    if constexpr (N_sites != 1) {
      static_cast<void>(args);
      static_cast<void>(density_orbitals);
      return T(0.0);
    } else {
      // Closed form, τ-independent. At U=∞ the projected Hilbert space is
      // {|0⟩, |↓⟩, |↑⟩}; the density n_σ(0) projects the trace onto start
      // state |σ⟩. Because the operator chain is particle-number balanced at
      // each vertex, the e^{μ(Στ−Στ′)} factors cancel across the diagram —
      // so the per-leaf trace is either 0 (operator string can't close on
      // |σ⟩) or e^{βμ}/Z_∞ × permutation_sign (the |σ⟩ Boltzmann weight).
      //
      // Cases on density_orbitals.size():
      //   0   — no density; defer to the vacuum closed form.
      //   1   — single n_σ(0) at the marked vertex.
      //   2   — coincident pair. Same spin: n²=n, collapses to one density;
      //         opposite spin: double-occupancy projector, identically zero
      //         at U=∞ (no |↑↓⟩ state).
      if (density_orbitals.empty()) return this->G0n_infinite_U(args);
      if (density_orbitals.size() == 2 && density_orbitals[0] != density_orbitals[1]) return T(0.0);
      int sigma = density_orbitals[0];

      if (!args.operator_sequence_is_valid_infinite_U()) return T(0.0);

      // Order-0 with no hopping operators: just the |σ⟩ Boltzmann weight.
      if (args.order == 0) {
        return exp(this->params.beta * this->params.mu) * (T(1.0) / this->Z_infinite_U) * args.permutation_sign;
      }

      // ops are stored with ops[0]=largest τ, ops[n-1]=smallest τ (matches
      // operator_sequence_is_valid_infinite_U's iteration convention). The
      // trace closes on |σ⟩ iff the largest-τ operator is c†_σ (creation of
      // σ): the time-ordered product is applied to |σ⟩ right-to-left
      // (smallest τ first); the LAST op applied is ops[0], and it must end
      // the chain at |σ⟩. c†_σ acting on |0⟩ → |σ⟩ closes the trace; the
      // mirror "c_σ at largest τ" would end at |0⟩, not |σ⟩.
      auto last_op = args.ops[0];
      if (last_op.get_action() != 1 || last_op.get_orbital_index() != (uint8_t)sigma) {
        return T(0.0);
      }

      return exp(this->params.beta * this->params.mu) * (T(1.0) / this->Z_infinite_U) * args.permutation_sign;
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
    uint8_t c_op_id                    = static_cast<uint8_t>(orbital);
    uint8_t c_dag_op_id                = static_cast<uint8_t>(ACTION_BIT_LOCAL | orbital);

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
