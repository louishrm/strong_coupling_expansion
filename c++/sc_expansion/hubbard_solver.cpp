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

  template <int N_sites, typename T> HubbardSolver<N_sites, T>::HubbardSolver(Parameters<T> const &params_) : params(params_) {
    for (int i = 0; i < N_OPS; ++i) { this->operators[i] = FermionOperator<N_sites, T>(static_cast<uint8_t>(i)); }

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
      // Even-parity basis: |e+> = (|0101>+|1010>)/√2, |e-> = (|1001>+|0110>)/√2
      // 2×2 block: H_even = [[U-2mu, -2t], [-2t, -2mu]]
      // Eigenvector for eigenvalue E: (U-2mu-E)*a = 2t*b → b/a = (U-2mu-E)/(2t)
      // Fock-basis coefficients: norm = a/√2, component = b/√2
      T Ep = Eplus(t, U, mu);
      T Em = Eminus(t, U, mu);
      using std::sqrt;

      // At t=0 the even-parity block is diagonal: eigenstates are pure |e+> and |e->.
      // The ratio = (U-2mu-E)/(2t) diverges, so we handle this limit explicitly.
      T norm_plus, norm_minus, component_plus, component_minus;
      auto t_val = [](auto const &x) -> double {
        if constexpr (std::is_same_v<std::decay_t<decltype(x)>, Dual>) {
          return x.value;
        } else {
          return x;
        }
      };
      if (std::abs(t_val(t)) < 1e-15) {
        // t=0 limit: E+ = U-2mu (eigenvector |e+>), E- = -2mu (eigenvector |e->)
        norm_plus       = T(SQRT2_INV);
        component_plus  = T(0.0);
        norm_minus      = T(0.0);
        component_minus = T(SQRT2_INV);
      } else {
        T ratio_plus    = (U - T(2.0) * mu - Ep) / (T(2.0) * t);
        T ratio_minus   = (U - T(2.0) * mu - Em) / (T(2.0) * t);
        norm_plus       = T(SQRT2_INV) / sqrt(T(1.0) + ratio_plus * ratio_plus);
        norm_minus      = T(SQRT2_INV) / sqrt(T(1.0) + ratio_minus * ratio_minus);
        component_plus  = ratio_plus * norm_plus;
        component_minus = ratio_minus * norm_minus;
      }

      this->all_eigenstates[6] = Eigenstate<T>{{{5, T(norm_plus)}, {10, T(norm_plus)}, {9, T(component_plus)}, {6, T(component_plus)}}, Ep};
      this->all_eigenstates[7] = Eigenstate<T>{{{5, T(norm_minus)}, {10, T(norm_minus)}, {9, T(component_minus)}, {6, T(component_minus)}}, Em};

      // N=2, Sz=0, parity = odd
      this->all_eigenstates[8] = Eigenstate<T>{{{5, T(SQRT2_INV)}, {10, -T(SQRT2_INV)}}, U - 2.0 * mu};
      this->all_eigenstates[9] = Eigenstate<T>{{{9, T(SQRT2_INV)}, {6, -T(SQRT2_INV)}}, -2.0 * mu};

      // N=2, Sz=+1 |up,up>
      this->all_eigenstates[10] = Eigenstate<T>{{{12, T(1.0)}}, -2.0 * mu};

      // N=3, Sz = -1/2, |down up, down> ± |down, down up>
      // H = [[U-3mu, -t], [-t, U-3mu]] → symmetric (1,1)/√2 has eigenvalue U-t-3mu
      this->all_eigenstates[11] = Eigenstate<T>{{{7, T(SQRT2_INV)}, {11, T(SQRT2_INV)}}, U - t - 3.0 * mu};
      this->all_eigenstates[12] = Eigenstate<T>{{{7, T(SQRT2_INV)}, {11, -T(SQRT2_INV)}}, U + t - 3.0 * mu};

      // N=3, Sz = +1/2, |down up, up> ± |up, down up>
      this->all_eigenstates[13] = Eigenstate<T>{{{13, T(SQRT2_INV)}, {14, T(SQRT2_INV)}}, U - t - 3.0 * mu};
      this->all_eigenstates[14] = Eigenstate<T>{{{13, T(SQRT2_INV)}, {14, -T(SQRT2_INV)}}, U + t - 3.0 * mu};

      // N=4, Sz=0 |down up, down up>
      this->all_eigenstates[15] = Eigenstate<T>{{{15, T(1.0)}}, 2.0 * U - 4.0 * mu};
    }
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
    T result = T(0.0);

    if (!args.operator_sequence_is_valid()) { return result; }

    int n = args.order;

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
