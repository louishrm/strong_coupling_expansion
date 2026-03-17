#include "./hubbard_atom.hpp"
#include "./combinatorics.hpp"
#include "./dual.hpp"

namespace sc_expansion {

  Args::Args(std::vector<double> taus_, std::vector<int> spins_) : taus(std::move(taus_)), spins(std::move(spins_)) {
    if (taus.size() != spins.size()) { throw std::runtime_error("Error in Args constructor: taus and spins must have the same size"); }

    this->order            = taus.size();
    this->permutation_sign = 1.0; // Placeholder, will be set after sorting
    this->sort_args();
  }

  void Args::sort_args() {
    int n = this->taus.size();
    std::vector<int> argsort(n);
    std::iota(argsort.begin(), argsort.end(), 0);

    std::vector<int> ops_local(n); // 0 = cup, 1 = cdn, 2 = cdag_up, 3 = cdag_dn

    for (size_t i = 0; i < n; ++i) {

      if (i % 2 == 0) {
        ops_local[i] = 2 + this->spins[i];
      } // Even index: creation operator (cdag_up or cdag_dn)
      else {
        ops_local[i] = this->spins[i];
      } // Odd index: annihilation operator (cup or cdn)
    }

    std::stable_sort(argsort.begin(), argsort.end(), [&](int i, int j) { return this->taus[i] > this->taus[j]; });

    this->permutation_sign = compute_permutation_sign(argsort);

    std::vector<double> sorted_times(n);
    std::vector<int> sorted_ops(n);

    for (size_t i = 0; i < n; ++i) {
      sorted_times[i] = this->taus[argsort[i]];
      sorted_ops[i]   = ops_local[argsort[i]];
    }

    this->taus = std::move(sorted_times);
    this->ops  = std::move(sorted_ops);
  }

  bool Args::verify_consecutive_terms_infinite_U() const {
    int n = this->order;

    // 2. Internal Sequence Checks
    for (int i = 0; i + 1 < n; ++i) {

      int op_late  = this->ops[i];     // Applied second
      int op_early = this->ops[i + 1]; // Applied first

      bool early_is_create = (op_early >= 2);
      bool late_is_create  = (op_late >= 2);

      // Rule A: Strict Alternation
      if (early_is_create == late_is_create) return false;

      // Rule B: Spin Conservation while occupied
      if (early_is_create) {
        if (op_late != op_early - 2) return false;
      }
    }

    int op_beta = this->ops[0];     // The very last operator applied (latest time)
    int op_zero = this->ops[n - 1]; // The very first operator applied (earliest time)

    bool beta_is_create = (op_beta >= 2);

    if (beta_is_create) {
      if (op_zero != op_beta - 2) { return false; }
    }

    return true;
  }

  template <typename T>
  const std::array<Transition, 16> HubbardAtom<T>::lookup_table = []() {
    std::array<Transition, 16> table;

    // Iterate through all 4 states (00, 01, 10, 11)
    for (int state = 0; state < 4; ++state) {

      // Iterate through all 4 operators (0=c_up, 1=c_down, 2=cdag_up, 3=cdag_down)
      for (int op = 0; op < 4; ++op) {

        // The 1D array index: (State * 4) + Operator
        int index = (state << 2) | op;

        int next_state = 0;
        double mel     = 0.0;

        // Decode the operator type
        bool is_create = (op >= 2);     // Ops 2 and 3 are creation
        bool is_down   = (op % 2 != 0); // Ops 1 and 3 are down spin

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

  template <typename T>
  HubbardAtom<T>::HubbardAtom(Parameters<T> const &params_) : params(params_) {
    this->Z            = T(1.0) + T(2.0) * exp(params.beta * params.mu) + exp(params.beta * (T(2.0) * params.mu - params.U)); //Partition function at finite U
    this->Z_infinite_U = T(1.0) + T(2.0) * exp(params.beta * params.mu);                                 //Partition function at infinite U
    this->E            = {T(0.0), -params.mu, -params.mu, params.U - T(2.0) * params.mu};
  }

  template <typename T>
  T HubbardAtom<T>::G0(std::vector<double> const &taus, std::vector<int> const &spins) const {
    Args args(taus, spins);
    T G0_value = T(0.0);

    //get the only two states that the list op in the sorted list can act on
    int last_op = args.ops.back();

    for (int initial_state : this->valid_start_states[last_op]) {

      int current_state        = initial_state;
      T current_trace_val = T(1.0);

      for (int i = args.order - 1; i >= 0; --i) {

        int table_index = current_state * 4 + args.ops[i];
        Transition t    = this->lookup_table[table_index];

        if (t.matrix_element == 0.0) {
          current_trace_val = T(0.0);
          break;
        }
        int next_state     = t.connected_state;
        T energy_diff = this->E[next_state] - this->E[current_state];
        current_trace_val = current_trace_val * T(t.matrix_element) * exp(T(args.taus[i]) * energy_diff);
        current_state = next_state;
      }

      // Add to value for Dual properly
      if (current_state == initial_state) {
        current_trace_val = current_trace_val * exp(-this->params.beta * this->E[current_state]);
        G0_value = G0_value + current_trace_val;
      }
    }

    return T(1.0) / this->Z * T(args.permutation_sign) * G0_value;
  }

  template <typename T>
  T HubbardAtom<T>::G0_infinite_U(std::vector<double> const &taus, std::vector<int> const &spins) const {
    Args args(taus, spins);
    T G0_value = T(0.0);

    if (!args.verify_consecutive_terms_infinite_U()) { return T(0.0); }

    // If the sequence is valid, we can directly compute the contribution from the first operator.
    // The contribution from the rest of the sequence is guaranteed to be 1 due to the strict alternation and spin conservation rules.

    int first_op = args.ops[0];
    if (first_op >= 2) {
      // First operator is a creation operator
      G0_value = exp(this->params.beta * this->params.mu) / this->Z_infinite_U;
    } else {
      // First operator is a destruction operator
      G0_value = T(1.0) / this->Z_infinite_U;
    }

    return T(args.permutation_sign) * G0_value;
  }

  template class HubbardAtom<double>;
  template class HubbardAtom<Dual>;

} // namespace sc_expansion

// namespace sc_expansion {

//   using namespace triqs::operators;

//   HubbardAtom::HubbardAtom(double U, double beta, double mu) : U(U), beta(beta), mu(mu) {
//     this->Z_infinite_U = 1 + 2 * std::exp(beta * mu); //Partition function at infinite U
//     auto fops          = make_fops();
//     this->H            = make_H0(U, mu);
//     this->ad           = triqs::atom_diag::atom_diag<false>(this->H, fops, {});
//     this->rho0         = triqs::atom_diag::atomic_density_matrix(this->ad, beta)[0];
//     this->c_dn         = this->ad.c_matrix(0, 0);
//     this->c_up         = this->ad.c_matrix(1, 0);
//     this->cdag_dn      = this->ad.cdag_matrix(0, 0);
//     this->cdag_up      = this->ad.cdag_matrix(1, 0);
//   }

//   triqs::hilbert_space::fundamental_operator_set HubbardAtom::make_fops() {
//     triqs::hilbert_space::fundamental_operator_set fops = {};
//     for (int o : triqs::arrays::range(1)) fops.insert("dn", o);
//     for (int o : triqs::arrays::range(1)) fops.insert("up", o);
//     return fops;
//   }

//   triqs::operators::many_body_operator_generic<double> HubbardAtom::make_H0(double U, double mu) {
//     triqs::operators::many_body_operator_generic<double> h = U * n("up", 0) * n("dn", 0);
//     h += -mu * (n("up", 0) + n("dn", 0));
//     return h;
//   }

//   int HubbardAtom::calculate_permutation_sign(std::vector<int> const &perm) {

//     int inversions = 0;
//     for (size_t i = 0; i < perm.size(); ++i) {
//       for (size_t j = i + 1; j < perm.size(); ++j) {
//         if (perm[i] > perm[j]) { ++inversions; }
//       }
//     }
//     return (inversions % 2 == 0) ? 1 : -1;
//   }

//   std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, int>
//   HubbardAtom::sort_operators(const std::vector<double> &times, const std::vector<int> &spins, const std::vector<int> &flags) {

//     int n = times.size();
//     std::vector<int> argsort(n);
//     std::iota(argsort.begin(), argsort.end(), 0);

//     std::stable_sort(argsort.begin(), argsort.end(), [&times](int i, int j) { return times[i] > times[j]; });

//     int sign = calculate_permutation_sign(argsort);

//     std::vector<double> sorted_times(n);
//     std::vector<int> sorted_spins(n);
//     std::vector<int> sorted_flags(n);

//     for (int i = 0; i < n; ++i) {
//       sorted_times[i] = times[argsort[i]];
//       sorted_spins[i] = spins[argsort[i]];
//       sorted_flags[i] = flags[argsort[i]];
//     }

//     return {sorted_times, sorted_spins, sorted_flags, sign};
//   }

//   bool HubbardAtom::verify_consecutive_terms(const std::vector<int> &sorted_spins, const std::vector<int> &sorted_flags) {
//     // Vectors to store flags for each spin sector
//     std::vector<int> flags_up;
//     std::vector<int> flags_down;

//     // 1. Separate the stream by spin
//     for (size_t i = 0; i < sorted_spins.size(); ++i) {
//       if (sorted_spins[i] == 1) { // Assuming 1 is UP
//         flags_up.push_back(sorted_flags[i]);
//       } else { // Assuming 0 or -1 is DOWN
//         flags_down.push_back(sorted_flags[i]);
//       }
//     }

//     // 2. Helper lambda to check for strict alternation in a single sector
//     auto is_invalid_sector = [](const std::vector<int> &flags) {
//       for (size_t i = 0; i + 1 < flags.size(); ++i) {
//         // If two consecutive flags are identical (1,1 or 0,0), it's invalid
//         if (flags[i] == flags[i + 1]) { return true; }
//       }
//       return false;
//     };

//     // 3. Verify both sectors
//     if (is_invalid_sector(flags_up)) return false;
//     if (is_invalid_sector(flags_down)) return false;

//     return true;
//   }

//   nda::matrix<double> HubbardAtom::make_interaction_picture_destroy_op(double tau, int state_index) const {
//     double Z01 = triqs::atom_diag::partition_function(this->ad, tau);
//     double Z02 = triqs::atom_diag::partition_function(this->ad, -tau);

//     auto const &cmat = (state_index == 0 ? this->c_dn : this->c_up);
//     auto time_evol_1 = triqs::atom_diag::atomic_density_matrix(this->ad, -tau)[0];
//     auto time_evol_2 = triqs::atom_diag::atomic_density_matrix(this->ad, tau)[0];
//     auto ctau        = Z01 * Z02 * time_evol_1 * cmat * time_evol_2;

//     return ctau;
//   }

//   nda::matrix<double> HubbardAtom::make_interaction_picture_create_op(double tau, int state_index) const {
//     double Z01 = triqs::atom_diag::partition_function(this->ad, tau);
//     double Z02 = triqs::atom_diag::partition_function(this->ad, -tau);

//     auto const &cmat = (state_index == 0 ? this->cdag_dn : this->cdag_up);
//     auto time_evol_1 = triqs::atom_diag::atomic_density_matrix(this->ad, -tau)[0];
//     auto time_evol_2 = triqs::atom_diag::atomic_density_matrix(this->ad, tau)[0];
//     auto cdagtau     = Z01 * Z02 * time_evol_1 * cmat * time_evol_2;

//     return cdagtau;
//   }

//   double HubbardAtom::G0(cumul_args const &unprimed_args, cumul_args const &primed_args) const {
//     int order = unprimed_args.size();
//     if (order != primed_args.size()) { throw std::runtime_error("Error in G0: unprimed_args and primed_args must have the same size"); }

//     int n_ops = 2 * order;
//     std::vector<double> times;
//     std::vector<int> spins;
//     std::vector<int> flags;
//     times.reserve(n_ops);
//     spins.reserve(n_ops);
//     flags.reserve(n_ops);

//     for (int i = 0; i < order; ++i) {
//       // Push the Creation operator first
//       times.push_back(primed_args[i].first);
//       spins.push_back(primed_args[i].second);
//       flags.push_back(1);

//       // Push the Destruction operator second
//       times.push_back(unprimed_args[i].first);
//       spins.push_back(unprimed_args[i].second);
//       flags.push_back(0);
//     }

//     auto [sorted_times, sorted_spins, sorted_flags, sign] = sort_operators(times, spins, flags);
//     if (!verify_consecutive_terms(sorted_spins, sorted_flags)) { return 0.0; }

//     nda::matrix<double> op = this->rho0;

//     for (int i = 0; i < n_ops; ++i) {
//       if (sorted_flags[i] == 1) {
//         op *= this->make_interaction_picture_create_op(sorted_times[i], sorted_spins[i]);
//       } else {
//         op *= this->make_interaction_picture_destroy_op(sorted_times[i], sorted_spins[i]);
//       }
//     }
//     double G0_value = sign * trace(op);

//     return G0_value;
//   }

//   bool HubbardAtom::verify_consecutive_terms_infinite_U(const std::vector<int> &sorted_spins, const std::vector<int> &sorted_flags, int n_ops) {
//     int total_ops = 2 * n_ops;
//     for (int i = 0; i < total_ops - 1; ++i) {
//       if (sorted_flags[i] == sorted_flags[i + 1]) { return false; } //Two consecutive creates or destroys (violates alternation/no double occupancy)

//       if (sorted_flags[i] == 0 && sorted_flags[i + 1] == 1) {
//         if (sorted_spins[i] != sorted_spins[i + 1]) return false;
//       } //Destroying a particle with different spin than the one created at the previous (later) time
//     }

//     // Check boundary condition for trace: if the trace ends with an occupied state (first op is create),
//     // then it must have started with the same occupied state (last op is destroy).
//     if (sorted_flags[0] == 1 && sorted_spins[0] != sorted_spins[total_ops - 1]) return false;

//     return true;
//   }

//   double HubbardAtom::G0_infinite_U(cumul_args const &unprimed_args, cumul_args const &primed_args) const {

//     int order = unprimed_args.size();
//     if (order != primed_args.size()) { throw std::runtime_error("Error in G0: unprimed_args and primed_args must have the same size"); }

//     int n_ops = 2 * order;
//     std::vector<double> times;
//     std::vector<int> spins;
//     std::vector<int> flags;
//     times.reserve(n_ops);
//     spins.reserve(n_ops);
//     flags.reserve(n_ops);

//     for (auto [t, s] : unprimed_args) {
//       times.push_back(t);
//       spins.push_back(s);
//       flags.push_back(0); // destroy
//     }

//     for (auto [t, s] : primed_args) {
//       times.push_back(t);
//       spins.push_back(s);
//       flags.push_back(1); // create
//     }

//     auto [sorted_times, sorted_spins, sorted_flags, sign] = sort_operators(times, spins, flags);

//     if (!verify_consecutive_terms_infinite_U(sorted_spins, sorted_flags, order)) { return 0.0; }

//     if (sorted_flags[0] == 1) {
//       return sign * std::exp(this->beta * this->mu) / this->Z_infinite_U;
//     } //First operator is create
//     else {
//       return sign * 1.0 / this->Z_infinite_U;
//     } //First operator is destroy
//   }
// } // namespace sc_expansion