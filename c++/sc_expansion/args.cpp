#include "args.hpp"
#include "dual.hpp"
#include <numeric>
#include <algorithm>
#include <stdexcept>

namespace sc_expansion {

  template <int N_sites, typename T>
  Args<N_sites, T>::Args(std::vector<double> taus_, std::vector<FermionOperator<N_sites, T>> ops_)
     : taus(std::move(taus_)), ops(std::move(ops_)) {
    if (taus.size() != ops.size()) { throw std::runtime_error("Error in Args constructor: taus and ops must have the same size"); }

    this->order            = taus.size();
    this->permutation_sign = 1.0; // overwritten by sort_args
    this->sort_args();
  }

  template <int N_sites, typename T> void Args<N_sites, T>::sort_args() {

    std::vector<int> argsort_local(this->order);
    std::iota(argsort_local.begin(), argsort_local.end(), 0);
    std::stable_sort(argsort_local.begin(), argsort_local.end(), [&](int i, int j) { return this->taus[i] > this->taus[j]; });

    std::vector<double> sorted_taus;
    sorted_taus.reserve(this->order);
    std::vector<FermionOperator<N_sites, T>> sorted_ops;
    sorted_ops.reserve(this->order);
    for (int i : argsort_local) {
      sorted_taus.push_back(this->taus[i]);
      sorted_ops.push_back(this->ops[i]);
    }
    this->permutation_sign = (double)compute_permutation_sign(argsort_local);
    this->argsort          = std::move(argsort_local);
    this->taus             = std::move(sorted_taus);
    this->ops              = std::move(sorted_ops);
  }

  template <int N_sites, typename T> bool Args<N_sites, T>::operator_sequence_is_valid() const {
    // Both checks below — Pauli adjacency and even-touch trace closure — assume H_0
    // is diagonal in the (site×spin) orbital basis, so c†_orbital commutes with
    // exp(τ H_0). That holds only for atomic (N_sites=1). For clusters with
    // intra-cluster hopping (N_sites>=2), the time-ordered product
    // c†(τ_2) exp((τ_2−τ_1)H_0) c†(τ_1) is generically non-zero because H_0 mixes
    // orbitals between the two operator applications. In that regime, defer to
    // G0n's amplitude propagation, which is exact.
    if constexpr (N_sites != 1) { return true; }

    int last_action[N_ORBITALS];
    std::fill(std::begin(last_action), std::end(last_action), -1);
    int ops_count[N_ORBITALS] = {0};

    for (auto const &f_op : ops) {
      int orbital_index = f_op.get_orbital_index();
      int action        = f_op.get_action();
      if (last_action[orbital_index] == action) { return false; }
      last_action[orbital_index] = action;
      ops_count[orbital_index]++;
    }

    for (int i = 0; i < N_ORBITALS; ++i) {
      if (ops_count[i] % 2 != 0) { return false; }
    }
    return true;
  }

  template <int N_sites, typename T> bool Args<N_sites, T>::operator_sequence_is_valid_infinite_U() const {
    if constexpr (N_sites != 1) {
      // Infinite-U projection is single-site only. The dimer MCMC uses a uniform
      // reference weight and never invokes this path; signal "no" defensively.
      return false;
    } else {
      if (!this->operator_sequence_is_valid()) return false;

      // Allowed start states for single-site infinite-U: |0>, |down>, |up> (state 0,1,2).
      // Double occupancy (state 11 = 3) is forbidden.
      // Iterate ops from smallest tau (ops[n-1]) to largest tau (ops[0]) so the
      // check matches the physical matrix product <s| O_1 ... O_n |s>.
      for (int start_state : {0, 1, 2}) {
        int current_state = start_state;
        bool path_valid   = true;

        for (int i = (int)this->ops.size() - 1; i >= 0; --i) {
          auto const &f_op = this->ops[i];
          int orbital      = f_op.get_orbital_index();
          int action       = f_op.get_action();

          if (action == 1) { // Create
            if (current_state & (1 << orbital)) {
              path_valid = false;
              break;
            }
            current_state |= (1 << orbital);
          } else { // Destroy
            if (!(current_state & (1 << orbital))) {
              path_valid = false;
              break;
            }
            current_state &= ~(1 << orbital);
          }

          if (current_state == 3) { // Double occupancy reached
            path_valid = false;
            break;
          }
        }

        if (path_valid && current_state == start_state) return true;
      }

      return false;
    }
  }

  template <int N_sites, typename T>
  std::pair<Args<N_sites, T>, Args<N_sites, T>> Args<N_sites, T>::split_from_raw(std::vector<double> const &taus,
                                                                                 std::vector<uint8_t> const &op_ids) {
    if (taus.size() != op_ids.size()) { throw std::runtime_error("Error in Args::split_from_raw: taus and op_ids must have the same size"); }

    std::vector<double> t_u, t_p;
    std::vector<FermionOperator<N_sites, T>> o_u, o_p;

    for (size_t i = 0; i < taus.size(); ++i) {
      FermionOperator<N_sites, T> op(op_ids[i]);
      if (op.get_action() == 0) { // annihilator
        t_u.push_back(taus[i]);
        o_u.push_back(op);
      } else { // creator
        t_p.push_back(taus[i]);
        o_p.push_back(op);
      }
    }
    return {Args<N_sites, T>(std::move(t_u), std::move(o_u)), Args<N_sites, T>(std::move(t_p), std::move(o_p))};
  }

  template struct Args<1, double>;
  template struct Args<1, Dual>;
  template struct Args<2, double>;
  template struct Args<2, Dual>;

} // namespace sc_expansion
