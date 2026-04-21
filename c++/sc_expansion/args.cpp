#include "args.hpp"
#include "dual.hpp"
#include <numeric>
#include <algorithm>
#include <stdexcept>

namespace sc_expansion {

  template <typename T>
  Args<T>::Args(std::vector<double> taus_, std::vector<FermionOperator<T>> ops_) : taus(std::move(taus_)), ops(std::move(ops_)) {
    if (taus.size() != ops.size()) { throw std::runtime_error("Error in Args constructor: taus and ops must have the same size"); }

    this->order            = taus.size();
    this->permutation_sign = 1.0; // Placeholder, will be set after sorting
    this->sort_args();
  }

  template <typename T> void Args<T>::sort_args() {

    std::vector<int> argsort_local(this->order);
    std::iota(argsort_local.begin(), argsort_local.end(), 0);
    std::stable_sort(argsort_local.begin(), argsort_local.end(), [&](int i, int j) { return this->taus[i] > this->taus[j]; });

    // Rearrange ops and taus according to the sorted order
    std::vector<double> sorted_taus;
    sorted_taus.reserve(this->order);
    std::vector<FermionOperator<T>> sorted_ops;
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

  template <typename T> bool Args<T>::operator_sequence_is_valid() const {
    int last_action[Args<T>::N_ORBITALS];
    std::fill(std::begin(last_action), std::end(last_action), -1);
    int ops_count[Args<T>::N_ORBITALS] = {0};

    //1. ensure Pauli exclusion is respected
    for (auto const &f_op : ops) {
      int orbital_index = f_op.get_orbital_index();
      int action        = f_op.get_action();
      if (last_action[orbital_index] == action) { return false; } // Two consecutive creates or destroys on the same orbital
      last_action[orbital_index] = action;
      ops_count[orbital_index]++;
    }

    //ensure that trace loop closes: each orbital must be acted on an even number of times (created and destroyed the same number of times)
    for (int i = 0; i < Args<T>::N_ORBITALS; ++i) {
      if (ops_count[i] % 2 != 0) { return false; }
    }
    return true;
  }

  template <typename T> bool Args<T>::operator_sequence_is_valid_infinite_U() const {
    if (!this->operator_sequence_is_valid()) return false;

    // Start states: |0>, |down>, |up> (bit representation: 00, 01, 10)
    // Double occupancy (state 11 = 3) is forbidden.
    // Iterate ops from smallest tau (ops[n-1]) to largest tau (ops[0]) so the
    // check matches the physical matrix product <s| O_1 ... O_n |s> (rightmost
    // factor applied first). The forward order happens to agree with the
    // reverse order when every tau is distinct, but breaks for self-loops where
    // a (c, c^dag) pair sits at the same tau: forward admits mixed-spin start
    // states that are physically blocked, giving a spurious e^{beta*mu}/Z_inf.
    for (int start_state : {0, 1, 2}) {
      int current_state = start_state;
      bool path_valid  = true;

      for (int i = (int)this->ops.size() - 1; i >= 0; --i) {
        auto const &f_op = this->ops[i];
        int orbital = f_op.get_orbital_index();
        int action  = f_op.get_action();

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

  template <typename T>
  std::pair<Args<T>, Args<T>> Args<T>::split_from_raw(std::vector<double> const &taus, std::vector<uint8_t> const &op_ids) {
    if (taus.size() != op_ids.size()) { throw std::runtime_error("Error in Args::split_from_raw: taus and op_ids must have the same size"); }

    std::vector<double> t_u, t_p;
    std::vector<FermionOperator<T>> o_u, o_p;

    for (size_t i = 0; i < taus.size(); ++i) {
      FermionOperator<T> op(op_ids[i]);
      if (op.get_action() == 0) { // Destruction
        t_u.push_back(taus[i]);
        o_u.push_back(op);
      } else { // Creation
        t_p.push_back(taus[i]);
        o_p.push_back(op);
      }
    }
    // Note: The Args constructor automatically performs the canonical time-sorting and computes the permutation_sign.
    return {Args<T>(std::move(t_u), std::move(o_u)), Args<T>(std::move(t_p), std::move(o_p))};
  }

  template struct Args<double>;
  template struct Args<Dual>;

} // namespace sc_expansion
