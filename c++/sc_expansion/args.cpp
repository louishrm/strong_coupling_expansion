#include "args.hpp"
#include "dual.hpp"
#include <numeric>
#include <algorithm>

namespace sc_expansion {

  template <int N_sites, typename T>
  Args<N_sites, T>::Args(std::vector<double> taus_, std::vector<FermionOperator<N_sites, T>> ops_) : taus(std::move(taus_)), ops(std::move(ops_)) {
    if (taus.size() != ops.size()) { throw std::runtime_error("Error in Args constructor: taus and ops must have the same size"); }

    this->order            = taus.size();
    this->permutation_sign = 1.0; // Placeholder, will be set after sorting
    this->sort_args();
  }

  template <int N_sites, typename T> void Args<N_sites, T>::sort_args() {

    std::vector<int> argsort(this->order);
    std::iota(argsort.begin(), argsort.end(), 0);
    std::stable_sort(argsort.begin(), argsort.end(), [&](int i, int j) { return this->taus[i] > this->taus[j]; });

    // Rearrange ops and taus according to the sorted order
    std::vector<double> sorted_taus;
    sorted_taus.reserve(this->order);
    std::vector<FermionOperator<N_sites, T>> sorted_ops;
    sorted_ops.reserve(this->order);
    for (int i : argsort) {
      sorted_taus.push_back(this->taus[i]);
      sorted_ops.push_back(this->ops[i]);
    }
    this->permutation_sign = (double)compute_permutation_sign(argsort);
    this->taus             = std::move(sorted_taus);
    this->ops              = std::move(sorted_ops);
  }

  template <int N_sites, typename T> bool Args<N_sites, T>::operator_sequence_is_valid() const {
    int last_action[Args<N_sites, T>::N_ORBITALS];
    std::fill(std::begin(last_action), std::end(last_action), -1);
    int ops_count[Args<N_sites, T>::N_ORBITALS] = {0};

    //1. ensure Pauli exclusion is respected
    for (auto const &f_op : ops) {
      int orbital_index = f_op.get_orbital_index();
      int action        = f_op.get_action();
      if (last_action[orbital_index] == action) { return false; } // Two consecutive creates or destroys on the same orbital
      last_action[orbital_index] = action;
      ops_count[orbital_index]++;
    }

    //ensure that trace loop closes: each orbital must be acted on an even number of times (created and destroyed the same number of times)
    for (int i = 0; i < Args<N_sites, T>::N_ORBITALS; ++i) {
      if (ops_count[i] % 2 != 0) { return false; }
    }
    return true;
  }

  template <int N_sites, typename T> bool Args<N_sites, T>::operator_sequence_is_valid_infinite_U() const {
    if constexpr (N_sites != 1) return true; // Only defined for Hubbard Atom currently

    if (!this->operator_sequence_is_valid()) return false;

    // Start states: |0>, |down>, |up> (bit representation: 00, 01, 10)
    // Double occupancy (state 11 = 3) is forbidden
    for (int start_state : {0, 1, 2}) {
      int current_state = start_state;
      bool path_valid  = true;

      for (auto const &f_op : ops) {
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

  template struct Args<1, double>;
  template struct Args<1, Dual>;
  template struct Args<2, double>;
  template struct Args<2, Dual>;

} // namespace sc_expansion
