#include "operator_sequence.hpp"
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <map>

// Declaration of compute_permutation_sign (assuming it's in combinatorics.hpp but since we just need it here)
namespace sc_expansion {
    double compute_permutation_sign(const std::vector<int>& p);
}

namespace sc_expansion {

  OperatorSequence::OperatorSequence(std::vector<double> taus_, std::vector<int> ops_)
     : taus(std::move(taus_)), ops(std::move(ops_)) {
    if (this->taus.size() != this->ops.size()) {
      throw std::runtime_error("Error in OperatorSequence constructor: taus and ops must have the same size");
    }

    this->order            = this->taus.size();
    this->permutation_sign = 1.0;
    this->sort_sequence();
  }

  void OperatorSequence::sort_sequence() {
    int n = this->order;
    std::vector<int> argsort(n);
    std::iota(argsort.begin(), argsort.end(), 0);

    // stable_sort by imaginary time descending
    std::stable_sort(argsort.begin(), argsort.end(), [&](int i, int j) {
      return this->taus[i] > this->taus[j];
    });

    this->permutation_sign = compute_permutation_sign(argsort);

    std::vector<double> sorted_times(n);
    std::vector<int> sorted_ops(n);

    for (int i = 0; i < n; ++i) {
      sorted_times[i] = this->taus[argsort[i]];
      sorted_ops[i]   = this->ops[argsort[i]];
    }

    this->taus = std::move(sorted_times);
    this->ops  = std::move(sorted_ops);
  }

  bool OperatorSequence::verify_flavor_alternation(int action_bit_index) const {
    // Dictionary to track the expected next action for each flavor
    std::map<int, int> expected_action_for_flavor;

    for (int op : this->ops) {
      int action = (op >> action_bit_index) & 1;
      // The flavor is all bits except the action bit
      int flavor = op & ~(1 << action_bit_index);

      if (expected_action_for_flavor.find(flavor) != expected_action_for_flavor.end()) {
        if (expected_action_for_flavor[flavor] != action) {
          return false; // Two creations or two destructions in a row for this flavor
        }
      }

      // If we see action A, the next expected action is !A
      expected_action_for_flavor[flavor] = !action;
    }

    // Since total operators per flavor must be even (equal creates and destroys)
    // and they alternate, they will naturally be valid globally.
    // We check that the last expected action matches the first action of that flavor.
    std::map<int, int> first_action_for_flavor;
    for (int op : this->ops) {
      int action = (op >> action_bit_index) & 1;
      int flavor = op & ~(1 << action_bit_index);
      if (first_action_for_flavor.find(flavor) == first_action_for_flavor.end()) {
        first_action_for_flavor[flavor] = action;
      }
    }

    for (const auto& kv : expected_action_for_flavor) {
      if (kv.second != first_action_for_flavor[kv.first]) {
        return false;
      }
    }

    return true;
  }

  bool OperatorSequence::verify_infinite_U_constraints() const {
    int n = this->order;

    for (int i = 0; i + 1 < n; ++i) {
      int op_late  = this->ops[i];     // Applied second
      int op_early = this->ops[i + 1]; // Applied first

      // For HubbardAtom, Action is Bit 0 (LSB), Spin is Bit 1 (MSB)
      int late_action = op_late & 1;
      int early_action = op_early & 1;

      // Rule A: Strict Alternation globally
      if (early_action == late_action) return false;

      // Rule B: Spin Conservation while occupied
      if (early_action == 1) { // early is create
        int late_spin = (op_late >> 1) & 1;
        int early_spin = (op_early >> 1) & 1;
        if (late_spin != early_spin) return false;
      }
    }

    // Check beta/zero boundary
    if (n > 0) {
      int op_beta = this->ops[0];     // The very last operator applied (latest time)
      int op_zero = this->ops[n - 1]; // The very first operator applied (earliest time)
      int beta_action = op_beta & 1;
      int zero_action = op_zero & 1;

      if (beta_action == zero_action) return false;

      if (beta_action == 1) { // beta is create
        int beta_spin = (op_beta >> 1) & 1;
        int zero_spin = (op_zero >> 1) & 1;
        if (beta_spin != zero_spin) return false;
      }
    }

    return true;
  }

} // namespace sc_expansion
