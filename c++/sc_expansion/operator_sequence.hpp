#pragma once

#include <vector>

namespace sc_expansion {

  struct OperatorSequence {
    std::vector<double> taus;
    std::vector<int> ops;
    double permutation_sign;
    int order;

    OperatorSequence(std::vector<double> taus, std::vector<int> ops);

    bool verify_flavor_alternation(int action_bit_index) const;
    bool verify_infinite_U_constraints() const;

    private:
    void sort_sequence();
  };

} // namespace sc_expansion