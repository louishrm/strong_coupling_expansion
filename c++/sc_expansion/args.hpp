#pragma once
#include "fock_space.hpp"
#include "combinatorics.hpp"
#include <vector>

namespace sc_expansion {

  template <typename T> struct Args {
    std::vector<double> taus;               // imaginary times
    std::vector<FermionOperator<T>> ops;    // fermion operators, should always come as c^\dag, c, c^\dag, c...

    Args(std::vector<double> taus_, std::vector<FermionOperator<T>> ops_);

    int order;
    static constexpr int N_ORBITALS = 2; // Number of orbitals (spin up and down)

    double permutation_sign; // sign from sorting operators by imaginary time

    // Permutation applied during sort_args():
    //   taus[sorted_pos] = (input_taus_before_sort)[argsort[sorted_pos]]
    // i.e. argsort[sorted_pos] = the original (pre-sort) position of the element now at sorted_pos.
    std::vector<int> argsort;

    void sort_args();
    bool operator_sequence_is_valid() const;
    bool operator_sequence_is_valid_infinite_U() const;

    static std::pair<Args<T>, Args<T>> split_from_raw(std::vector<double> const &taus, std::vector<uint8_t> const &op_ids);
  };
} // namespace sc_expansion
