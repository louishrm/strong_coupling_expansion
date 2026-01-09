#pragma once

#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/atom_diag/functions.hpp>
#include <vector>
#include <nda/nda.hpp>
#include <utility>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <tuple>

namespace sc_expansion {

  class HubbardAtom {
    public:
    using cumul_args = std::vector<std::pair<double, int>>;

    // Attributes
    double U;
    double beta;
    double mu;
    triqs::operators::many_body_operator_generic<double> H;
    triqs::atom_diag::atom_diag<false> ad;

    // Constructor
    HubbardAtom(double U, double beta, double mu);

    // Methods
    static int calculate_permutation_sign(const std::vector<int> &p);

    static std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, int>
    sort_operators(const std::vector<double> &times, const std::vector<int> &spins, const std::vector<int> &flags);

    nda::matrix<double> make_interaction_picture_destroy_op(double tau, int state_index) const;

    nda::matrix<double> make_interaction_picture_create_op(double tau, int state_index) const;

    double G0(cumul_args const &unprimed_args, cumul_args const &primed_args) const;

    private:
    static triqs::hilbert_space::fundamental_operator_set make_fops();
    static triqs::operators::many_body_operator_generic<double> make_H0(double U, double mu);
  };

} // namespace sc_expansion