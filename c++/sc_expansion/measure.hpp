#pragma once
#include "configuration.hpp"
#include <triqs/arrays.hpp>
#include <vector>
#include <algorithm>

struct measure {

  Configuration *config;
  double average_sign;
  int count;

  // Storage for accumulation (references to external)
  std::vector<double> &signs_vec;
  std::vector<double> &refs_vec;

  // Storage for results (references to external, rank 0)
  nda::array<double, 1> &signs;
  nda::array<double, 1> &reference_vals;

  measure(Configuration *config_, std::vector<double> &s_vec, std::vector<double> &r_vec, nda::array<double, 1> &s, nda::array<double, 1> &r)
     : config(config_), average_sign(0.0), count(0), signs_vec(s_vec), refs_vec(r_vec), signs(s), reference_vals(r) {}

  void accumulate(double) {
    average_sign += config->sign;
    signs_vec.push_back(config->sign);
    refs_vec.push_back(config->get_reference_weight() / config->weight);
    count++;
  }

  void collect_results(mpi::communicator c) {

    double sum_signs = mpi::reduce(average_sign, c);
    int total_count  = mpi::reduce(count, c);

    // Convert local vectors to nda::array for reduction
    nda::array<double, 1> local_signs(signs_vec.size());
    std::copy(signs_vec.begin(), signs_vec.end(), local_signs.begin());

    nda::array<double, 1> local_refs(refs_vec.size());
    std::copy(refs_vec.begin(), refs_vec.end(), local_refs.begin());

    // Reduce to rank 0 and store in shared result arrays
    signs          = mpi::reduce(local_signs, c);
    reference_vals = mpi::reduce(local_refs, c);

    if (c.rank() == 0) {
      double final_average_sign = sum_signs / total_count;
      std::cout << "Average Sign: " << final_average_sign << std::endl;
    }
  }
};