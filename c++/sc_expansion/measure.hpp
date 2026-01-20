#pragma once
#include "configuration.hpp"
#include <triqs/arrays.hpp>
#include <vector>
#include <algorithm>

struct measure {

  Configuration *config;
  double average_sign;
  double average_ref;
  int count;

  // Storage for accumulation (references to external)
  std::vector<double> &signs_vec;
  std::vector<double> &refs_vec;

  // Storage for results (references to external, rank 0)
  nda::array<double, 1> &signs;
  nda::array<double, 1> &reference_vals;

  measure(Configuration *config_, std::vector<double> &s_vec, std::vector<double> &r_vec, nda::array<double, 1> &s, nda::array<double, 1> &r)
     : config(config_), average_sign(0.0), average_ref(0.0), count(0), signs_vec(s_vec), refs_vec(r_vec), signs(s), reference_vals(r) {}

  void accumulate(double) {
    average_sign += config->sign;
    average_ref += config->get_reference_weight() / config->weight;
    signs_vec.push_back(config->sign);
    refs_vec.push_back(config->get_reference_weight() / config->weight);
    count++;
  }

  void collect_results(mpi::communicator c) {

    double sum_signs = mpi::reduce(average_sign, c);
    double sum_refs  = mpi::reduce(average_ref, c);
    int total_count  = mpi::reduce(count, c);

    // Convert local vectors to nda::array for gathering
    nda::array<double, 1> local_signs(signs_vec.size());
    std::copy(signs_vec.begin(), signs_vec.end(), local_signs.begin());

    nda::array<double, 1> local_refs(refs_vec.size());
    std::copy(refs_vec.begin(), refs_vec.end(), local_refs.begin());

    // Gather to rank 0 and store in shared result arrays
    auto gathered_signs = mpi::gather(local_signs, c, 0);
    auto gathered_refs  = mpi::gather(local_refs, c, 0);

    nda::array<double, 1> all_signs;
    nda::array<double, 1> all_refs;

    all_signs = gathered_signs;
    all_refs  = gathered_refs;

    if (c.rank() == 0) {

      this->signs          = all_signs;
      this->reference_vals = all_refs;

      std::cout << "Size " << all_signs.size() << " " << all_refs.size() << std::endl;

      double final_average_sign = sum_signs / total_count;
      double final_average_ref  = sum_refs / total_count;
      std::cout << "Average Sign: " << final_average_sign << std::endl;
      std::cout << "Average Ref: " << final_average_ref << std::endl;
    }
  }
};