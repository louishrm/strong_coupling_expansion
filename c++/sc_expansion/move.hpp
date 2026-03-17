#pragma once
#include "configuration.hpp"
#include <triqs/mc_tools/random_generator.hpp>

// Optimized Move Struct
template <typename T> struct move {
  double current_weight;

  // Proposed values
  double proposed_weight;
  double proposed_integrand;
  double proposed_ref_integrand;

  // Store only what is needed to revert the move, not the whole state
  int changed_index;
  double old_tau;
  double new_tau;

  Configuration<T> *config;
  triqs::mc_tools::random_generator &RNG;

  move(Configuration<T> *config_, triqs::mc_tools::random_generator &RNG_) : config(config_), RNG(RNG_) {
    // Cache initial weight
    current_weight = config->metropolis_weight;
  }

  double attempt() {
    // 1. Pick random change
    changed_index = RNG(config->state.size());
    new_tau       = RNG(config->beta);

    // 2. Save old value for potential revert
    old_tau = config->state[changed_index];

    // 3. APPLY CHANGE DIRECTLY (No vector copy)
    config->state[changed_index] = new_tau;

    // 4. Calculate new integrands
    auto [new_finite, new_infinite] = config->get_integrands();
    proposed_integrand              = new_finite;
    proposed_ref_integrand          = new_infinite;

    // 5. Calculate new weight
    proposed_weight = config->calculate_weight(proposed_integrand, proposed_ref_integrand);

    // 6. Return acceptance ratio
    if (current_weight == 0.0) {
      // If current weight is 0, we accept any proposal with non-zero weight
      return (proposed_weight > 0.0) ? 1.0e100 : 1.0;
    }

    return proposed_weight / current_weight;
  }

  double accept() {
    // State is already updated in attempt(), just commit the values
    config->commit_update(proposed_integrand, proposed_ref_integrand);
    current_weight = proposed_weight;
    return 1.0;
  }

  void reject() {
    // REVERT the state
    config->state[changed_index] = old_tau;
    // No need to revert weight in config as we didn't commit it
  }
};
