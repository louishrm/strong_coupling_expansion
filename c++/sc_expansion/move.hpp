#pragma once
#include "configuration.hpp"
#include <triqs/mc_tools/random_generator.hpp>

// Optimized Move Struct
struct move {
  double current_weight;
  double current_sign;

  // Store only what is needed to revert the move, not the whole state
  int changed_index;
  double old_tau;
  double new_tau;

  Configuration *config;
  triqs::mc_tools::random_generator &RNG;

  move(Configuration *config_, triqs::mc_tools::random_generator &RNG_) : config(config_), RNG(RNG_) {
    // Cache initial weight
    current_weight = config->weight;
    current_sign   = config->sign;
  }

  double attempt() {
    // 1. Pick random change
    changed_index = RNG(config->order);
    new_tau       = RNG(0.0, config->beta);

    // 2. Save old value for potential revert
    old_tau = config->state[changed_index];

    // 3. APPLY CHANGE DIRECTLY (No vector copy)
    config->state[changed_index] = new_tau;

    // 4. Calculate new weight
    // Optimization: If possible, create a function that computes
    // weight ratio based on the change, rather than full recompute.
    auto [new_weight, new_sign] = config->weight_and_sign();

    // Store potential new values
    // We don't overwrite config->weight yet, only on accept
    double acceptance_ratio = new_weight / current_weight;

    // Update temp values for accept() to use
    current_weight = new_weight;
    current_sign   = new_sign;

    return acceptance_ratio;
  }

  double accept() {
    // State is already updated in attempt(), just commit the weight/sign
    config->commit_update(current_weight, current_sign);
    return 1.0;
  }

  void reject() {
    // REVERT the state
    config->state[changed_index] = old_tau;
    // No need to revert weight/sign as we didn't commit them to config
  }
};