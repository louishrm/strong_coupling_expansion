#pragma once
#include "configuration_base.hpp"
#include <triqs/mc_tools/random_generator.hpp>

template <typename T> struct move {
  double current_weight;
  double proposed_weight;

  // Store only what is needed to revert the move
  int changed_index;
  double old_tau;

  ConfigurationBase<T> *config;
  triqs::mc_tools::random_generator &RNG;

  move(ConfigurationBase<T> *config_, triqs::mc_tools::random_generator &RNG_) : config(config_), RNG(RNG_) {
    current_weight = config->metropolis_weight;
  }

  double attempt() {
    // 1. Pick random tau index and new tau value
    changed_index = RNG(config->state.size());
    new_tau       = RNG(config->beta);

    // 2. Save old value for potential revert
    old_tau = config->state[changed_index];

    // 3. Apply change directly
    config->state[changed_index] = new_tau;

    // 4. Evaluate proposed weight (configuration stashes integrands internally)
    proposed_weight = config->evaluate_proposed();

    // 5. Return acceptance ratio
    if (current_weight == 0.0) { return (proposed_weight > 0.0) ? 1.0e100 : 1.0; }

    return proposed_weight / current_weight;
  }

  double accept() {
    config->commit_proposal();
    current_weight = proposed_weight;
    return 1.0;
  }

  void reject() {
    // Revert the state; config's committed cache is untouched
    config->state[changed_index] = old_tau;
  }

  private:
  double new_tau;
};
