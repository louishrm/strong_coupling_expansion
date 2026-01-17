#pragma once
#include "configuration.hpp"
#include <triqs/mc_tools/random_generator.hpp>

struct move {

  double weight, sign;
  std::vector<double> new_state;
  Configuration *config;
  triqs::mc_tools::random_generator &RNG;

  move(Configuration *config_, triqs::mc_tools::random_generator &RNG_) : config(config_), RNG(RNG_) {}

  double attempt() {

    int random_time_index        = RNG(config->order);
    double new_tau               = RNG(0.0, config->beta);
    new_state                    = config->state;
    new_state[random_time_index] = new_tau;

    auto [new_weight, new_sign] = config->weight_and_sign(new_state);
    weight                      = new_weight;
    sign                        = new_sign;

    return weight / config->weight;
  }

  double accept() {
    config->update_state(new_state, weight, sign);
    return 1.0;
  }

  void reject() {}
};
