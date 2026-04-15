#pragma once
#include "diagmc_configuration.hpp"
#include <triqs/mc_tools/random_generator.hpp>
#include <cmath>

template <typename T> struct DiagMCMove {

  DiagMCConfiguration<T> *config;
  triqs::mc_tools::random_generator &RNG;

  // Probability of tau move (vs diagram move). Default: 10 tau moves per 1 diagram move.
  double tau_move_prob;

  // Per-diagram defensive constant c = alpha / N
  double c;

  // Internal state for revert
  bool is_diagram_move;
  int changed_tau_index;
  double old_tau;
  int proposed_diagram;
  double proposed_abs;
  double proposed_signed;

  DiagMCMove(DiagMCConfiguration<T> *config_, triqs::mc_tools::random_generator &RNG_, double tau_move_prob_ = 10.0 / 11.0)
     : config(config_), RNG(RNG_), tau_move_prob(tau_move_prob_), c(config_->get_defensive_per_diagram()) {}

  double attempt() {
    double r = RNG(1.0);

    if (r < this->tau_move_prob) {
      // === TAU MOVE ===
      this->is_diagram_move  = false;
      this->changed_tau_index = RNG(config->state.size());
      double new_tau         = RNG(config->beta);
      this->old_tau          = config->state[this->changed_tau_index];
      config->state[this->changed_tau_index] = new_tau;

      // Mark dirty and evaluate current diagram
      config->mark_diagram_tau_dirty(config->current_diagram, this->changed_tau_index);
      this->proposed_signed = config->evaluate_diagram(config->current_diagram);
      this->proposed_abs    = std::abs(this->proposed_signed) + this->c;

    } else {
      // === DIAGRAM MOVE ===
      this->is_diagram_move  = true;
      // Propose uniformly from {1, ..., N}
      this->proposed_diagram = 1 + RNG(config->get_n_diagrams());

      if (this->proposed_diagram == config->current_diagram) {
        // Self-transition: accept trivially, no evaluation needed
        this->proposed_abs    = config->current_abs;
        this->proposed_signed = config->current_signed;
      } else {
        // Evaluate proposed diagram from scratch
        config->mark_diagram_all_dirty(this->proposed_diagram);
        this->proposed_signed = config->evaluate_diagram(this->proposed_diagram);
        this->proposed_abs    = std::abs(this->proposed_signed) + this->c;
      }
    }

    // Metropolis acceptance ratio (current_abs > 0 always due to defensive floor)
    return this->proposed_abs / config->current_abs;
  }

  double accept() {
    if (this->is_diagram_move) { config->current_diagram = this->proposed_diagram; }
    config->current_abs    = this->proposed_abs;
    config->current_signed = this->proposed_signed;
    return 1.0;
  }

  void reject() {
    if (!this->is_diagram_move) {
      // Revert the tau and re-dirty so the stale local cache is not reused
      config->state[this->changed_tau_index] = this->old_tau;
      config->mark_diagram_tau_dirty(config->current_diagram, this->changed_tau_index);
    }
    // Diagram moves don't modify state on proposal, nothing to revert.
  }
};
