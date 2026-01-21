#pragma once

#include "sc_expansion/diagram.hpp"
#include <./triqs/mc_tools/random_generator.hpp>
#include <functional>

class Configuration {

  public:
  double beta;
  int order;
  double weight;
  double sign;
  std::vector<double> state;

  Configuration(double U_, double beta_, double mu_, std::vector<sc_expansion::adjmat> diagram_mats_, int order_)
     : beta(beta_), order(order_), state(order_), U(U_), mu(mu_), diagram_mats(diagram_mats_) {

    triqs::mc_tools::random_generator RNG("mt19937", 23432);
    for (int i = 0; i < this->order; i++) { this->state[i] = RNG(this->beta); }

    this->diagrams.reserve(this->diagram_mats.size());
    for (auto const &mat : this->diagram_mats) { this->diagrams.push_back(sc_expansion::Diagram(mat, this->U, this->beta, this->mu)); }

    recompute_weight_and_sign();
  }

  std::pair<double, double> weight_and_sign(std::vector<double> const &st) const {
    double res = 0.0;

    for (auto const &diagram : this->diagrams) { res -= diagram.evaluate_at_taus(st); }
    double w = std::abs(res);
    double s = (res < 0) ? -1.0 : 1.0;
    return {w, s};
  }

  // Fast commit method
  void commit_update(double new_weight, double new_sign) {
    this->weight            = new_weight;
    this->sign              = new_sign;
    this->ref_weight_calced = false;
  }

  std::pair<double, double> weight_and_sign() const { return weight_and_sign(this->state); }
  void set_reference_weight_function(std::function<double(std::vector<double> const &)> f) { reference_weight_func = f; }

  double get_reference_weight() const {
    if (!reference_weight_func) return 1.0;
    if (!ref_weight_calced) {
      cached_ref_weight = reference_weight_func(this->state);
      ref_weight_calced = true;
    }
    return cached_ref_weight;
  }

  private:
  std::function<double(std::vector<double> const &)> reference_weight_func;
  mutable double cached_ref_weight;
  mutable bool ref_weight_calced = false;

  void recompute_weight_and_sign() {
    auto [w, s]             = this->weight_and_sign();
    this->weight            = w;
    this->sign              = s;
    this->ref_weight_calced = false;
  }

  double U;
  double mu;
  std::vector<sc_expansion::adjmat> diagram_mats;
  std::vector<sc_expansion::Diagram> diagrams;
};