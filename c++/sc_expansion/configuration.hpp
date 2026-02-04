#pragma once

#include <vector>
#include <cmath>
#include <utility>
#include "diagram.hpp"
#include <./triqs/mc_tools/random_generator.hpp>

class Configuration {

  public:
  double beta;
  std::vector<double> state; //set of imaginary times
  double metropolis_weight;  //|alpha Uinf+ (1-alpha) Ufin|

  double integrand;           //Omega(configuration)
  double reference_integrand; //Omega_infinite_U(configuration)

  Configuration(double U_, double beta_, double mu_, int order_, double alpha_, std::vector<sc_expansion::adjmat> const &diagram_mats_)
     : beta(beta_), U(U_), mu(mu_), order(order_), alpha(alpha_), diagram_mats(diagram_mats_) {

    this->state.resize(this->order);
    triqs::mc_tools::random_generator RNG("mt19937", 23432);
    for (int i = 0; i < this->order; i++) { this->state[i] = RNG(this->beta); }

    for (auto const &mat : diagram_mats) { diagrams.push_back(sc_expansion::Diagram(mat, U, beta, mu)); }

    this->recompute_integrands();
    this->metropolis_weight = this->get_metropolis_weight();
  }

  std::pair<double, double> get_integrands() {

    double finite_U   = 0.0;
    double infinite_U = 0.0;
    for (auto const &diagram : this->diagrams) {
      finite_U += diagram.evaluate_at_taus(this->state, false);
      infinite_U += diagram.evaluate_at_taus(this->state, true);
    }
    return {finite_U, infinite_U};
  }

  double calculate_weight(double finite_U, double infinite_U) const { return std::abs(alpha * infinite_U + (1 - alpha) * finite_U); }

  double get_metropolis_weight() { return std::abs(alpha * this->reference_integrand + (1 - alpha) * this->integrand); }

  void commit_update(double new_integrand, double new_reference_integrand) {
    this->integrand           = new_integrand;
    this->reference_integrand = new_reference_integrand;
    this->metropolis_weight   = this->get_metropolis_weight();
  }

  private:
  double U;
  double mu;
  int order;
  double alpha;
  std::vector<sc_expansion::adjmat> diagram_mats;
  std::vector<sc_expansion::Diagram> diagrams;

  void recompute_integrands() {
    auto [finite_U, infinite_U] = this->get_integrands();
    this->integrand             = finite_U;
    this->reference_integrand   = infinite_U;
  }
};