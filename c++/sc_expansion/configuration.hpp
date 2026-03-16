#pragma once

#include <vector>
#include <deque>
#include <cmath>
#include <utility>
#include "diagram.hpp"
#include "generate_diagrams.hpp"
#include "dual.hpp"
#include <./triqs/mc_tools/random_generator.hpp>

class Configuration {

  public:
  double beta;
  std::vector<double> state; //set of imaginary times
  double metropolis_weight;  //|alpha Uinf+ (1-alpha) Ufin|
  bool bipartite;

  double integrand;           //Omega(configuration)
  double reference_integrand; //Omega_infinite_U(configuration)

  Configuration(sc_expansion::Parameters<double> const &params, int order_, double alpha_)
     : beta(params.beta), bipartite(params.bipartite), U(params.U), mu(params.mu), order(order_), alpha(alpha_) {

    this->state.resize(this->order);
    triqs::mc_tools::random_generator RNG("mt19937", 23432);
    for (int i = 0; i < this->order; i++) { this->state[i] = RNG(this->beta); }

    // Generate diagrams
    sc_expansion::VacuumDiagramGenerator gen(this->order, params.bipartite);
    gen.generate();
    const auto &unique_graphs = gen.get_unique_graphs();

    this->diagrams.reserve(unique_graphs.size());
    this->evaluators.reserve(unique_graphs.size());

    // Construct Diagram objects from the unique graphs
    for (auto const &g : unique_graphs) { this->diagrams.emplace_back(sc_expansion::Diagram(g)); }

    // Construct DiagramEvaluators for each diagram
    for (auto const &diag : this->diagrams) { this->evaluators.emplace_back(diag, params); }

    this->recompute_integrands();
    this->metropolis_weight = this->get_metropolis_weight();
  }

  Configuration(sc_expansion::Parameters<Dual> const &params, int order_, double alpha_)
     : beta(params.beta.value), bipartite(params.bipartite), U(params.U.value), mu(params.mu.value), order(order_), alpha(alpha_) {

    this->state.resize(this->order);
    triqs::mc_tools::random_generator RNG("mt19937", 23432);
    for (int i = 0; i < this->order; i++) { this->state[i] = RNG(this->beta); }

    // Generate diagrams
    sc_expansion::VacuumDiagramGenerator gen(this->order, params.bipartite);
    gen.generate();
    const auto &unique_graphs = gen.get_unique_graphs();

    this->diagrams.reserve(unique_graphs.size());
    this->dual_evaluators.reserve(unique_graphs.size());

    // Construct Diagram objects from the unique graphs
    for (auto const &g : unique_graphs) { this->diagrams.emplace_back(sc_expansion::Diagram(g)); }

    // Construct DiagramEvaluators for each diagram
    for (auto const &diag : this->diagrams) { this->dual_evaluators.emplace_back(diag, params); }

    this->use_dual = true;
    this->recompute_integrands();
    this->metropolis_weight = this->get_metropolis_weight();
  }

  int get_order() const { return this->order; }
  double get_U() const { return this->U; }

  std::pair<double, double> get_integrands() {
    double finite_U   = 0.0;
    double infinite_U = 0.0;

    if (this->use_dual) {
      for (auto const &evaluator : this->dual_evaluators) {
        finite_U += evaluator.evaluate_at_taus(this->state, false, true).derivative;
        infinite_U += evaluator.evaluate_at_taus(this->state, true, true).derivative;
      }
    } else {
      for (auto const &evaluator : this->evaluators) {
        finite_U += evaluator.evaluate_at_taus(this->state, false, true);
        infinite_U += evaluator.evaluate_at_taus(this->state, true, true);
      }
    }
    return {finite_U, infinite_U};
  }

  double calculate_weight(double finite_U, double infinite_U) const { return std::abs(finite_U - infinite_U) + alpha * std::abs(infinite_U); }

  double get_metropolis_weight() { return std::abs(this->integrand - this->reference_integrand) + alpha * std::abs(this->reference_integrand); }

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
  bool use_dual = false;
  std::vector<sc_expansion::Diagram> diagrams;
  std::vector<sc_expansion::DiagramEvaluator<double>> evaluators;
  std::vector<sc_expansion::DiagramEvaluator<Dual>> dual_evaluators;

  void recompute_integrands() {
    auto [finite_U, infinite_U] = this->get_integrands();
    this->integrand             = finite_U;
    this->reference_integrand   = infinite_U;
  }
};