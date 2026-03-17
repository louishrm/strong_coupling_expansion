#pragma once

#include <vector>
#include <deque>
#include <cmath>
#include <utility>
#include "diagram.hpp"
#include "generate_diagrams.hpp"
#include "dual.hpp"
#include <./triqs/mc_tools/random_generator.hpp>

namespace sc_expansion {
  template <typename T> double get_val(const T &x) { return x; }
  template <> inline double get_val<Dual>(const Dual &x) { return x.value; }
} // namespace sc_expansion

template <typename T> class Configuration {

  public:
  sc_expansion::Parameters<T> const &params;
  double beta;
  bool bipartite;
  std::vector<double> state; //set of imaginary times
  double metropolis_weight;  //|alpha Uinf+ (1-alpha) Ufin|

  double integrand;           //Omega(configuration)
  double reference_integrand; //Omega_infinite_U(configuration)

  Configuration(sc_expansion::Parameters<T> const &params_, int order_, double alpha_)
     : params(params_), beta(sc_expansion::get_val(params_.beta)), bipartite(params_.bipartite), order(order_), alpha(alpha_) {

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

  int get_order() const { return this->order; }
  double get_U() const { return sc_expansion::get_val(this->params.U); }

  std::pair<double, double> get_integrands() {
    double finite_U   = 0.0;
    double infinite_U = 0.0;

    for (auto const &evaluator : this->evaluators) {
      if constexpr (std::is_same_v<T, Dual>) {
        finite_U += evaluator.evaluate_at_taus(this->state, false, true).derivative;
        infinite_U += evaluator.evaluate_at_taus(this->state, true, true).derivative;
      } else {
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
  int order;
  double alpha;
  std::vector<sc_expansion::Diagram> diagrams;
  std::vector<sc_expansion::DiagramEvaluator<T>> evaluators;

  void recompute_integrands() {
    auto [finite_U, infinite_U] = this->get_integrands();
    this->integrand             = finite_U;
    this->reference_integrand   = infinite_U;
  }
};
