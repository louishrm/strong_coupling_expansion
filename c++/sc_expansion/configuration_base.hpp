#pragma once

#include <vector>
#include <cmath>
#include "hubbard_solver.hpp"
#include "dual.hpp"

namespace sc_expansion {
  template <typename T> double get_val(const T &x) { return x; }
  template <> inline double get_val<Dual>(const Dual &x) { return x.value; }
} // namespace sc_expansion

template <typename T> class ConfigurationBase {

  public:
  sc_expansion::Parameters<T> const &params;
  double beta;
  bool bipartite;
  std::vector<double> state; // set of imaginary times
  double metropolis_weight;  // |W| for acceptance ratio

  // Called by move: evaluate all diagrams at current state, return proposed weight.
  // Internally stashes whatever the derived class needs for commit.
  virtual double evaluate_proposed() = 0;

  // Commit the internally-stashed proposal into the committed state.
  virtual void commit_proposal() = 0;

  // Called by measure: expose the committed integrand values.
  virtual double get_integrand() const           = 0;
  virtual double get_reference_integrand() const = 0; // returns 0 if unused

  // Called by move to propagate dirty marking to diagrams.
  // Default no-op; subclasses with a FreeEnergyCalculator override this.
  virtual void mark_tau_dirty(int /*tau_index*/) {}

  int get_order() const { return this->order; }
  double get_U() const { return sc_expansion::get_val(this->params.U); }

  virtual ~ConfigurationBase() = default;

  protected:
  int order;

  ConfigurationBase(sc_expansion::Parameters<T> const &params_, int order_)
     : params(params_), beta(sc_expansion::get_val(params_.beta)), bipartite(params_.bipartite), metropolis_weight(0.0), order(order_) {}
};
