#pragma once

#include <vector>
#include "free_energy_calculator.hpp"
#include "dual.hpp"

template <typename T> class DiagMCConfiguration {

  public:
  DiagMCConfiguration(sc_expansion::Parameters<T> const &params, int order,
                      sc_expansion::FreeEnergyCalculator<2, T> &calculator, double alpha = 0.001);

  // Evaluate a single physical diagram (1-indexed: d=1..N maps to diagrams[d-1]).
  // Returns the signed value as double (extracts derivative for Dual).
  double evaluate_diagram(int d);

  // Mark dirty for a specific physical diagram (1-indexed)
  void mark_diagram_tau_dirty(int d, int tau_index);
  void mark_diagram_all_dirty(int d);

  // Clear global VertexType caches (call after tau moves, not after diagram moves)
  void clear_global_caches();

  // Accessors
  double get_alpha() const;
  double get_defensive_per_diagram() const;
  int get_n_diagrams() const;
  int get_order() const;
  double get_U() const;
  sc_expansion::FreeEnergyCalculator<2, T> const &get_calculator() const;

  // Public state (accessible by move and measure)
  sc_expansion::Parameters<T> const &params;
  double beta;
  std::vector<double> state;   // tau values (size = order)
  int current_diagram;         // 1..N (physical diagrams only)
  double current_abs;          // |D_d(tau)| + c = Metropolis weight (defensive)
  double current_signed;       // D_d(tau) (signed value for measurement)

  private:
  int order;
  int n_diagrams;              // N = number of physical diagrams
  double alpha;                // total defensive weight
  double defensive_per_diagram; // c = alpha / N
  sc_expansion::FreeEnergyCalculator<2, T> &calculator;
};
