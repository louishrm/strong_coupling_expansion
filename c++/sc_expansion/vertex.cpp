#include "vertex.hpp"

namespace sc_expansion {

  template <int N_sites, typename T> VertexType<N_sites, T>::VertexType(int n_legs) : n_legs(n_legs){};

  template <int N_sites, typename T>
  T VertexType<N_sites, T>::evaluate(std::vector<double> const &taus, std::vector<int> const &op_indices, HubbardSolver<N_sites, T> const &solver,
                                     bool infinite_U) const {
    // Check cache first
    VertexCacheKey key{op_indices, taus, solver.is_infinite_U()};
    auto it = global_cache.find(key);
    if (it != global_cache.end()) { return it->second; }

    // Map args to lists
    auto [unprimed, primed] = map_args_to_list(taus, op_indices);

    // Compute the vertex value using the solver and the mapped arguments
    T vertex_value = compute_cumulant_decomposition(unprimed, primed, solver, infinite_U);

    // Store in cache before returning
    global_cache[key] = vertex_value;
    return vertex_value;
  }

  template <int N_sites, typename T>
  std::pair<ArgList, ArgList> VertexType<N_sites, T>::map_args_to_list(std::vector<double> const &taus, std::vector<int> const &op_indices) const {
    ArgList unprimed, primed;
    for (size_t i = 0; i < op_indices.size(); ++i) {
      int op     = op_indices[i];
      double tau = taus[i];
      if (op % 2 == 0) { // Even op index -> unprimed
        unprimed.emplace_back(tau, op);
      } else { // Odd op index -> primed
        primed.emplace_back(tau, op);
      }
    }
    return {unprimed, primed};
  }
} // namespace sc_expansion