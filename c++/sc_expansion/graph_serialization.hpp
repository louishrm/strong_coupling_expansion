#pragma once
#include <string>
#include <vector>
#include "graph.hpp"

namespace sc_expansion {

  // On-disk cache for generated vacuum diagrams. Stored as a binary blob in the
  // project-root "diagrams/" directory so MCMC apps can skip regeneration at
  // high orders where diagram enumeration is expensive.

  std::string bipartite_diagrams_path(int order);

  void save_bipartite_graphs(int order, const std::vector<Graph> &graphs);

  // Returns true and fills out_graphs on success; false if the file is missing.
  bool load_bipartite_graphs(int order, std::vector<Graph> &out_graphs);

} // namespace sc_expansion
