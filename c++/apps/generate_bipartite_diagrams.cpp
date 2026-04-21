// Standalone utility: enumerate bipartite vacuum diagrams at a given order.
//
// Usage:
//   generate_bipartite_diagrams <order>

#include "sc_expansion/generate_diagrams.hpp"
#include <cstdlib>
#include <iomanip>
#include <iostream>

using namespace sc_expansion;

static void print_bipartite_diagrams(int order) {
  VacuumDiagramGenerator gen(order, true);
  gen.generate();
  auto const &graphs = gen.get_unique_graphs();

  for (size_t d = 0; d < graphs.size(); d++) {
    auto const &g = graphs[d];
    int V         = g.get_V();
    auto adj      = g.get_canonical_form();
    std::cout << "# diagram " << d << "  V=" << V << "  aut=" << g.get_automorphism_count()
              << "  free_mult=" << g.get_free_multiplicity() << std::endl;
    for (int i = 0; i < V; i++) {
      for (int j = 0; j < V; j++) { std::cout << std::setw(3) << (int)adj[i * V + j]; }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << "Total bipartite diagrams at order " << order << ": " << graphs.size() << std::endl;

  save_bipartite_graphs(order, graphs);
  std::cout << "Saved to " << bipartite_diagrams_path(order) << std::endl;
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <order>" << std::endl;
    return 1;
  }
  int order = std::atoi(argv[1]);
  print_bipartite_diagrams(order);
  return 0;
}
