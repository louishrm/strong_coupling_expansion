// Standalone utility: enumerate bipartite vacuum diagrams at a given order.
//
// Usage:
//   generate_bipartite_diagrams <order>

#include "sc_expansion/generate_diagrams.hpp"
#include "sc_expansion/graph.hpp"
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>

using namespace sc_expansion;

static void print_bipartite_diagrams(int order) {
  VacuumDiagramGenerator gen(order, true);
  reset_embedding_timer();
  auto t0 = std::chrono::steady_clock::now();
  gen.generate();
  auto t1 = std::chrono::steady_clock::now();
  double gen_seconds = std::chrono::duration<double>(t1 - t0).count();
  double embed_seconds = get_embedding_time_seconds();
  auto const &graphs = gen.get_unique_graphs();
  std::cout << "Total bipartite diagrams at order " << order << ": " << graphs.size() << std::endl;
  std::cout << "Generation time: " << std::fixed << std::setprecision(3) << gen_seconds << " s"
            << "  (embedding: " << embed_seconds << " s"
            << ", rest: " << (gen_seconds - embed_seconds) << " s)" << std::endl;

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
