// Standalone utility: enumerate the full (non-bipartite) vacuum-diagram set used
// by the staggered-dimer free-energy expansion and write it to the on-disk cache.
//
// Run this ONCE per order (a single process) before launching a large mcmc_dimer
// scan. Every rank of every scan job then loads the cache instead of repeating
// the (CPU-heavy) enumeration. mcmc_dimer also self-warms the cache on a miss,
// so this prep step is an optimization, not a correctness requirement — but it
// removes redundant generation across the array and keeps job startup cheap.
//
// Usage:
//   generate_dimer_diagrams <order>

#include "sc_expansion/generate_diagrams.hpp"
#include "sc_expansion/graph.hpp"
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>

using namespace sc_expansion;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <order>" << std::endl;
    return 1;
  }
  int order = std::atoi(argv[1]);

  // Idempotent: if the cache already exists, do nothing (so re-running the prep
  // step, or a SLURM requeue, is free).
  std::vector<Graph> existing;
  if (load_vacuum_graphs(order, /*bipartite_only=*/false, existing)) {
    std::cout << "Cache already present at " << vacuum_diagrams_path(order, false) << " (" << existing.size()
              << " diagrams); nothing to do." << std::endl;
    return 0;
  }

  VacuumDiagramGenerator gen(order, /*bipartite_only=*/false);
  reset_embedding_timer();
  auto t0 = std::chrono::steady_clock::now();
  gen.generate();
  auto t1              = std::chrono::steady_clock::now();
  double gen_seconds   = std::chrono::duration<double>(t1 - t0).count();
  double embed_seconds = get_embedding_time_seconds();
  auto const &graphs   = gen.get_unique_graphs();

  std::cout << "Total non-bipartite vacuum diagrams at order " << order << ": " << graphs.size() << std::endl;
  std::cout << "Generation time: " << std::fixed << std::setprecision(3) << gen_seconds << " s"
            << "  (embedding: " << embed_seconds << " s, rest: " << (gen_seconds - embed_seconds) << " s)" << std::endl;

  save_vacuum_graphs(order, /*bipartite_only=*/false, graphs);
  std::cout << "Saved to " << vacuum_diagrams_path(order, false) << std::endl;
  return 0;
}
