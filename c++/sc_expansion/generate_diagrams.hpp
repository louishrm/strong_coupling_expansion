#pragma once
#include <string>
#include <vector>
#include <unordered_set>
#include "graph.hpp"
#include "combinatorics.hpp"

namespace sc_expansion {

  // On-disk cache for generated vacuum diagrams. Stored as a binary blob in the
  // project-root "diagrams/" directory so MCMC apps can skip regeneration at
  // high orders where diagram enumeration is expensive.
  std::string bipartite_diagrams_path(int order, bool allow_self_loops = false);
  void save_bipartite_graphs(int order, const std::vector<Graph> &graphs, bool allow_self_loops = false);
  // Returns true and fills out_graphs on success; false if the file is missing.
  bool load_bipartite_graphs(int order, std::vector<Graph> &out_graphs, bool allow_self_loops = false);


  struct VectorHasher {
    size_t operator()(const std::vector<uint8_t> &v) const {
      size_t seed = v.size();
      for (auto x : v) { seed ^= std::hash<uint8_t>{}(x) + 0x9e3779b9 + (seed << 6) + (seed >> 2); }
      return seed;
    }
  };

  std::vector<uint8_t> generate_n_cycle_adjacency_matrix(int n);

  int calculate_n_cycle_free_multiplicity(int n, bool bipartite);

  class VacuumDiagramGenerator {

    public:
    VacuumDiagramGenerator(int order, bool bipartite_only = true, bool allow_self_loops = false);

    void generate();
    const std::vector<Graph> &get_unique_graphs() const { return graphs; }

    private:
    int order;
    bool bipartite_only;
    bool allow_self_loops;
    PartitionGenerator partitions;
    std::unordered_set<std::vector<uint8_t>, VectorHasher> unique_adjmats;
    std::vector<Graph> graphs;

    void propose_create_process(const std::vector<int> &partition);

    void process_graph(const Graph &graph);

    std::vector<uint8_t> fill_matrix(const std::vector<int> &source, const std::vector<int> &target, int V) const;
  };
} // namespace sc_expansion
