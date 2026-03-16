#pragma once
#include <vector>
#include <unordered_set>
#include "graph.hpp"
#include "combinatorics.hpp"

namespace sc_expansion {

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
    VacuumDiagramGenerator(int order, bool bipartite_only = true);

    void generate();
    const std::vector<Graph> &get_unique_graphs() const { return graphs; }

    private:
    int order;
    bool bipartite_only;
    PartitionGenerator partitions;
    std::unordered_set<std::vector<uint8_t>, VectorHasher> unique_adjmats;
    std::vector<Graph> graphs;

    void propose_create_process(const std::vector<int> &partition);

    void process_graph(const Graph &graph);

    std::vector<uint8_t> fill_matrix(const std::vector<int> &source, const std::vector<int> &target, int V) const;
  };
} // namespace sc_expansion
