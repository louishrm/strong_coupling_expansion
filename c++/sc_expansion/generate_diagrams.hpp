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

  class VacuumDiagramGenerator {

    public:
    VacuumDiagramGenerator(int order);

    void generate();
    const std::unordered_set<std::vector<uint8_t>, VectorHasher> &get_unique_graphs() const { return unique_graphs; }

    private:
    int order;
    PartitionGenerator partitions;
    std::unordered_set<std::vector<uint8_t>, VectorHasher> unique_graphs;

    void propose_create_process(const std::vector<int> &partition);

    void process_graph(const Graph &graph);

    std::vector<uint8_t> fill_matrix(const std::vector<int> &source, const std::vector<int> &target, int V) const;
  };
} // namespace sc_expansion