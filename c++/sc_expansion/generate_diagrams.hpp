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
  std::string bipartite_diagrams_path(int order);
  void save_bipartite_graphs(int order, const std::vector<Graph> &graphs);
  // Returns true and fills out_graphs on success; false if the file is missing.
  bool load_bipartite_graphs(int order, std::vector<Graph> &out_graphs);


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

  // Rooted diagram generator. Consumes the vacuum catalog at the given order
  // (loads from disk if present, otherwise runs VacuumDiagramGenerator) and
  // for each vacuum graph G enumerates ways to place n_marks marks on its
  // vertices, dedupes by colored canonical form, and finally computes the
  // shell-indexed embedding count for each unique rooted graph.
  //
  // n_marks must be 1 or 2. The two marks are interchangeable (same color)
  // so the unordered pair {i, j} is canonicalized once.
  class RootedDiagramGenerator {

    public:
    RootedDiagramGenerator(int order, int n_marks, bool bipartite_only = true);

    void generate();
    std::vector<RootedGraph> const &get_unique_rooted_graphs() const { return rooted_graphs; }

    private:
    int order;
    int n_marks;
    bool bipartite_only;

    // Dedup key: (canonical adjacency, canonical marks).
    struct RootedKey {
      std::vector<uint8_t> canonical_adj;
      std::vector<int> canonical_marks;
      bool operator==(RootedKey const &o) const {
        return this->canonical_adj == o.canonical_adj && this->canonical_marks == o.canonical_marks;
      }
    };
    struct RootedKeyHasher {
      size_t operator()(RootedKey const &k) const {
        size_t seed = VectorHasher{}(k.canonical_adj);
        for (int m : k.canonical_marks) seed ^= std::hash<int>{}(m) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
      }
    };

    std::unordered_set<RootedKey, RootedKeyHasher> unique_keys;
    std::vector<RootedGraph> rooted_graphs;

    void try_emit(Graph const &G, std::vector<int> marks);
  };
} // namespace sc_expansion
