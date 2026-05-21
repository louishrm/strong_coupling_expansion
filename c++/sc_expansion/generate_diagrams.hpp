#pragma once
#include <map>
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

    // Generate vacuum graphs at the constructor-given expansion order but
    // restricted to exactly V vertices. Clears any previously generated
    // graphs first; the caller reads them via get_unique_graphs() after each
    // call. When V == order, short-circuits to the n-cycle directly (the only
    // connected vacuum topology with V == order, and the most expensive to
    // canonicalize otherwise).
    void generate_at_vertex_count(int V);

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

    // Emit the n-cycle (V == order fast path). Used by both generate() and
    // generate_at_vertex_count(order). Pushes onto graphs and unique_adjmats;
    // no canonicalization sweep, no embedding sweep.
    void emit_n_cycle();
  };

  // Dedup key for rooted diagrams, shared by RootedDiagramGenerator and
  // DistanceRootedDiagramGenerator.
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

    std::unordered_set<RootedKey, RootedKeyHasher> unique_keys;
    std::vector<RootedGraph> rooted_graphs;

    void try_emit(Graph const &G, std::vector<int> marks);
  };

  // Distance-targeted rooted diagram generator. Parameterized by a lattice
  // displacement vector r and an expansion order n; emits rooted topologies
  // whose two marks can embed at displacement r on the lattice, restricted to
  // graphs with V in [2*|r|_1, n]. See distance_rooted_diagrams.md for the
  // rationale (high-d, fixed-n queries where the full unconstrained catalog
  // is overwhelmingly noise).
  //
  // Embedding count is NOT computed here: per atomic/diagram.hpp, the
  // r-indexed lattice multiplier is the measurement layer's job. The output
  // is a map V -> unique rooted topologies (graph + canonical marks + rooted
  // symmetry factor), nothing more.
  class DistanceRootedDiagramGenerator {

    public:
    // r: lattice displacement (e.g. {1,0} 1st nn, {1,1} 2nd nn, {2,0} 3rd nn).
    // n: truncation expansion order for this call. Throws if n < 2 * |r|_1.
    DistanceRootedDiagramGenerator(std::vector<int> r, int n, bool bipartite_only = true);

    void generate();

    // Rooted graphs keyed by vertex count V (only V in [2*|r|_1, n] populated).
    std::map<int, std::vector<RootedGraph>> const &get_rooted_graphs() const { return rooted_graphs_by_V; }

    private:
    std::vector<int> r;
    int n;
    int d; // |r|_1
    bool bipartite_only;

    std::map<int, std::vector<RootedGraph>> rooted_graphs_by_V;

    // All-pairs shortest path on a vacuum graph (BFS from each vertex,
    // treating any nonzero adjacency entry as a unit edge). Returns a V*V
    // row-major distance matrix; unreachable entries are -1 (won't occur on
    // connected vacuum graphs, but the check is cheap).
    static std::vector<int> bfs_all_pairs(Graph const &G);

    // Try to canonicalize-and-emit a (G, {i,j}) candidate into bucket V.
    void try_emit(Graph const &G, std::vector<int> marks, int V,
                  std::unordered_set<RootedKey, RootedKeyHasher> &unique_keys);
  };
} // namespace sc_expansion
