#pragma once
#include <vector>
#include <queue>
#include <unordered_map>
#include "combinatorics.hpp"

namespace sc_expansion {

  class Graph {

    public:
    Graph(std::vector<uint8_t> adjacency_matrix, int V, bool bipartite_only = true);
    Graph(std::vector<uint8_t> adjacency_matrix, int V, int automorphism_count, int symmetry_factor, int free_multiplicity,
          bool bipartite_only = true);
    uint8_t operator()(int i, int j) const;

    int get_V() const { return this->V; }
    int get_order() const { return this->order; }

    double get_symmetry_factor() const { return (double)this->symmetry_factor; }
    bool get_connectivity() const { return this->connected; }
    std::vector<uint8_t> get_canonical_form() const { return this->canonical_matrix; }
    double get_free_multiplicity() const { return (double)this->free_multiplicity; }
    bool get_bipartite() const { return this->bipartite; }
    bool get_bipartite_only() const { return this->bipartite_only; }
    int get_automorphism_count() const { return this->automorphism_count; }

    // Compute and store the lattice embedding count for this graph. The
    // constructor leaves free_multiplicity = 0 by default because embedding
    // is the most expensive per-graph step; callers (e.g. the diagram
    // generator) should defer it until after deduplication so it runs once
    // per unique canonical form rather than once per candidate.
    void compute_free_multiplicity();

    private:
    std::vector<uint8_t> adjacency_matrix;
    int V; //number of vertices
    int order;
    std::vector<uint8_t> canonical_matrix;
    std::vector<int> degrees;
    bool connected;
    bool bipartite;
    bool bipartite_only;
    int automorphism_count;
    int symmetry_factor;
    int free_multiplicity;

    int get_degree_of_vertex(int vertex) const;
    void check_connectivity();
    bool check_bipartite_dfs(int vertex, std::vector<int> &colors) const;
    void check_if_bipartite();
    void compute_canonical_form();
  };

  // Profiling: cumulative time spent inside compute_lattice_free_multiplicity
  // across all Graph constructions on the current thread. Reset before the
  // workload of interest, read after.
  void reset_embedding_timer();
  double get_embedding_time_seconds();
} // namespace sc_expansion
