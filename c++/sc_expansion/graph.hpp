#pragma once
#include <vector>
#include <queue>
#include <unordered_map>
#include "combinatorics.hpp"

namespace sc_expansion {

  class Graph {

    public:
    Graph(std::vector<uint8_t> adjacency_matrix, int V, bool bipartite_only = true);
    Graph(std::vector<uint8_t> adjacency_matrix, int V, int automorphism_count, int symmetry_factor, int free_multiplicity, bool bipartite_only = true);
    uint8_t operator()(int i, int j) const;

    int get_V() const { return this->V; }
    int get_order() const { return this->order; }

    double get_symmetry_factor() const { return (double)this->symmetry_factor; }
    bool get_connectivity() const { return this->connected; }
    std::vector<uint8_t> get_canonical_form() const { return this->canonical_matrix; }
    double get_free_multiplicity() const { return (double)this->free_multiplicity; }
    bool get_bipartite() const { return this->bipartite; }

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
    void compute_free_multiplicity();

    // Helper structs and methods for free multiplicity
    struct Point {
      int x;
      int y;
      Point() : x(0), y(0) {}
      Point(int x_, int y_) : x(x_), y(y_) {}
    };

    long solve_embedding_recursive(int placed_count, std::vector<bool> &placed, std::vector<Point> &coords);
  };
} // namespace sc_expansion
