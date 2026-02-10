#pragma once
#include <vector>
#include <queue>
#include "combinatorics.hpp"

namespace sc_expansion {

  class Graph {

    public:
    Graph(std::vector<uint8_t> adjacency_matrix, int V);
    uint8_t operator()(int i, int j) const;

    double get_symmetry_factor() const { return (double)this->symmetry_factor; }
    bool get_connectivity() const { return this->connected; }
    std::vector<uint8_t> get_canonical_form() const { return this->canonical_matrix; }

    private:
    std::vector<uint8_t> adjacency_matrix;
    int V; //number of vertices
    std::vector<uint8_t> canonical_matrix;
    bool connected;
    int automorphism_count;
    int symmetry_factor;

    void check_connectivity();
    void compute_canonical_form();
  };
} // namespace sc_expansion