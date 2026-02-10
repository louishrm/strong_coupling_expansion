#include "graph.hpp"

namespace sc_expansion {

  Graph::Graph(std::vector<uint8_t> adjacency_matrix_, int V_) : adjacency_matrix(adjacency_matrix_), V(V_), canonical_matrix(adjacency_matrix_) {
    this->check_connectivity();

    if (this->connected) { this->compute_canonical_form(); }
  }

  uint8_t Graph::operator()(int i, int j) const { return this->adjacency_matrix[i * this->V + j]; }

  void Graph::check_connectivity() {

    if (this->V == 0) {
      this->connected = false;
      return;
    }

    //Use BFS to check connectivity
    std::vector<bool> visited(this->V, false);
    std::queue<int> q;
    q.push(0);
    visited[0] = true;

    int visited_count = 1;
    while (!q.empty()) {
      int vertex = q.front(); //get the next vertex to visit
      q.pop();                //remove it from the queue
      //check neighbors
      for (int neighbor = 0; neighbor < this->V; neighbor++) {
        //check for connection in either direction
        bool is_neighbor_connected = ((*this)(vertex, neighbor) > 0) || ((*this)(neighbor, vertex) > 0);
        if (is_neighbor_connected && !visited[neighbor]) { //if there is a connection and neighbor not visited
          visited[neighbor] = true;                        //mark as visited
          q.push(neighbor);                                //add to queue for further exploration
          visited_count++;                                 //increment visited count
        }
      }
    }
    this->connected = (visited_count == this->V);
  }

  void Graph::compute_canonical_form() {
    // 1. Initialize "Best" as the current matrix (assuming it's the best so far)
    this->canonical_matrix = this->adjacency_matrix;

    // Initialize symmetry count (starts at 0, or 1 if we count Identity immediately)
    // Better strategy: Count strictly inside the loop.
    int auto_count = 0;

    auto permutations = sc_expansion::generate_permutations(this->V);

    for (const auto &perm : permutations) {
      std::vector<uint8_t> candidate(this->V * this->V);

      for (int i = 0; i < this->V; i++) {
        for (int j = 0; j < this->V; j++) {
          uint8_t val                = (*this)(perm[i], perm[j]);
          candidate[i * this->V + j] = val;
        }
      }
      // If candidate < canonical, we found a "smaller" representation
      if (candidate < this->canonical_matrix) {
        this->canonical_matrix = candidate; // Update unique ID
        auto_count             = 1;         // Reset count (new best shape found)
      } else if (candidate == this->canonical_matrix) {
        auto_count++; // Another permutation produces the best shape
      }
    }
    this->automorphism_count = auto_count;
    int factorial_product    = 1;
    for (const auto &entry : this->adjacency_matrix) { factorial_product *= sc_expansion::factorial(entry); }

    this->symmetry_factor = auto_count * factorial_product;
  }
} // namespace sc_expansion
