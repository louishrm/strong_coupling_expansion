#include "graph.hpp"
#include <numeric>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace sc_expansion {

  Graph::Graph(std::vector<uint8_t> adjacency_matrix_, int V_) : adjacency_matrix(adjacency_matrix_), V(V_), canonical_matrix(adjacency_matrix_) {

    // Calculate Order (Total number of lines)
    this->order = 0;
    for (auto val : this->adjacency_matrix) { this->order += val; }

    this->check_connectivity();

    if (this->connected) {
      this->compute_canonical_form();
      this->compute_free_multiplicity();
    } else {
      this->symmetry_factor    = 0;
      this->free_multiplicity  = 0;
      this->automorphism_count = 0;
    }
  }

  uint8_t Graph::operator()(int i, int j) const { return this->adjacency_matrix[i * this->V + j]; }

  // --- 1. Connectivity Check (BFS) ---
  void Graph::check_connectivity() {
    if (this->V == 0) {
      this->connected = false;
      return;
    }

    std::vector<bool> visited(this->V, false);
    std::queue<int> q;

    q.push(0);
    visited[0]        = true;
    int visited_count = 0; // Count pops to be safe

    while (!q.empty()) {
      int vertex = q.front();
      q.pop();
      visited_count++;

      for (int neighbor = 0; neighbor < this->V; neighbor++) {
        // Check undirected connection
        bool is_connected = ((*this)(vertex, neighbor) > 0) || ((*this)(neighbor, vertex) > 0);

        if (is_connected && !visited[neighbor]) {
          visited[neighbor] = true;
          q.push(neighbor);
        }
      }
    }
    this->connected = (visited_count == this->V);
  }

  // --- 2. Canonicalization (Min-Lex + Symmetry) ---
  void Graph::compute_canonical_form() {
    this->canonical_matrix = this->adjacency_matrix;
    int auto_count         = 0;

    std::vector<int> p(this->V);
    std::iota(p.begin(), p.end(), 0); // {0, 1, 2, ...}

    std::vector<uint8_t> candidate(this->V * this->V); // Allocate once

    do {
      // Apply permutation: M'[i][j] = M[p[i]][p[j]]
      // We map the row p[i] to position i
      for (int i = 0; i < this->V; i++) {
        for (int j = 0; j < this->V; j++) { candidate[i * this->V + j] = (*this)(p[i], p[j]); }
      }

      // Min-Lex Comparison
      if (candidate < this->canonical_matrix) {
        this->canonical_matrix = candidate;
        auto_count             = 1;
      } else if (candidate == this->canonical_matrix) {
        auto_count++;
      }
    } while (std::next_permutation(p.begin(), p.end()));

    this->automorphism_count = auto_count;

    // Compute total symmetry factor g(D)
    int factorial_product = 1;
    for (const auto &entry : this->adjacency_matrix) {
      if (entry > 1) factorial_product *= sc_expansion::factorial(entry);
    }

    this->symmetry_factor = auto_count * factorial_product;
  }

  void Graph::compute_free_multiplicity() {

    std::vector<Point> coords(this->V);
    std::vector<bool> placed(this->V, false);

    // Fix Vertex 0 at origin
    coords[0] = Point(0, 0);
    placed[0] = true;
    // Start recursion with 1 vertex placed
    this->free_multiplicity = (int)solve_embedding_recursive(1, placed, coords);
  }

  long Graph::solve_embedding_recursive(int placed_count, std::vector<bool> &placed, std::vector<Point> &coords) {
    // Base Case: All vertices placed
    if (placed_count == this->V) return 1;

    // A. Find a target node and an anchor
    // We look for any unplaced node connected to a placed node
    int anchor      = -1;
    int target_node = -1;

    for (int candidate = 0; candidate < this->V; ++candidate) {
      if (!placed[candidate]) {
        for (int p = 0; p < this->V; ++p) {
          if (placed[p]) {
            uint8_t val = this->canonical_matrix[p * this->V + candidate] + this->canonical_matrix[candidate * this->V + p];
            if (val > 0) {
              target_node = candidate;
              anchor      = p;
              goto found_target;
            }
          }
        }
      }
    }
  found_target:;

    if (target_node == -1) return 0; // Should not happen for connected graphs

    long count = 0;

    // B. Try 4 Lattice Directions
    const int dx[] = {1, -1, 0, 0};
    const int dy[] = {0, 0, 1, -1};

    Point anchor_pos = coords[anchor];

    for (int dir = 0; dir < 4; ++dir) {
      Point candidate_pos = Point(anchor_pos.x + dx[dir], anchor_pos.y + dy[dir]);

      // C. Check Consistency with ALL placed neighbors
      bool valid = true;
      for (int i = 0; i < this->V; ++i) {
        if (placed[i]) {
          // Check if target_node is connected to i
          uint8_t links = this->canonical_matrix[target_node * this->V + i] + this->canonical_matrix[i * this->V + target_node];

          if (links > 0) {
            // Must be exactly dist 1
            int dist = std::abs(candidate_pos.x - coords[i].x) + std::abs(candidate_pos.y - coords[i].y);
            if (dist != 1) {
              valid = false;
              break;
            }
          }
        }
      }

      if (valid) {
        coords[target_node] = candidate_pos;
        placed[target_node] = true;

        count += solve_embedding_recursive(placed_count + 1, placed, coords);

        placed[target_node] = false;
      }
    }
    return count;
  }

} // namespace sc_expansion