#include "graph.hpp"
#include <numeric>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace sc_expansion {

  Graph::Graph(std::vector<uint8_t> adjacency_matrix_, int V_, bool bipartite_only_) : adjacency_matrix(adjacency_matrix_), V(V_), canonical_matrix(adjacency_matrix_), bipartite_only(bipartite_only_) {

    // Calculate Order (Total number of lines)
    this->order = 0;
    for (auto val : this->adjacency_matrix) { this->order += val; }

    // Pre-calculate Degrees
    this->degrees.resize(this->V, 0);
    for (int i = 0; i < this->V; i++) { this->degrees[i] = this->get_degree_of_vertex(i); }

    this->check_connectivity();

    this->check_if_bipartite();

    if ((this->connected) && (!this->bipartite_only || this->bipartite)) {
      this->compute_canonical_form();
      this->compute_free_multiplicity();
    } else {
      this->symmetry_factor    = 0;
      this->free_multiplicity  = 0;
      this->automorphism_count = 0;
    }
  }

  Graph::Graph(std::vector<uint8_t> adjacency_matrix_, int V_, int automorphism_count_, int symmetry_factor_, int free_multiplicity_, bool bipartite_only_)
    : adjacency_matrix(adjacency_matrix_), V(V_), canonical_matrix(adjacency_matrix_), connected(true), bipartite(true), bipartite_only(bipartite_only_),
      automorphism_count(automorphism_count_), symmetry_factor(symmetry_factor_), free_multiplicity(free_multiplicity_) {

    // Calculate Order (Total number of lines)
    this->order = 0;
    for (auto val : this->adjacency_matrix) { this->order += val; }

    // Pre-calculate Degrees
    this->degrees.resize(this->V, 0);
    for (int i = 0; i < this->V; i++) { this->degrees[i] = this->get_degree_of_vertex(i); }
  }

  int Graph::get_degree_of_vertex(int vertex) const {
    //degree = number of outgoing+incoming lines
    int degree = 0;
    for (int j = 0; j < this->V; j++) { degree += (*this)(vertex, j) + this->adjacency_matrix[j * this->V + vertex]; }
    return degree;
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
        bool is_connected = ((*this)(vertex, neighbor) > 0) || (this->adjacency_matrix[neighbor * this->V + vertex] > 0);

        if (is_connected && !visited[neighbor]) {
          visited[neighbor] = true;
          q.push(neighbor);
        }
      }
    }
    this->connected = (visited_count == this->V);
  }

  bool Graph::check_bipartite_dfs(int vertex, std::vector<int> &colors) const {

    for (int neighbor = 0; neighbor < this->V; neighbor++) {

      bool is_connected = ((*this)(vertex, neighbor) > 0) || (this->adjacency_matrix[neighbor * this->V + vertex] > 0);

      if (is_connected) {
        if (colors[neighbor] == 0) {
          colors[neighbor] = -colors[vertex];
          if (!check_bipartite_dfs(neighbor, colors)) return false;
        }

        else if (colors[neighbor] == colors[vertex]) {
          return false;
        }
      }
    }
    return true;
  }

  void Graph::check_if_bipartite() {
    std::vector<int> colors(this->V, 0); // 0 = uncolored, 1 = color A, -1 = color B
    for (int vertex = 0; vertex < this->V; vertex++) {
      if (colors[vertex] == 0) {
        colors[vertex] = 1;
        if (!check_bipartite_dfs(vertex, colors)) {
          this->bipartite = false;
          return;
        }
      }
    }
    this->bipartite = true;
  }

  // --- 2. Canonicalization (Min-Lex + Symmetry) ---
  void Graph::compute_canonical_form() {

    // 1. Establish a canonical starting labeling by sorting vertices by degree (non-decreasing)
    std::vector<int> p_sort(this->V);
    std::iota(p_sort.begin(), p_sort.end(), 0);
    std::stable_sort(p_sort.begin(), p_sort.end(), [this](int a, int b) { return this->degrees[a] < this->degrees[b]; });

    std::vector<uint8_t> sorted_matrix(this->V * this->V);
    std::vector<int> sorted_degrees(this->V);
    for (int i = 0; i < this->V; i++) {
      sorted_degrees[i] = this->degrees[p_sort[i]];
      for (int j = 0; j < this->V; j++) { sorted_matrix[i * this->V + j] = (*this)(p_sort[i], p_sort[j]); }
    }

    this->canonical_matrix = sorted_matrix;
    int auto_count         = 0;

    std::vector<int> p(this->V);
    std::iota(p.begin(), p.end(), 0);

    std::vector<uint8_t> candidate(this->V * this->V);

    do {
      // 2. Only apply permutations that preserve the sorted degree profile
      bool can_swap = true;
      for (int i = 0; i < this->V; i++) {
        if (sorted_degrees[p[i]] != sorted_degrees[i]) {
          can_swap = false;
          break;
        }
      }
      if (!can_swap) continue;

      for (int i = 0; i < this->V; i++) {
        for (int j = 0; j < this->V; j++) { candidate[i * this->V + j] = sorted_matrix[p[i] * this->V + p[j]]; }
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

    // B. Lattice Directions
    // Square lattice neighbors: (1,0), (-1,0), (0,1), (0,-1)
    // Triangular lattice neighbors: (1,0), (-1,0), (-1,1), (0,1), (0,-1), (1,-1)
    std::vector<int> dx, dy;
    if (this->bipartite_only) {
      dx = {1, -1, 0, 0};
      dy = {0, 0, 1, -1};
    } else {
      dx = {1, -1, -1, 0, 0, 1};
      dy = {0, 0, 1, 1, -1, -1};
    }

    Point anchor_pos = coords[anchor];

    auto is_neighbor = [this](Point p1, Point p2) {
      int dx = p1.x - p2.x;
      int dy = p1.y - p2.y;
      if (this->bipartite_only) {
        return std::abs(dx) + std::abs(dy) == 1;
      } else {
        // Triangular lattice neighbors: dist 1 in square sense OR (-1,1) or (1,-1)
        if (std::abs(dx) + std::abs(dy) == 1) return true;
        if (dx == -1 && dy == 1) return true;
        if (dx == 1 && dy == -1) return true;
        return false;
      }
    };

    for (size_t dir = 0; dir < dx.size(); ++dir) {
      Point candidate_pos = Point(anchor_pos.x + dx[dir], anchor_pos.y + dy[dir]);

      // C. Check Consistency with ALL placed neighbors
      bool valid = true;
      for (int i = 0; i < this->V; ++i) {
        if (placed[i]) {
          // Check if target_node is connected to i
          uint8_t links = this->canonical_matrix[target_node * this->V + i] + this->canonical_matrix[i * this->V + target_node];

          if (links > 0) {
            // Must be exactly dist 1
            if (!is_neighbor(candidate_pos, coords[i])) {
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
