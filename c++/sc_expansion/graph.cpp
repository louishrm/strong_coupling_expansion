#include "graph.hpp"
#include "canonicalize.hpp"
#include "diagram_common.hpp"
#include <algorithm>
#include <chrono>
#include <queue>
#include <stdexcept>

namespace sc_expansion {

  Graph::Graph(std::vector<uint8_t> adjacency_matrix_, int V_, bool bipartite_only_, bool defer_embedding)
     : adjacency_matrix(adjacency_matrix_), V(V_), canonical_matrix(adjacency_matrix_), bipartite_only(bipartite_only_) {

    // Calculate Order (Total number of lines)
    this->order = 0;
    for (auto val : this->adjacency_matrix) { this->order += val; }

    // Pre-calculate Degrees
    this->degrees.resize(this->V, 0);
    for (int i = 0; i < this->V; i++) { this->degrees[i] = this->get_degree_of_vertex(i); }

    this->check_connectivity();

    this->check_if_bipartite();

    this->free_multiplicity = 0;
    if ((this->connected) && (!this->bipartite_only || this->bipartite)) {
      this->compute_canonical_form();
      // Embedding is the most expensive per-graph step. Generators that
      // produce many candidates per unique graph should pass
      // defer_embedding=true and call compute_free_multiplicity() once per
      // unique graph after deduplication.
      if (!defer_embedding) this->compute_free_multiplicity();
    } else {
      this->symmetry_factor    = 0;
      this->automorphism_count = 0;
    }
  }

  Graph::Graph(std::vector<uint8_t> adjacency_matrix_, int V_, int automorphism_count_, int symmetry_factor_, int free_multiplicity_,
               bool bipartite_only_)
     : adjacency_matrix(adjacency_matrix_),
       V(V_),
       canonical_matrix(adjacency_matrix_),
       connected(true),
       bipartite(V_ % 2 == 0),
       bipartite_only(bipartite_only_),
       automorphism_count(automorphism_count_),
       symmetry_factor(symmetry_factor_),
       free_multiplicity(free_multiplicity_) {

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

      // Self-loops are on-site density insertions — they don't break the
      // lattice-bipartite property, so ignore diagonal entries here.
      if (neighbor == vertex) continue;

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

  // --- 2. Canonicalization (bliss) ---
  void Graph::compute_canonical_form() {
    // Uncolored canonicalization: all vertices share color 0.
    std::vector<int> color(this->V, 0);
    auto result              = canonicalize(this->adjacency_matrix, this->V, color);
    this->canonical_matrix   = std::move(result.canonical_matrix);
    this->automorphism_count = static_cast<int>(result.automorphism_count);

    // Multi-edge correction to the symmetry factor: each (i,j) entry with
    // multiplicity m > 1 contributes m! permutations of parallel lines that
    // bliss does NOT count (it sees only the multiset of edges).
    int factorial_product = 1;
    for (const auto &entry : this->adjacency_matrix) {
      if (entry > 1) factorial_product *= sc_expansion::factorial(entry);
    }

    this->symmetry_factor = this->automorphism_count * factorial_product;
  }

  // Forward declaration — implementation lives in diagram.cpp
  int compute_lattice_free_multiplicity(Graph const &graph);

  namespace {
    int64_t embedding_total_ns = 0;
  }

  void reset_embedding_timer() { embedding_total_ns = 0; }
  double get_embedding_time_seconds() { return static_cast<double>(embedding_total_ns) * 1e-9; }

  void Graph::compute_free_multiplicity() {
    auto t0                 = std::chrono::steady_clock::now();
    this->free_multiplicity = compute_lattice_free_multiplicity(*this);
    auto t1                 = std::chrono::steady_clock::now();
    embedding_total_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
  }

  // ---------------------------------------------------------------------------
  // RootedGraph
  // ---------------------------------------------------------------------------

  RootedGraph::RootedGraph(Graph const &parent, std::vector<int> marks_in, bool defer_embedding)
     : V(parent.get_V()), order(parent.get_order()), bipartite_only(parent.get_bipartite_only()) {

    if (marks_in.empty() || marks_in.size() > 2)
      throw std::invalid_argument("RootedGraph: marks must have size 1 or 2");
    for (int m : marks_in) {
      if (m < 0 || m >= this->V) throw std::invalid_argument("RootedGraph: mark index out of range");
    }

    // Colored canonicalization: marks share color 1, unmarked share color 0.
    // The two marks are interchangeable (same color), so bliss's automorphism
    // count is exactly |Stab_{Aut(G)}({m0, m1})| as an unordered set.
    std::vector<int> color(this->V, 0);
    for (int m : marks_in) color[m] = 1;
    auto adj    = parent.get_adjacency_matrix();
    auto result = canonicalize(adj, this->V, color);

    this->canonical_matrix = result.canonical_matrix;
    this->adjacency_matrix = result.canonical_matrix;

    // Map marks into canonical labeling and sort so the pair is order-free.
    this->marks.reserve(marks_in.size());
    for (int m : marks_in) this->marks.push_back(static_cast<int>(result.canonical_permutation[m]));
    std::sort(this->marks.begin(), this->marks.end());

    this->rooted_automorphism_count = static_cast<int>(result.automorphism_count);

    // Multi-edge correction to the rooted symmetry factor: parallel edges
    // between fixed vertex pairs can still be permuted, and that factor is
    // not seen by bliss (each parallel edge already lives on its own aux
    // vertex in the encoded graph; aux permutations are divided out inside
    // canonicalize). Apply the m! correction once per (i,j) entry.
    int factorial_product = 1;
    for (auto e : this->adjacency_matrix)
      if (e > 1) factorial_product *= sc_expansion::factorial(e);
    this->rooted_symmetry_factor = this->rooted_automorphism_count * factorial_product;

    if (!defer_embedding) this->compute_shell_multiplicity();
  }

  void RootedGraph::compute_shell_multiplicity() {
    // Wrap the (already canonical) adjacency as a Graph so we can reuse the
    // embedding recursion. The override constructor below skips the connect-
    // ivity / bipartite / canonicalize work that we don't need to redo.
    Graph wrapper(this->adjacency_matrix, this->V, /*aut*/ 0, /*sym*/ 0, /*fm*/ 0, this->bipartite_only);
    int m0 = this->marks[0];
    int m1 = this->marks.size() == 2 ? this->marks[1] : -1;

    auto t0                  = std::chrono::steady_clock::now();
    this->shell_multiplicity = compute_rooted_shell_multiplicity(wrapper, m0, m1);
    auto t1                  = std::chrono::steady_clock::now();
    embedding_total_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
  }

} // namespace sc_expansion
