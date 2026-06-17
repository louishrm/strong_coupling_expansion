#include "graph.hpp"
#include "canonicalize.hpp"
#include "diagram_common.hpp"
#include <algorithm>
#include <chrono>
#include <queue>
#include <stdexcept>
#include <tuple>

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

    // Colored canonicalization: the color of vertex v equals the *number* of
    // marks placed on v (0 unmarked, 1 single mark, 2 if both marks coincide
    // on v). Distinct color values keep single-mark and double-mark vertices
    // in separate bliss equivalence classes, so {v,v} same-vertex placements
    // do not collapse with {v} single-mark forms — and bliss's automorphism
    // count is |Stab_{Aut(G)}(M)| where M is the unordered mark multiset.
    std::vector<int> color(this->V, 0);
    for (int m : marks_in) color[m] += 1;
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

  // ---------------------------------------------------------------------------
  // DimerRootedGraph
  // ---------------------------------------------------------------------------

  DimerRootedGraph::DimerRootedGraph(Graph const &parent, std::vector<int> marks_in, std::vector<int> sites_in)
     : V(parent.get_V()), order(parent.get_order()) {

    if (marks_in.size() != sites_in.size())
      throw std::invalid_argument("DimerRootedGraph: marks and sites size mismatch");
    if (marks_in.empty() || marks_in.size() > 2)
      throw std::invalid_argument("DimerRootedGraph: marks must have size 1 or 2");
    for (int m : marks_in) {
      if (m < 0 || m >= this->V) throw std::invalid_argument("DimerRootedGraph: mark index out of range");
    }
    for (int s : sites_in) {
      if (s != 0 && s != 1) throw std::invalid_argument("DimerRootedGraph: site must be 0 or 1");
    }

    // Colored canonicalisation. Color of vertex v encodes which marks (if any)
    // sit on v together with their within-dimer-site assignments:
    //   color 0           : unmarked
    //   color 1           : single mark at site 0
    //   color 2           : single mark at site 1
    //   color 3           : two coincident marks, multiset {0,0}
    //   color 4           : two coincident marks, multiset {0,1}
    //   color 5           : two coincident marks, multiset {1,1}
    // The coincident-mark encoding uses 3 + (lo + hi) so the multiset {ms_a,
    // ms_b} canonically maps to a single color regardless of input order.
    //
    // Dimer-inversion symmetry (180° rotation about a dimer center, which
    // swaps sites 0 ↔ 1 globally) is a symmetry of the staggered tiling.
    // We fold it into the canonical key by canonicalising both the input
    // (G, marks, sites) and its site-flipped counterpart (G, marks, 1-sites),
    // and picking the lex-min (canonical_adj, marks, sites) tuple. Sectors
    // (0,0) ↔ (1,1) and (0,1) ↔ (1,0) collapse to a single catalog entry.
    auto adj = parent.get_adjacency_matrix();

    auto canonicalise_with_sites = [&](std::vector<int> const &sites_v) {
      std::vector<int> color(this->V, 0);
      std::vector<std::vector<int>> site_buckets(this->V);
      for (size_t k = 0; k < marks_in.size(); ++k) site_buckets[marks_in[k]].push_back(sites_v[k]);
      for (int v = 0; v < this->V; ++v) {
        auto &sb = site_buckets[v];
        if (sb.empty()) {
          color[v] = 0;
        } else if (sb.size() == 1) {
          color[v] = 1 + sb[0];
        } else {
          int lo   = std::min(sb[0], sb[1]);
          int hi   = std::max(sb[0], sb[1]);
          color[v] = 3 + lo + hi;
        }
      }
      auto result = canonicalize(adj, this->V, color);
      std::vector<std::pair<int, int>> ms_pairs;
      ms_pairs.reserve(marks_in.size());
      for (size_t k = 0; k < marks_in.size(); ++k) {
        ms_pairs.emplace_back(static_cast<int>(result.canonical_permutation[marks_in[k]]), sites_v[k]);
      }
      std::sort(ms_pairs.begin(), ms_pairs.end());
      std::vector<int> out_marks, out_sites;
      out_marks.reserve(ms_pairs.size());
      out_sites.reserve(ms_pairs.size());
      for (auto const &p : ms_pairs) {
        out_marks.push_back(p.first);
        out_sites.push_back(p.second);
      }
      return std::make_tuple(result.canonical_matrix, std::move(out_marks), std::move(out_sites),
                             static_cast<int>(result.automorphism_count));
    };

    std::vector<int> sites_flipped(sites_in.size());
    for (size_t k = 0; k < sites_in.size(); ++k) sites_flipped[k] = 1 - sites_in[k];

    auto cand_a = canonicalise_with_sites(sites_in);
    auto cand_b = canonicalise_with_sites(sites_flipped);

    // Lex-min on (canonical_adj, marks, sites).
    auto key_a = std::tie(std::get<0>(cand_a), std::get<1>(cand_a), std::get<2>(cand_a));
    auto key_b = std::tie(std::get<0>(cand_b), std::get<1>(cand_b), std::get<2>(cand_b));
    auto &chosen = (key_a <= key_b) ? cand_a : cand_b;

    this->canonical_matrix          = std::get<0>(chosen);
    this->adjacency_matrix          = std::get<0>(chosen);
    this->marks                     = std::get<1>(chosen);
    this->sites                     = std::get<2>(chosen);
    this->rooted_automorphism_count = std::get<3>(chosen);

    // Multi-edge factorial correction (same as RootedGraph): parallel edges
    // between fixed vertex pairs are permutable but invisible to bliss.
    int factorial_product = 1;
    for (auto e : this->adjacency_matrix)
      if (e > 1) factorial_product *= sc_expansion::factorial(e);
    this->rooted_symmetry_factor = this->rooted_automorphism_count * factorial_product;
  }

} // namespace sc_expansion
