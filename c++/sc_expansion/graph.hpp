#pragma once
#include <map>
#include <queue>
#include <unordered_map>
#include <vector>
#include "combinatorics.hpp"

namespace sc_expansion {

  class Graph {

    public:
    // Standard constructor: canonicalizes and computes the lattice embedding
    // count. Set defer_embedding=true to skip the embedding step (useful for
    // generators that want to dedup first and only embed the unique graphs);
    // the caller is then responsible for invoking compute_free_multiplicity().
    Graph(std::vector<uint8_t> adjacency_matrix, int V, bool bipartite_only = true, bool defer_embedding = false);
    Graph(std::vector<uint8_t> adjacency_matrix, int V, int automorphism_count, int symmetry_factor, int free_multiplicity,
          bool bipartite_only = true);
    uint8_t operator()(int i, int j) const;

    int get_V() const { return this->V; }
    int get_order() const { return this->order; }

    double get_symmetry_factor() const { return (double)this->symmetry_factor; }
    bool get_connectivity() const { return this->connected; }
    std::vector<uint8_t> get_canonical_form() const { return this->canonical_matrix; }
    std::vector<uint8_t> get_adjacency_matrix() const { return this->adjacency_matrix; }
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

  // ---------------------------------------------------------------------------
  // Rooted graph: a Graph with 1 or 2 "marked" vertices (density insertions
  // for the ⟨S^z(r) S^z(0)⟩ correlator). Topology lives alongside the vacuum
  // Graph because the underlying adjacency, bipartite check, and embedding
  // recursion are shared; only the canonicalization coloring, the symmetry
  // factor semantics, and the embedding output differ:
  //
  //   - canonicalization: marked vertices share a color distinct from the
  //     unmarked ones, so bliss returns |Stab_{Aut(G)}(M)| directly (no
  //     orbit computation needed).
  //   - embedding: mark[0] anchors at the origin; for n_marks==2, the output
  //     is shell-indexed (map d² → count). For n_marks==1 the map has a
  //     single entry at d²=0.
  //
  // RootedGraph never re-runs connectivity / bipartite checks — it is built
  // from a Graph that has already passed them.
  class RootedGraph {

    public:
    // Construct from a parent (vacuum) Graph and the indices of the marked
    // vertices in the parent's labeling. If defer_embedding=true, the shell
    // multiplicity is left empty and compute_shell_multiplicity() must be
    // called explicitly (used by the generator to embed only the unique
    // rooted forms).
    RootedGraph(Graph const &parent, std::vector<int> marks, bool defer_embedding = false);

    int get_V() const { return this->V; }
    int get_order() const { return this->order; }
    int get_n_marks() const { return static_cast<int>(this->marks.size()); }

    // Marks expressed in the *canonical* labeling of this rooted graph
    // (sorted ascending; for n_marks==2 the pair is interchangeable so this
    // is a canonical representative of the unordered pair).
    std::vector<int> const &get_marks() const { return this->marks; }

    std::vector<uint8_t> get_canonical_form() const { return this->canonical_matrix; }
    double get_rooted_symmetry_factor() const { return (double)this->rooted_symmetry_factor; }
    int get_rooted_automorphism_count() const { return this->rooted_automorphism_count; }
    bool get_bipartite_only() const { return this->bipartite_only; }

    // d² → embedding count. Empty until compute_shell_multiplicity() runs.
    std::map<int, int> const &get_shell_multiplicity() const { return this->shell_multiplicity; }

    void compute_shell_multiplicity();

    private:
    std::vector<uint8_t> adjacency_matrix; // canonical labeling
    std::vector<uint8_t> canonical_matrix; // == adjacency_matrix here
    int V;
    int order;
    std::vector<int> marks; // in canonical labeling
    bool bipartite_only;
    int rooted_automorphism_count;
    int rooted_symmetry_factor;
    std::map<int, int> shell_multiplicity;
  };

  // ---------------------------------------------------------------------------
  // Rooted graph for the staggered-dimer expansion. Like RootedGraph, but each
  // mark also carries a within-dimer-site index (0 or 1). The dimer cumulant
  // at each vertex is a 2-site Hubbard cumulant, and a density-density
  // insertion's physical position is (dimer position, within-dimer site) — so
  // the catalog must distinguish (G, marks, sites) tuples, not just (G, marks).
  //
  // Canonicalisation uses vertex colors that encode the within-dimer-site
  // multiset at each marked vertex, so physically distinct site assignments
  // produce distinct canonical forms. See dimer-rooted-constraints memory.
  //
  // Embedding is not performed here (no shell-multiplicity analog yet): the
  // staggered-tiling embedding count at a target physical displacement r is
  // the consumer's job, mirroring the atomic side's count_lattice_embeddings.
  class DimerRootedGraph {

    public:
    // marks_in: vertex indices (parent labeling) hosting marks. n_marks = 1 or 2.
    //   Coincident marks (both on one vertex) are allowed.
    // sites_in: within-dimer site (0 or 1) per mark, parallel to marks_in.
    DimerRootedGraph(Graph const &parent, std::vector<int> marks_in, std::vector<int> sites_in);

    int get_V() const { return this->V; }
    int get_order() const { return this->order; }
    int get_n_marks() const { return static_cast<int>(this->marks.size()); }

    // Marks in this graph's canonical labeling. For coincident marks both
    // entries are the same vertex index.
    std::vector<int> const &get_marks() const { return this->marks; }

    // Within-dimer sites, one per mark, in the same order as get_marks().
    // For non-coincident marks the order matches the canonical mark order
    // (ascending vertex index); for coincident marks the sites pair is
    // lex-sorted so {0,1} and {1,0} both serialise as (0,1).
    std::vector<int> const &get_sites() const { return this->sites; }

    std::vector<uint8_t> get_canonical_form() const { return this->canonical_matrix; }
    double get_rooted_symmetry_factor() const { return (double)this->rooted_symmetry_factor; }
    int get_rooted_automorphism_count() const { return this->rooted_automorphism_count; }
    // Dimer expansion always works with the full (non-bipartite) graph set —
    // the staggered superlattice is triangular. Exposed for symmetry with
    // RootedGraph; always false here.
    bool get_bipartite_only() const { return false; }

    private:
    std::vector<uint8_t> adjacency_matrix; // canonical labeling
    std::vector<uint8_t> canonical_matrix; // == adjacency_matrix here
    int V;
    int order;
    std::vector<int> marks; // canonical labeling
    std::vector<int> sites; // parallel to marks
    int rooted_automorphism_count;
    int rooted_symmetry_factor;
  };

} // namespace sc_expansion
