#pragma once
#include <vector>
#include <deque>
#include <map>
#include <utility>
#include "../graph.hpp"
#include "../hubbard_solver.hpp"
#include "diagram.hpp"
#include "vertex.hpp"

namespace sc_expansion::atomic {

  // Build the rooted catalog for a density-density measurement at lattice
  // displacement r. Runs DistanceRootedDiagramGenerator and flattens its
  // (V → unique rooted graphs) output into parallel vectors of (Graph, marks)
  // — each graph already carries the rooted symmetry factor and
  // free_multiplicity = 1 via the Graph override constructor. Marks are in
  // the graph's canonical labeling. Intended as the rank-0-only step of an
  // MPI driver, with (graphs_out, marks_out) broadcast to the other ranks.
  void build_rooted_catalog(int order, bool bipartite, std::vector<int> const &r,
                            std::vector<Graph> &graphs_out,
                            std::vector<std::vector<int>> &marks_out);

  // Diagram-list owner + evaluator. Handles the vacuum/free-energy series and
  // the rooted density-density series with the same diagram/solver/vertex-type
  // infrastructure. The rooted constructor enables density-density mode; in
  // that mode is_density_density_mode() returns true and the target d² shell
  // is exposed via get_target_d_sq().
  template <typename T> class SumDiagrams {
    public:
    // Vacuum / free-energy constructors.
    SumDiagrams(Parameters<T> const &params, int order, int override_fm = -1);
    SumDiagrams(Parameters<T> const &params, int order, std::vector<Graph> const &prebuilt_graphs, int override_fm = -1);

    // Rooted density-density constructor (Flavor 2): builds the catalog of
    // rooted topologies whose two marks can embed at displacement r on the
    // bipartite hypercubic lattice, forcing free_multiplicity = 1 and
    // lattice_multiplier = {|r|²: 1} per diagram.
    SumDiagrams(Parameters<T> const &params, int order, std::vector<int> r, int s1, int s2);

    // Rooted density-density constructor with prebuilt catalog. Mirrors the
    // vacuum prebuilt-graphs path so MPI drivers can have rank 0 enumerate
    // the rooted catalog once and broadcast (graphs, marks) to all ranks.
    // graphs[i] must already carry the rooted symmetry factor and
    // free_multiplicity = 1 (i.e. constructed via Graph's override ctor);
    // marks[i] gives the two mark indices in graphs[i]'s labeling.
    SumDiagrams(Parameters<T> const &params, int order,
                std::vector<Graph> const &prebuilt_rooted_graphs,
                std::vector<std::vector<int>> const &prebuilt_marks,
                std::vector<int> const &r, int s1, int s2);

    // Scalar accumulator over the vacuum/free-energy diagram list.
    T free_energy(std::vector<double> const &taus, bool infinite_U) const;

    // Per-d² accumulator over the rooted diagram list. Under the rooted
    // constructor (Flavor 2) the map has exactly one entry keyed by |r|².
    std::map<int, T> density_density(std::vector<double> const &taus, bool infinite_U) const;

    // Returns {abs sum, signed sum} × β^n/n!.
    std::pair<double, double> free_energy_infinite_U_coefficient();

    // One {abs, signed} pair per d² entry, both × β^n/n!.
    std::map<int, std::pair<double, double>> density_density_infinite_U_coefficient();

    void mark_tau_dirty(int tau_index);
    void mark_all_dirty();

    std::deque<Graph> const &get_graphs() const { return this->graphs; }
    std::deque<VertexType<T>> const &get_vertex_types() const { return this->vertex_types; }
    std::deque<Diagram<T>> const &get_diagrams() const { return this->diagrams; }
    HubbardSolver<1, T> const &get_solver() const { return this->solver; }

    int get_n_diagrams() const { return (int)this->diagrams.size(); }

    bool is_density_density_mode() const { return this->target_d_sq >= 0; }
    int get_target_d_sq() const { return this->target_d_sq; }

    private:
    void init_from_graphs(std::vector<Graph> const &source_graphs, int override_fm);
    void init_from_rooted_catalog(std::vector<Graph> const &rooted_graphs,
                                  std::vector<std::vector<int>> const &marks,
                                  std::vector<int> const &r, int s1, int s2);

    // Shared n!-permutation SJT sweep used by both *_infinite_U_coefficient
    // implementations. `per_perm(taus, accum)` is invoked once per permutation
    // after mark_all_dirty(); `accum` is the caller-owned accumulator.
    template <class Accum, class PerPerm> void sjt_sweep(Accum &accum, PerPerm &&per_perm);

    Parameters<T> const &params;
    int order;
    std::deque<Graph> graphs;
    std::deque<Diagram<T>> diagrams;
    std::deque<VertexType<T>> vertex_types; // index = cumulant_order - 1
    HubbardSolver<1, T> solver;

    // Density-density state. target_d_sq < 0 ⇒ vacuum/free-energy mode.
    int target_d_sq = -1;
    std::vector<std::map<int, int>> lattice_multiplier; // per-diagram
  };

} // namespace sc_expansion::atomic
