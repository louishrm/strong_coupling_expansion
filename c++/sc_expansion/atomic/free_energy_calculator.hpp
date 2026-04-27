#pragma once
#include <vector>
#include <deque>
#include "../graph.hpp"
#include "../hubbard_solver.hpp"
#include "diagram.hpp"
#include "vertex.hpp"

namespace sc_expansion::atomic {

  template <typename T> class FreeEnergyCalculator {
    public:
    FreeEnergyCalculator(Parameters<T> const &params, int order, int override_fm = -1);

    // Construct from pre-built graphs (avoids redundant diagram generation across MPI ranks).
    FreeEnergyCalculator(Parameters<T> const &params, int order, std::vector<Graph> const &prebuilt_graphs, int override_fm = -1);

    T compute_sum_diagrams(std::vector<double> const &taus, bool infinite_U) const;

    std::pair<double, double> compute_infinite_U_coefficient();

    void mark_tau_dirty(int tau_index);
    void mark_all_dirty();

    std::deque<Graph> const &get_graphs() const { return this->graphs; }
    std::deque<VertexType<T>> const &get_vertex_types() const { return this->vertex_types; }
    std::deque<Diagram<T>> const &get_diagrams() const { return this->diagrams; }
    HubbardSolver<1, T> const &get_solver() const { return this->solver; }

    int get_n_diagrams() const { return (int)this->diagrams.size(); }

    private:
    void init_from_graphs(std::vector<Graph> const &source_graphs, int override_fm);

    Parameters<T> const &params;
    int order;
    std::deque<Graph> graphs;
    std::deque<Diagram<T>> diagrams;
    std::deque<VertexType<T>> vertex_types; // index = cumulant_order - 1
    HubbardSolver<1, T> solver;
  };

} // namespace sc_expansion::atomic
