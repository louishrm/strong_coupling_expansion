#pragma once
#include <vector>
#include <deque>
#include "diagram.hpp"
#include "hubbard_solver.hpp"

namespace sc_expansion {

  template <typename T>
  class FreeEnergyCalculator {
    public:
    FreeEnergyCalculator(Parameters<T> const &params, int order, int override_fm = -1);

    // Construct from pre-built graphs (avoids redundant diagram generation across MPI ranks)
    FreeEnergyCalculator(Parameters<T> const &params, int order, std::vector<Graph> const &prebuilt_graphs);

    T compute_sum_diagrams(std::vector<double> const &taus, bool infinite_U) const;

    std::pair<double, double> compute_infinite_U_coefficient() const;

    void clear_all_caches() const;
    void mark_tau_dirty(int tau_index);
    void mark_all_dirty();

    std::deque<Graph> const &get_graphs() const { return this->graphs; }
    std::deque<VertexType<T>> const &get_vertex_types() const { return this->vertex_types; }
    std::deque<Diagram<T>> const &get_diagrams() const { return this->diagrams; }

    int get_n_diagrams() const { return (int)this->diagrams.size(); }

    private:
    Parameters<T> const &params;
    int order;
    std::deque<Graph> graphs;
    std::deque<Diagram<T>> diagrams;
    std::deque<VertexType<T>> vertex_types; // One per cumulant order (index = cumulant_order - 1)
    HubbardSolver<T> solver;
  };

} // namespace sc_expansion
