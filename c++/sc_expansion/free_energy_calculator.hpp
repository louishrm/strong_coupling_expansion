#pragma once
#include <vector>
#include <deque>
#include "diagram.hpp"
#include "hubbard_solver.hpp"

namespace sc_expansion {

  template <typename T>
  class FreeEnergyCalculator {
    public:
    FreeEnergyCalculator(Parameters<T> const &params, int order, int override_fm = -1, bool allow_self_loops = false);

    // Construct from pre-built graphs (avoids redundant diagram generation across MPI ranks)
    FreeEnergyCalculator(Parameters<T> const &params, int order, std::vector<Graph> const &prebuilt_graphs, int override_fm = -1);

    T compute_sum_diagrams(std::vector<double> const &taus, bool infinite_U) const;

    std::pair<double, double> compute_infinite_U_coefficient();

    void mark_tau_dirty(int tau_index);
    void mark_all_dirty();

    std::deque<Graph> const &get_graphs() const { return this->graphs; }
    std::deque<VertexType<T>> const &get_vertex_types() const { return this->vertex_types; }
    std::deque<Diagram<T>> const &get_diagrams() const { return this->diagrams; }
    HubbardSolver<T> const &get_solver() const { return this->solver; }

    int get_n_diagrams() const { return (int)this->diagrams.size(); }

    // Self-loop diagram (V=1, n self-loops = order): tau-independent, precomputed exactly.
    bool has_self_loop_diagram() const { return this->self_loop_present; }
    T get_self_loop_contribution_finite() const;
    T get_self_loop_contribution_infinite() const;

    private:
    void init_from_graphs(std::vector<Graph> const &source_graphs, int override_fm);
    void precompute_self_loop_diagram();

    Parameters<T> const &params;
    int order;
    std::deque<Graph> graphs;
    std::deque<Diagram<T>> diagrams;
    std::deque<VertexType<T>> vertex_types; // One per cumulant order (index = cumulant_order - 1)
    HubbardSolver<T> solver;

    bool self_loop_present     = false;
    T    self_loop_val_finite  = T(0.0);
    T    self_loop_val_infinite = T(0.0);
  };

} // namespace sc_expansion
