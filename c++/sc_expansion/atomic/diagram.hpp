#pragma once

#include <vector>
#include "../graph.hpp"
#include "../hubbard_solver.hpp"
#include "../diagram_common.hpp"
#include "vertex.hpp"

namespace sc_expansion::atomic {

  template <typename T> class Diagram {
    public:
    explicit Diagram(Graph const &graph, std::vector<VertexType<T> *> const &vertex_types);

    T evaluate(std::vector<double> const &taus, HubbardSolver<1, T> const &solver, bool infinite_U);

    void mark_tau_dirty(int tau_index);
    void mark_all_dirty();

    double get_free_multiplicity() const { return this->graph.get_free_multiplicity(); }
    std::vector<ValidGlobalConfig> const &get_valid_configurations() const { return this->valid_configurations; }
    double get_diagram_sign() const { return (double)this->diagram_sign; }
    Graph const &get_graph() const { return this->graph; }

    private:
    Graph const &graph;
    std::vector<VertexType<T> *> vertex_type_ptrs;
    std::vector<ValidGlobalConfig> valid_configurations;
    std::vector<std::vector<LegInfo>> legs_per_vertex;
    Lines hopping_lines;
    int diagram_sign = 1;

    // [config_idx][vertex_idx] — local cache per (config, vertex) pair
    std::vector<std::vector<VertexInstance<T>>> vertex_instances;

    // tau_index → vertex indices that depend on it.
    std::vector<std::vector<int>> tau_to_vertices;

    void setup_vertices(std::vector<VertexType<T> *> const &vertex_types);
    void compute_valid_configurations();
    void build_vertex_instances();
  };

} // namespace sc_expansion::atomic
