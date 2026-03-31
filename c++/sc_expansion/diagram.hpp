#pragma once

#include <vector>
#include <numeric>   // For std::accumulate
#include <algorithm> // For std::next_permutation
#include <queue>
#include "./hubbard_solver.hpp"
#include "./cumulant.hpp"
#include "./graph.hpp"
#include "./dual.hpp"

namespace sc_expansion {

  struct Line {
    int from_vertex;
    int to_vertex;
  };

  struct Vertex {
    std::vector<int> outgoing_lines; //annihilation indices
    std::vector<int> incoming_lines; //creation indices

    std::vector<int> local_spin_configs;

    int degree() const { return outgoing_lines.size() + incoming_lines.size(); }
  };

  class Diagram {

    public:
    explicit Diagram(Graph const &graph);

    std::vector<Line> const &get_hopping_lines() const { return this->hopping_lines; }
    std::vector<Vertex> const &get_vertices() const { return this->vertices; }
    double get_diagram_sign() const { return (double)this->diagram_sign; }
    Graph const &get_graph() const { return this->graph; }
    const std::vector<uint64_t> &get_global_configs() const { return this->global_spin_configurations; }
    const std::vector<double> &get_config_weights() const { return this->config_weights; }

    private:
    Graph graph;
    std::vector<Line> hopping_lines;
    std::vector<Vertex> vertices;
    std::vector<uint64_t> global_spin_configurations;
    std::vector<double> config_weights;
    int diagram_sign;

    void compute_hopping_lines_and_vertex_structures();
    void compute_global_spin_configurations();
    void compute_local_spin_configurations();
    void compute_diagram_sign();
  };

  template <int N_sites, typename T> class DiagramEvaluator {

    public:
    explicit DiagramEvaluator(Diagram const &diagram, Parameters<T> const &params);

    T evaluate_at_taus(std::vector<double> const &taus, bool infinite_U, bool use_cache) const;
    T evaluate_at_taus_dimer(std::vector<double> const &taus, bool infinite_U, bool use_cache) const;
    Diagram const &get_diagram() const { return this->diagram; }

    private:
    const Diagram &diagram;
    HubbardSolver<N_sites, T> solver;
    mutable std::vector<double> current_taus;
    mutable std::vector<std::vector<T>> cache_finite;
    mutable std::vector<std::vector<T>> cache_infinite;

    void check_vertex(int v_idx, std::vector<double> const &taus) const;
    void recompute_vertex(int v_idx, std::vector<double> const &taus) const;
    std::pair<Args<N_sites, T>, Args<N_sites, T>> get_local_cumul_args(int v_idx, std::vector<double> const &taus, uint32_t local_mask) const;
  };

} // namespace sc_expansion
