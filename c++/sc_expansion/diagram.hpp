#pragma once

#include <vector>
#include <numeric>   // For std::accumulate
#include <algorithm> // For std::next_permutation
#include <queue>
#include "./hubbard_atom.hpp"
#include "./cumulant.hpp"
#include "./graph.hpp" // Include Graph for compute_free_multiplicity logic

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

  template <typename T>
  class DiagramEvaluator {

    public:
    explicit DiagramEvaluator(Diagram const &diagram, Parameters<T> const &params);

    T evaluate_at_taus(std::vector<double> const &taus, bool infinite_U, bool use_cache) const;
    T evaluate_at_taus_dimer(std::vector<double> const &taus, bool infinite_U, bool use_cache) const;
    Diagram const &get_diagram() const { return this->diagram; }

    private:
    const Diagram &diagram;
    HubbardAtom<T> atom;
    mutable std::vector<double> current_taus;
    mutable std::vector<std::vector<T>> cache_finite;
    mutable std::vector<std::vector<T>> cache_infinite;

    void check_vertex(int v_idx, std::vector<double> const &taus) const;
    void recompute_vertex(int v_idx, std::vector<double> const &taus) const;
    std::pair<ArgList, ArgList> get_local_cumul_args(int v_idx, std::vector<double> const &taus,
                                                                                     uint32_t local_mask) const;
  };

  // class Diagram {

  //   public:
  //   explicit Diagram(Graph const &graph);

  //   double get_diagram_sign() const { return (double)this->diagram_sign; }
  //   std::vector<long> get_valid_spin_configurations() const { return this->valid_spin_configurations; }
  //   std::vector<Line> get_hopping_lines() const { return this->hopping_lines; }
  //   Graph const &get_graph() const { return this->graph; }

  //   const std::vector<int> &get_unprimed_indices(int vertex) const { return this->unprimed_line_indices_per_vertex[vertex]; }
  //   const std::vector<int> &get_primed_indices(int vertex) const { return this->primed_line_indices_per_vertex[vertex]; }

  //   private:
  //   void compute_hopping_lines();
  //   void compute_valid_spin_configurations();
  //   void compute_diagram_sign();
  //   void compute_vertex_structures();

  //   Graph graph;
  //   std::vector<Line> hopping_lines;
  //   std::vector<std::vector<int>> unprimed_line_indices_per_vertex;
  //   std::vector<std::vector<int>> primed_line_indices_per_vertex;
  //   std::vector<long> valid_spin_configurations;
  //   int diagram_sign;
  // };

  // class DiagramEvaluator {

  //   public:
  //   explicit DiagramEvaluator(Diagram const &diagram, Parameters const &params);

  //   double evaluate_at_taus(std::vector<double> const &taus, bool infinite_U) const;

  //   private:
  //   double evaluate_at_points(HubbardAtom::cumul_args const &args, bool infinite_U) const;

  //   mutable std::vector<double> current_taus;

  //   void recompute_vertex(int vertex, std::vector<double> const &taus, bool infinite_U);
  //   void solve_vertex(int vertex, std::vector<double> const &taus, bool infinite_U);

  //   VertexCache vertex_cache;

  //   Diagram const &diagram;
  //   HubbardAtom atom;
  // };

  // struct VertexCache {

  //   std::unordered_map<int, std::unordered_map<long, double>> cache;

  //   std::vector<double> current_times;

  //   void mark_as_corrupted(int vertex);

  //   void update_cache(int vertex, double new_value);

  //   void assign_times_to_vertices(std::vector<double> const &taus);
  // };

} // namespace sc_expansion
