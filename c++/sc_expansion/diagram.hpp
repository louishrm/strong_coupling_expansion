#pragma once

#include <vector>
#include <map>
#include <set>
#include <numeric>
#include <algorithm>
#include <queue>
#include "./hubbard_solver.hpp"
#include "./cumulant.hpp"
#include "./graph.hpp"
#include "./dual.hpp"
#include "./fock_space.hpp"
#include "./vertex.hpp"

namespace sc_expansion {

  struct HoppingLine {
    int from_vertex;
    int to_vertex;
  };

  struct Lines {
    std::vector<HoppingLine> lines;
  };

  // Helper to store valid physics states during generation
  struct GlobalConfiguration {
    std::vector<uint8_t> config;

    bool operator<(const GlobalConfiguration &other) const { return this->config < other.config; }
    bool operator==(const GlobalConfiguration &other) const { return this->config == other.config; }
  };

  // Helper for the cumulant symmetry group (spin flip).
  template <typename T> struct SymmetryGroup {
    using Config = GlobalConfiguration;

    static Config apply_spin_flip(Config const &c) {
      Config res = c;
      for (auto &op_id : res.config) { op_id = FermionOperator<T>(op_id).apply_spin_flip().op; }
      return res;
    }

    static std::vector<Config> get_orbit(Config const &c) {
      std::vector<Config> orbit;
      orbit.push_back(c);
      orbit.push_back(apply_spin_flip(c));
      std::sort(orbit.begin(), orbit.end());
      orbit.erase(std::unique(orbit.begin(), orbit.end()), orbit.end());
      return orbit;
    }

    static Config get_canonical(Config const &c) {
      auto orbit = get_orbit(c);
      return orbit[0];
    }
  };

  // Per-leg info at a vertex: which hopping line it belongs to and its role
  struct LegInfo {
    int line_index; // Index into hopping_lines
    bool is_source; // true = annihilation (from_vertex), false = creation (to_vertex)
  };

  // A symmetry-reduced global configuration with its total weight.
  struct ValidGlobalConfig {
    std::vector<uint8_t> config;
    double weight;
  };

  // Lattice embedding count for a graph on the square (bipartite) or triangular (non-bipartite) lattice.
  int compute_lattice_free_multiplicity(Graph const &graph);

  template <typename T> class Diagram {

    public:
    explicit Diagram(Graph const &graph, std::vector<VertexType<T> *> const &vertex_types);

    T evaluate(std::vector<double> const &taus, HubbardSolver<T> const &solver, bool infinite_U);

    void mark_tau_dirty(int tau_index);
    void mark_all_dirty();

    double get_free_multiplicity() const { return this->graph.get_free_multiplicity(); }

    const std::vector<ValidGlobalConfig> &get_valid_configurations() const { return this->valid_configurations; }
    double get_diagram_sign() const { return (double)this->diagram_sign; }
    Graph const &get_graph() const { return this->graph; }

    private:
    Graph const &graph;
    std::vector<VertexType<T> *> vertex_type_ptrs;
    std::vector<ValidGlobalConfig> valid_configurations;
    std::vector<std::vector<LegInfo>> legs_per_vertex;
    Lines hopping_lines;
    int diagram_sign = 1;

    // VertexInstances: [config_idx][vertex_idx]
    std::vector<std::vector<VertexInstance<T>>> vertex_instances;

    // Precomputed inverse map: tau_index → list of vertex indices that depend on it.
    std::vector<std::vector<int>> tau_to_vertices;

    void compute_hopping_lines();
    void setup_vertices(std::vector<VertexType<T> *> const &vertex_types);
    void compute_valid_configurations();
    void compute_diagram_sign();
    void build_vertex_instances();
  };
}; // namespace sc_expansion
