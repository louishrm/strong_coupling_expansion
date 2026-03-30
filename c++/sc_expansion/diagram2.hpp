#pragma once

#include <vector>
#include <numeric>   // For std::accumulate
#include <algorithm> // For std::next_permutation
#include <queue>
#include "./hubbard_solver.hpp"
#include "./cumulant.hpp"
#include "./graph.hpp" // Include Graph for compute_free_multiplicity logic
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
  template <int N_sites> struct GlobalConfiguration {
    std::vector<uint8_t> config; // Size = Total number of legs in the diagram

    bool operator<(const GlobalConfiguration &other) const { return this->config < other.config; }
    bool operator==(const GlobalConfiguration &other) const { return this->config == other.config; }
  };

  // Helper for symmetry operations (Spin Flip and Lattice Reflection)
  template <int N_sites, typename T> struct SymmetryGroup {
    using Config = GlobalConfiguration<N_sites>;

    static Config apply_spin_flip(Config const &c) {
      Config res = c;
      for (auto &op_id : res.config) { op_id = FermionOperator<N_sites, T>(op_id).apply_spin_flip().op; }
      return res;
    }

    static Config apply_reflection(Config const &c) {
      Config res = c;
      if constexpr (N_sites == 2) {
        for (auto &op_id : res.config) { op_id = FermionOperator<N_sites, T>(op_id).apply_reflection().op; }
      }
      return res;
    }

    // Computes all unique configurations in the orbit of the given configuration
    static std::vector<Config> get_orbit(Config const &c) {
      std::vector<Config> orbit;
      orbit.push_back(c);
      orbit.push_back(apply_spin_flip(c));
      if constexpr (N_sites == 2) {
        Config refl = apply_reflection(c);
        orbit.push_back(refl);
        orbit.push_back(apply_spin_flip(refl));
      }
      std::sort(orbit.begin(), orbit.end());
      orbit.erase(std::unique(orbit.begin(), orbit.end()), orbit.end());
      return orbit;
    }

    static Config get_canonical(Config const &c) {
      auto orbit = get_orbit(c);
      return orbit[0]; // Lexicographically smallest
    }
  };

  // Storage for the orbit representative and its multiplicity
  struct OrbitalConfiguration {
    std::vector<uint16_t> vertex_config_ids; // Indices into vertices[v].unique_configs
    double weight;                           // Number of configurations in the symmetry orbit
  };

  template <int N_sites, typename T> class Diagram2 {

    public:
    explicit Diagram2(Graph const &graph, std::vector<VertexType<N_sites, T> *> const &vertex_types);

    T evaluate(std::vector<double> const &taus, HubbardSolver<N_sites, T> const &solver, bool infinite_U);

    private:
    Graph const &graph;
    std::vector<VertexInstance<N_sites, T>> vertices;
    std::vector<OrbitalConfiguration> valid_configurations;
    Lines hopping_lines;

    void compute_hopping_lines();
    void setup_vertices();
    void compute_valid_configurations();
    void compute_diagram_sign();
  };
}; // namespace sc_expansion
