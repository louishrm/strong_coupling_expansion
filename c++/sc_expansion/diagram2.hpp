#pragma once

#include <vector>
#include <map>
#include <numeric>   // For std::accumulate
#include <algorithm> // For std::next_permutation
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

  // Stores a unique spatial embedding pattern for the dimer (N_sites=2) case.
  // Each entry in `directions` is 0 (left/dx<0) or 1 (right/dx>0) for the corresponding hopping line.
  // `weight` is the total number of lattice embeddings producing this pattern.
  struct SpatialConfiguration {
    std::vector<uint8_t> directions; // Per hopping line: 0=left, 1=right
    double weight;
  };

  // Lattice embedding count for a graph on the square (bipartite) or triangular (non-bipartite) lattice.
  // Moved here from Graph — this is the N_sites=1 (single-site) free multiplicity.
  int compute_lattice_free_multiplicity(Graph const &graph);

  template <int N_sites, typename T> class Diagram2 {

    public:
    explicit Diagram2(Graph const &graph, std::vector<VertexType<N_sites, T> *> const &vertex_types);

    T evaluate(std::vector<double> const &taus, HubbardSolver<N_sites, T> const &solver, bool infinite_U);

    const std::vector<SpatialConfiguration> &get_spatial_configurations() const { return this->spatial_configurations; }

    double get_free_multiplicity() const {
      double total = 0.0;
      for (auto const &sc : this->spatial_configurations) { total += sc.weight; }
      return total;
    }

    private:
    Graph const &graph;
    std::vector<VertexInstance<N_sites, T>> vertices;
    std::vector<OrbitalConfiguration> valid_configurations;
    std::vector<SpatialConfiguration> spatial_configurations;
    Lines hopping_lines;

    void compute_hopping_lines();
    void setup_vertices(std::vector<VertexType<N_sites, T> *> const &vertex_types);
    void compute_spatial_configurations();
    void compute_valid_configurations();
    void compute_diagram_sign();

    // Helpers for dimer spatial embedding on the triangular superlattice
    void solve_dimer_embedding(int placed_count, std::vector<bool> &placed, std::vector<std::pair<int, int>> &coords,
                               std::map<std::vector<uint8_t>, int> &config_counts) const;

    std::vector<uint8_t> canonicalize_directions(std::vector<uint8_t> const &dirs,
                                                 std::vector<std::vector<int>> const &automorphisms) const;

    std::vector<uint8_t> apply_automorphism_to_directions(std::vector<uint8_t> const &dirs, std::vector<int> const &perm) const;
  };
}; // namespace sc_expansion
