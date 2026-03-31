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
      // For spatial embeddings, the site-swap (reflection) is already accounted for
      // by the spatial symmetry group (C2 inversion on the lattice).
      // Applying it here would double count the symmetry if we multiply by orbit size.
      // So we only consider spin-flip for internal physics symmetry.
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

  struct Point {
    int x;
    int y;
    Point() : x(0), y(0) {}
    Point(int x_, int y_) : x(x_), y(y_) {}
    bool operator==(const Point &other) const { return x == other.x && y == other.y; }
    bool operator<(const Point &other) const { return (x == other.x) ? (y < other.y) : (x < other.x); }
  };

  struct SpatialOrbit {
    std::vector<Point> representative;
    double spatial_weight;
    std::vector<int> signature; // for N_sites=2
  };

  struct LatticeSymmetry {
    static int get_internal_site(int dx, int dy) {
      // Staggered tiling rule: a consistent map from neighbor displacement to {0, 1}
      if (dx == 1 && dy == 0) return 1;
      if (dx == -1 && dy == 0) return 0;
      if (dx == 0 && dy == 1) return 1;
      if (dx == 0 && dy == -1) return 0;
      if (dx == 1 && dy == 1) return 1;
      if (dx == -1 && dy == -1) return 0;
      return 0;
    }

    static std::vector<std::vector<Point>> apply_D4(const std::vector<Point> &coords) {
      std::vector<std::vector<Point>> orbit;
      std::vector<std::pair<int, int>> transforms = {
        {1, 1}, {-1, 1}, {1, -1}, {-1, -1}
      };
      for (auto [sx, sy] : transforms) {
        std::vector<Point> p1, p2;
        for (const auto &pt : coords) {
          p1.push_back(Point(sx * pt.x, sy * pt.y));
          p2.push_back(Point(sx * pt.y, sy * pt.x));
        }
        orbit.push_back(p1);
        orbit.push_back(p2);
      }
      return orbit;
    }

    static std::vector<std::vector<Point>> apply_C2(const std::vector<Point> &coords) {
      std::vector<std::vector<Point>> orbit;
      orbit.push_back(coords);
      std::vector<Point> inv;
      for (const auto &pt : coords) {
        inv.push_back(Point(-pt.x, -pt.y));
      }
      orbit.push_back(inv);
      return orbit;
    }

    template <int N_sites>
    static std::vector<std::vector<Point>> get_orbit(const std::vector<Point> &coords, bool is_square) {
      std::vector<std::vector<Point>> orbit;
      if constexpr (N_sites == 1) {
        if (is_square) orbit = apply_D4(coords);
        else orbit = apply_C2(coords); // default for triangular atomic
      } else {
        orbit = apply_C2(coords);
      }
      std::sort(orbit.begin(), orbit.end());
      orbit.erase(std::unique(orbit.begin(), orbit.end()), orbit.end());
      return orbit;
    }
  };

  template <int N_sites, typename T> class Diagram2 {

    public:
    explicit Diagram2(Graph const &graph, std::vector<VertexType<N_sites, T> *> const &vertex_types);

    T evaluate(std::vector<double> const &taus, HubbardSolver<N_sites, T> const &solver, bool infinite_U);

    std::vector<std::vector<Point>> find_spatial_embeddings() const;
    std::vector<SpatialOrbit> group_spatial_embeddings(const std::vector<std::vector<Point>> &embeddings) const;

    std::vector<OrbitalConfiguration> const& get_valid_configs() const { return this->valid_configurations; }

    private:
    Graph const &graph;
    std::vector<VertexInstance<N_sites, T>> vertices;
    std::vector<OrbitalConfiguration> valid_configurations;
    Lines hopping_lines;

    void compute_hopping_lines();
    void setup_vertices(std::vector<VertexType<N_sites, T> *> const &vertex_types);
    void compute_valid_configurations();
    void compute_diagram_sign();

    void solve_embedding_recursive(int placed_count, std::vector<bool> &placed, std::vector<Point> &coords, std::vector<std::vector<Point>> &all_embeddings) const;
  };
}; // namespace sc_expansion
