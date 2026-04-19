#include "diagram.hpp"
#include <numeric>
#include <cmath>

namespace sc_expansion {

  // =====================================================================
  // compute_lattice_free_multiplicity: count the number of distinct lattice
  // embeddings of a graph (single-site free multiplicity).
  // Bipartite graphs are embedded on the square lattice (4 NN),
  // non-bipartite graphs on the triangular lattice (6 NN).
  // =====================================================================

  namespace {

    struct Point {
      int x, y;
    };

    long solve_embedding_recursive(Graph const &graph, bool bipartite_only, int V, int placed_count, std::vector<bool> &placed,
                                   std::vector<Point> &coords) {

      if (placed_count == V) return 1;

      int anchor = -1, target_node = -1;
      for (int candidate = 0; candidate < V; ++candidate) {
        if (!placed[candidate]) {
          for (int p = 0; p < V; ++p) {
            if (placed[p]) {
              uint8_t val = graph(candidate, p) + graph(p, candidate);
              if (val > 0) {
                target_node = candidate;
                anchor      = p;
                goto found_target;
              }
            }
          }
        }
      }
    found_target:;
      if (target_node == -1) return 0;

      long count = 0;

      static constexpr int sq_dx[4]  = {1, -1, 0, 0};
      static constexpr int sq_dy[4]  = {0, 0, 1, -1};
      static constexpr int tri_dx[6] = {1, -1, -1, 0, 0, 1};
      static constexpr int tri_dy[6] = {0, 0, 1, 1, -1, -1};

      const int *dx    = bipartite_only ? sq_dx : tri_dx;
      const int *dy    = bipartite_only ? sq_dy : tri_dy;
      int n_dirs       = bipartite_only ? 4 : 6;
      Point anchor_pos = coords[anchor];

      auto is_neighbor = [bipartite_only](Point p1, Point p2) -> bool {
        int ddx = p1.x - p2.x;
        int ddy = p1.y - p2.y;
        if (bipartite_only) { return std::abs(ddx) + std::abs(ddy) == 1; }
        if (std::abs(ddx) + std::abs(ddy) == 1) return true;
        if (ddx == -1 && ddy == 1) return true;
        if (ddx == 1 && ddy == -1) return true;
        return false;
      };

      for (int dir = 0; dir < n_dirs; ++dir) {
        Point candidate_pos = {anchor_pos.x + dx[dir], anchor_pos.y + dy[dir]};

        bool valid = true;
        for (int i = 0; i < V; ++i) {
          if (placed[i]) {
            uint8_t links = graph(target_node, i) + graph(i, target_node);
            if (links > 0) {
              if (!is_neighbor(candidate_pos, coords[i])) {
                valid = false;
                break;
              }
            }
          }
        }

        if (valid) {
          coords[target_node] = candidate_pos;
          placed[target_node] = true;
          count += solve_embedding_recursive(graph, bipartite_only, V, placed_count + 1, placed, coords);
          placed[target_node] = false;
        }
      }
      return count;
    }

  } // anonymous namespace

  int compute_lattice_free_multiplicity(Graph const &graph) {
    int V = graph.get_V();
    std::vector<Point> coords(V, {0, 0});
    std::vector<bool> placed(V, false);
    coords[0] = {0, 0};
    placed[0] = true;

    bool bipartite_only = graph.get_bipartite_only();
    return (int)solve_embedding_recursive(graph, bipartite_only, V, 1, placed, coords);
  }

  // =====================================================================
  // Constructor
  // =====================================================================
  template <typename T>
  Diagram<T>::Diagram(Graph const &graph_, std::vector<VertexType<T> *> const &vertex_types) : graph(graph_) {
    this->compute_hopping_lines();
    this->setup_vertices(vertex_types);
    this->compute_valid_configurations();
    this->compute_diagram_sign();
    this->build_vertex_instances();
  }

  // =====================================================================
  // compute_hopping_lines: enumerate all directed hopping lines from the adjacency matrix
  // =====================================================================
  template <typename T> void Diagram<T>::compute_hopping_lines() {

    for (int i = 0; i < this->graph.get_V(); i++) {
      for (int j = 0; j < this->graph.get_V(); j++) {

        int lines_ij = this->graph(i, j);
        if (lines_ij != 0) {

          for (int k = 0; k < lines_ij; k++) {
            HoppingLine line;
            line.from_vertex = i;
            line.to_vertex   = j;
            this->hopping_lines.lines.push_back(line);
          }
        }
      }
    }
  }

  // =====================================================================
  // setup_vertices: map each graph vertex to its VertexType for caching.
  // =====================================================================
  template <typename T>
  void Diagram<T>::setup_vertices(std::vector<VertexType<T> *> const &vertex_types) {
    int V = this->graph.get_V();
    this->vertex_type_ptrs.resize(V, nullptr);

    for (int v = 0; v < V; v++) {
      int degree = 0;
      for (int j = 0; j < V; j++) { degree += this->graph(v, j) + this->graph(j, v); }
      int cumulant_order = degree / 2;
      int vt_idx         = cumulant_order - 1;
      if (vt_idx >= 0 && vt_idx < (int)vertex_types.size()) { this->vertex_type_ptrs[v] = vertex_types[vt_idx]; }
    }
  }

  // =====================================================================
  // compute_valid_configurations: enumerate all 2^n_lines spin assignments,
  // check spin conservation at each vertex, canonicalize under spin-flip,
  // and store the canonical representative with its weight.
  // =====================================================================
  template <typename T> void Diagram<T>::compute_valid_configurations() {

    int V       = this->graph.get_V();
    int n_lines = (int)this->hopping_lines.lines.size();
    double sym_factor = this->graph.get_symmetry_factor();
    double free_mult  = this->graph.get_free_multiplicity();

    constexpr uint8_t ACTION_BIT = FermionOperator<T>::ACTION_BIT;

    // --- Step 1: Build per-vertex leg info ---
    this->legs_per_vertex.resize(V);
    for (int k = 0; k < n_lines; k++) {
      int from = this->hopping_lines.lines[k].from_vertex;
      int to   = this->hopping_lines.lines[k].to_vertex;
      this->legs_per_vertex[from].push_back({k, true});  // annihilation (source)
      this->legs_per_vertex[to].push_back({k, false});    // creation (dest)
    }

    // --- Step 2: Enumerate spin assignments ---
    using Config = GlobalConfiguration;
    std::map<std::vector<uint8_t>, double> canonical_weights;
    int n_spin_configs = 1 << n_lines;

    std::set<std::vector<uint8_t>> seen;

    for (int spin_mask = 0; spin_mask < n_spin_configs; spin_mask++) {

      Config global;
      global.config.reserve(2 * n_lines);
      bool valid = true;

      for (int v = 0; v < V && valid; v++) {
        int up_cre = 0, up_ann = 0, dn_cre = 0, dn_ann = 0;

        for (auto const &leg : this->legs_per_vertex[v]) {
          int spin = (spin_mask >> leg.line_index) & 1; // 0 = down, 1 = up

          uint8_t orbital = spin; // single site: orbital == spin bit
          uint8_t op_id   = leg.is_source ? orbital : (ACTION_BIT | orbital);
          global.config.push_back(op_id);

          if (leg.is_source) {
            if (spin == 1) up_ann++;
            else dn_ann++;
          } else {
            if (spin == 1) up_cre++;
            else dn_cre++;
          }
        }

        if (up_cre != up_ann || dn_cre != dn_ann) { valid = false; }
      }

      if (!valid) continue;

      // --- Step 3: Canonicalize under the spin-flip symmetry group ---
      auto orbit      = SymmetryGroup<T>::get_orbit(global);
      auto &canonical = orbit[0];
      int orbit_size  = (int)orbit.size();

      if (!seen.insert(canonical.config).second) continue;

      double weight = free_mult * orbit_size / sym_factor;
      canonical_weights[canonical.config] += weight;
    }

    // --- Step 4: Store results ---
    for (auto &[config, weight] : canonical_weights) {
      this->valid_configurations.push_back({config, weight});
    }
  }

  // =====================================================================
  // compute_diagram_sign: count the number of fermion loops.
  // =====================================================================
  template <typename T> void Diagram<T>::compute_diagram_sign() {

    int V       = this->graph.get_V();
    int n_lines = (int)this->hopping_lines.lines.size();

    std::vector<int> successor_map(n_lines);

    for (int v = 0; v < V; v++) {
      std::vector<int> incoming, outgoing;
      for (auto const &leg : this->legs_per_vertex[v]) {
        if (leg.is_source) {
          outgoing.push_back(leg.line_index);
        } else {
          incoming.push_back(leg.line_index);
        }
      }
      for (size_t i = 0; i < incoming.size(); i++) { successor_map[incoming[i]] = outgoing[i]; }
    }

    int num_loops = 0;
    std::vector<bool> visited(n_lines, false);
    for (int i = 0; i < n_lines; i++) {
      if (!visited[i]) {
        num_loops++;
        int current = i;
        while (!visited[current]) {
          visited[current] = true;
          current          = successor_map[current];
        }
      }
    }

    this->diagram_sign = (num_loops % 2 == 0) ? 1 : -1;
  }

  // =====================================================================
  // build_vertex_instances: create VertexInstances for each (config, vertex) pair.
  // =====================================================================
  template <typename T> void Diagram<T>::build_vertex_instances() {

    int V       = this->graph.get_V();
    int n_lines = (int)this->hopping_lines.lines.size();

    bool has_any_type = false;
    for (int v = 0; v < V; v++) {
      if (this->vertex_type_ptrs[v] != nullptr) { has_any_type = true; break; }
    }
    if (!has_any_type) return;

    this->tau_to_vertices.resize(n_lines);
    for (int v = 0; v < V; v++) {
      for (auto const &leg : this->legs_per_vertex[v]) { this->tau_to_vertices[leg.line_index].push_back(v); }
    }
    for (auto &vlist : this->tau_to_vertices) {
      std::sort(vlist.begin(), vlist.end());
      vlist.erase(std::unique(vlist.begin(), vlist.end()), vlist.end());
    }

    this->vertex_instances.resize(this->valid_configurations.size());
    for (int gc_idx = 0; gc_idx < (int)this->valid_configurations.size(); gc_idx++) {
      auto const &gc = this->valid_configurations[gc_idx];
      int cfg_offset = 0;

      for (int v = 0; v < V; v++) {
        int n_legs = (int)this->legs_per_vertex[v].size();

        std::vector<int> tau_indices(n_legs);
        for (int i = 0; i < n_legs; i++) { tau_indices[i] = this->legs_per_vertex[v][i].line_index; }

        std::vector<uint8_t> op_ids(gc.config.begin() + cfg_offset, gc.config.begin() + cfg_offset + n_legs);
        cfg_offset += n_legs;

        this->vertex_instances[gc_idx].emplace_back(this->vertex_type_ptrs[v], std::move(tau_indices), std::move(op_ids));
      }
    }
  }

  // =====================================================================
  // mark_tau_dirty
  // =====================================================================
  template <typename T> void Diagram<T>::mark_tau_dirty(int tau_index) {
    for (int v : this->tau_to_vertices[tau_index]) {
      for (auto &config_instances : this->vertex_instances) { config_instances[v].mark_dirty(); }
    }
  }

  // =====================================================================
  // mark_all_dirty
  // =====================================================================
  template <typename T> void Diagram<T>::mark_all_dirty() {
    for (auto &config_instances : this->vertex_instances) {
      for (auto &vi : config_instances) { vi.mark_dirty(); }
    }
  }

  // =====================================================================
  // evaluate
  // =====================================================================
  template <typename T>
  T Diagram<T>::evaluate(std::vector<double> const &taus, HubbardSolver<T> const &solver, bool infinite_U) {

    int V = this->graph.get_V();
    T sum = T(0.0);

    for (int gc_idx = 0; gc_idx < (int)this->valid_configurations.size(); gc_idx++) {
      T product = T(1.0);
      for (int v = 0; v < V; v++) { product = product * this->vertex_instances[gc_idx][v].get_value(taus, solver, infinite_U); }
      sum = sum + T(this->valid_configurations[gc_idx].weight) * product;
    }

    T prefactor = (T(-1.0) / solver.params.beta) * T(this->diagram_sign);
    return prefactor * sum;
  }

  // =====================================================================
  // Explicit template instantiations
  // =====================================================================
  template class Diagram<double>;
  template class Diagram<Dual>;

} // namespace sc_expansion
