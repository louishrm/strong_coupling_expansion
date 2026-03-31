#include "diagram2.hpp"
#include <numeric>
#include <cmath>

namespace sc_expansion {

  // =====================================================================
  // compute_lattice_free_multiplicity: count the number of distinct lattice
  // embeddings of a graph (single-site / N_sites=1 free multiplicity).
  // Bipartite graphs are embedded on the square lattice (4 NN),
  // non-bipartite graphs on the triangular lattice (6 NN).
  // Moved here from Graph — this is a lattice-specific computation.
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

      // Square lattice: (1,0), (-1,0), (0,1), (0,-1)
      // Triangular lattice: (1,0), (-1,0), (-1,1), (0,1), (0,-1), (1,-1)
      static constexpr int sq_dx[4]  = {1, -1, 0, 0};
      static constexpr int sq_dy[4]  = {0, 0, 1, -1};
      static constexpr int tri_dx[6] = {1, -1, -1, 0, 0, 1};
      static constexpr int tri_dy[6] = {0, 0, 1, 1, -1, -1};

      const int *dx     = bipartite_only ? sq_dx : tri_dx;
      const int *dy     = bipartite_only ? sq_dy : tri_dy;
      int n_dirs        = bipartite_only ? 4 : 6;
      Point anchor_pos  = coords[anchor];

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

    bool bipartite_only = graph.get_bipartite_only(); // square lattice if true, triangular otherwise
    return (int)solve_embedding_recursive(graph, bipartite_only, V, 1, placed, coords);
  }

  // =====================================================================
  // Constructor
  // =====================================================================
  template <int N_sites, typename T>
  Diagram2<N_sites, T>::Diagram2(Graph const &graph_, std::vector<VertexType<N_sites, T> *> const &vertex_types) : graph(graph_) {
    this->compute_hopping_lines();
    this->compute_spatial_configurations();
    this->setup_vertices(vertex_types);
    this->compute_valid_configurations();
    this->compute_diagram_sign();
  }

  // =====================================================================
  // compute_hopping_lines: enumerate all directed hopping lines from the adjacency matrix
  // =====================================================================
  template <int N_sites, typename T> void Diagram2<N_sites, T>::compute_hopping_lines() {

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
  // compute_spatial_configurations: for N_sites=2, embed the graph on the
  // triangular dimer superlattice and group embeddings by hopping direction
  // pattern. For N_sites=1, a single trivial config with the free multiplicity.
  // =====================================================================
  template <int N_sites, typename T> void Diagram2<N_sites, T>::compute_spatial_configurations() {

    if constexpr (N_sites == 1) {
      // Single-site expansion: one trivial spatial config
      SpatialConfiguration config;
      config.directions.resize(this->hopping_lines.lines.size(), 0);
      config.weight = compute_lattice_free_multiplicity(this->graph);
      this->spatial_configurations.push_back(config);
      return;
    }

    if constexpr (N_sites == 2) {
      int V       = this->graph.get_V();
      int n_lines = (int)this->hopping_lines.lines.size();

      // --- Step 1: Recursive embedding on the triangular superlattice ---
      // Collect raw direction vectors with their embedding counts.
      std::map<std::vector<uint8_t>, int> raw_counts;

      std::vector<std::pair<int, int>> coords(V, {0, 0});
      std::vector<bool> placed(V, false);
      coords[0] = {0, 0};
      placed[0] = true;

      this->solve_dimer_embedding(1, placed, coords, raw_counts);

      // --- Step 2: Compute graph automorphisms ---
      // An automorphism is a vertex permutation pi such that
      //   adj[v][w] == adj[pi(v)][pi(w)] for all v, w.

      // Precompute degrees for pruning
      std::vector<int> degrees(V, 0);
      for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) { degrees[i] += this->graph(i, j) + this->graph(j, i); }
      }

      std::vector<std::vector<int>> automorphisms;
      std::vector<int> perm(V);
      std::iota(perm.begin(), perm.end(), 0);

      do {
        // Prune: check degree preservation
        bool degree_ok = true;
        for (int i = 0; i < V; i++) {
          if (degrees[perm[i]] != degrees[i]) {
            degree_ok = false;
            break;
          }
        }
        if (!degree_ok) continue;

        // Check full adjacency preservation
        bool is_auto = true;
        for (int i = 0; i < V && is_auto; i++) {
          for (int j = 0; j < V && is_auto; j++) {
            if (this->graph(i, j) != this->graph(perm[i], perm[j])) { is_auto = false; }
          }
        }

        if (is_auto) { automorphisms.push_back(perm); }

      } while (std::next_permutation(perm.begin(), perm.end()));

      // --- Step 3: Canonicalize each raw config and merge ---
      std::map<std::vector<uint8_t>, double> merged;
      for (auto &[dirs, count] : raw_counts) {
        auto canonical = this->canonicalize_directions(dirs, automorphisms);
        merged[canonical] += count;
      }

      // --- Step 4: Store results ---
      for (auto &[dirs, weight] : merged) {
        SpatialConfiguration config;
        config.directions = dirs;
        config.weight     = weight;
        this->spatial_configurations.push_back(config);
      }
    }
  }

  // =====================================================================
  // solve_dimer_embedding: recursive backtracking to place vertices on the
  // triangular superlattice. For each complete embedding, compute the
  // direction vector (per hopping line: 0=left, 1=right) and accumulate.
  //
  // The staggered dimer tiling maps to a triangular superlattice with
  // 6 nearest neighbours. In abstract coords (n1, n2), the real x-displacement
  // is dx_real = 2*dn1 + dn2. This is always nonzero for NN, giving a
  // clean 3-left / 3-right split:
  //   Right (dx_real > 0): (1,0), (0,1), (1,-1)
  //   Left  (dx_real < 0): (-1,0), (0,-1), (-1,1)
  // =====================================================================
  template <int N_sites, typename T>
  void Diagram2<N_sites, T>::solve_dimer_embedding(int placed_count, std::vector<bool> &placed,
                                                   std::vector<std::pair<int, int>> &coords,
                                                   std::map<std::vector<uint8_t>, int> &config_counts) const {

    int V = this->graph.get_V();

    // Base case: all vertices placed
    if (placed_count == V) {
      // Compute the direction vector for this embedding
      int n_lines = (int)this->hopping_lines.lines.size();
      std::vector<uint8_t> dirs(n_lines);

      for (int k = 0; k < n_lines; k++) {
        int from = this->hopping_lines.lines[k].from_vertex;
        int to   = this->hopping_lines.lines[k].to_vertex;
        int dn1  = coords[to].first - coords[from].first;
        int dn2  = coords[to].second - coords[from].second;
        // Real x-displacement on the staggered dimer lattice
        int dx_real = 2 * dn1 + dn2;
        dirs[k]     = (dx_real > 0) ? 1 : 0;
      }

      config_counts[dirs]++;
      return;
    }

    // Find an unplaced vertex (target) connected to a placed vertex (anchor)
    int anchor = -1, target = -1;
    for (int c = 0; c < V; ++c) {
      if (!placed[c]) {
        for (int p = 0; p < V; ++p) {
          if (placed[p]) {
            uint8_t links = this->graph(c, p) + this->graph(p, c);
            if (links > 0) {
              target = c;
              anchor = p;
              goto found_target;
            }
          }
        }
      }
    }
  found_target:;
    if (target == -1) return; // Should not happen for connected graphs

    // Triangular lattice: 6 nearest-neighbour directions
    static constexpr int dx[6] = {1, -1, -1, 0, 0, 1};
    static constexpr int dy[6] = {0, 0, 1, 1, -1, -1};

    auto is_triangular_neighbor = [](int x1, int y1, int x2, int y2) -> bool {
      int ddx = x1 - x2;
      int ddy = y1 - y2;
      if (std::abs(ddx) + std::abs(ddy) == 1) return true;
      if (ddx == -1 && ddy == 1) return true;
      if (ddx == 1 && ddy == -1) return true;
      return false;
    };

    int ax = coords[anchor].first;
    int ay = coords[anchor].second;

    for (int dir = 0; dir < 6; ++dir) {
      int cx = ax + dx[dir];
      int cy = ay + dy[dir];

      // Check consistency with ALL placed vertices
      bool valid = true;
      for (int i = 0; i < V; ++i) {
        if (placed[i]) {
          uint8_t links = this->graph(target, i) + this->graph(i, target);
          if (links > 0) {
            if (!is_triangular_neighbor(cx, cy, coords[i].first, coords[i].second)) {
              valid = false;
              break;
            }
          }
        }
      }

      if (valid) {
        coords[target] = {cx, cy};
        placed[target]  = true;
        this->solve_dimer_embedding(placed_count + 1, placed, coords, config_counts);
        placed[target] = false;
      }
    }
  }

  // =====================================================================
  // apply_automorphism_to_directions: given a direction vector and a vertex
  // permutation (graph automorphism), compute the permuted direction vector.
  //
  // Under automorphism pi, hopping line (v -> w, mult m) maps to
  // (pi(v) -> pi(w), mult m). The direction of the mapped line is the same
  // as the original.
  // =====================================================================
  template <int N_sites, typename T>
  std::vector<uint8_t> Diagram2<N_sites, T>::apply_automorphism_to_directions(std::vector<uint8_t> const &dirs,
                                                                              std::vector<int> const &perm) const {

    int V       = this->graph.get_V();
    int n_lines = (int)dirs.size();

    // Build (from, to, mult_idx) for each hopping line
    // and a reverse lookup: start_idx[from][to] = first line index for pair (from, to)
    std::vector<int> line_from(n_lines), line_to(n_lines), line_mult(n_lines);
    std::vector<std::vector<int>> start_idx(V, std::vector<int>(V, -1));
    {
      int idx = 0;
      for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
          int count = this->graph(i, j);
          if (count > 0 && start_idx[i][j] == -1) { start_idx[i][j] = idx; }
          for (int k = 0; k < count; k++) {
            line_from[idx] = i;
            line_to[idx]   = j;
            line_mult[idx] = k;
            idx++;
          }
        }
      }
    }

    // Apply the permutation
    std::vector<uint8_t> result(n_lines);
    for (int l = 0; l < n_lines; l++) {
      int new_from = perm[line_from[l]];
      int new_to   = perm[line_to[l]];
      int new_idx  = start_idx[new_from][new_to] + line_mult[l];
      result[new_idx] = dirs[l];
    }

    return result;
  }

  // =====================================================================
  // canonicalize_directions: compute the lexicographic minimum of the
  // direction vector under all graph automorphisms x {identity, lattice inversion}.
  //
  // Lattice inversion (x -> -x on the superlattice) flips all directions
  // (0 <-> 1), corresponding to the dimer site swap (Reflect symmetry).
  // =====================================================================
  template <int N_sites, typename T>
  std::vector<uint8_t> Diagram2<N_sites, T>::canonicalize_directions(std::vector<uint8_t> const &dirs,
                                                                     std::vector<std::vector<int>> const &automorphisms) const {

    auto min_dirs = dirs;

    for (auto const &perm : automorphisms) {
      // Apply automorphism
      auto permuted = this->apply_automorphism_to_directions(dirs, perm);
      if (permuted < min_dirs) { min_dirs = permuted; }

      // Apply automorphism + lattice inversion (flip all directions)
      auto inverted = permuted;
      for (auto &d : inverted) { d = 1 - d; }
      if (inverted < min_dirs) { min_dirs = inverted; }
    }

    return min_dirs;
  }

  // =====================================================================
  // Stubs for methods not yet implemented
  // =====================================================================
  template <int N_sites, typename T>
  void Diagram2<N_sites, T>::setup_vertices(std::vector<VertexType<N_sites, T> *> const & /*vertex_types*/) {
    // TODO: Create VertexInstance objects for each graph vertex,
    // assigning tau_indices and linking to the appropriate VertexType.
  }

  template <int N_sites, typename T> void Diagram2<N_sites, T>::compute_valid_configurations() {
    // TODO: For each spatial config, enumerate valid spin assignments,
    // build global configurations, and reduce by symmetry group.
  }

  template <int N_sites, typename T> void Diagram2<N_sites, T>::compute_diagram_sign() {
    // TODO: Compute the overall fermion sign of the diagram.
  }

  template <int N_sites, typename T>
  T Diagram2<N_sites, T>::evaluate(std::vector<double> const & /*taus*/, HubbardSolver<N_sites, T> const & /*solver*/, bool /*infinite_U*/) {
    // TODO: Sum vertex contributions over all valid configurations.
    return T{};
  }

  // =====================================================================
  // Explicit template instantiations
  // =====================================================================
  template class Diagram2<1, double>;
  template class Diagram2<1, Dual>;
  template class Diagram2<2, double>;
  template class Diagram2<2, Dual>;

} // namespace sc_expansion
