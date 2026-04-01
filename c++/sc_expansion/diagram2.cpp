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
  void Diagram2<N_sites, T>::solve_dimer_embedding(int placed_count, std::vector<bool> &placed, std::vector<std::pair<int, int>> &coords,
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
        placed[target] = true;
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
  std::vector<uint8_t> Diagram2<N_sites, T>::apply_automorphism_to_directions(std::vector<uint8_t> const &dirs, std::vector<int> const &perm) const {

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
      int new_from    = perm[line_from[l]];
      int new_to      = perm[line_to[l]];
      int new_idx     = start_idx[new_from][new_to] + line_mult[l];
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
  // setup_vertices: map each graph vertex to its VertexType for caching.
  // vertex_types is indexed by (cumulant_order - 1), where
  // cumulant_order = degree / 2 for a vertex.
  // =====================================================================
  template <int N_sites, typename T>
  void Diagram2<N_sites, T>::setup_vertices(std::vector<VertexType<N_sites, T> *> const &vertex_types) {
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
  // compute_valid_configurations: for each spatial config, enumerate all
  // 2^n_lines spin assignments, check spin conservation at each vertex,
  // build the global op_id config, canonicalize under the symmetry group
  // {Id, SpinFlip, Reflect, SpinFlip+Reflect}, and store the canonical
  // representative with weight = spatial_weight * orbit_size / automorphism_count.
  // =====================================================================
  template <int N_sites, typename T> void Diagram2<N_sites, T>::compute_valid_configurations() {

    int V       = this->graph.get_V();
    int n_lines = (int)this->hopping_lines.lines.size();
    int auto_count = this->graph.get_automorphism_count();

    constexpr uint8_t ACTION_BIT = FermionOperator<N_sites, T>::ACTION_BIT;

    // --- Step 1: Build per-vertex leg info ---
    this->legs_per_vertex.resize(V);
    for (int k = 0; k < n_lines; k++) {
      int from = this->hopping_lines.lines[k].from_vertex;
      int to   = this->hopping_lines.lines[k].to_vertex;
      this->legs_per_vertex[from].push_back({k, true});  // annihilation (source)
      this->legs_per_vertex[to].push_back({k, false});    // creation (dest)
    }

    // --- Step 2: Enumerate spin assignments for each spatial config ---
    using Config = GlobalConfiguration<N_sites>;
    std::map<std::vector<uint8_t>, double> canonical_weights;
    int n_spin_configs = 1 << n_lines;

    for (auto const &spatial : this->spatial_configurations) {

      // Track which canonicals we've already seen from THIS spatial config
      // to avoid double-counting orbit members enumerated from different spins.
      std::set<std::vector<uint8_t>> seen_this_spatial;

      for (int spin_mask = 0; spin_mask < n_spin_configs; spin_mask++) {

        Config global;
        global.config.reserve(2 * n_lines); // total legs = 2 * n_lines
        bool valid = true;

        for (int v = 0; v < V && valid; v++) {
          int up_cre = 0, up_ann = 0, dn_cre = 0, dn_ann = 0;

          for (auto const &leg : this->legs_per_vertex[v]) {
            int spin    = (spin_mask >> leg.line_index) & 1; // 0 = down, 1 = up
            uint8_t dir = spatial.directions[leg.line_index];

            // Site: for N_sites=1, always 0. For N_sites=2, determined by direction.
            uint8_t site;
            if constexpr (N_sites == 1) {
              site = 0;
            } else {
              // Right (dir=1): source site = 1, dest site = 0
              // Left  (dir=0): source site = 0, dest site = 1
              site = leg.is_source ? dir : (1 - dir);
            }

            uint8_t orbital = site + spin * N_sites;
            uint8_t op_id   = leg.is_source ? orbital : (ACTION_BIT | orbital);
            global.config.push_back(op_id);

            // Accumulate spin counts for conservation check
            if (leg.is_source) {
              if (spin == 1) up_ann++;
              else dn_ann++;
            } else {
              if (spin == 1) up_cre++;
              else dn_cre++;
            }
          }

          // Spin conservation: #up_creation == #up_annihilation (and same for down)
          if (up_cre != up_ann || dn_cre != dn_ann) { valid = false; }
        }

        if (!valid) continue;

        // --- Step 3: Canonicalize under the cumulant symmetry group ---
        auto orbit      = SymmetryGroup<N_sites, T>::get_orbit(global);
        auto &canonical = orbit[0]; // lex smallest
        int orbit_size  = (int)orbit.size();

        // Only count each canonical once per spatial config
        if (!seen_this_spatial.insert(canonical.config).second) continue;

        double weight = spatial.weight * orbit_size / auto_count;
        canonical_weights[canonical.config] += weight;
      }
    }

    // --- Step 4: Store results ---
    for (auto &[config, weight] : canonical_weights) {
      this->valid_configurations.push_back({config, weight});
    }
  }

  // =====================================================================
  // compute_diagram_sign: count the number of fermion loops in the diagram.
  // At each vertex, pair the i-th incoming line with the i-th outgoing line
  // (in hopping-line-index order) to build a successor map. Then count
  // cycles. Sign = (-1)^num_loops.
  //
  // This is stored as a member so evaluate() can use it.
  // =====================================================================
  template <int N_sites, typename T> void Diagram2<N_sites, T>::compute_diagram_sign() {

    int V       = this->graph.get_V();
    int n_lines = (int)this->hopping_lines.lines.size();

    // Build successor map: for each line, which line follows it through a vertex?
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
      // Pair i-th incoming with i-th outgoing
      for (size_t i = 0; i < incoming.size(); i++) { successor_map[incoming[i]] = outgoing[i]; }
    }

    // Count loops
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
  // evaluate: sum over all valid global configurations.
  // For each config, multiply vertex cumulants together, weight by the
  // config weight, and sum. Uses VertexType global cache when available.
  //
  // taus: one entry per hopping line (size = order of the diagram).
  //       taus[k] = imaginary time of the k-th hopping event.
  //
  // Returns: (-1/beta) * sum_configs weight * product_vertices C_v(taus, ops)
  // =====================================================================
  template <int N_sites, typename T>
  T Diagram2<N_sites, T>::evaluate(std::vector<double> const &taus, HubbardSolver<N_sites, T> const &solver, bool infinite_U) {

    int V = this->graph.get_V();
    T sum = T(0.0);

    for (auto const &gc : this->valid_configurations) {
      T product      = T(1.0);
      int cfg_offset = 0;

      for (int v = 0; v < V; v++) {
        int n_legs = (int)this->legs_per_vertex[v].size();

        // Extract local taus and op_ids for this vertex from the global config
        std::vector<double> local_taus(n_legs);
        std::vector<uint8_t> local_ops(n_legs);

        for (int i = 0; i < n_legs; i++) {
          local_taus[i] = taus[this->legs_per_vertex[v][i].line_index];
          local_ops[i]  = gc.config[cfg_offset + i];
        }
        cfg_offset += n_legs;

        // Split into annihilation (unprimed) and creation (primed), sort by time
        auto [unprimed, primed] = Args<N_sites, T>::split_from_raw(local_taus, local_ops);

        // Evaluate the cumulant — use VertexType cache if available
        T vertex_val;
        if (this->vertex_type_ptrs[v] != nullptr) {
          vertex_val = this->vertex_type_ptrs[v]->evaluate_canonical(unprimed, primed, solver, infinite_U);
        } else {
          vertex_val = compute_cumulant_decomposition(unprimed, primed, solver, infinite_U);
        }

        // Restore the diagram's leg convention: the canonical value was computed
        // with time-sorted operators; the permutation signs undo that sorting.
        vertex_val = vertex_val * T(unprimed.permutation_sign) * T(primed.permutation_sign);

        product = product * vertex_val;
      }

      sum = sum + T(gc.weight) * product;
    }

    // Free-energy prefactor: -1/beta * diagram_sign
    T prefactor = (T(-1.0) / solver.params.beta) * T(this->diagram_sign);
    return prefactor * sum;
  }

  // =====================================================================
  // Explicit template instantiations
  // =====================================================================
  template class Diagram2<1, double>;
  template class Diagram2<1, Dual>;
  template class Diagram2<2, double>;
  template class Diagram2<2, Dual>;

} // namespace sc_expansion
