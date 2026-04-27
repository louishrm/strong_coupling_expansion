#include "diagram.hpp"
#include "../dual.hpp"
#include "../fock_space.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <numeric>
#include <set>
#include <utility>

namespace sc_expansion::dimer {

  template <typename T>
  Diagram<T>::Diagram(Graph const &graph_, std::vector<VertexType<T> *> const &vertex_types) : graph(graph_) {
    this->hopping_lines = compute_hopping_lines(this->graph);
    this->compute_spatial_configurations();
    this->setup_vertices(vertex_types);
    this->compute_valid_configurations();
    this->diagram_sign = compute_diagram_sign(this->graph.get_V(), this->hopping_lines, this->legs_per_vertex);

    bool has_any_type = false;
    for (auto *p : this->vertex_type_ptrs) {
      if (p != nullptr) {
        has_any_type = true;
        break;
      }
    }
    if (has_any_type) this->build_local_state_tables();

    this->build_vertex_instances();
  }

  template <typename T>
  Diagram<T>::Diagram(Graph const &graph_, std::vector<VertexType<T> *> const &vertex_types,
                       std::vector<std::pair<int, int>> const &cluster_positions, int n_cluster_sites)
     : graph(graph_) {
    this->hopping_lines = compute_hopping_lines(this->graph);
    this->compute_spatial_configurations_cluster(cluster_positions, n_cluster_sites);
    this->setup_vertices(vertex_types);
    this->compute_valid_configurations();
    this->diagram_sign = compute_diagram_sign(this->graph.get_V(), this->hopping_lines, this->legs_per_vertex);

    bool has_any_type = false;
    for (auto *p : this->vertex_type_ptrs) {
      if (p != nullptr) {
        has_any_type = true;
        break;
      }
    }
    if (has_any_type) this->build_local_state_tables();

    this->build_vertex_instances();
  }

  // ---------------------------------------------------------------------------
  // Spatial embedding on the staggered (triangular) dimer superlattice.
  // ---------------------------------------------------------------------------

  // The staggered tiling alternates row offsets: dimer (u, v) covers
  // (2u + v%2, v) and (2u + v%2 + 1, v). Consequently, the set of 6 NN
  // dimer-coord offsets — and the (src_site, dst_site) carried by each bond —
  // depends on whether the *source* row is even or odd.
  //
  //   v even: { (+1,0), (0,+1), (0,-1) }  →  src=1, dst=0  (label 0)
  //           { (-1,0), (-1,+1), (-1,-1) } →  src=0, dst=1  (label 1)
  //   v odd:  { (+1,0), (+1,+1), (+1,-1) } →  src=1, dst=0  (label 0)
  //           { (-1,0), (0,+1), (0,-1) }   →  src=0, dst=1  (label 1)

  static inline bool v_is_odd(int v) { return (v & 1) != 0; }

  // 6 NN offsets in dimer-coord space, listed with the three src=1 directions
  // first (so callers that iterate this set can rely on the first 3 being
  // label 0 if needed).
  static inline std::array<std::pair<int, int>, 6> tri_offsets(int v) {
    if (v_is_odd(v)) {
      return {{{+1, 0}, {+1, +1}, {+1, -1}, {-1, 0}, {0, +1}, {0, -1}}};
    }
    return {{{+1, 0}, {0, +1}, {0, -1}, {-1, 0}, {-1, +1}, {-1, -1}}};
  }

  // Returns the binary bond label (0 or 1) for an offset (du, dv) from a
  // source at row v_source, or -1 if not a valid NN.
  static inline int tri_bond_label(int du, int dv, int v_source) {
    auto offsets = tri_offsets(v_source);
    bool ok      = false;
    for (auto const &o : offsets) {
      if (o.first == du && o.second == dv) {
        ok = true;
        break;
      }
    }
    if (!ok) return -1;

    if (du == 0) return v_is_odd(v_source) ? 1 : 0;
    return du > 0 ? 0 : 1;
  }

  static inline bool is_tri_neighbor(int x1, int y1, int x2, int y2) {
    int du = x2 - x1, dv = y2 - y1;
    auto offsets = tri_offsets(y1);
    for (auto const &o : offsets) {
      if (o.first == du && o.second == dv) return true;
    }
    return false;
  }

  template <typename T> void Diagram<T>::compute_spatial_configurations() {
    int V = this->graph.get_V();

    // Step 1: enumerate raw embeddings on the rectangular superlattice.
    std::map<std::vector<uint8_t>, int> raw_counts;

    std::vector<std::pair<int, int>> coords(V, {0, 0});
    std::vector<bool> placed(V, false);
    coords[0] = {0, 0};
    placed[0] = true;

    this->solve_dimer_embedding(1, placed, coords, raw_counts);

    // Step 2: graph automorphisms.
    std::vector<int> degrees(V, 0);
    for (int i = 0; i < V; ++i) {
      for (int j = 0; j < V; ++j) degrees[i] += this->graph(i, j) + this->graph(j, i);
    }

    std::vector<std::vector<int>> automorphisms;
    std::vector<int> perm(V);
    std::iota(perm.begin(), perm.end(), 0);

    do {
      bool degree_ok = true;
      for (int i = 0; i < V; ++i) {
        if (degrees[perm[i]] != degrees[i]) {
          degree_ok = false;
          break;
        }
      }
      if (!degree_ok) continue;

      bool is_auto = true;
      for (int i = 0; i < V && is_auto; ++i) {
        for (int j = 0; j < V && is_auto; ++j) {
          if (this->graph(i, j) != this->graph(perm[i], perm[j])) is_auto = false;
        }
      }

      if (is_auto) automorphisms.push_back(perm);
    } while (std::next_permutation(perm.begin(), perm.end()));

    // Step 3: canonicalise raw configs and merge.
    std::map<std::vector<uint8_t>, double> merged;
    for (auto &[dirs, count] : raw_counts) {
      auto canonical = this->canonicalize_directions(dirs, automorphisms);
      merged[canonical] += count;
    }

    for (auto &[dirs, weight] : merged) {
      this->spatial_configurations.push_back({dirs, weight});
    }
  }

  template <typename T>
  void Diagram<T>::compute_spatial_configurations_cluster(std::vector<std::pair<int, int>> const &cluster_positions, int n_cluster_sites) {
    int V = this->graph.get_V();

    std::map<std::vector<uint8_t>, int> raw_counts;
    for (auto const &start_pos : cluster_positions) {
      std::vector<std::pair<int, int>> coords(V, {0, 0});
      std::vector<bool> placed(V, false);
      coords[0] = start_pos;
      placed[0] = true;
      this->solve_cluster_embedding(1, placed, coords, raw_counts, cluster_positions);
    }

    std::vector<int> degrees(V, 0);
    for (int i = 0; i < V; ++i) {
      for (int j = 0; j < V; ++j) degrees[i] += this->graph(i, j) + this->graph(j, i);
    }

    std::vector<std::vector<int>> automorphisms;
    std::vector<int> perm(V);
    std::iota(perm.begin(), perm.end(), 0);

    do {
      bool degree_ok = true;
      for (int i = 0; i < V; ++i) {
        if (degrees[perm[i]] != degrees[i]) {
          degree_ok = false;
          break;
        }
      }
      if (!degree_ok) continue;

      bool is_auto = true;
      for (int i = 0; i < V && is_auto; ++i) {
        for (int j = 0; j < V && is_auto; ++j) {
          if (this->graph(i, j) != this->graph(perm[i], perm[j])) is_auto = false;
        }
      }
      if (is_auto) automorphisms.push_back(perm);
    } while (std::next_permutation(perm.begin(), perm.end()));

    std::map<std::vector<uint8_t>, double> merged;
    for (auto &[dirs, count] : raw_counts) {
      auto canonical = this->canonicalize_directions(dirs, automorphisms);
      merged[canonical] += count;
    }

    for (auto &[dirs, weight] : merged) {
      this->spatial_configurations.push_back({dirs, weight / (double)n_cluster_sites});
    }
  }

  template <typename T>
  void Diagram<T>::solve_cluster_embedding(int placed_count, std::vector<bool> &placed, std::vector<std::pair<int, int>> &coords,
                                            std::map<std::vector<uint8_t>, int> &config_counts,
                                            std::vector<std::pair<int, int>> const &cluster_positions) const {
    int V = this->graph.get_V();

    if (placed_count == V) {
      int n_lines = (int)this->hopping_lines.lines.size();
      std::vector<uint8_t> dirs(n_lines);

      for (int k = 0; k < n_lines; ++k) {
        int from = this->hopping_lines.lines[k].from_vertex;
        int to   = this->hopping_lines.lines[k].to_vertex;
        int ddx  = coords[to].first - coords[from].first;
        int ddy  = coords[to].second - coords[from].second;
        int lbl  = tri_bond_label(ddx, ddy, coords[from].second);
        if (lbl < 0) return;
        dirs[k] = (uint8_t)lbl;
      }
      config_counts[dirs]++;
      return;
    }

    int target = -1;
    for (int c = 0; c < V; ++c) {
      if (!placed[c]) {
        for (int p = 0; p < V; ++p) {
          if (placed[p]) {
            uint8_t links = this->graph(c, p) + this->graph(p, c);
            if (links > 0) {
              target = c;
              goto found_cluster_target;
            }
          }
        }
      }
    }
  found_cluster_target:;
    if (target == -1) return;

    for (auto const &pos : cluster_positions) {
      int cx = pos.first;
      int cy = pos.second;

      bool valid = true;
      for (int i = 0; i < V; ++i) {
        if (placed[i]) {
          uint8_t links = this->graph(target, i) + this->graph(i, target);
          if (links > 0 && !is_tri_neighbor(cx, cy, coords[i].first, coords[i].second)) {
            valid = false;
            break;
          }
        }
      }

      if (valid) {
        coords[target] = {cx, cy};
        placed[target] = true;
        this->solve_cluster_embedding(placed_count + 1, placed, coords, config_counts, cluster_positions);
        placed[target] = false;
      }
    }
  }

  // Staggered dimer tiling: dimer at superlattice (u, v) covers physical sites
  // (2u + v%2, v) and (2u + v%2 + 1, v). The superlattice is triangular with
  // 6 NN directions, each mapping uniquely to a (source_site, dest_site)
  // pair (see SpatialConfiguration in diagram.hpp for the bond-label table).
  template <typename T>
  void Diagram<T>::solve_dimer_embedding(int placed_count, std::vector<bool> &placed, std::vector<std::pair<int, int>> &coords,
                                          std::map<std::vector<uint8_t>, int> &config_counts) const {
    int V = this->graph.get_V();

    if (placed_count == V) {
      int n_lines = (int)this->hopping_lines.lines.size();
      std::vector<uint8_t> dirs(n_lines);

      for (int k = 0; k < n_lines; ++k) {
        int from = this->hopping_lines.lines[k].from_vertex;
        int to   = this->hopping_lines.lines[k].to_vertex;
        int ddx  = coords[to].first - coords[from].first;
        int ddy  = coords[to].second - coords[from].second;
        int lbl  = tri_bond_label(ddx, ddy, coords[from].second);
        if (lbl < 0) return; // not a valid NN pair; should not happen if placement is correct
        dirs[k] = (uint8_t)lbl;
      }
      config_counts[dirs]++;
      return;
    }

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
    if (target == -1) return;

    int ax       = coords[anchor].first;
    int ay       = coords[anchor].second;
    auto offsets = tri_offsets(ay);

    for (auto const &o : offsets) {
      int cx = ax + o.first;
      int cy = ay + o.second;

      bool valid = true;
      for (int i = 0; i < V; ++i) {
        if (placed[i]) {
          uint8_t links = this->graph(target, i) + this->graph(i, target);
          if (links > 0 && !is_tri_neighbor(cx, cy, coords[i].first, coords[i].second)) {
            valid = false;
            break;
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

  template <typename T>
  std::vector<uint8_t> Diagram<T>::apply_automorphism_to_directions(std::vector<uint8_t> const &dirs, std::vector<int> const &perm) const {
    int V       = this->graph.get_V();
    int n_lines = (int)dirs.size();

    std::vector<int> line_from(n_lines), line_to(n_lines), line_mult(n_lines);
    std::vector<std::vector<int>> start_idx(V, std::vector<int>(V, -1));
    {
      int idx = 0;
      for (int i = 0; i < V; ++i) {
        for (int j = 0; j < V; ++j) {
          int count = this->graph(i, j);
          if (count > 0 && start_idx[i][j] == -1) start_idx[i][j] = idx;
          for (int k = 0; k < count; ++k) {
            line_from[idx] = i;
            line_to[idx]   = j;
            line_mult[idx] = k;
            idx++;
          }
        }
      }
    }

    std::vector<uint8_t> result(n_lines);
    for (int l = 0; l < n_lines; ++l) {
      int new_from    = perm[line_from[l]];
      int new_to      = perm[line_to[l]];
      int new_idx     = start_idx[new_from][new_to] + line_mult[l];
      result[new_idx] = dirs[l];
    }
    return result;
  }

  // Lex-min over all graph automorphisms × {identity, lattice inversion}.
  // Lattice inversion (2-fold rotation about a dimer center) swaps site 0
  // <-> site 1 within each dimer, mapping bond labels 0 <-> 1.
  template <typename T>
  std::vector<uint8_t> Diagram<T>::canonicalize_directions(std::vector<uint8_t> const &dirs,
                                                            std::vector<std::vector<int>> const &automorphisms) const {
    static constexpr uint8_t inversion_map[2] = {1, 0};

    auto min_dirs = dirs;
    for (auto const &perm : automorphisms) {
      auto permuted = this->apply_automorphism_to_directions(dirs, perm);
      if (permuted < min_dirs) min_dirs = permuted;

      auto inverted = permuted;
      for (auto &d : inverted) d = inversion_map[d];
      if (inverted < min_dirs) min_dirs = inverted;
    }
    return min_dirs;
  }

  // ---------------------------------------------------------------------------
  // Vertex setup, valid-config enumeration, instances.
  // ---------------------------------------------------------------------------

  template <typename T> void Diagram<T>::setup_vertices(std::vector<VertexType<T> *> const &vertex_types) {
    int V = this->graph.get_V();
    this->vertex_type_ptrs.resize(V, nullptr);

    for (int v = 0; v < V; ++v) {
      int degree = 0;
      for (int j = 0; j < V; ++j) degree += this->graph(v, j) + this->graph(j, v);
      int cumulant_order = degree / 2;
      int vt_idx         = cumulant_order - 1;
      if (vt_idx >= 0 && vt_idx < (int)vertex_types.size()) this->vertex_type_ptrs[v] = vertex_types[vt_idx];
    }
  }

  template <typename T> void Diagram<T>::compute_valid_configurations() {
    int V             = this->graph.get_V();
    int n_lines       = (int)this->hopping_lines.lines.size();
    double sym_factor = this->graph.get_symmetry_factor();

    constexpr uint8_t ACTION_BIT = FermionOperator<2, T>::ACTION_BIT;

    this->legs_per_vertex = compute_legs_per_vertex(this->graph, this->hopping_lines);

    using Config = GlobalConfiguration<2>;
    std::map<std::vector<uint8_t>, double> canonical_weights;
    int n_spin_configs = 1 << n_lines;

    for (auto const &spatial : this->spatial_configurations) {
      std::set<std::vector<uint8_t>> seen_this_spatial;

      for (int spin_mask = 0; spin_mask < n_spin_configs; ++spin_mask) {
        Config global;
        global.config.reserve(2 * n_lines);
        bool valid = true;

        for (int v = 0; v < V && valid; ++v) {
          int up_cre = 0, up_ann = 0, dn_cre = 0, dn_ann = 0;

          for (auto const &leg : this->legs_per_vertex[v]) {
            int spin    = (spin_mask >> leg.line_index) & 1; // 0 = down, 1 = up
            uint8_t dir = spatial.directions[leg.line_index];

            // Binary bond label → site assignment for the staggered (triangular)
            // dimer tiling. All NN bonds cross the A/B sublattice, so the local
            // cumulant depends only on (src, dst), not on the geometric direction.
            static constexpr uint8_t source_site[2] = {1, 0};
            static constexpr uint8_t dest_site[2]   = {0, 1};
            uint8_t site                            = leg.is_source ? source_site[dir] : dest_site[dir];

            // Orbital index: low bit = site, high bit = spin (matches FermionOperator
            // convention "low N_sites bits = down block, next N_sites bits = up block").
            uint8_t orbital = site + spin * 2;
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

          if (up_cre != up_ann || dn_cre != dn_ann) valid = false;
        }

        if (!valid) continue;

        auto orbit      = SymmetryGroup<2, T>::get_orbit(global);
        auto &canonical = orbit[0];
        int orbit_size  = (int)orbit.size();

        if (!seen_this_spatial.insert(canonical.config).second) continue;

        double weight = spatial.weight * orbit_size / sym_factor;
        canonical_weights[canonical.config] += weight;
      }
    }

    for (auto &[config, weight] : canonical_weights) this->valid_configurations.push_back({config, weight});
  }

  template <typename T> void Diagram<T>::build_vertex_instances() {
    int V       = this->graph.get_V();
    int n_lines = (int)this->hopping_lines.lines.size();

    bool has_any_type = false;
    for (int v = 0; v < V; ++v) {
      if (this->vertex_type_ptrs[v] != nullptr) {
        has_any_type = true;
        break;
      }
    }
    if (!has_any_type) return;

    this->tau_to_vertices.resize(n_lines);
    for (int v = 0; v < V; ++v) {
      for (auto const &leg : this->legs_per_vertex[v]) this->tau_to_vertices[leg.line_index].push_back(v);
    }
    for (auto &vlist : this->tau_to_vertices) {
      std::sort(vlist.begin(), vlist.end());
      vlist.erase(std::unique(vlist.begin(), vlist.end()), vlist.end());
    }

    this->vertex_instances.resize(this->valid_configurations.size());
    for (int gc_idx = 0; gc_idx < (int)this->valid_configurations.size(); ++gc_idx) {
      auto const &gc = this->valid_configurations[gc_idx];
      int cfg_offset = 0;

      for (int v = 0; v < V; ++v) {
        int n_legs = (int)this->legs_per_vertex[v].size();

        std::vector<int> tau_indices(n_legs);
        for (int i = 0; i < n_legs; ++i) tau_indices[i] = this->legs_per_vertex[v][i].line_index;

        std::vector<uint8_t> op_ids(gc.config.begin() + cfg_offset, gc.config.begin() + cfg_offset + n_legs);
        cfg_offset += n_legs;

        this->vertex_instances[gc_idx].emplace_back(this->vertex_type_ptrs[v], std::move(tau_indices), std::move(op_ids));
      }
    }
  }

  template <typename T> void Diagram<T>::build_local_state_tables() {
    int V         = this->graph.get_V();
    int n_configs = (int)this->valid_configurations.size();

    this->local_states.resize(V);
    this->config_to_local.resize(n_configs, std::vector<int>(V));

    std::vector<int> offsets(V);
    for (int v = 1; v < V; ++v) offsets[v] = offsets[v - 1] + (int)this->legs_per_vertex[v - 1].size();

    for (int v = 0; v < V; ++v) {
      std::map<std::vector<uint8_t>, int> state_map;
      int n_legs = (int)this->legs_per_vertex[v].size();

      for (int gc_idx = 0; gc_idx < n_configs; ++gc_idx) {
        std::vector<uint8_t> ops(this->valid_configurations[gc_idx].config.begin() + offsets[v],
                                 this->valid_configurations[gc_idx].config.begin() + offsets[v] + n_legs);

        auto [it, inserted] = state_map.try_emplace(ops, (int)this->local_states[v].size());
        if (inserted) this->local_states[v].push_back(ops);
        this->config_to_local[gc_idx][v] = it->second;
      }
    }

    this->local_values.resize(V);
    for (int v = 0; v < V; ++v) {
      this->local_values[v].resize(this->local_states[v].size(), T(0.0));
    }

    this->vertex_dirty_finite.assign(V, true);
  }

  template <typename T> void Diagram<T>::mark_tau_dirty(int tau_index) {
    if (!this->local_states.empty()) {
      for (int v : this->tau_to_vertices[tau_index]) {
        this->vertex_dirty_finite[v] = true;
      }
      return;
    }
    for (int v : this->tau_to_vertices[tau_index]) {
      for (auto &config_instances : this->vertex_instances) config_instances[v].mark_dirty();
    }
  }

  template <typename T> void Diagram<T>::mark_all_dirty() {
    if (!this->local_states.empty()) {
      std::fill(this->vertex_dirty_finite.begin(), this->vertex_dirty_finite.end(), true);
      return;
    }
    for (auto &config_instances : this->vertex_instances) {
      for (auto &vi : config_instances) vi.mark_dirty();
    }
  }

  // ---------------------------------------------------------------------------
  // Plan building (factored path).
  //
  // CumulantPlan is τ-independent and U-independent.
  // ---------------------------------------------------------------------------

  template <typename T> void Diagram<T>::build_local_plans(HubbardSolver<2, T> const & /*solver*/) {
    int V = this->graph.get_V();
    this->local_plans_finite.assign(V, {});

    for (int v = 0; v < V; ++v) {
      int n_legs = (int)this->legs_per_vertex[v].size();
      std::vector<double> dummy_taus(n_legs, 0.5);

      int ns = (int)this->local_states[v].size();
      this->local_plans_finite[v].reserve(ns);

      for (int s = 0; s < ns; ++s) {
        auto [u, p] = Args<2, T>::split_from_raw(dummy_taus, this->local_states[v][s]);

        CumulantPlan plan;
        CumulantSolver<2, T> builder(u, p);
        builder.record_plan(plan);
        this->local_plans_finite[v].push_back(std::move(plan));
      }
    }

    this->local_plans_built = true;
  }

  template <typename T>
  T Diagram<T>::evaluate_factored(std::vector<double> const &taus, HubbardSolver<2, T> const &solver) {
    if (!this->local_plans_built) this->build_local_plans(solver);

    int V        = this->graph.get_V();
    auto &dirty  = this->vertex_dirty_finite;
    auto &values = this->local_values;
    auto &plans  = this->local_plans_finite;

    auto t_p1 = std::chrono::high_resolution_clock::now();
    for (int v = 0; v < V; ++v) {
      if (!dirty[v]) {
        this->local_cache_hits++;
        continue;
      }
      this->local_cache_misses++;

      int n_legs = (int)this->legs_per_vertex[v].size();
      std::vector<double> local_taus(n_legs);
      for (int i = 0; i < n_legs; ++i) local_taus[i] = taus[this->legs_per_vertex[v][i].line_index];

      for (int s = 0; s < (int)this->local_states[v].size(); ++s) {
        auto [unprimed, primed] = Args<2, T>::split_from_raw(local_taus, this->local_states[v][s]);
        T val                   = evaluate_plan(plans[v][s], unprimed, primed, solver, /*infinite_U=*/false);
        values[v][s]            = val * T(unprimed.permutation_sign) * T(primed.permutation_sign);
      }

      dirty[v] = false;
    }
    auto t_p2 = std::chrono::high_resolution_clock::now();

    T sum = T(0.0);
    for (int gc_idx = 0; gc_idx < (int)this->valid_configurations.size(); ++gc_idx) {
      T product = T(1.0);
      for (int v = 0; v < V; ++v) product = product * values[v][this->config_to_local[gc_idx][v]];
      sum = sum + T(this->valid_configurations[gc_idx].weight) * product;
    }
    auto t_end = std::chrono::high_resolution_clock::now();

    this->phase1_seconds += std::chrono::duration<double>(t_p2 - t_p1).count();
    this->phase2_seconds += std::chrono::duration<double>(t_end - t_p2).count();

    // Absorb the (-t)^n factor that each hopping line contributes, so the
    // returned coefficient is the n-th order term in a series in t (not -t).
    int n_lines    = (int)this->hopping_lines.lines.size();
    double t_sign  = (n_lines % 2 == 0) ? 1.0 : -1.0;
    T prefactor    = (T(-1.0) / solver.params.beta) * T(this->diagram_sign * (int)t_sign);
    return prefactor * sum;
  }

  template <typename T>
  T Diagram<T>::evaluate(std::vector<double> const &taus, HubbardSolver<2, T> const &solver) {
    if (!this->local_states.empty()) return this->evaluate_factored(taus, solver);

    int V = this->graph.get_V();
    T sum = T(0.0);

    if (this->vertex_instances.empty()) {
      // No factored tables and no VertexInstance fallback — caller forgot to
      // pass VertexTypes. Return zero rather than silently miscomputing.
      return sum;
    }

    for (int gc_idx = 0; gc_idx < (int)this->valid_configurations.size(); ++gc_idx) {
      T product = T(1.0);
      for (int v = 0; v < V; ++v) product = product * this->vertex_instances[gc_idx][v].get_value(taus, solver);
      sum = sum + T(this->valid_configurations[gc_idx].weight) * product;
    }

    // Absorb the (-t)^n factor that each hopping line contributes, so the
    // returned coefficient is the n-th order term in a series in t (not -t).
    int n_lines    = (int)this->hopping_lines.lines.size();
    double t_sign  = (n_lines % 2 == 0) ? 1.0 : -1.0;
    T prefactor    = (T(-1.0) / solver.params.beta) * T(this->diagram_sign * (int)t_sign);
    return prefactor * sum;
  }

  template <typename T>
  std::vector<T> Diagram<T>::evaluate_per_config(std::vector<double> const &taus, HubbardSolver<2, T> const &solver) {
    if (!this->local_plans_built) this->build_local_plans(solver);

    int V = this->graph.get_V();
    this->mark_all_dirty();
    auto &values = this->local_values;
    auto &plans  = this->local_plans_finite;
    auto &dirty  = this->vertex_dirty_finite;

    for (int v = 0; v < V; ++v) {
      int n_legs = (int)this->legs_per_vertex[v].size();
      std::vector<double> local_taus(n_legs);
      for (int i = 0; i < n_legs; ++i) local_taus[i] = taus[this->legs_per_vertex[v][i].line_index];
      for (int s = 0; s < (int)this->local_states[v].size(); ++s) {
        auto [unprimed, primed] = Args<2, T>::split_from_raw(local_taus, this->local_states[v][s]);
        T val                   = evaluate_plan(plans[v][s], unprimed, primed, solver, /*infinite_U=*/false);
        values[v][s]            = val * T(unprimed.permutation_sign) * T(primed.permutation_sign);
      }
      dirty[v] = false;
    }

    std::vector<T> contribs;
    contribs.reserve(this->valid_configurations.size());
    T sign = T((double)this->diagram_sign);
    for (int gc_idx = 0; gc_idx < (int)this->valid_configurations.size(); ++gc_idx) {
      T product = T(1.0);
      for (int v = 0; v < V; ++v) product = product * values[v][this->config_to_local[gc_idx][v]];
      contribs.push_back(sign * T(this->valid_configurations[gc_idx].weight) * product);
    }
    return contribs;
  }

  template class Diagram<double>;
  template class Diagram<Dual>;

} // namespace sc_expansion::dimer
