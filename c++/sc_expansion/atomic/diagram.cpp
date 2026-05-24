#include "diagram.hpp"
#include "../dual.hpp"
#include "../fock_space.hpp"
#include <algorithm>
#include <cstdlib>
#include <functional>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>

namespace sc_expansion::atomic {

  template <typename T>
  Diagram<T>::Diagram(Graph const &graph_, std::vector<VertexType<T> *> const &vertex_types) : graph(graph_) {
    this->hopping_lines    = compute_hopping_lines(this->graph);
    this->legs_per_vertex  = compute_legs_per_vertex(this->graph, this->hopping_lines);
    this->setup_vertices(vertex_types);
    this->compute_valid_configurations();
    this->diagram_sign = compute_diagram_sign(this->graph.get_V(), this->hopping_lines, this->legs_per_vertex);
    this->build_vertex_instances();
  }

  template <typename T>
  Diagram<T>::Diagram(Graph const &graph_, std::vector<VertexType<T> *> const &vertex_types,
                     std::vector<int> marks_, std::vector<int> mark_spins_,
                     bool flip_mark_order_,
                     MarkEncoding mark_encoding_)
     : graph(graph_), is_rooted(true), flip_mark_order(flip_mark_order_),
       mark_encoding(mark_encoding_),
       marks(std::move(marks_)), mark_spins(std::move(mark_spins_)) {

    if (this->marks.size() != this->mark_spins.size())
      throw std::invalid_argument("Diagram(rooted): marks and mark_spins size mismatch");
    if (this->marks.empty() || this->marks.size() > 2)
      throw std::invalid_argument("Diagram(rooted): marks size must be 1 or 2");
    for (int s : this->mark_spins)
      if (s != 0 && s != 1) throw std::invalid_argument("Diagram(rooted): mark_spins entries must be 0 or 1");
    for (int m : this->marks)
      if (m < 0 || m >= this->graph.get_V()) throw std::invalid_argument("Diagram(rooted): mark index out of range");

    this->hopping_lines   = compute_hopping_lines(this->graph);
    this->legs_per_vertex = compute_legs_per_vertex(this->graph, this->hopping_lines);
    this->setup_vertices(vertex_types);
    this->compute_valid_configurations();
    this->diagram_sign = compute_diagram_sign(this->graph.get_V(), this->hopping_lines, this->legs_per_vertex);
    this->build_vertex_instances();
  }

  template <typename T> void Diagram<T>::setup_vertices(std::vector<VertexType<T> *> const &vertex_types) {
    int V = this->graph.get_V();
    this->vertex_type_ptrs.resize(V, nullptr);

    // For rooted diagrams, each mark on a vertex contributes one extra
    // (c†_σ, c_σ) pair — i.e. one extra cumulant order. A vertex with two
    // coincident marks (same-site density-density insertion) gets +2.
    std::vector<int> mark_bonus(V, 0);
    if (this->is_rooted)
      for (int m : this->marks) mark_bonus[m] += 1;

    for (int v = 0; v < V; ++v) {
      int degree = 0;
      for (int j = 0; j < V; ++j) { degree += this->graph(v, j) + this->graph(j, v); }
      int cumulant_order = degree / 2 + mark_bonus[v];
      int vt_idx         = cumulant_order - 1;
      if (vt_idx >= 0 && vt_idx < (int)vertex_types.size()) { this->vertex_type_ptrs[v] = vertex_types[vt_idx]; }
    }
  }

  template <typename T> void Diagram<T>::compute_valid_configurations() {
    int V             = this->graph.get_V();
    int n_lines       = (int)this->hopping_lines.lines.size();
    double sym_factor = this->graph.get_symmetry_factor();
    double free_mult  = this->graph.get_free_multiplicity();

    constexpr uint8_t ACTION_BIT = FermionOperator<1, T>::ACTION_BIT;

    // For rooted diagrams, each mark on a vertex contributes one extra
    // (c†_σ, c_σ) pair appended after the hopping legs (order matches
    // build_vertex_instances below). A vertex with two coincident marks
    // carries TWO such pairs (one per mark, each with its own spin). The
    // mark spin σ is fixed by mark_spins; appending each pair gives +1 to
    // both creation and annihilation of σ, conserving balance trivially.
    std::vector<std::vector<int>> mark_idxs_at(V);
    if (this->is_rooted)
      for (size_t i = 0; i < this->marks.size(); ++i) mark_idxs_at[this->marks[i]].push_back((int)i);

    using Config = GlobalConfiguration<1>;
    std::map<std::vector<uint8_t>, double> canonical_weights;
    int n_spin_configs = 1 << n_lines;
    std::set<std::vector<uint8_t>> seen;

    for (int spin_mask = 0; spin_mask < n_spin_configs; ++spin_mask) {
      Config global;
      global.config.reserve(2 * n_lines + 2 * this->marks.size());
      bool valid = true;

      for (int v = 0; v < V && valid; ++v) {
        int up_cre = 0, up_ann = 0, dn_cre = 0, dn_ann = 0;

        for (auto const &leg : this->legs_per_vertex[v]) {
          int spin        = (spin_mask >> leg.line_index) & 1; // 0 = down, 1 = up
          uint8_t orbital = spin;                              // single site: orbital == spin bit
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

        // BlockConstraint encoding only: append one mark insertion (c†_σ, c_σ)
        // per mark on this vertex. With descending-τ stable sort and both
        // ops at τ=0, the input order chosen here pins the equal-time
        // convention: default (c†, c) puts c†_σ at the later effective time
        // → reads as n_σ; flipped (c, c†) puts c at the later effective
        // time → reads as 1 − n_σ.
        //
        // StaticDensity encoding skips this block entirely — the density
        // is attached to the per-vertex cumulant leaf via add_static_density
        // in build_vertex_instances, not as an operator in op_ids. The
        // creation/annihilation balance is unaffected either way (n_σ
        // contributes +1 to both, so it never breaks the constraint).
        if (this->mark_encoding == MarkEncoding::BlockConstraint) {
          for (int mi : mark_idxs_at[v]) {
            int mark_spin   = this->mark_spins[mi];
            uint8_t orbital = (uint8_t)mark_spin;
            if (this->flip_mark_order) {
              global.config.push_back(orbital);                          // c_σ
              global.config.push_back((uint8_t)(ACTION_BIT | orbital)); // c†_σ
            } else {
              global.config.push_back((uint8_t)(ACTION_BIT | orbital)); // c†_σ
              global.config.push_back(orbital);                          // c_σ
            }
            if (mark_spin == 1) { up_cre++; up_ann++; }
            else { dn_cre++; dn_ann++; }
          }
        }

        if (up_cre != up_ann || dn_cre != dn_ann) valid = false;
      }

      if (!valid) continue;

      if (this->is_rooted) {
        // Marks fix the spin assignment of the external operators, breaking
        // the global spin-flip symmetry used for orbit folding in the
        // free-energy case. Each internal-spin assignment contributes once.
        if (!seen.insert(global.config).second) continue;
        // Rooted weight = free_mult · n_mark_orbit / sym_factor.
        // n_mark_orbit is the size of the S_m orbit of the mark→vertex
        // assignment: m! divided by the per-vertex coincidence stabilizer
        // ∏_v n_marks_at[v]! (marks on the same vertex are physically
        // identical, so permuting their labels gives the same labeled
        // embedding and must not be counted).
        double n_mark_orbit = 1.0;
        for (int k = 2; k <= (int)this->marks.size(); ++k) n_mark_orbit *= (double)k;
        std::map<int, int> n_marks_at_count;
        for (int mv : this->marks) ++n_marks_at_count[mv];
        for (auto const &[v, cnt] : n_marks_at_count) {
          for (int k = 2; k <= cnt; ++k) n_mark_orbit /= (double)k;
        }
        double weight = free_mult * n_mark_orbit / sym_factor;
        canonical_weights[global.config] += weight;
      } else {
        auto orbit      = SymmetryGroup<1, T>::get_orbit(global);
        auto &canonical = orbit[0];
        int orbit_size  = (int)orbit.size();

        if (!seen.insert(canonical.config).second) continue;

        double weight = free_mult * orbit_size / sym_factor;
        canonical_weights[canonical.config] += weight;
      }
    }

    for (auto &[config, weight] : canonical_weights) { this->valid_configurations.push_back({config, weight}); }
  }

  template <typename T> void Diagram<T>::build_vertex_instances() {
    int V       = this->graph.get_V();
    int n_lines = (int)this->hopping_lines.lines.size();
    int pinned_tau_idx = n_lines; // reserved slot for τ=0 mark insertions

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
      for (auto const &leg : this->legs_per_vertex[v]) { this->tau_to_vertices[leg.line_index].push_back(v); }
    }
    for (auto &vlist : this->tau_to_vertices) {
      std::sort(vlist.begin(), vlist.end());
      vlist.erase(std::unique(vlist.begin(), vlist.end()), vlist.end());
    }

    std::vector<int> n_marks_at(V, 0);
    if (this->is_rooted)
      for (int m : this->marks) ++n_marks_at[m];

    constexpr uint8_t ACTION_BIT = FermionOperator<1, T>::ACTION_BIT;

    this->vertex_instances.resize(this->valid_configurations.size());
    for (int gc_idx = 0; gc_idx < (int)this->valid_configurations.size(); ++gc_idx) {
      auto const &gc = this->valid_configurations[gc_idx];
      int cfg_offset = 0;

      for (int v = 0; v < V; ++v) {
        int n_hop_legs = (int)this->legs_per_vertex[v].size();
        // StaticDensity encoding skips the appended mark pairs entirely;
        // BlockConstraint reserves 2 op slots per mark on this vertex.
        int n_mark_op_slots = (this->mark_encoding == MarkEncoding::BlockConstraint) ? 2 * n_marks_at[v] : 0;
        int n_legs          = n_hop_legs + n_mark_op_slots;

        std::vector<int> tau_indices(n_legs);
        for (int i = 0; i < n_hop_legs; ++i) tau_indices[i] = this->legs_per_vertex[v][i].line_index;
        for (int i = 0; i < n_mark_op_slots; ++i) tau_indices[n_hop_legs + i] = pinned_tau_idx;

        std::vector<uint8_t> op_ids(gc.config.begin() + cfg_offset, gc.config.begin() + cfg_offset + n_legs);
        cfg_offset += n_legs;

        this->vertex_instances[gc_idx].emplace_back(this->vertex_type_ptrs[v], std::move(tau_indices), std::move(op_ids));

        if (n_marks_at[v] > 0) {
          if (this->mark_encoding == MarkEncoding::BlockConstraint) {
            // Pin each mark's (c†_σ, c_σ) pair into its own partition block.
            // The mark ops occupy the last 2·n_marks_at[v] entries of op_ids
            // in (c†, c) order (flip_mark_order swaps inside each pair).
            // Input positions in the unprimed / primed lists (per
            // Args::split_from_raw) follow the relative order of ann/cre ops
            // in op_ids — so the first mark's ann lands at n_u_before and
            // its cre at n_p_before, the next mark's at n_u_before+1 /
            // n_p_before+1, etc.
            int n_u_before = 0, n_p_before = 0;
            int op_base = cfg_offset - n_legs;
            for (int k = 0; k < n_hop_legs; ++k) {
              bool is_creator = (gc.config[op_base + k] & ACTION_BIT) != 0;
              if (is_creator) ++n_p_before; else ++n_u_before;
            }
            for (int mi = 0; mi < n_marks_at[v]; ++mi) {
              this->vertex_instances[gc_idx].back().add_block_constraint(n_u_before + mi, n_p_before + mi);
            }

            // Coincidence groups: marks on this vertex sharing the same
            // spin collapse via n_σ^k = n_σ. Emit one group per spin with
            // >= 2 marks. Block indices passed are the per-vertex-mark
            // positions [0..n_marks_at[v]); map through mark_idxs_at[v]
            // to recover the global mark index → spin.
            std::map<int, std::vector<int>> by_spin;
            for (int mi = 0; mi < n_marks_at[v]; ++mi) {
              int global_mark_idx = -1;
              int seen            = 0;
              for (int gi = 0; gi < (int)this->marks.size(); ++gi) {
                if (this->marks[gi] == v) {
                  if (seen == mi) {
                    global_mark_idx = gi;
                    break;
                  }
                  ++seen;
                }
              }
              int spin = this->mark_spins[global_mark_idx];
              by_spin[spin].push_back(mi);
            }
            for (auto const &[spin, blk_idxs] : by_spin) {
              if (blk_idxs.size() < 2) continue;
              this->vertex_instances[gc_idx].back().add_coincidence_group(spin, blk_idxs);
            }
          } else {
            // StaticDensity encoding: register one add_static_density(spin)
            // per mark on this vertex; the (c†, c) pair never enters
            // op_ids. The cumulant plan treats each density as a bosonic
            // atom routed across partition blocks independently.
            for (size_t mi = 0; mi < this->marks.size(); ++mi) {
              if (this->marks[mi] == v) { this->vertex_instances[gc_idx].back().add_static_density(this->mark_spins[mi]); }
            }
          }
        }
      }
    }
  }

  template <typename T> void Diagram<T>::mark_tau_dirty(int tau_index) {
    for (int v : this->tau_to_vertices[tau_index]) {
      for (auto &config_instances : this->vertex_instances) { config_instances[v].mark_dirty(); }
    }
  }

  template <typename T> void Diagram<T>::mark_all_dirty() {
    for (auto &config_instances : this->vertex_instances) {
      for (auto &vi : config_instances) vi.mark_dirty();
    }
  }

  template <typename T>
  T Diagram<T>::evaluate(std::vector<double> const &taus, HubbardSolver<1, T> const &solver, bool infinite_U) {
    int V = this->graph.get_V();
    T sum = T(0.0);

    for (int gc_idx = 0; gc_idx < (int)this->valid_configurations.size(); ++gc_idx) {
      T product = T(1.0);
      for (int v = 0; v < V; ++v) { product = product * this->vertex_instances[gc_idx][v].get_value(taus, solver, infinite_U); }
      sum = sum + T(this->valid_configurations[gc_idx].weight) * product;
    }

    // Free-energy prefactor: -1/β · sign (comes from F = -log Z / β).
    // Correlator (rooted) prefactor: just the diagram sign — the -1/β factor
    // is free-energy-specific. Any additional convention (e.g. an overall
    // sign from external operator ordering) will be calibrated against the
    // order-2 ED reference; for now we emit the bare sign × Σ.
    // Rooted (correlator) prefactor: diagram_sign. Density marks n_σ(0) =
    // c†_σ c_σ at coincident time/site don't add new closed fermion loops
    // on their own — they splice into a loop of the underlying hopping
    // skeleton — so the (−1)^L parity of the rooted diagram equals that of
    // the corresponding vacuum graph. At order 2 the digon has 1 loop
    // (diagram_sign = −1), recovering the previously hardcoded −1. At
    // order 4 different topologies carry different loop counts (e.g. the
    // (2,2) multi-digon has 2 loops → +1), and the relative signs are
    // needed for the cumulant cancellation.
    // Free-energy prefactor stays −sign/β.
    T prefactor = this->is_rooted ? T(this->diagram_sign) : (T(-1.0) / solver.params.beta) * T(this->diagram_sign);
    return prefactor * sum;
  }

  template <typename T>
  std::vector<T> Diagram<T>::evaluate_per_config(std::vector<double> const &taus, HubbardSolver<1, T> const &solver, bool infinite_U) {
    int V = this->graph.get_V();
    int n_cfg = (int)this->valid_configurations.size();
    std::vector<T> out(n_cfg, T(0.0));

    // Rooted (correlator) prefactor: diagram_sign. Density marks n_σ(0) =
    // c†_σ c_σ at coincident time/site don't add new closed fermion loops
    // on their own — they splice into a loop of the underlying hopping
    // skeleton — so the (−1)^L parity of the rooted diagram equals that of
    // the corresponding vacuum graph. At order 2 the digon has 1 loop
    // (diagram_sign = −1), recovering the previously hardcoded −1. At
    // order 4 different topologies carry different loop counts (e.g. the
    // (2,2) multi-digon has 2 loops → +1), and the relative signs are
    // needed for the cumulant cancellation.
    // Free-energy prefactor stays −sign/β.
    T prefactor = this->is_rooted ? T(this->diagram_sign) : (T(-1.0) / solver.params.beta) * T(this->diagram_sign);

    for (int gc_idx = 0; gc_idx < n_cfg; ++gc_idx) {
      T product = T(1.0);
      for (int v = 0; v < V; ++v) {
        product = product * this->vertex_instances[gc_idx][v].get_value(taus, solver, infinite_U);
      }
      out[gc_idx] = prefactor * T(this->valid_configurations[gc_idx].weight) * product;
    }
    return out;
  }

  template class Diagram<double>;
  template class Diagram<Dual>;

  int count_lattice_embeddings(Graph const &graph, std::vector<int> const &marks, std::vector<int> const &r) {
    if (marks.size() != 2 || r.size() != 2) return 0;
    int V     = graph.get_V();
    int mark0 = marks[0];
    int mark1 = marks[1];
    if (mark0 < 0 || mark0 >= V || mark1 < 0 || mark1 >= V) return 0;

    // Simple-graph adjacency: multi-edges collapse to one constraint (L¹=1
    // between the two endpoint sites), self-loops are skipped (they would
    // impose L¹=0 on a vertex with itself, which can never satisfy L¹=1).
    std::vector<std::vector<int>> adj(V);
    for (int i = 0; i < V; ++i) {
      for (int j = i + 1; j < V; ++j) {
        if (graph(i, j) != 0 || graph(j, i) != 0) {
          adj[i].push_back(j);
          adj[j].push_back(i);
        }
      }
    }

    // BFS from mark0: guarantees every later-placed vertex has ≥1 already-
    // placed neighbor we can use as a lattice anchor (≤4 candidate sites).
    std::vector<int> order;
    order.reserve(V);
    std::vector<char> in_queue(V, 0);
    std::queue<int> q;
    q.push(mark0);
    in_queue[mark0] = 1;
    while (!q.empty()) {
      int u = q.front();
      q.pop();
      order.push_back(u);
      for (int v : adj[u])
        if (!in_queue[v]) {
          in_queue[v] = 1;
          q.push(v);
        }
    }
    if ((int)order.size() != V) return 0; // disconnected — vacuum diagrams shouldn't be.

    std::vector<std::pair<int, int>> pos(V, {0, 0});
    std::vector<char> placed(V, 0);
    pos[mark0]    = {0, 0};
    placed[mark0] = 1;
    if (mark1 != mark0) {
      pos[mark1]    = {r[0], r[1]};
      placed[mark1] = 1;
    } else if (r[0] != 0 || r[1] != 0) {
      return 0; // coincident marks force r=(0,0).
    }

    auto verify = [&](int u, int x, int y) {
      for (int v : adj[u]) {
        if (placed[v]) {
          int d = std::abs(pos[v].first - x) + std::abs(pos[v].second - y);
          if (d != 1) return false;
        }
      }
      return true;
    };

    constexpr int dx[4] = {1, -1, 0, 0};
    constexpr int dy[4] = {0, 0, 1, -1};

    std::function<int(int)> dfs = [&](int idx) -> int {
      // Skip past pre-placed vertices (mark0, mark1), verifying constraints
      // against whatever is currently placed.
      while (idx < (int)order.size() && placed[order[idx]]) {
        int u = order[idx];
        if (!verify(u, pos[u].first, pos[u].second)) return 0;
        ++idx;
      }
      if (idx == (int)order.size()) return 1;

      int u      = order[idx];
      int anchor = -1;
      for (int v : adj[u])
        if (placed[v]) {
          anchor = v;
          break;
        }
      if (anchor == -1) return 0; // unreachable in a connected graph + BFS order.

      int total = 0;
      for (int d = 0; d < 4; ++d) {
        int nx = pos[anchor].first + dx[d];
        int ny = pos[anchor].second + dy[d];
        if (!verify(u, nx, ny)) continue;
        pos[u]    = {nx, ny};
        placed[u] = 1;
        total += dfs(idx + 1);
        placed[u] = 0;
      }
      return total;
    };

    return dfs(0);
  }

} // namespace sc_expansion::atomic
