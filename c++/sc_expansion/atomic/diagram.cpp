#include "diagram.hpp"
#include "../dual.hpp"
#include "../fock_space.hpp"
#include <algorithm>
#include <map>
#include <set>

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

  template <typename T> void Diagram<T>::setup_vertices(std::vector<VertexType<T> *> const &vertex_types) {
    int V = this->graph.get_V();
    this->vertex_type_ptrs.resize(V, nullptr);

    for (int v = 0; v < V; ++v) {
      int degree = 0;
      for (int j = 0; j < V; ++j) { degree += this->graph(v, j) + this->graph(j, v); }
      int cumulant_order = degree / 2;
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

    using Config = GlobalConfiguration<1>;
    std::map<std::vector<uint8_t>, double> canonical_weights;
    int n_spin_configs = 1 << n_lines;
    std::set<std::vector<uint8_t>> seen;

    for (int spin_mask = 0; spin_mask < n_spin_configs; ++spin_mask) {
      Config global;
      global.config.reserve(2 * n_lines);
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

        if (up_cre != up_ann || dn_cre != dn_ann) valid = false;
      }

      if (!valid) continue;

      auto orbit      = SymmetryGroup<1, T>::get_orbit(global);
      auto &canonical = orbit[0];
      int orbit_size  = (int)orbit.size();

      if (!seen.insert(canonical.config).second) continue;

      double weight = free_mult * orbit_size / sym_factor;
      canonical_weights[canonical.config] += weight;
    }

    for (auto &[config, weight] : canonical_weights) { this->valid_configurations.push_back({config, weight}); }
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
      for (auto const &leg : this->legs_per_vertex[v]) { this->tau_to_vertices[leg.line_index].push_back(v); }
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

    T prefactor = (T(-1.0) / solver.params.beta) * T(this->diagram_sign);
    return prefactor * sum;
  }

  template class Diagram<double>;
  template class Diagram<Dual>;

} // namespace sc_expansion::atomic
