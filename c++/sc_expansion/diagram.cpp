#include "diagram.hpp"

namespace sc_expansion {

  Diagram::Diagram(Graph const &graph_) : graph(graph_) {
    this->compute_hopping_lines_and_vertex_structures();
    this->compute_global_spin_configurations();
    this->compute_local_spin_configurations();
    this->compute_diagram_sign();
  }

  void Diagram::compute_hopping_lines_and_vertex_structures() {
    std::vector<Line> hopping_lines;
    int V = this->graph.get_V();

    this->vertices.resize(V);

    for (int i = 0; i < V; i++) {
      for (int j = 0; j < V; j++) {
        int line_count = this->graph(i, j);

        for (int k = 0; k < line_count; k++) {
          Line line;
          line.from_vertex = i;
          line.to_vertex   = j;
          hopping_lines.push_back(line);
        }
      }
    }
    this->hopping_lines = hopping_lines;

    for (size_t line_idx = 0; line_idx < this->hopping_lines.size(); line_idx++) {
      const Line &line = this->hopping_lines[line_idx];
      this->vertices[line.from_vertex].outgoing_lines.push_back(line_idx);
      this->vertices[line.to_vertex].incoming_lines.push_back(line_idx);
    }
  }

  void Diagram::compute_global_spin_configurations() {
    std::vector<uint64_t> valid_configs;
    this->config_weights.clear();
    uint64_t num_configs = 1ULL << this->graph.get_order();

    for (uint64_t config = 0; config < num_configs; ++config) {
      std::vector<int> spin_balance(this->graph.get_V(), 0);
      for (int i = 0; i < this->graph.get_order(); ++i) {
        if ((config >> i) & 1) {
          spin_balance[this->hopping_lines[i].to_vertex]++;
          spin_balance[this->hopping_lines[i].from_vertex]--;
        }
      }

      bool valid = true;
      for (int b : spin_balance) {
        if (b != 0) {
          valid = false;
          break;
        }
      }

      if (valid) {
        uint64_t inverted = (~config) & (num_configs - 1);
        if (config <= inverted) {
          valid_configs.push_back(config);
          this->config_weights.push_back(config == inverted ? 1.0 : 2.0);
        }
      }
    }
    this->global_spin_configurations = valid_configs;
  }

  void Diagram::compute_local_spin_configurations() {
    for (uint64_t config : this->global_spin_configurations) {
      for (size_t i = 0; i < this->vertices.size(); i++) {
        int local_mask = 0;
        int bit_pos    = 0;
        Vertex &v      = this->vertices[i];
        for (int line_idx : v.outgoing_lines) {
          if ((config >> line_idx) & 1) local_mask |= (1 << bit_pos);
          bit_pos++;
        }
        for (int line_idx : v.incoming_lines) {
          if ((config >> line_idx) & 1) local_mask |= (1 << bit_pos);
          bit_pos++;
        }
        v.local_spin_configs.push_back(local_mask);
      }
    }
  }

  void Diagram::compute_diagram_sign() {
    int num_loops = 0;
    std::vector<bool> visited_lines(this->graph.get_order(), false);
    auto lines = this->hopping_lines;
    std::vector<int> successor_map(this->graph.get_order());
    for (int v = 0; v < this->graph.get_V(); ++v) {
      std::vector<int> incoming_indices;
      std::vector<int> outgoing_indices;
      for (int i = 0; i < this->graph.get_order(); ++i) {
        if (lines[i].to_vertex == v) incoming_indices.push_back(i);
        if (lines[i].from_vertex == v) outgoing_indices.push_back(i);
      }
      for (size_t i = 0; i < incoming_indices.size(); ++i) { successor_map[incoming_indices[i]] = outgoing_indices[i]; }
    }

    for (int i = 0; i < this->graph.get_order(); ++i) {
      if (!visited_lines[i]) {
        num_loops++;
        int current_line = i;
        while (!visited_lines[current_line]) {
          visited_lines[current_line] = true;
          current_line                = successor_map[current_line];
        }
      }
    }
    this->diagram_sign = (num_loops % 2 == 0) ? 1 : -1;
  }

  template <typename T>
  DiagramEvaluator<T>::DiagramEvaluator(Diagram const &diagram_, Parameters<T> const &params)
     : diagram(diagram_), atom(params.U, params.beta, params.mu) {
    int order = this->diagram.get_graph().get_order();
    this->current_taus.assign(order, -1.0);
    for (const Vertex &v : this->diagram.get_vertices()) {
      this->cache_finite.emplace_back(1 << v.degree(), T(0.0));
      this->cache_infinite.emplace_back(1 << v.degree(), T(0.0));
    }
  }

  template <typename T> void DiagramEvaluator<T>::check_vertex(int v_idx, std::vector<double> const &taus) const {
    const auto &v     = diagram.get_vertices()[v_idx];
    bool is_corrupted = false;

    for (int idx : v.outgoing_lines) {
      if (taus[idx] != this->current_taus[idx]) {
        is_corrupted = true;
        break;
      }
    }
    if (!is_corrupted) {
      for (int idx : v.incoming_lines) {
        if (taus[idx] != this->current_taus[idx]) {
          is_corrupted = true;
          break;
        }
      }
    }

    if (is_corrupted) { this->recompute_vertex(v_idx, taus); }
  }

  template <typename T>
  std::pair<ArgList, ArgList> DiagramEvaluator<T>::get_local_cumul_args(int v_idx, std::vector<double> const &taus, uint32_t local_mask) const {
    const auto &v = this->diagram.get_vertices()[v_idx];
    ArgList unprimed_args, primed_args;

    int bit_pos = 0;
    for (int idx : v.outgoing_lines) {
      int spin = (local_mask >> bit_pos) & 1;
      unprimed_args.push_back({taus[idx], spin});
      bit_pos++;
    }

    for (int idx : v.incoming_lines) {
      int spin = (local_mask >> bit_pos) & 1;
      primed_args.push_back({taus[idx], spin});
      bit_pos++;
    }

    return {unprimed_args, primed_args};
  }

  template <typename T> void DiagramEvaluator<T>::recompute_vertex(int v_idx, std::vector<double> const &taus) const {
    const auto &v = diagram.get_vertices()[v_idx];
    std::vector<bool> already_done(1 << v.degree(), false);

    for (uint32_t mask : v.local_spin_configs) {
      if (already_done[mask]) continue;

      auto args                         = this->get_local_cumul_args(v_idx, taus, mask);
      this->cache_finite[v_idx][mask]   = compute_cumulant_decomposition(args.first, args.second, this->atom, false);
      this->cache_infinite[v_idx][mask] = compute_cumulant_decomposition(args.first, args.second, this->atom, true);

      already_done[mask] = true;
    }
  }

  template <typename T> T DiagramEvaluator<T>::evaluate_at_taus(std::vector<double> const &taus, bool infinite_U, bool use_cache) const {
    int V = this->diagram.get_vertices().size();
    for (int v = 0; v < V; ++v) { this->check_vertex(v, taus); }
    this->current_taus = taus;

    T sum                      = T(0.0);
    const auto &global_configs = this->diagram.get_global_configs();
    const auto &weights        = this->diagram.get_config_weights();
    const auto &vertices       = this->diagram.get_vertices();
    const auto &cache          = infinite_U ? cache_infinite : cache_finite;

    for (size_t g_idx = 0; g_idx < global_configs.size(); ++g_idx) {
      T product = T(1.0);
      for (int v_idx = 0; v_idx < V; ++v_idx) {
        int mask = vertices[v_idx].local_spin_configs[g_idx];

        if (use_cache) {
          product = product * cache[v_idx][mask];
        } else {
          auto args = this->get_local_cumul_args(v_idx, taus, mask);
          product   = product * compute_cumulant_decomposition(args.first, args.second, this->atom, infinite_U);
        }
      }
      sum = sum + T(weights[g_idx]) * product;
    }
    T sign            = T(this->diagram.get_diagram_sign());
    T symmetry_factor = T(this->diagram.get_graph().get_symmetry_factor());
    T fm              = T(this->diagram.get_graph().get_free_multiplicity());
    T prefactor       = (T(-1.0) / this->atom.beta) * sign / symmetry_factor; // * fm;

    return prefactor * sum;
  }

  template <typename T> T DiagramEvaluator<T>::evaluate_at_taus_dimer(std::vector<double> const &taus, bool infinite_U, bool use_cache) const {
    int V = this->diagram.get_vertices().size();
    for (int v = 0; v < V; ++v) { this->check_vertex(v, taus); }
    this->current_taus = taus;

    T sum                      = T(0.0);
    const auto &global_configs = this->diagram.get_global_configs();
    const auto &weights        = this->diagram.get_config_weights();
    const auto &vertices       = this->diagram.get_vertices();
    const auto &cache          = infinite_U ? cache_infinite : cache_finite;

    for (size_t g_idx = 0; g_idx < global_configs.size(); ++g_idx) {
      T product = T(1.0);
      for (int v_idx = 0; v_idx < V; ++v_idx) {
        int mask = vertices[v_idx].local_spin_configs[g_idx];

        if (use_cache) {
          product = product * cache[v_idx][mask];
        } else {
          auto args = this->get_local_cumul_args(v_idx, taus, mask);
          product   = product * compute_cumulant_decomposition(args.first, args.second, this->atom, infinite_U);
        }
      }
      sum = sum + T(weights[g_idx]) * product;
    }
    T sign            = T(this->diagram.get_diagram_sign());
    T symmetry_factor = T(this->diagram.get_graph().get_symmetry_factor());
    T fm              = T(1.0); // Fixed for dimer calculation
    T prefactor       = (T(-1.0) / this->atom.beta) * sign / symmetry_factor * fm;

    return prefactor * sum;
  }

  template class DiagramEvaluator<double>;
  template class DiagramEvaluator<Dual>;
} // namespace sc_expansion
