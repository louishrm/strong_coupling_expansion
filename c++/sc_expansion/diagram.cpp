#include "diagram.hpp"

namespace sc_expansion {

  Diagram::Diagram(Graph const &graph_) : graph(graph_) {

    this->compute_hopping_lines();
    this->compute_diagram_sign();
    this->compute_valid_spin_configurations();
  }

  void Diagram::compute_hopping_lines() {

    std::vector<Line> hopping_lines;
    int V = this->graph.get_V();
    this->unprimed_line_indices_per_vertex.assign(V, {});
    this->primed_line_indices_per_vertex.assign(V, {});

    int line_idx = 0;
    for (int i = 0; i < V; ++i) {
      for (int j = 0; j < V; ++j) {
        int line_count = this->graph(i, j);
        for (int k = 0; k < line_count; ++k) {
          Line line;
          line.from_vertex = i;
          line.to_vertex   = j;
          hopping_lines.push_back(line);
          this->unprimed_line_indices_per_vertex[i].push_back(line_idx);
          this->primed_line_indices_per_vertex[j].push_back(line_idx);
          line_idx++;
        }
      }
    }
    this->hopping_lines = hopping_lines;
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
      // Pair the i-th incoming line with the i-th outgoing line
      for (size_t i = 0; i < incoming_indices.size(); ++i) { successor_map[incoming_indices[i]] = outgoing_indices[i]; }
    }

    // Step 2: Traverse and count
    for (int i = 0; i < this->graph.get_order(); ++i) {
      if (!visited_lines[i]) {
        num_loops++;
        int current_line = i;
        // Follow the loop until we've marked all its lines
        while (!visited_lines[current_line]) {
          visited_lines[current_line] = true;
          current_line                = successor_map[current_line];
        }
      }
    }

    // Step 3: Determine the sign
    this->diagram_sign = (num_loops % 2 == 0) ? 1 : -1;
  }

  void Diagram::compute_valid_spin_configurations() {
    std::vector<long> valid_configs;
    long num_configs = 1L << this->graph.get_order(); // 2^(number of lines)

    for (long config = 0; config < num_configs; ++config) {

      // Check spin conservation
      // We only track the flow of UP spins. Since total particle number is conserved (checked elsewhere),
      // if UP spin flow is conserved, DOWN spin flow is also conserved.
      std::vector<int> spin_balance(this->graph.get_V(), 0);
      for (int i = 0; i < this->graph.get_order(); ++i) {
        if ((config >> i) & 1) { // if spin is UP
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
        // Spin Inversion Symmetry:
        // The Hubbard model (with B=0) is symmetric under spin inversion (UP <-> DOWN).
        // If 'config' is valid, then 'inverted' (bitwise NOT of config) is also valid
        // and yields the same diagram value.
        // We store only the representative (lexicographically smaller) to save computation.
        long inverted = (~config) & (num_configs - 1);
        if (config <= inverted) { valid_configs.push_back(config); }
      }
    }
    this->valid_spin_configurations = valid_configs;
  }

  // --- DiagramEvaluator Implementation ---

  DiagramEvaluator::DiagramEvaluator(Diagram const &diagram_, Parameters const &params) : diagram(diagram_), atom(params.U, params.beta, params.mu) {}

  double DiagramEvaluator::evaluate_at_points(HubbardAtom::cumul_args const &args, bool infinite_U) const {

    int V = this->diagram.get_graph().get_V();

    double prod = 1.0;
    for (int vertex = 0; vertex < V; vertex++) {
      const auto &u_indices = this->diagram.get_unprimed_indices(vertex);
      const auto &p_indices = this->diagram.get_primed_indices(vertex);

      HubbardAtom::cumul_args unprimed_args;
      HubbardAtom::cumul_args primed_args;
      unprimed_args.reserve(u_indices.size());
      primed_args.reserve(p_indices.size());

      for (int idx : u_indices) unprimed_args.push_back(args[idx]);
      for (int idx : p_indices) primed_args.push_back(args[idx]);

      double new_factor = compute_cumulant_decomposition(unprimed_args, primed_args, this->atom, infinite_U);

      if (new_factor == 0.0) { return 0.0; }
      prod *= new_factor;
    }

    return prod;
  }

  double DiagramEvaluator::evaluate_at_taus(std::vector<double> const &taus, bool infinite_U) const {

    double spin_sum     = 0.0;
    int order           = this->diagram.get_graph().get_order();
    const auto &configs = this->diagram.get_valid_spin_configurations();

    HubbardAtom::cumul_args args(order);
    for (long config : configs) {
      for (int i = 0; i < order; i++) {
        int spin = (config & (1L << i)) ? 1 : 0;
        args[i]  = {taus[i], spin};
      }

      double val = this->evaluate_at_points(args, infinite_U);

      long inverted = (~config) & ((1L << order) - 1);
      if (config != inverted) { val *= 2.0; }
      spin_sum += val;
    }

    double sign            = this->diagram.get_diagram_sign();
    double symmetry_factor = this->diagram.get_graph().get_symmetry_factor();
    double fm              = this->diagram.get_graph().get_free_multiplicity();

    return (-1.0 / this->atom.beta) * sign * spin_sum / symmetry_factor * fm;
  }

} // namespace sc_expansion