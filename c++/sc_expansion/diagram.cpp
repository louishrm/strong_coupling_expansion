#include "diagram.hpp"

namespace sc_expansion {

  Diagram::Diagram(adjmat adjacency_matrix, double U, double beta, double mu) : atom(U, beta, mu) {

    this->adjacency_matrix = adjacency_matrix;

    this->V = adjacency_matrix.size();

    int order = 0;
    for (const auto &row : adjacency_matrix) {
      for (const auto &entry : row) { order += entry; }
    }
    this->n = order;
  }

  bool Diagram::is_particle_number_conserving() const {

    for (int i = 0; i < this->V; i++) {
      int total_in  = 0;
      int total_out = 0;
      for (int j = 0; j < this->V; j++) {
        total_in += this->adjacency_matrix[j][i];
        total_out += this->adjacency_matrix[i][j];
      }
      if (total_in != total_out) { return false; }
    }
    return true;
  }

  bool Diagram::is_connected() const {

    //Use BFS to check connectivity
    std::vector<bool> visited(this->V, false);
    std::queue<int> q;
    q.push(0);

    int visited_count = 0;
    while (!q.empty()) {
      int vertex = q.front(); //get the next vertex to visit
      q.pop();                //remove it from the queue
      //check neighbors
      for (int neighbor = 0; neighbor < this->V; neighbor++) {
        //check for connection in either direction
        bool connected = (this->adjacency_matrix[vertex][neighbor] > 0) || (this->adjacency_matrix[neighbor][vertex] > 0);
        if (connected && !visited[neighbor]) { //if there is a connection and neighbor not visited
          visited[neighbor] = true;            //mark as visited
          q.push(neighbor);                    //add to queue for further exploration
          visited_count++;                     //increment visited count
        }
      }
    }
    return visited_count == this->V;
  }

  std::vector<std::vector<int>> generate_permutations(std::vector<int> vertices) {

    std::vector<std::vector<int>> permutations;
    std::sort(vertices.begin(), vertices.end());
    do { permutations.push_back(vertices); } while (std::next_permutation(vertices.begin(), vertices.end()));
    return permutations;
  }

  int factorial(int k) {
    if (k <= 1) return 1;
    return k * factorial(k - 1);
  }

  int Diagram::diagram_sign() const {

    int num_loops = 0;
    std::vector<bool> visited_lines(this->n, false);
    auto lines = this->get_hopping_lines();

    // Step 1: Build the pairing map (successor map)
    // For each line 'i', this map will tell you the index of the next line in its loop.
    std::vector<int> successor_map(this->n);
    for (int v = 0; v < this->V; ++v) {
      std::vector<int> incoming_indices;
      std::vector<int> outgoing_indices;
      for (int i = 0; i < this->n; ++i) {
        if (lines[i].to_vertex == v) incoming_indices.push_back(i);
        if (lines[i].from_vertex == v) outgoing_indices.push_back(i);
      }
      // Pair the i-th incoming line with the i-th outgoing line
      for (size_t i = 0; i < incoming_indices.size(); ++i) { successor_map[incoming_indices[i]] = outgoing_indices[i]; }
    }

    // Step 2: Traverse and count
    for (int i = 0; i < this->n; ++i) {
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
    return (num_loops % 2 == 0) ? 1 : -1;
  }

  int Diagram::get_symmetry_factor() const {

    //get all V! permutations of vertices.
    std::vector<int> vertices(this->V);
    std::iota(vertices.begin(), vertices.end(), 0); // Fill with 0, 1, ..., V-1
    auto permutations = generate_permutations(vertices);

    //For each permutation, check if the permuted adjacency matrix matches the original.
    int symmetry_count = 0;
    for (const auto &perm : permutations) {

      //Create a new adjacency matrix based on the permutation
      adjmat permuted_matrix(this->V, std::vector<int>(this->V, 0));
      for (int i = 0; i < this->V; i++) {
        for (int j = 0; j < this->V; j++) { permuted_matrix[i][j] = this->adjacency_matrix[perm[i]][perm[j]]; }
        //Check if the permuted matrix matches the original
        if (permuted_matrix == this->adjacency_matrix) { symmetry_count++; }
      }
    }
    //product of the factorials of the entries
    int factorial_product = 1;
    for (int i = 0; i < this->V; i++) {
      for (int j = 0; j < this->V; j++) { factorial_product *= factorial(this->adjacency_matrix[i][j]); }
    }
    return symmetry_count * factorial_product;
  }

  std::vector<Diagram::Line> Diagram::get_hopping_lines() const {

    std::vector<Line> hopping_lines;
    for (int i = 0; i < this->V; i++) {
      for (int j = 0; j < this->V; j++) {
        int line_count = this->adjacency_matrix[i][j];
        for (int k = 0; k < line_count; k++) {
          Line line;
          line.from_vertex = i;
          line.to_vertex   = j;
          hopping_lines.push_back(line);
        }
      }
    }
    return hopping_lines;
  }

  double Diagram::evaluate_at_points(HubbardAtom::cumul_args const &args) const {

    //evaluates a diagram at a given set of time-spin args
    auto hopping_lines = this->get_hopping_lines();
    std::vector<HubbardAtom::cumul_args> unprimed_args_per_vertex(this->V);
    std::vector<HubbardAtom::cumul_args> primed_args_per_vertex(this->V);

    // Monotonic shift to enforce strict ordering and uniqueness
    // delta, 2*delta, 3*delta, ...
    //constexpr double DELTA = 1.0e-12;

    for (size_t line_idx = 0; line_idx < hopping_lines.size(); line_idx++) { //loop through each hopping line

      auto line = hopping_lines[line_idx]; //get the destroy and create vertices for this line

      // Unprimed (Annihilation): Shift by 2*i * delta
      auto arg_unprimed = args[line_idx];
      // arg_unprimed.first += (2 * line_idx) * DELTA;
      unprimed_args_per_vertex[line.from_vertex].push_back(arg_unprimed);

      // Primed (Creation): Shift by (2*i + 1) * delta
      // Ensures Creation > Annihilation (Normal Ordering)
      auto arg_primed = args[line_idx];
      // arg_primed.first += (2 * line_idx + 1) * DELTA;
      primed_args_per_vertex[line.to_vertex].push_back(arg_primed);
    }
    double prod = 1.0;
    for (int vertex = 0; vertex < this->V; vertex++) {

      HubbardAtom::cumul_args unprimed_args = unprimed_args_per_vertex[vertex];

      HubbardAtom::cumul_args primed_args = primed_args_per_vertex[vertex];
      prod *= compute_cumulant_decomposition(unprimed_args, primed_args, this->atom, false);
    }
    return prod;
  }

  double Diagram::evaluate_at_taus(std::vector<double> const &taus) const {

    //evaluates a diagram at a given set of times, summing over all spin indices
    double spin_sum  = 0.0;
    long num_configs = 1 << this->n;                        //2^n spin configurations
    for (long config = 0; config < num_configs; config++) { //TODO: only consider spin-conserving configs

      HubbardAtom::cumul_args args;
      for (int i = 0; i < this->n; i++) {
        int spin = (config & (1 << i)) ? 1 : 0; //extract spin from bit representation
        args.push_back({taus[i], spin});
      }
      spin_sum += this->evaluate_at_points(args);
    }

    double symmetry_factor   = (double)this->get_symmetry_factor();
    double diagram_sign      = (double)this->diagram_sign();
    double free_multiplicity = 1.0; //assume 1 for now, can be modified later if needed
    return diagram_sign * spin_sum * free_multiplicity / symmetry_factor;
  }

} // namespace sc_expansion