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

    this->hopping_lines   = this->compute_hopping_lines();
    this->sign            = (double)this->compute_diagram_sign();
    this->symmetry_factor = (double)this->compute_symmetry_factor();
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

  int Diagram::compute_diagram_sign() const {

    int num_loops = 0;
    std::vector<bool> visited_lines(this->n, false);
    auto lines = this->hopping_lines;
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

  // Public accessor returning cached value
  int Diagram::diagram_sign() const { return (int)this->sign; }

  int Diagram::compute_symmetry_factor() const {

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

  // Public accessor returning cached value
  int Diagram::get_symmetry_factor() const { return (int)this->symmetry_factor; }

  std::vector<Diagram::Line> Diagram::compute_hopping_lines() const {

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

  // Public accessor returning cached value
  std::vector<Diagram::Line> Diagram::get_hopping_lines() const { return this->hopping_lines; }

  double Diagram::evaluate_at_points(HubbardAtom::cumul_args const &args) const {

    //evaluates a diagram at a given set of time-spin args
    // Use cached hopping_lines
    // auto hopping_lines = this->get_hopping_lines(); // This is now cheap, but we can access member directly.
    // Use reference to avoid copy if possible, though Line is small.
    const auto &lines = this->hopping_lines;

    std::vector<HubbardAtom::cumul_args> unprimed_args_per_vertex(this->V);
    std::vector<HubbardAtom::cumul_args> primed_args_per_vertex(this->V);

    // Monotonic shift to enforce strict ordering and uniqueness
    // delta, 2*delta, 3*delta, ...
    //constexpr double DELTA = 1.0e-12;

    for (size_t line_idx = 0; line_idx < lines.size(); line_idx++) { //loop through each hopping line

      auto line = lines[line_idx]; //get the destroy and create vertices for this line

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

    //evaluates a diagram at a given set of time-spin args
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

    double free_multiplicity = 1.0; //assume 1 for now, can be modified later if needed
    // Use cached values
    return this->sign * spin_sum * free_multiplicity / this->symmetry_factor;
  }

  void next_step(adjmat &A, std::vector<int> &sequence, int vertex, int V, int order) {

    if (sequence.size() == order + 1) { return; }
    for (int j = 0; j < V; j++) {
      if (A[vertex][j] > 0) {
        sequence.push_back(j);
        A[vertex][j]--;
        next_step(A, sequence, j, V, order);
      }
    }
  }

  std::vector<int> generate_walk_sequence(adjmat A, int order) {

    int V = A.size();
    std::vector<int> sequence;
    sequence.push_back(0); // Start from vertex 0
    next_step(A, sequence, 0, V, order);
    return sequence;
  }

  Point::Point() : x(0), y(0) {}

  Point::Point(int x_, int y_) : x(x_), y(y_) {}

  SquareLattice::SquareLattice() {}

  std::vector<Point> SquareLattice::get_neighbors(Point const &r) const {
    return {Point(r.x + 1, r.y), Point(r.x - 1, r.y), Point(r.x, r.y + 1), Point(r.x, r.y - 1)};
  }
  int SquareLattice::manhattan_distance(Point const &r) const { return std::abs(r.x) + std::abs(r.y); }

  bool SquareLattice::prune(Point const &r, int current_distance, int order) const {

    int remaining_steps = order - current_distance;
    int dist            = this->manhattan_distance(r);
    return dist > remaining_steps;
  }

  bool SquareLattice::is_neighbor(Point const &r1, Point const &r2) const { return (std::abs(r1.x - r2.x) + std::abs(r1.y - r2.y)) == 1; }

  void place_next_vertex(std::unordered_map<int, Point> &placed_vertices, SquareLattice const &lattice, std::vector<int> const &sequence,
                         int current_vertex_index, int order, int &free_multiplicity, int hopping_count) {

    /* Place the next vertex in the non self-avoiding lattice walk. 

    Args: place_vertices: a map whose keys are vertex index and values are their position. 
          lattice: the lattice object to get neighbors and prune positions.
          sequence: the sequence of vertex indices in the walk.
          current_vertex_index: the index of the last placed vertex in the sequence.
          order: total number of hopping lines (length of the walk).
          free_multiplicity: reference to the count of valid placements found so far.
          hopping_count: reference to the count of hopping lines placed so far.

    Start on vertex 0.

    1. Base case: check if the hopping count is equal to the order. If yes, increment free multiplicity and stop. 

    2. Check if the next vertex in the sequence is already placed. If yes, make sure it is a neighbor of the last placed vertex. If yes, recurse to place the next vertex stop. 

    3. If not placed, get the neighbors of the last placed vertex. For each neighbor, check if placing the next vertex there is valid (not pruned). If valid, place the vertex and recurse. After recursion, remove the vertex to backtrack.
    */

    if (hopping_count == order) {
      free_multiplicity++;
      return;
    }

    Point last_position          = placed_vertices[current_vertex_index];
    std::vector<Point> neighbors = lattice.get_neighbors(last_position);

    int next_vertex_index = sequence[hopping_count + 1];
    if (placed_vertices.find(next_vertex_index) != placed_vertices.end()) {

      Point neighbor = placed_vertices[next_vertex_index];
      if (lattice.is_neighbor(last_position, neighbor)) {
        place_next_vertex(placed_vertices, lattice, sequence, next_vertex_index, order, free_multiplicity, hopping_count + 1);
      }
      return;
    }

    for (auto const &neighbor : neighbors) {
      if (lattice.prune(neighbor, hopping_count + 1, order)) { continue; }
      placed_vertices[next_vertex_index] = neighbor;
      place_next_vertex(placed_vertices, lattice, sequence, next_vertex_index, order, free_multiplicity, hopping_count + 1);
      placed_vertices.erase(next_vertex_index);
    }
  }

  int compute_free_multiplicity(adjmat A, int order) {

    SquareLattice lattice;
    std::vector<int> sequence = generate_walk_sequence(A, order);

    std::unordered_map<int, Point> placed_vertices;
    placed_vertices[0] = Point(0, 0); //place the first vertex at origin

    int free_multiplicity = 0;
    place_next_vertex(placed_vertices, lattice, sequence, 0, order, free_multiplicity, 0);
    return free_multiplicity;
  }
} // namespace sc_expansion