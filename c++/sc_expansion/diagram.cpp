#include "diagram.hpp"

Diagram::Diagram(const adjmat &adjacency_matrix) {

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

double Diagram::evaluate_at_points(triqs::atom_diag::atom_diag<false> ad, double beta, hubbard_atom::cumul_args args) const {

  //evaluates a diagram at a given set of time-spin args

  auto hopping_lines = this->get_hopping_lines();
  std::vector<hubbard_atom::cumul_args> unprimed_args_per_vertex(this->V);
  std::vector<hubbard_atom::cumul_args> primed_args_per_vertex(this->V);

  for (size_t line_idx = 0; line_idx < hopping_lines.size(); line_idx++) { //loop through each hopping line

    auto line = hopping_lines[line_idx];                                  //get the destroy and create vertices for this line
    unprimed_args_per_vertex[line.from_vertex].push_back(args[line_idx]); //assign the unprimed args to the from vertex
    primed_args_per_vertex[line.to_vertex].push_back(args[line_idx]);     //assign the primed args to the to vertex
  }

  double prod = 1.0;
  for (int vertex = 0; vertex < this->V; vertex++) {

    hubbard_atom::cumul_args unprimed_args = unprimed_args_per_vertex[vertex];
    hubbard_atom::cumul_args primed_args   = primed_args_per_vertex[vertex];
    prod *= compute_cumulant_decomposition(unprimed_args, primed_args, ad, beta);
  }
  return prod;
}
