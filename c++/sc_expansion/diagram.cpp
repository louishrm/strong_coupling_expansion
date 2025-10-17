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