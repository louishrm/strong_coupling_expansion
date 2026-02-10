#include "generate_diagrams.hpp"

namespace sc_expansion {

  VacuumDiagramGenerator::VacuumDiagramGenerator(int order_) : order(order_), partitions(order_, order_ / 2) {};

  void VacuumDiagramGenerator::generate() {
    this->partitions.reset();
    while (this->partitions.is_valid()) {
      this->propose_create_process(this->partitions.current());
      if (!this->partitions.next()) break;
    }
  }

  void VacuumDiagramGenerator::propose_create_process(const std::vector<int> &partition) {

    /*Given a partition of number of lines, creates a source vector whose entries are the vertex index of the outgoing lines
    Then, copies this vector and iterates over permutations of the target vector, creating a new graph for each valid permutation
    For each graph, check if it is connected, then canonicalize and check it hasn't been visited
    If it hasn't, add it to the hash table*/

    std::vector<int> source(this->order);
    int current_vertex = 0;
    int V              = partition.size(); //number of vertices for this candidate

    //fill the source vector
    int k = 0;
    for (auto const &entry : partition) {
      for (int i = 0; i < entry; ++i) { source[k++] = current_vertex; }
      current_vertex++;
    }

    std::vector<int> target = source;

    do {
      bool valid_connections = true;
      for (size_t i = 0; i < this->order; ++i) {
        if (source[i] == target[i]) {
          valid_connections = false; //no self-connections allowed
          break;
        }
      }

      if (!valid_connections) continue;
      std::vector<uint8_t> adjmat = this->fill_matrix(source, target, V);
      Graph graph(adjmat, V);
      this->process_graph(graph);
    } while (std::next_permutation(target.begin(), target.end()));
  }

  std::vector<uint8_t> VacuumDiagramGenerator::fill_matrix(const std::vector<int> &source, const std::vector<int> &target, int V) const {

    std::vector<uint8_t> adjmat(V * V, 0);
    for (size_t i = 0; i < source.size(); ++i) { adjmat[source[i] * V + target[i]]++; }
    return adjmat;
  }

  void VacuumDiagramGenerator::process_graph(const Graph &graph) {

    /*Check if a graph is connected
    if it is, canonicalize and check it hasn't been visited
    if it hasn't, add it to the hash table*/
    if (!graph.get_connectivity()) return; //discard disconnected graphs

    std::vector<uint8_t> canonical = graph.get_canonical_form();
    if (this->unique_graphs.find(canonical) == this->unique_graphs.end()) { this->unique_graphs.insert(canonical); }
  }
} // namespace sc_expansion