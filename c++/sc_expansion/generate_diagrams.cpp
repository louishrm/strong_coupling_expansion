#include "generate_diagrams.hpp"

namespace sc_expansion {

  std::vector<uint8_t> generate_n_cycle_adjacency_matrix(int n) {
    if (n <= 1) {
      std::vector<uint8_t> adjmat(n * n, 0);
      if (n == 1) adjmat[0] = 1; // Handle 1-cycle
      return adjmat;
    }

    std::vector<int> dest(n, -1);
    std::vector<bool> in_degree_taken(n, false);

    // Greedily build the lexicographically largest sequence of destinations
    for (int u = 0; u < n; ++u) {
      for (int v = n - 1; v >= 0; --v) {
        if (in_degree_taken[v]) continue;

        int curr   = v;
        int length = 1;
        while (dest[curr] != -1) {
          curr = dest[curr];
          length++;
        }

        // If the path loops back to 'u', it's a cycle.
        if (curr == u) {
          if (length == n) {
            dest[u]            = v;
            in_degree_taken[v] = true;
            break;
          }
        } else {
          // No cycle detected, safe to add this edge
          dest[u]            = v;
          in_degree_taken[v] = true;
          break;
        }
      }
    }

    // Construct the final min-lex flattened matrix
    std::vector<uint8_t> canonical_adjmat(n * n, 0);
    for (int u = 0; u < n; ++u) { canonical_adjmat[u * n + dest[u]] = 1; }

    return canonical_adjmat;
  }

  int calculate_n_cycle_free_multiplicity(int n, bool bipartite) {
    if (bipartite) {
      // Guard against odd lengths on a bipartite (square) lattice
      if (n % 2 != 0) return 0;

      int nCn2 = binomial_coefficient(n, n / 2);
      return nCn2 * nCn2;
    } else {
      int result = 0;
      for (int k = 0; k <= n; ++k) {

        int sign       = ((n - k) % 2 == 0) ? 1 : -1;
        int power_term = sign * (1 << (n - k));

        int current_term = binomial_coefficient(n, k) * power_term;

        int inner_sum = 0;
        for (int j = 0; j <= k; ++j) {
          int binom = binomial_coefficient(k, j);
          inner_sum += binom * binom * binom;
        }
        result += current_term * inner_sum;
      }
      return result;
    }
  }

  VacuumDiagramGenerator::VacuumDiagramGenerator(int order_, bool bipartite_only_)
     : order(order_), bipartite_only(bipartite_only_), partitions(order_, order_ / 2) {};

  //manually fill the n-cycle by hand to save time (implies n!n! redundant checks)

  void VacuumDiagramGenerator::generate() {
    //manually fill the n-cycle by hand to save time (implies n!n! redundant checks)
    // Only add n-cycle if it matches the bipartite criteria
    if (!this->bipartite_only || this->order % 2 == 0) {
      std::vector<uint8_t> n_cycle = generate_n_cycle_adjacency_matrix(this->order);

      if (this->bipartite_only && this->order % 2 == 0) {
        // Optimization: Override constructor for n-cycle on square lattice: symmetry factor is n, free multiplicity is [nC(n/2)]^2
        int fm = calculate_n_cycle_free_multiplicity(this->order, true);
        // The specialized constructor bypasses canonicalization and sets canonical_matrix = n_cycle
        if (fm > 0) {
          this->unique_adjmats.insert(n_cycle);
          this->graphs.emplace_back(sc_expansion::Graph(n_cycle, this->order, this->order, this->order, fm, this->bipartite_only));
        }
      } else {
        // Non-bipartite or other cases: use standard constructor to compute FM accurately on the target lattice
        int fm = calculate_n_cycle_free_multiplicity(this->order, false);
        if (fm > 0) {
          this->unique_adjmats.insert(n_cycle);
          this->graphs.push_back(sc_expansion::Graph(n_cycle, this->order, this->order, this->order, fm, this->bipartite_only));
        }
      }
    }

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

    if (std::all_of(partition.begin(), partition.end(), [](int i) { return i == 1; })) { return; }

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
      Graph graph(adjmat, V, this->bipartite_only);
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
    if (!graph.get_connectivity()) return;                      //discard disconnected graphs
    if (this->bipartite_only && !graph.get_bipartite()) return; //discard non-bipartite graphs only if requested
    if (graph.get_free_multiplicity() <= 0) return;             //discard graphs that cannot be embedded on the lattice

    std::vector<uint8_t> canonical = graph.get_canonical_form();
    if (this->unique_adjmats.find(canonical) == this->unique_adjmats.end()) {
      this->unique_adjmats.insert(canonical);
      this->graphs.push_back(graph);
    }
  }
} // namespace sc_expansion
