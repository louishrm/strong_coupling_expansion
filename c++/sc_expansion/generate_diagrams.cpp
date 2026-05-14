#include "generate_diagrams.hpp"
#include "canonicalize.hpp"
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace sc_expansion {

  // Binary layout:
  //   int32 n_graphs
  //   per graph:
  //     int32 V
  //     int32 automorphism_count
  //     int32 symmetry_factor
  //     int32 free_multiplicity
  //     uint8 bipartite_only
  //     uint8[V*V] adjacency (canonical form)
  // SC_EXPANSION_PROJECT_ROOT is injected at configure time (see CMakeLists.txt)
  // so the cache lives at <project>/diagrams/ regardless of the executable's CWD.
  static const std::string DIAGRAMS_DIR = std::string(SC_EXPANSION_PROJECT_ROOT) + "/diagrams";

  std::string bipartite_diagrams_path(int order) {
    return DIAGRAMS_DIR + "/bipartite_order_" + std::to_string(order);
  }

  void save_bipartite_graphs(int order, const std::vector<Graph> &graphs) {
    std::filesystem::create_directories(DIAGRAMS_DIR);
    std::ofstream out(bipartite_diagrams_path(order), std::ios::binary);
    if (!out) { throw std::runtime_error("Failed to open " + bipartite_diagrams_path(order) + " for writing"); }

    int32_t n = static_cast<int32_t>(graphs.size());
    out.write(reinterpret_cast<const char *>(&n), sizeof(n));

    for (auto const &g : graphs) {
      int32_t V     = g.get_V();
      int32_t aut   = g.get_automorphism_count();
      int32_t sym   = static_cast<int32_t>(g.get_symmetry_factor());
      int32_t fm    = static_cast<int32_t>(g.get_free_multiplicity());
      uint8_t bonly = g.get_bipartite_only() ? 1 : 0;
      auto adj      = g.get_canonical_form();

      out.write(reinterpret_cast<const char *>(&V), sizeof(V));
      out.write(reinterpret_cast<const char *>(&aut), sizeof(aut));
      out.write(reinterpret_cast<const char *>(&sym), sizeof(sym));
      out.write(reinterpret_cast<const char *>(&fm), sizeof(fm));
      out.write(reinterpret_cast<const char *>(&bonly), sizeof(bonly));
      out.write(reinterpret_cast<const char *>(adj.data()), static_cast<std::streamsize>(adj.size()));
    }
  }

  bool load_bipartite_graphs(int order, std::vector<Graph> &out_graphs) {
    std::string path = bipartite_diagrams_path(order);
    if (!std::filesystem::exists(path)) return false;

    std::ifstream in(path, std::ios::binary);
    if (!in) return false;

    int32_t n = 0;
    in.read(reinterpret_cast<char *>(&n), sizeof(n));
    if (!in) throw std::runtime_error("Corrupt diagram cache: " + path);

    out_graphs.clear();
    out_graphs.reserve(n);

    for (int32_t i = 0; i < n; i++) {
      int32_t V, aut, sym, fm;
      uint8_t bonly;
      in.read(reinterpret_cast<char *>(&V), sizeof(V));
      in.read(reinterpret_cast<char *>(&aut), sizeof(aut));
      in.read(reinterpret_cast<char *>(&sym), sizeof(sym));
      in.read(reinterpret_cast<char *>(&fm), sizeof(fm));
      in.read(reinterpret_cast<char *>(&bonly), sizeof(bonly));

      std::vector<uint8_t> adj(static_cast<size_t>(V) * V);
      in.read(reinterpret_cast<char *>(adj.data()), static_cast<std::streamsize>(adj.size()));
      if (!in) throw std::runtime_error("Corrupt diagram cache: " + path);

      out_graphs.emplace_back(adj, V, aut, sym, fm, bonly != 0);
    }
    return true;
  }

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

    // Order 1 has no valid vacuum diagram without self-loops.
    if (this->order == 1) return;

    //manually fill the n-cycle by hand to save time (implies n!n! redundant checks)
    // Only add n-cycle if it matches the bipartite criteria
    if (!this->bipartite_only || this->order % 2 == 0) {
      std::vector<uint8_t> n_cycle = generate_n_cycle_adjacency_matrix(this->order);
      // Canonicalize so the stored canonical_matrix matches what the standard
      // pipeline would produce; the override constructor below stores its
      // input as both adjacency_matrix and canonical_matrix verbatim.
      auto canon = canonicalize(n_cycle, this->order, std::vector<int>(this->order, 0));

      int fm = calculate_n_cycle_free_multiplicity(this->order, this->bipartite_only && this->order % 2 == 0);
      if (fm > 0) {
        this->unique_adjmats.insert(canon.canonical_matrix);
        this->graphs.emplace_back(sc_expansion::Graph(canon.canonical_matrix, this->order, this->order, this->order, fm, this->bipartite_only));
      }
    }

    this->partitions.reset();

    while (this->partitions.is_valid()) {
      this->propose_create_process(this->partitions.current());
      if (!this->partitions.next()) break;
    }

    // Compute lattice embeddings now that the unique set is final. The
    // n-cycle fast-path graphs already carry a precomputed free_multiplicity
    // (set via the override constructor) and so are skipped.
    for (auto &g : this->graphs) {
      if (g.get_free_multiplicity() == 0) g.compute_free_multiplicity();
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
      for (size_t i = 0; i < (size_t)this->order; ++i) {
        if (source[i] == target[i]) {
          valid_connections = false; //no self-connections allowed
          break;
        }
      }
      if (!valid_connections) continue;

      std::vector<uint8_t> adjmat = this->fill_matrix(source, target, V);
      Graph graph(adjmat, V, this->bipartite_only, /*defer_embedding=*/true);
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

    std::vector<uint8_t> canonical = graph.get_canonical_form();
    if (this->unique_adjmats.find(canonical) == this->unique_adjmats.end()) {
      this->unique_adjmats.insert(canonical);
      this->graphs.push_back(graph);
    }
  }

  // ---------------------------------------------------------------------------
  // RootedDiagramGenerator
  // ---------------------------------------------------------------------------

  RootedDiagramGenerator::RootedDiagramGenerator(int order_, int n_marks_, bool bipartite_only_)
     : order(order_), n_marks(n_marks_), bipartite_only(bipartite_only_) {
    if (n_marks_ != 1 && n_marks_ != 2) throw std::invalid_argument("RootedDiagramGenerator: n_marks must be 1 or 2");
  }

  void RootedDiagramGenerator::generate() {
    // Obtain the vacuum catalog: prefer the on-disk cache; otherwise generate
    // fresh. Only the canonical adjacency and V matter for rooting — the
    // vacuum symmetry / multiplicity fields are not consumed here.
    std::vector<Graph> vacuum;
    if (!(this->bipartite_only && load_bipartite_graphs(this->order, vacuum))) {
      VacuumDiagramGenerator gen(this->order, this->bipartite_only);
      gen.generate();
      vacuum = gen.get_unique_graphs();
    }

    for (auto const &G : vacuum) {
      int V = G.get_V();
      if (this->n_marks == 1) {
        for (int v = 0; v < V; ++v) this->try_emit(G, {v});
      } else {
        for (int i = 0; i < V; ++i)
          for (int j = i + 1; j < V; ++j) this->try_emit(G, {i, j});
      }
    }

    // Embed only the unique rooted forms.
    for (auto &rg : this->rooted_graphs) rg.compute_shell_multiplicity();
  }

  void RootedDiagramGenerator::try_emit(Graph const &G, std::vector<int> marks) {
    RootedGraph rg(G, marks, /*defer_embedding=*/true);
    RootedKey key{rg.get_canonical_form(), rg.get_marks()};
    if (this->unique_keys.insert(key).second) this->rooted_graphs.push_back(std::move(rg));
  }

} // namespace sc_expansion
