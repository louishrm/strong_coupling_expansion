#include "generate_diagrams.hpp"
#include "canonicalize.hpp"
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <queue>
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

  std::string bipartite_diagrams_path(int order) { return DIAGRAMS_DIR + "/bipartite_order_" + std::to_string(order); }

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

  void VacuumDiagramGenerator::emit_n_cycle() {
    //manually fill the n-cycle by hand to save time (implies n!n! redundant checks)
    // Only add n-cycle if it matches the bipartite criteria
    if (this->bipartite_only && this->order % 2 != 0) return;

    std::vector<uint8_t> n_cycle = generate_n_cycle_adjacency_matrix(this->order);
    // Canonicalize so the stored canonical_matrix matches what the standard
    // pipeline would produce; the override constructor below stores its
    // input as both adjacency_matrix and canonical_matrix verbatim.
    auto canon = canonicalize(n_cycle, this->order, std::vector<int>(this->order, 0));

    int fm = calculate_n_cycle_free_multiplicity(this->order, this->bipartite_only && this->order % 2 == 0);
    if (fm > 0) {
      this->unique_adjmats.insert(canon.canonical_matrix);
      this->graphs.emplace_back(
         sc_expansion::Graph(canon.canonical_matrix, this->order, this->order, this->order, fm, this->bipartite_only));
    }
  }

  void VacuumDiagramGenerator::generate() {

    // Order 1 has no valid vacuum diagram without self-loops.
    if (this->order == 1) return;

    this->emit_n_cycle();

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

  void VacuumDiagramGenerator::generate_at_vertex_count(int V) {
    // Clear any previous results — this entry point is designed to be called
    // multiple times (once per V) by DistanceRootedDiagramGenerator.
    this->graphs.clear();
    this->unique_adjmats.clear();

    if (V <= 0 || V > this->order) return;
    if (this->order == 1) return; // no valid vacuum diagram at order 1

    // V == order fast path: the only connected vacuum topology with every
    // vertex at degree exactly 2 is the n-cycle. Skip the proposal sweep and
    // the (most-expensive) canonicalization of the densely-symmetric n-cycle.
    if (V == this->order) {
      this->emit_n_cycle();
      return;
    }

    // Generic case: sweep partitions of `order` whose length is exactly V.
    // PartitionGenerator yields all partitions of `order` with parts <= order/2;
    // we filter by length here.
    //
    // Vacuum free_multiplicity is intentionally NOT computed here: this entry
    // point exists for the distance-rooted pipeline, which consumes only the
    // adjacency / V of each vacuum graph. Callers that need vacuum embeddings
    // should use generate() instead.
    this->partitions.reset();
    while (this->partitions.is_valid()) {
      auto const &p = this->partitions.current();
      if (static_cast<int>(p.size()) == V) this->propose_create_process(p);
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
        // Include same-vertex pairs (i == j): both density operators on one
        // hopping vertex. These contribute to the connected correlator with
        // the unordered-pair convention; bliss colors handle dedup (a vertex
        // with two marks gets color 2, distinct from single-mark color 1).
        for (int i = 0; i < V; ++i)
          for (int j = i; j < V; ++j) this->try_emit(G, {i, j});
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

  // ---------------------------------------------------------------------------
  // DistanceRootedDiagramGenerator
  // ---------------------------------------------------------------------------

  DistanceRootedDiagramGenerator::DistanceRootedDiagramGenerator(std::vector<int> r_, int n_, bool bipartite_only_)
     : r(std::move(r_)), n(n_), bipartite_only(bipartite_only_) {
    int dd = 0;
    for (int c : this->r) dd += std::abs(c);
    this->d = dd;
    if (this->n < 2 * this->d) {
      throw std::invalid_argument("DistanceRootedDiagramGenerator: n (" + std::to_string(this->n)
                                  + ") must be >= 2 * |r|_1 (" + std::to_string(2 * this->d) + ")");
    }
  }

  std::vector<int> DistanceRootedDiagramGenerator::bfs_all_pairs(Graph const &G) {
    int V = G.get_V();
    std::vector<int> dist(static_cast<size_t>(V) * V, -1);
    for (int src = 0; src < V; ++src) {
      int *row = dist.data() + src * V;
      row[src] = 0;
      std::queue<int> q;
      q.push(src);
      while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (int v = 0; v < V; ++v) {
          if (v == u) continue;
          // Treat the multigraph as a simple graph for shortest-path purposes:
          // any nonzero adjacency entry is a unit edge.
          if ((G(u, v) != 0 || G(v, u) != 0) && row[v] < 0) {
            row[v] = row[u] + 1;
            q.push(v);
          }
        }
      }
    }
    return dist;
  }

  void DistanceRootedDiagramGenerator::try_emit(Graph const &G, std::vector<int> marks, int V,
                                                std::unordered_set<RootedKey, RootedKeyHasher> &unique_keys) {
    RootedGraph rg(G, std::move(marks), /*defer_embedding=*/true);
    RootedKey key{rg.get_canonical_form(), rg.get_marks()};
    if (unique_keys.insert(key).second) this->rooted_graphs_by_V[V].push_back(std::move(rg));
  }

  void DistanceRootedDiagramGenerator::generate() {
    VacuumDiagramGenerator gen(this->n, this->bipartite_only);

    for (int V = this->d + 1; V <= this->n; ++V) {
      gen.generate_at_vertex_count(V);

      // Per-V dedup set: different V's never collide (canonical adjacency
      // size differs), so scoping to V keeps the set small and avoids cross-V
      // bookkeeping.
      std::unordered_set<RootedKey, RootedKeyHasher> unique_keys;

      for (auto const &G : gen.get_unique_graphs()) {
        std::vector<int> dist_G = bfs_all_pairs(G);

        // Required diameter prune: V >= 2d is necessary but not sufficient
        // for diameter >= d (e.g. a 2d-cycle with an extra chord shortens
        // diameter below d). One BFS pass already done — just scan dist_G.
        int diameter = 0;
        for (int v : dist_G) diameter = std::max(diameter, v);
        if (diameter < this->d) continue;

        // Parity filter: when both the lattice AND the abstract graph G are
        // bipartite, every walk from i to j on the lattice has length parity
        // matching graph-distance(i, j) parity; the lattice's bipartite
        // coloring then forces Manhattan(phi(i), phi(j)) to share that
        // parity. If G has an odd cycle, walks of either parity exist
        // between some vertex pairs, so the filter would be unsound — we
        // drop it for those graphs. Non-bipartite lattices also drop it.
        bool use_parity_filter = this->bipartite_only && G.get_bipartite();

        // Enumerate unordered mark pairs {i, j} (i <= j) whose graph-distance
        // matches the target Manhattan distance up to parity. Same-vertex
        // pairs (i == j) are needed for the on-site case (d == 0): both
        // density operators on one hopping vertex. For d > 0 the dij < d
        // filter prunes (i, i) automatically. The n-cycle (V == n)
        // contributes multiple distance classes here; canonical dedup
        // collapses pairs equivalent under Aut(G), so e.g. the n-cycle
        // surfaces one rooted topology per distinct graph-distance >= d.
        for (int i = 0; i < V; ++i) {
          for (int j = i; j < V; ++j) {
            int dij = dist_G[i * V + j];
            if (dij < this->d) continue;
            if (use_parity_filter && ((dij - this->d) % 2 != 0)) continue;
            this->try_emit(G, {i, j}, V, unique_keys);
          }
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // DimerDistanceRootedDiagramGenerator
  // ---------------------------------------------------------------------------

  // Floor-divide-by-2 of a non-negative integer parity. Used only on values
  // derived from the parity-allowed sector check, so the numerator is always
  // even here — but std::abs would also work.
  static inline int positive_mod2(int x) { return ((x % 2) + 2) % 2; }

  int DimerDistanceRootedDiagramGenerator::d_super(int rx, int ry, int ms0, int ms1) {
    int du     = ms0 - ms1 + rx;
    int dv     = ry;
    int abs_du = std::abs(du);
    int abs_dv = std::abs(dv);
    // d_super = max(|Δv|, (|Δû| + |Δv|) / 2). Parity-allowed sectors have
    // (ms0 + ms1) ≡ (rx + ry) mod 2 ⟹ (Δû + Δv) even ⟹ (|Δû| + |Δv|) even,
    // so the integer division is exact.
    return std::max(abs_dv, (abs_du + abs_dv) / 2);
  }

  std::vector<int> DimerDistanceRootedDiagramGenerator::bfs_all_pairs(Graph const &G) {
    int V = G.get_V();
    std::vector<int> dist(static_cast<size_t>(V) * V, -1);
    for (int src = 0; src < V; ++src) {
      int *row = dist.data() + src * V;
      row[src] = 0;
      std::queue<int> q;
      q.push(src);
      while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (int v = 0; v < V; ++v) {
          if (v == u) continue;
          if ((G(u, v) != 0 || G(v, u) != 0) && row[v] < 0) {
            row[v] = row[u] + 1;
            q.push(v);
          }
        }
      }
    }
    return dist;
  }

  DimerDistanceRootedDiagramGenerator::DimerDistanceRootedDiagramGenerator(std::vector<int> r_, int n_)
     : r(std::move(r_)), n(n_) {
    if (this->r.size() != 2)
      throw std::invalid_argument("DimerDistanceRootedDiagramGenerator: r must have size 2");

    int rx       = this->r[0];
    int ry       = this->r[1];
    int r_parity = positive_mod2(rx + ry);

    // Enumerate the 2 parity-allowed sectors (out of 4 (ms₀, ms₁) ∈ {0,1}²)
    // with their d_super values. The other 2 sectors have non-integer u' for
    // mark₁'s dimer anchor and are off-lattice.
    for (int ms0 = 0; ms0 < 2; ++ms0) {
      for (int ms1 = 0; ms1 < 2; ++ms1) {
        if (positive_mod2(ms0 + ms1) != r_parity) continue;
        int d = d_super(rx, ry, ms0, ms1);
        this->allowed_sectors.emplace_back(std::make_pair(ms0, ms1), d);
      }
    }
    if (this->allowed_sectors.empty()) {
      // Cannot occur — every r-parity admits exactly 2 of 4 sectors. Guard
      // against future refactors that change the parity logic.
      throw std::runtime_error("DimerDistanceRootedDiagramGenerator: no parity-allowed sectors (logic bug?)");
    }

    // Catalog-level n_min = 2 · min over sectors of d_super. If n is below
    // every sector's n_min, no graph at this order can host both marks at
    // displacement r.
    int n_min = std::numeric_limits<int>::max();
    for (auto const &kv : this->allowed_sectors) n_min = std::min(n_min, 2 * kv.second);
    if (this->n < n_min) {
      throw std::invalid_argument("DimerDistanceRootedDiagramGenerator: n (" + std::to_string(this->n)
                                  + ") is below n_min = 2 * min(d_super) (" + std::to_string(n_min) + ") for r");
    }
  }

  void DimerDistanceRootedDiagramGenerator::try_emit(Graph const &G, std::vector<int> marks, std::vector<int> sites, int V,
                                                       std::unordered_set<DimerRootedKey, DimerRootedKeyHasher> &unique_keys) {
    DimerRootedGraph rg(G, std::move(marks), std::move(sites));
    DimerRootedKey key{rg.get_canonical_form(), rg.get_marks(), rg.get_sites()};
    if (unique_keys.insert(key).second) this->rooted_graphs_by_V[V].push_back(std::move(rg));
  }

  void DimerDistanceRootedDiagramGenerator::generate() {
    // Dimer expansion always uses the full (non-bipartite) graph set — the
    // staggered superlattice is triangular and admits odd cycles.
    VacuumDiagramGenerator gen(this->n, /*bipartite_only=*/false);

    // Smallest V_min across allowed sectors. VacuumDiagramGenerator skips
    // V == 1 (would require self-loops, not generated), so the V = 1 case
    // (intra-dimer / on-site correlator at zeroth order) is not produced by
    // this generator — that contribution is computed separately.
    int v_min_overall = std::numeric_limits<int>::max();
    for (auto const &kv : this->allowed_sectors) v_min_overall = std::min(v_min_overall, kv.second + 1);
    if (v_min_overall < 2) v_min_overall = 2;

    for (int V = v_min_overall; V <= this->n; ++V) {
      gen.generate_at_vertex_count(V);

      // Per-V dedup set: canonical adjacency size differs across V's, so
      // scoping the set to V keeps it small and avoids cross-V bookkeeping.
      std::unordered_set<DimerRootedKey, DimerRootedKeyHasher> unique_keys;

      for (auto const &G : gen.get_unique_graphs()) {
        std::vector<int> dist_G = bfs_all_pairs(G);

        // Required diameter prune: skip the graph if no allowed sector's
        // d_super is reachable given the graph diameter.
        int diameter = 0;
        for (int v : dist_G) diameter = std::max(diameter, v);

        bool graph_ok = false;
        for (auto const &kv : this->allowed_sectors) {
          if (V >= kv.second + 1 && diameter >= kv.second) {
            graph_ok = true;
            break;
          }
        }
        if (!graph_ok) continue;

        // Enumerate unordered mark pairs {i, j} (i <= j). Same-vertex pairs
        // (i == j) are needed for the on-dimer / on-site correlator (small r,
        // d_super == 0) — the dij < d_super filter prunes them when d_super
        // > 0. No graph-bipartite parity filter (superlattice is triangular).
        for (int i = 0; i < V; ++i) {
          for (int j = i; j < V; ++j) {
            int dij = dist_G[i * V + j];
            if (dij < 0) continue; // shouldn't happen for a connected G

            // Each parity-allowed sector adds an emit candidate when its
            // d_super is satisfied. The DimerRootedGraph constructor's
            // colored canonicalisation handles dedup across input sectors
            // that map to the same canonical form (e.g. sector (0,1) and
            // (1,0) on a graph whose auto swaps the two marks both produce
            // the same canonical form — bliss merges them via colors).
            for (auto const &kv : this->allowed_sectors) {
              int ms0 = kv.first.first;
              int ms1 = kv.first.second;
              int d_s = kv.second;
              if (V < d_s + 1) continue;
              if (dij < d_s) continue;
              this->try_emit(G, {i, j}, {ms0, ms1}, V, unique_keys);
            }
          }
        }
      }
    }
  }

} // namespace sc_expansion
