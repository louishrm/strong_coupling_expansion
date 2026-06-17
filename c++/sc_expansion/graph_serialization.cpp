#include "graph_serialization.hpp"
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
  // SC_EXPANSION_PROJECT_ROOT is injected at configure time (see src/c++/sc_expansion/CMakeLists.txt)
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
      int32_t V    = g.get_V();
      int32_t aut  = g.get_automorphism_count();
      int32_t sym  = static_cast<int32_t>(g.get_symmetry_factor());
      int32_t fm   = static_cast<int32_t>(g.get_free_multiplicity());
      uint8_t bonly = g.get_bipartite_only() ? 1 : 0;
      auto adj     = g.get_canonical_form();

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

} // namespace sc_expansion
