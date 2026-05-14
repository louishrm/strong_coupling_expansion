#include "diagram_common.hpp"
#include <cmath>
#include <cstdlib>
#include <map>

namespace sc_expansion {

  // ---------------------------------------------------------------------------
  // Lattice embedding count (single-site convention).
  // Bipartite graphs embed on the square lattice (4 NN); non-bipartite on the
  // triangular lattice (6 NN).
  //
  // The recursive body is shared between vacuum (single-site) and rooted
  // (marked-vertex) embeddings via an optional "track_vertex"/histogram pair:
  //   - vacuum: histogram == nullptr, returns total count.
  //   - rooted 2-mark: track_vertex = canonical index of mark[1], histogram
  //     accumulates d^2 = Δx²+Δy² of mark[1] relative to the anchored mark[0].
  // ---------------------------------------------------------------------------

  namespace {

    struct Point {
      int x, y;
    };

    long solve_embedding_recursive(Graph const &graph, bool bipartite_only, int V, int placed_count, std::vector<bool> &placed,
                                   std::vector<Point> &coords, int track_vertex, std::map<int, int> *histogram) {
      if (placed_count == V) {
        if (histogram) {
          Point p = coords[track_vertex];
          (*histogram)[p.x * p.x + p.y * p.y]++;
        }
        return 1;
      }

      int anchor = -1, target_node = -1;
      for (int candidate = 0; candidate < V; ++candidate) {
        if (!placed[candidate]) {
          for (int p = 0; p < V; ++p) {
            if (placed[p]) {
              uint8_t val = graph(candidate, p) + graph(p, candidate);
              if (val > 0) {
                target_node = candidate;
                anchor      = p;
                goto found_target;
              }
            }
          }
        }
      }
    found_target:;
      if (target_node == -1) return 0;

      long count = 0;

      static constexpr int sq_dx[4]  = {1, -1, 0, 0};
      static constexpr int sq_dy[4]  = {0, 0, 1, -1};
      static constexpr int tri_dx[6] = {1, -1, -1, 0, 0, 1};
      static constexpr int tri_dy[6] = {0, 0, 1, 1, -1, -1};

      const int *dx    = bipartite_only ? sq_dx : tri_dx;
      const int *dy    = bipartite_only ? sq_dy : tri_dy;
      int n_dirs       = bipartite_only ? 4 : 6;
      Point anchor_pos = coords[anchor];

      auto is_neighbor = [bipartite_only](Point p1, Point p2) -> bool {
        int ddx = p1.x - p2.x;
        int ddy = p1.y - p2.y;
        if (bipartite_only) { return std::abs(ddx) + std::abs(ddy) == 1; }
        if (std::abs(ddx) + std::abs(ddy) == 1) return true;
        if (ddx == -1 && ddy == 1) return true;
        if (ddx == 1 && ddy == -1) return true;
        return false;
      };

      for (int dir = 0; dir < n_dirs; ++dir) {
        Point candidate_pos = {anchor_pos.x + dx[dir], anchor_pos.y + dy[dir]};

        bool valid = true;
        for (int i = 0; i < V; ++i) {
          if (placed[i]) {
            uint8_t links = graph(target_node, i) + graph(i, target_node);
            if (links > 0) {
              if (!is_neighbor(candidate_pos, coords[i])) {
                valid = false;
                break;
              }
            }
          }
        }

        if (valid) {
          coords[target_node] = candidate_pos;
          placed[target_node] = true;
          count += solve_embedding_recursive(graph, bipartite_only, V, placed_count + 1, placed, coords, track_vertex, histogram);
          placed[target_node] = false;
        }
      }
      return count;
    }

  } // namespace

  int compute_lattice_free_multiplicity(Graph const &graph) {
    int V = graph.get_V();
    std::vector<Point> coords(V, {0, 0});
    std::vector<bool> placed(V, false);
    coords[0] = {0, 0};
    placed[0] = true;

    bool bipartite_only = graph.get_bipartite_only();
    return (int)solve_embedding_recursive(graph, bipartite_only, V, 1, placed, coords, /*track*/ -1, /*hist*/ nullptr);
  }

  // Rooted embeddings: anchor mark0_vertex at the origin and recurse. If
  // mark1_vertex >= 0, the returned map buckets final positions of mark1 by
  // squared distance d² = Δx² + Δy²; sum(values) == total embedding count.
  // For the single-mark case (mark1_vertex < 0), the map has one entry at
  // d² = 0 equal to the total count.
  std::map<int, int> compute_rooted_shell_multiplicity(Graph const &graph, int mark0_vertex, int mark1_vertex) {
    int V = graph.get_V();
    std::vector<Point> coords(V, {0, 0});
    std::vector<bool> placed(V, false);
    coords[mark0_vertex] = {0, 0};
    placed[mark0_vertex] = true;

    std::map<int, int> histogram;
    bool bipartite_only = graph.get_bipartite_only();
    long total          = solve_embedding_recursive(graph, bipartite_only, V, 1, placed, coords, mark1_vertex,
                                           mark1_vertex >= 0 ? &histogram : nullptr);
    if (mark1_vertex < 0) histogram[0] = static_cast<int>(total);
    return histogram;
  }

  // ---------------------------------------------------------------------------
  // Hopping-line and leg bookkeeping.
  // ---------------------------------------------------------------------------

  Lines compute_hopping_lines(Graph const &graph) {
    Lines out;
    int V = graph.get_V();
    for (int i = 0; i < V; ++i) {
      for (int j = 0; j < V; ++j) {
        int n_ij = graph(i, j);
        for (int k = 0; k < n_ij; ++k) { out.lines.push_back({i, j}); }
      }
    }
    return out;
  }

  std::vector<std::vector<LegInfo>> compute_legs_per_vertex(Graph const &graph, Lines const &hopping_lines) {
    int V = graph.get_V();
    std::vector<std::vector<LegInfo>> out(V);
    for (int k = 0; k < (int)hopping_lines.lines.size(); ++k) {
      out[hopping_lines.lines[k].from_vertex].push_back({k, true});
      out[hopping_lines.lines[k].to_vertex].push_back({k, false});
    }
    return out;
  }

  int compute_diagram_sign(int V, Lines const &hopping_lines, std::vector<std::vector<LegInfo>> const &legs_per_vertex) {
    int n_lines = (int)hopping_lines.lines.size();
    std::vector<int> successor_map(n_lines);

    for (int v = 0; v < V; ++v) {
      std::vector<int> incoming, outgoing;
      for (auto const &leg : legs_per_vertex[v]) {
        if (leg.is_source) outgoing.push_back(leg.line_index);
        else incoming.push_back(leg.line_index);
      }
      for (size_t i = 0; i < incoming.size(); ++i) successor_map[incoming[i]] = outgoing[i];
    }

    int num_loops = 0;
    std::vector<bool> visited(n_lines, false);
    for (int i = 0; i < n_lines; ++i) {
      if (!visited[i]) {
        ++num_loops;
        int current = i;
        while (!visited[current]) {
          visited[current] = true;
          current          = successor_map[current];
        }
      }
    }
    return (num_loops % 2 == 0) ? 1 : -1;
  }

} // namespace sc_expansion
