#pragma once

#include <vector>
#include <numeric>   // For std::accumulate
#include <algorithm> // For std::next_permutation
#include <queue>
#include "./hubbard_atom.hpp"
#include "./cumulant.hpp"

namespace sc_expansion {

  using adjmat = std::vector<std::vector<int>>; //adjacency matrix is a VxV matrix where V is the number of vertices.
  //Each entry Aij is the number of directed lines from i to j.

  class Diagram {

    public:
    Diagram(adjmat adjacency_matrix, double U, double beta, double mu);

    bool is_connected() const;
    bool is_particle_number_conserving() const;
    Diagram get_canonical_form() const;

    struct Line {
      int from_vertex;
      int to_vertex;
    };
    std::vector<Line> get_hopping_lines() const;

    double evaluate_at_points(HubbardAtom::cumul_args const &args) const;

    double evaluate_at_taus(std::vector<double> const &taus) const;

    int diagram_sign() const;
    int get_symmetry_factor() const;
    int get_free_multiplicity() const;

    private:
    int compute_symmetry_factor() const;
    int compute_diagram_sign() const;
    int compute_free_multiplicity() const;
    std::vector<Line> compute_hopping_lines() const;
    std::vector<std::vector<int>> compute_valid_spin_configurations() const;

    int n; //order= number of hopping lines
    int V; //number of vertices
    adjmat adjacency_matrix;
    HubbardAtom atom;

    double sign;
    double symmetry_factor;
    double fm; //free multiplicity on inifnite 2d square lattice.
    std::vector<Line> hopping_lines;
    std::vector<std::vector<int>> valid_spin_configurations;
  };

  struct Point {
    int x;
    int y;

    Point();
    Point(int x, int y);
  };

  struct SquareLattice {

    std::vector<Point> get_neighbors(Point const &r) const;

    int manhattan_distance(Point const &r) const;
    bool prune(Point const &r, int current_distance, int order) const;
    bool is_neighbor(Point const &r1, Point const &r2) const;

    SquareLattice();
  };

  void next_step(adjmat &A, std::vector<int> &sequence, int vertex, int V, int order);
  std::vector<int> generate_walk_sequence(adjmat A, int order);

  void place_next_vertex(std::unordered_map<int, Point> &placed_vertices, SquareLattice const &lattice, std::vector<int> const &sequence,
                         int current_vertex_index, int order, int &free_multiplicity, int hopping_count);

  int compute_free_multiplicity(adjmat A, int order);

} // namespace sc_expansion
