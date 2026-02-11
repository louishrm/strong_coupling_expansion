#pragma once

#include <vector>
#include <numeric>   // For std::accumulate
#include <algorithm> // For std::next_permutation
#include <queue>
#include "./hubbard_atom.hpp"
#include "./cumulant.hpp"
#include "./graph.hpp" // Include Graph for compute_free_multiplicity logic

namespace sc_expansion {

  struct Line {
    int from_vertex;
    int to_vertex;
  };

  struct Parameters {
    double U;
    double beta;
    double mu;
  };

  class Diagram {

    public:
    explicit Diagram(Graph const &graph);

    double get_diagram_sign() const { return (double)this->diagram_sign; }
    std::vector<long> get_valid_spin_configurations() const { return this->valid_spin_configurations; }
    std::vector<Line> get_hopping_lines() const { return this->hopping_lines; }
    Graph const &get_graph() const { return this->graph; }

    private:
    void compute_hopping_lines();
    void compute_valid_spin_configurations();
    void compute_diagram_sign();

    Graph graph;
    std::vector<Line> hopping_lines;
    std::vector<long> valid_spin_configurations;
    int diagram_sign;
  };

  class DiagramEvaluator {

    public:
    explicit DiagramEvaluator(Diagram const &diagram, Parameters const &params);

    double evaluate_at_taus(std::vector<double> const &taus, bool infinite_U) const;

    private:
    double evaluate_at_points(HubbardAtom::cumul_args const &args, bool infinite_U) const;

    Diagram const &diagram;
    HubbardAtom atom;
  };

} // namespace sc_expansion
