#include "diagram2.hpp"

namespace sc_expansion {

  template <int N_sites, typename T>
  Diagram2<N_sites, T>::Diagram2(Graph const &graph_, std::vector<VertexType<N_sites, T> *> const &vertex_types) : graph(graph_) {
    this->compute_hopping_lines();
    this->setup_vertices(vertex_types);
    this->compute_valid_configurations();
    this->compute_diagram_sign();
  }

  template <int N_sites, typename T> void Diagram2<N_sites, T>::compute_hopping_lines() {

    for (int i = 0; i < this->graph.get_V(); i++) {
      for (int j = 0; j < this->graph.get_V(); j++) {

        int lines_ij = this->graph(i, j);
        if (lines_ij != 0) {

          for (int k = 0; k < lines_ij; k++) {
            HoppingLine line;
            line.from_vertex = i;
            line.to_vertex   = j;
            this->hopping_lines.lines.push_back(line);
          }
        }
      }
    }
  }

} // namespace sc_expansion