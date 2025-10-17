#include <vector>
#include <numeric>   // For std::accumulate
#include <algorithm> // For std::next_permutation
#include <queue>

using adjmat = std::vector<std::vector<int>>; //adjacency matrix is a VxV matrix where V is the number of vertices.
//Each entry Aij is the number of directed lines from i to j.

class Diagram {

  public:
  Diagram(const adjmat &adjacency_matrix);

  bool is_connected() const;
  bool is_particle_number_conserving() const;
  Diagram get_canonical_form() const;
  int get_symmetry_factor() const;

  private:
  int n; //order= number of hopping lines
  int V; //number of vertices
  adjmat adjacency_matrix;
};