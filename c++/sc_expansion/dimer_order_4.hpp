#include "diagram.hpp"
#include "hubbard_atom.hpp"
#include "cumulant.hpp"

namespace sc_expansion {

  class order4 {

    public:
    adjmat D4a = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}; //4-cycle
    adjmat D4b = {{0, 1, 1}, {1, 0, 0}, {1, 0, 0}};                        //3-cycle with double lines
    adjmat D4c = {{0, 2}, {2, 0}};                                         //2-cycle with double lines

    std::vector<Diagram> diagrams;
    double beta;

    order4(double U, double mu, double beta_)
       : diagrams({Diagram(D4a, U, beta_, mu), Diagram(D4b, U, beta_, mu), Diagram(D4c, U, beta_, mu)}), beta(beta_) {}

    double compute_sum_diagrams(std::vector<double> taus) {

      double diagram_sum = 0.0;
      for (auto const &diagram : diagrams) {
        double val = diagram.evaluate_at_taus(taus);
        diagram_sum += val;
      }
      return diagram_sum;
    }
  };

  class order2 {

    public:
    adjmat D2a = {{0, 1}, {1, 0}}; //2-cycle

    std::vector<Diagram> diagrams;
    double beta;

    order2(double U, double mu, double beta_) : diagrams({Diagram(D2a, U, beta_, mu)}), beta(beta_) {}

    double compute_sum_diagrams(std::vector<double> taus) {

      double diagram_sum = 0.0;
      for (auto const &diagram : diagrams) {
        double val = diagram.evaluate_at_taus(taus);
        diagram_sum += val;
      }
      return diagram_sum;
    }
  }; // namespace sc_expansion
}; // namespace sc_expansion
