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
       : diagrams({Diagram(this->D4a, U, beta_, mu), Diagram(this->D4b, U, beta_, mu), Diagram(this->D4c, U, beta_, mu)}), beta(beta_) {}

    double compute_sum_diagrams(std::vector<double> taus) {

      double diagram_sum = 0.0;
      for (auto const &diagram : this->diagrams) {
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

    order2(double U, double mu, double beta_) : diagrams({Diagram(this->D2a, U, beta_, mu)}), beta(beta_) {}

    double compute_sum_diagrams(std::vector<double> taus) {

      double diagram_sum = 0.0;
      for (auto const &diagram : this->diagrams) {
        double val = diagram.evaluate_at_taus(taus);
        diagram_sum += val;
      }
      return diagram_sum;
    }
  };

  class order6 {

    public:
    adjmat D6a = {{0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}, {1, 0, 0, 0, 0, 0}}; //6-cycle
    adjmat D6b = {{0, 3}, {3, 0}};                                                                      //watermelon triple
    adjmat D6c = {{0, 1, 1, 1}, {1, 0, 0, 0}, {1, 0, 0, 0}, {1, 0, 0, 0}};                              //petal with 4 vertice
    adjmat D6d = {{0, 1, 1, 0, 0}, {1, 0, 0, 0, 0}, {0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}, {1, 0, 0, 0, 0}}; //square +digon
    adjmat D6e = {{0, 1, 1, 0}, {1, 0, 0, 1}, {1, 0, 0, 0}, {0, 1, 0, 0}};                              //crab diagram
    adjmat D6f = {{0, 2, 1}, {2, 0, 0}, {1, 0, 0}};                                                     //watermelon double + digon
    adjmat D6g = {{0, 2, 0, 0}, {1, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}; //square with one double line                           //chain of 4 digons

    std::vector<Diagram> diagrams;
    double beta;

    order6(double U, double mu, double beta_)
       : diagrams({Diagram(this->D6a, U, beta_, mu), Diagram(this->D6b, U, beta_, mu), Diagram(this->D6c, U, beta_, mu),
                   Diagram(this->D6d, U, beta_, mu), Diagram(this->D6e, U, beta_, mu), Diagram(this->D6f, U, beta_, mu),
                   Diagram(this->D6g, U, beta_, mu)}),
         beta(beta_) {}

    double compute_sum_diagrams(std::vector<double> taus) {

      double diagram_sum = 0.0;
      for (auto const &diagram : this->diagrams) {
        double val = diagram.evaluate_at_taus(taus);
        diagram_sum += val;
      }
      return diagram_sum;
    }
  };
} // namespace sc_expansion