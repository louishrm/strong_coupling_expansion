#include "diagram.hpp"
#include "hubbard_atom.hpp"
#include "cumulant.hpp"

namespace sc_expansion {

class order4 {

  public:

  adjmat D4a                    = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}; //4-cycle
  adjmat D4b                    = {{0, 1, 1}, {1, 0, 0}, {1, 0, 0}};                        //3-cycle with double lines
  adjmat D4c                    = {{0, 2}, {2, 0}};                                         //2-cycle with double lines

  // atom diag
  triqs::atom_diag::atom_diag<false> ad;
  std::vector<Diagram> diagrams;
  double beta;

  order4(double U, double mu, double beta_) :
    ad(hubbard_atom::make_H0(U, mu), hubbard_atom::make_fops(), {}),
    diagrams({Diagram(D4a), Diagram(D4b), Diagram(D4c)}),
    beta(beta_) {}

  double compute_sum_diagrams(std::vector<double> taus) {

    double diagram_sum = 0.0;
    for (auto const &diagram : diagrams) {
      double val = diagram.evaluate_at_taus(ad, beta, taus);
      diagram_sum += val;
    }
    return diagram_sum;

  }

};

}
