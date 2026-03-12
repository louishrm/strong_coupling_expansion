#include "./hubbard_dimer.hpp"

namespace sc_expansion {

  using namespace triqs::operators;

  template <typename T> HubbardDimer<T>::HubbardDimer(T t, T U, T beta, T mu) {

    //diagonalize with triqs
    triqs::hilbert_space::fundamental_operator_set fops = {};
    for (int o : triqs::arrays::range(2)) fops.insert("dn", o);
    for (int o : triqs::arrays::range(2)) fops.insert("up", o);

    triqs::operators::many_body_operator_generic<T> h;

    for (int i = 0; i < 2; ++i) {
      h += -mu * (n("up", i) + n("dn", i));
      h += U * n("up", i) * n("dn", i);
    }

    h += -t * (c_dag("up", 0) * c("up", 1) + c_dag("up", 1) * c("up", 0) + c_dag("dn", 0) * c("dn", 1) + c_dag("dn", 1) * c("dn", 0));

    triqs::atom_diag::atom_diag<false> ad(h, fops);
    this->Z = triqs::atom_diag::partition_function(ad, beta);
    std::cout << "Hubbard Dimer partition function Z: " << this->Z << std::endl;
  }

  template <typename T> T HubbardDimer<T>::G0(std::vector<double> const &taus, std::vector<int> const &spins, std::vector<int> const &sites) const {
    // Placeholder for future implementation
    return 0.0; // Placeholder for G0 value
  }

  template class HubbardDimer<double>;

  } // namespace sc_expansion