#include "./hubbard_dimer.hpp"
#include <triqs/atom_diag/functions.hpp>
using namespace triqs::operators;

namespace hubbard_dimer {

  triqs::hilbert_space::fundamental_operator_set make_fops() {
    triqs::hilbert_space::fundamental_operator_set fops;
    for (int o : triqs::arrays::range(2)) fops.insert("dn", o);
    for (int o : triqs::arrays::range(2)) fops.insert("up", o);
    return fops;
  }

  auto make_H0(double U, double mu) {

    auto h = U * (n("up", 0) * n("dn", 0) + n("up", 1) * n("dn", 1));
    h += -mu * (n("up", 0) + n("dn", 0) + n("up", 1) + n("dn", 1));

    return h;
  }

  auto boltz_weights_and_imag_time(auto energies, double tau1, double tau2) {

    //diagonal matrix e^{-(t1-t2) H0}
    //energies in the subspace
    std::vector<double> weights = std::vector<double>(energies.size());

    for (int i = 0; i < energies.size(); i++) { weights[i] = std::exp(-(tau1 - tau2) * energies[i]); }
    return nda::diag(weights);
  }

  double partition_function(auto energies, double beta) {

    //get partition function
    double Z = 0.0; //partition function

    for (int sp_index = 0; sp_index < energies.size(); sp_index++) {
      for (int i = 0; i < energies[sp_index].size(); i++) {
        Z += std::exp(-beta * energies[sp_index][i]); //sum over all states in the subspace
      }
    }
    return Z;
  }

  double occupation(double U, double mu, double beta) {

    auto fops = make_fops();                                  //hilbert space
    auto H0   = make_H0(U, mu);                               //unperturbed Hamiltonian
    auto ad   = triqs::atom_diag::atom_diag<false>(H0, fops); //atom_diag object

    auto energies = ad.get_energies();

    //get partition function
    double Z = partition_function(energies, beta); //partition function
    std::cout << "Partition function: " << Z << std::endl;

    int n_sub = energies.size(); //number of subspaces

    double occ = 0.0; //trace of the occupation operator

    for (int sp_index = 0; sp_index < n_sub; sp_index++) {

      auto Esp = energies[sp_index]; //energies in the subspace

      auto boltz = boltz_weights_and_imag_time(Esp, beta, 0.0); //boltzmann weights

      for (int i = 0; i < 4; i++) {
        int c_conn = ad.c_connection(i, sp_index);
        if (c_conn < 0) { continue; } //skip if no connection
        auto Nis = boltz * ad.cdag_matrix(i, c_conn) * ad.c_matrix(i, sp_index);
        occ += trace(Nis);
      }
    }

    return occ / Z;
  }

  double unperturbed_green_function(double U, double mu, double beta, std::vector<std::pair<int, int>> path, std::vector<double> tau) {

    auto fops = make_fops();                                  //hilbert space
    auto H0   = make_H0(U, mu);                               //unperturbed Hamiltonian
    auto ad   = triqs::atom_diag::atom_diag<false>(H0, fops); //atom_diag object

    auto energies = ad.get_energies();

    //get partition function
    double Z0 = partition_function(energies, beta); //partition function of H0
    int n_sub = energies.size();                    //number of subspaces

    double gf = 0.0;

    for (int sp_index = 0; sp_index < n_sub; sp_index++) {
      bool flag = false;
      auto Esp  = energies[sp_index];                          //energies in the subspace
      auto op   = boltz_weights_and_imag_time(Esp, beta, 0.0); //boltzmann weights

      for (size_t op_idx = 0; op_idx < path.size(); ++op_idx) { //loop through path operators
        int is1 = path[op_idx].first;                           //cdag(i, sigma1)
        int js2 = path[op_idx].second;                          //c(j,sigma2)

        int c_conn = ad.c_connection(js2, sp_index);
        if (c_conn < 0) {
          flag = true;
          break;
        } //contribution vanishes if no connection

        op *= ad.cdag_matrix(is1, c_conn) * ad.c_matrix(js2, sp_index);
      }

      //if (flag) { continue; }

      gf += trace(op);
    }

    return gf / Z0; //return the unperturbed Green's function
  };

} // namespace hubbard_dimer