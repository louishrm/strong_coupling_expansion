#include "./hubbard_atom.hpp"

namespace hubbard_atom {

  triqs::hilbert_space::fundamental_operator_set make_fops() {
    triqs::hilbert_space::fundamental_operator_set fops = {};
    for (int o : triqs::arrays::range(1)) fops.insert("dn", o);
    for (int o : triqs::arrays::range(1)) fops.insert("up", o);
    return fops;
  }

  triqs::operators::many_body_operator_generic<double> make_H0(double U, double mu) {

    triqs::operators::many_body_operator_generic<double> h = U * n("up", 0) * n("dn", 0);
    h += -mu * (n("up", 0) + n("dn", 0));

    return h;
  }

  double partition_function(triqs::atom_diag::atom_diag<false> ad, double beta) {

    double Z0        = triqs::atom_diag::partition_function(ad, beta); //Z0 from atom_diag
    double gs_energy = ad.get_gs_energy();
    Z0 *= std::exp(-beta * gs_energy); //Z0 = Z0 * exp(-beta * E_0)
    return Z0;
  }

  std::tuple<std::vector<double>, std::vector<int>, int> imag_time_sort_and_sign(std::vector<double> tau) {

    //sort imaginary times tau, return the sorted times, parity and the argsort
    std::vector<int> argsort = {};
    for (int i = 0; i < tau.size(); i++) argsort.push_back(i); //create an index vector

    int n      = tau.size(); //number of imaginary times
    int parity = 0;          //parity of the permutation

    for (int i = 0; i < n - 1; i++) {
      for (int j = i + 1; j < n; j++) {
        if (tau[i] < tau[j]) {
          std::swap(tau[i], tau[j]);         //swap the imaginary times
          std::swap(argsort[i], argsort[j]); //swap the index vector
          parity++;                          //increment the parity
        }
      }
    }
    return std::make_tuple(tau, argsort, std::pow(-1, parity));
  }

  bool isSpinConserved(const std::vector<int> &spins, const std::vector<int> &flags) {
    int netSpinUp   = 0;
    int netSpinDown = 0;

    for (size_t i = 0; i < spins.size(); ++i) {
      if (flags[i] == 0) { // Adding a spin
        if (spins[i] == 1) {
          netSpinUp++;
        } else {
          netSpinDown++;
        }
      } else if (flags[i] == 1) { // Removing a spin
        if (spins[i] == 1) {
          netSpinUp--;
        } else {
          netSpinDown--;
        }
      }
    }

    return netSpinUp == 0 && netSpinDown == 0;
  }

  nda::matrix<double> make_interaction_picture_destroy_op(triqs::atom_diag::atom_diag<false> ad, double tau, int state_index) {

    double Z01 = partition_function(ad, tau);  //Z0
    double Z02 = partition_function(ad, -tau); //Z0

    auto cmat        = ad.c_matrix(state_index, 0);
    auto time_evol_1 = triqs::atom_diag::atomic_density_matrix(ad, -tau)[0]; //time evolution operator
    // std::cout << "cmat: " << cmat << std::endl;

    auto time_evol_2 = triqs::atom_diag::atomic_density_matrix(ad, tau)[0]; //time evolution operator
    auto ctau        = Z01 * Z02 * time_evol_1 * cmat * time_evol_2;        //interaction picture destroy operator
    //std::cout << "Interaction picture destroy operator: " << ctau << std::endl;
    return ctau; //return the interaction picture destroy operator
  }

  nda::matrix<double> make_interaction_picture_create_op(triqs::atom_diag::atom_diag<false> ad, double tau, int state_index) {

    double Z01 = partition_function(ad, tau);  //Z0
    double Z02 = partition_function(ad, -tau); //Z0

    auto cmat        = ad.cdag_matrix(state_index, 0);
    auto time_evol_1 = triqs::atom_diag::atomic_density_matrix(ad, -tau)[0]; //time evolution operator

    // std::cout << "cdagmat: " << cmat << std::endl;
    // std::cout << "time evolution operator 1: " << time_evol_1 << std::endl;

    auto time_evol_2 = triqs::atom_diag::atomic_density_matrix(ad, tau)[0]; //time evolution operator

    auto cdagtau = Z01 * Z02 * time_evol_1 * cmat * time_evol_2; //interaction picture destroy operator

    // std::cout << "Interaction picture create operator: " << cdagtau << std::endl;

    return cdagtau; //return the interaction picture destroy operator
  }

  double G0(triqs::atom_diag::atom_diag<false> ad, double beta, std::vector<double> times, std::vector<int> spins, std::vector<int> flags) {

    //unperturbed n body local Green's function: < Tc(tau_1,s_1) ... c(tau_n,s_n) \ cdag(tau_1', s_1')... cdag(tau_n', s_n')>_0

    //times is a vector of imaginary times
    //every imaginary time has a corresponding spin and flag.
    //spins is a vector of spins, 0 for "dn" and 1 for "up".
    //flags is a vector of flags, 0 for create and 1 for destroy.

    //if (isSpinConserved(spins, flags) == false) { return 0.0; } //if the spin is not conserved, return 0

    nda::matrix<double> rho0 = triqs::atom_diag::atomic_density_matrix(ad, beta)[0]; //rho_0

    //first, sort the times, get the argsort and the sign of the UGF
    auto [sorted_times, argsort, sign] = imag_time_sort_and_sign(times);

    //now sort the spins according to the argsort
    std::vector<int> sorted_spins;
    std::vector<int> sorted_flags; 
    sorted_spins.reserve(spins.size());
    sorted_flags.reserve(flags.size());
    for (auto i : argsort) {
      sorted_spins.push_back(spins[i]);
      sorted_flags.push_back(flags[i]);
    }

    //now compute the Green's function
    nda::matrix<double> op = rho0;

    for (int i = 0; i < sorted_times.size(); i++) {
      if (sorted_flags[i] == 0) {
        op *= make_interaction_picture_create_op(ad, sorted_times[i], sorted_spins[i]);
      } else {
        op *= make_interaction_picture_destroy_op(ad, sorted_times[i], sorted_spins[i]);
      }
      //std::cout << "op after time " << sorted_times[i] << ": " << op << std::endl;
    }
    double G0_value = sign * trace(op);
    //std::cout << "G0 value: " << G0_value << std::endl;

    return G0_value;
  }

  double C02(triqs::atom_diag::atom_diag<false> ad, double beta, std::vector<double> times, std::vector<int> spins, std::vector<int> flags) {

    //Second order cumulant C(1,2|3,4) = G2(1,2|3,4) - G1(1|3) G1(2|4) + G1(1|4) G1(2|3)

    double G02 = G0(ad, beta, times, spins, flags); //G0(1,2|1',2')

    double G01_13 = G0(ad, beta, {times[0], times[2]}, {spins[0], spins[2]}, {flags[0], flags[2]}); //G(1|3)
    double G01_24 = G0(ad, beta, {times[1], times[3]}, {spins[1], spins[3]}, {flags[1], flags[3]}); //G(2|4)

    double G01_14 = G0(ad, beta, {times[0], times[3]}, {spins[0], spins[3]}, {flags[0], flags[3]}); //G(1|4)
    double G01_23 = G0(ad, beta, {times[1], times[2]}, {spins[1], spins[2]}, {flags[1], flags[2]}); //G(2|3)

    return G02 - G01_13 * G01_24 + G01_14 * G01_23; //return the second order cumulant
  }

} // namespace hubbard_atom