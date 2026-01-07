#include "./hubbard_atom.hpp"

namespace hubbard_atom {

  using cumul_args = std::vector<std::pair<double, int>>; //vector of pairs of imaginary time and spin

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

  double _partition_function(triqs::atom_diag::atom_diag<false> ad, double beta) {

    double Z0 = triqs::atom_diag::partition_function(ad, beta); //Z0 from atom_diag
    // double gs_energy = ad.get_gs_energy();
    // Z0 *= std::exp(-beta * gs_energy); //Z0 = Z0 * exp(-beta * E_0)
    return Z0;
  }

  int calculate_permutation_sign(const std::vector<int> &p) {
    int n = p.size();
    if (n == 0) return 1;

    std::vector<bool> visited(n, false);
    int num_cycles = 0;
    for (int i = 0; i < n; ++i) {
      if (!visited[i]) {
        num_cycles++;
        int j = i;
        while (!visited[j]) {
          visited[j] = true;
          j          = p[j];
        }
      }
    }
    // Sign of a permutation is (-1)^(n - number_of_cycles)
    return std::pow(-1, n - num_cycles);
  }

  std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, int>
  sort_operators(const std::vector<double> &times, const std::vector<int> &spins, const std::vector<int> &flags) {

    int n = times.size();
    std::vector<int> argsort(n);
    // 1. Create an index vector: 0, 1, 2, ...
    std::iota(argsort.begin(), argsort.end(), 0);

    // 2. Sort the index vector based on the corresponding time values (descending)
    std::stable_sort(argsort.begin(), argsort.end(), [&times](int i, int j) { return times[i] > times[j]; });

    // 3. Calculate the permutation sign
    int sign = calculate_permutation_sign(argsort);

    // 4. Create the sorted vectors using the argsort vector
    std::vector<double> sorted_times(n);
    std::vector<int> sorted_spins(n);
    std::vector<int> sorted_flags(n);

    for (int i = 0; i < n; ++i) {
      sorted_times[i] = times[argsort[i]];
      sorted_spins[i] = spins[argsort[i]];
      sorted_flags[i] = flags[argsort[i]];
    }

    return {sorted_times, sorted_spins, sorted_flags, sign};
  }

  nda::matrix<double> make_interaction_picture_destroy_op(triqs::atom_diag::atom_diag<false> ad, double tau, int state_index) {

    double Z01 = _partition_function(ad, tau);  //Z0
    double Z02 = _partition_function(ad, -tau); //Z0

    auto cmat        = ad.c_matrix(state_index, 0);
    auto time_evol_1 = triqs::atom_diag::atomic_density_matrix(ad, -tau)[0]; //time evolution operator
    auto time_evol_2 = triqs::atom_diag::atomic_density_matrix(ad, tau)[0];  //time evolution operator
    auto ctau        = Z01 * Z02 * time_evol_1 * cmat * time_evol_2;         //interaction picture destroy operator, watch -convention here

    return ctau; //return the interaction picture destroy operator
  }

  nda::matrix<double> make_interaction_picture_create_op(triqs::atom_diag::atom_diag<false> ad, double tau, int state_index) {

    double Z01 = _partition_function(ad, tau);  //Z0
    double Z02 = _partition_function(ad, -tau); //Z0

    auto cmat        = ad.cdag_matrix(state_index, 0);
    auto time_evol_1 = triqs::atom_diag::atomic_density_matrix(ad, -tau)[0]; //time evolution operator
    auto time_evol_2 = triqs::atom_diag::atomic_density_matrix(ad, tau)[0];  //time evolution operator
    auto cdagtau     = Z01 * Z02 * time_evol_1 * cmat * time_evol_2;         //interaction picture destroy operator

    return cdagtau; //return the interaction picture destroy operator
  }

  double G0(triqs::atom_diag::atom_diag<false> ad, double beta, cumul_args unprimed_args, cumul_args primed_args) {

    //unperturbed n body local Green's function: < Tc(tau_1,s_1) ... c(tau_n,s_n) \ cdag(tau_1', s_1')... cdag(tau_n', s_n')>_0

    int order = unprimed_args.size(); //order of the Green's function
    if (order != primed_args.size()) { throw std::runtime_error("Error in G0: unprimed_args and primed_args must have the same size"); }

    nda::matrix<double> rho0 = triqs::atom_diag::atomic_density_matrix(ad, beta)[0]; //rho_0 = exp(-beta H0)/Z0

    //create vectors for all the times, spins and flags
    int n_ops = 2 * order;
    std::vector<double> times;
    std::vector<int> spins;
    std::vector<int> flags;
    times.reserve(n_ops);
    spins.reserve(n_ops);
    flags.reserve(n_ops);

    for (auto [t, s] : unprimed_args) {
      times.push_back(t);
      spins.push_back(s);
      flags.push_back(0); //destroy
    }

    for (auto [t, s] : primed_args) {
      times.push_back(t);
      spins.push_back(s);
      flags.push_back(1); //create
    }

    //get the sorted times, spins and flags as well as the sign of the time-ordering permutation
    auto [sorted_times, sorted_spins, sorted_flags, sign] = sort_operators(times, spins, flags);

    //now compute the Green's function
    nda::matrix<double> op = rho0;

    for (int i = 0; i < n_ops; ++i) {
      if (sorted_flags[i] == 1) {
        op *= make_interaction_picture_create_op(ad, sorted_times[i], sorted_spins[i]);
      } else {
        op *= make_interaction_picture_destroy_op(ad, sorted_times[i], sorted_spins[i]);
      }
    }
    double G0_value = sign * trace(op);

    return G0_value;
  }

  double C02(triqs::atom_diag::atom_diag<false> ad, double beta, cumul_args unprimed_args, cumul_args primed_args) {

    double G02 = G0(ad, beta, unprimed_args, primed_args); //G0(1,2|1',2')

    double G0_11 = G0(ad, beta, {unprimed_args[0]}, {primed_args[0]}); //G(1|3)
    double G0_22 = G0(ad, beta, {unprimed_args[1]}, {primed_args[1]}); //G(2|4)

    double G0_12 = G0(ad, beta, {unprimed_args[0]}, {primed_args[1]}); //G(1|4)
    double G0_21 = G0(ad, beta, {unprimed_args[1]}, {primed_args[0]}); //G(2|3)

    return G02 - G0_11 * G0_22 + G0_12 * G0_21;
  } //return the second order cumulant

  //   //Second order cumulant C(1,2|3,4) = G2(1,2|3,4) - G1(1|3) G1(2|4) + G1(1|4) G1(2|3)

  //   double G02 = G0(ad, beta, times, spins, flags); //G0(1,2|1',2')

  //   double G01_13 = G0(ad, beta, {times[0], times[2]}, {spins[0], spins[2]}, {flags[0], flags[2]}); //G(1|3)
  //   double G01_24 = G0(ad, beta, {times[1], times[3]}, {spins[1], spins[3]}, {flags[1], flags[3]}); //G(2|4)

  //   double G01_14 = G0(ad, beta, {times[0], times[3]}, {spins[0], spins[3]}, {flags[0], flags[3]}); //G(1|4)
  //   double G01_23 = G0(ad, beta, {times[1], times[2]}, {spins[1], spins[2]}, {flags[1], flags[2]}); //G(2|3)

  //   return G02 - G01_13 * G01_24 + G01_14 * G01_23; //return the second order cumulant
  // }

} // namespace hubbard_atom