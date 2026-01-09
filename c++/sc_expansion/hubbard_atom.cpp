#include "./hubbard_atom.hpp"

namespace sc_expansion {

  using namespace triqs::operators;

  HubbardAtom::HubbardAtom(double U, double beta, double mu) : U(U), beta(beta), mu(mu) {
    auto fops = make_fops();
    H         = make_H0(U, mu);
    ad        = triqs::atom_diag::atom_diag<false>(H, fops, {});
  }

  triqs::hilbert_space::fundamental_operator_set HubbardAtom::make_fops() {
    triqs::hilbert_space::fundamental_operator_set fops = {};
    for (int o : triqs::arrays::range(1)) fops.insert("dn", o);
    for (int o : triqs::arrays::range(1)) fops.insert("up", o);
    return fops;
  }

  triqs::operators::many_body_operator_generic<double> HubbardAtom::make_H0(double U, double mu) {
    triqs::operators::many_body_operator_generic<double> h = U * n("up", 0) * n("dn", 0);
    h += -mu * (n("up", 0) + n("dn", 0));
    return h;
  }

  int HubbardAtom::calculate_permutation_sign(const std::vector<int> &p) {
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
    return std::pow(-1, n - num_cycles);
  }

  std::tuple<std::vector<double>, std::vector<int>, std::vector<int>, int>
  HubbardAtom::sort_operators(const std::vector<double> &times, const std::vector<int> &spins, const std::vector<int> &flags) {

    int n = times.size();
    std::vector<int> argsort(n);
    std::iota(argsort.begin(), argsort.end(), 0);

    std::stable_sort(argsort.begin(), argsort.end(), [&times](int i, int j) { return times[i] > times[j]; });

    int sign = calculate_permutation_sign(argsort);

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

  nda::matrix<double> HubbardAtom::make_interaction_picture_destroy_op(double tau, int state_index) const {

    double Z01 = triqs::atom_diag::partition_function(ad, tau);
    double Z02 = triqs::atom_diag::partition_function(ad, -tau);

    auto cmat        = ad.c_matrix(state_index, 0);
    auto time_evol_1 = triqs::atom_diag::atomic_density_matrix(ad, -tau)[0];
    auto time_evol_2 = triqs::atom_diag::atomic_density_matrix(ad, tau)[0];
    auto ctau        = Z01 * Z02 * time_evol_1 * cmat * time_evol_2;

    return ctau;
  }

  nda::matrix<double> HubbardAtom::make_interaction_picture_create_op(double tau, int state_index) const {

    double Z01 = triqs::atom_diag::partition_function(ad, tau);
    double Z02 = triqs::atom_diag::partition_function(ad, -tau);

    auto cmat        = ad.cdag_matrix(state_index, 0);
    auto time_evol_1 = triqs::atom_diag::atomic_density_matrix(ad, -tau)[0];
    auto time_evol_2 = triqs::atom_diag::atomic_density_matrix(ad, tau)[0];
    auto cdagtau     = Z01 * Z02 * time_evol_1 * cmat * time_evol_2;

    return cdagtau;
  }

  double HubbardAtom::G0(cumul_args const &unprimed_args, cumul_args const &primed_args) const {

    int order = unprimed_args.size();
    if (order != primed_args.size()) { throw std::runtime_error("Error in G0: unprimed_args and primed_args must have the same size"); }

    nda::matrix<double> rho0 = triqs::atom_diag::atomic_density_matrix(ad, beta)[0];

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
      flags.push_back(0); // destroy
    }

    for (auto [t, s] : primed_args) {
      times.push_back(t);
      spins.push_back(s);
      flags.push_back(1); // create
    }

    auto [sorted_times, sorted_spins, sorted_flags, sign] = sort_operators(times, spins, flags);

    nda::matrix<double> op = rho0;

    for (int i = 0; i < n_ops; ++i) {
      if (sorted_flags[i] == 1) {
        op *= make_interaction_picture_create_op(sorted_times[i], sorted_spins[i]);
      } else {
        op *= make_interaction_picture_destroy_op(sorted_times[i], sorted_spins[i]);
      }
    }
    double G0_value = sign * trace(op);

    return G0_value;
  }

} // namespace sc_expansion