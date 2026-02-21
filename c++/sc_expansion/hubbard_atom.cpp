#include "./hubbard_atom.hpp"

namespace sc_expansion {

  using namespace triqs::operators;

  HubbardAtom::HubbardAtom(double U, double beta, double mu) : U(U), beta(beta), mu(mu) {
    this->Z_infinite_U = 1 + 2 * std::exp(beta * mu); //Partition function at infinite U
    auto fops          = make_fops();
    this->H            = make_H0(U, mu);
    this->ad           = triqs::atom_diag::atom_diag<false>(this->H, fops, {});
    this->rho0         = triqs::atom_diag::atomic_density_matrix(this->ad, beta)[0];
    this->c_dn         = this->ad.c_matrix(0, 0);
    this->c_up         = this->ad.c_matrix(1, 0);
    this->cdag_dn      = this->ad.cdag_matrix(0, 0);
    this->cdag_up      = this->ad.cdag_matrix(1, 0);
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

  int HubbardAtom::calculate_permutation_sign(std::vector<int> const &perm) {

    int inversions = 0;
    for (size_t i = 0; i < perm.size(); ++i) {
      for (size_t j = i + 1; j < perm.size(); ++j) {
        if (perm[i] > perm[j]) { ++inversions; }
      }
    }
    return (inversions % 2 == 0) ? 1 : -1;
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

  bool HubbardAtom::verify_consecutive_terms(const std::vector<int> &sorted_spins, const std::vector<int> &sorted_flags) {
    // Vectors to store flags for each spin sector
    std::vector<int> flags_up;
    std::vector<int> flags_down;

    // 1. Separate the stream by spin
    for (size_t i = 0; i < sorted_spins.size(); ++i) {
      if (sorted_spins[i] == 1) { // Assuming 1 is UP
        flags_up.push_back(sorted_flags[i]);
      } else { // Assuming 0 or -1 is DOWN
        flags_down.push_back(sorted_flags[i]);
      }
    }

    // 2. Helper lambda to check for strict alternation in a single sector
    auto is_invalid_sector = [](const std::vector<int> &flags) {
      for (size_t i = 0; i + 1 < flags.size(); ++i) {
        // If two consecutive flags are identical (1,1 or 0,0), it's invalid
        if (flags[i] == flags[i + 1]) { return true; }
      }
      return false;
    };

    // 3. Verify both sectors
    if (is_invalid_sector(flags_up)) return false;
    if (is_invalid_sector(flags_down)) return false;

    return true;
  }

  nda::matrix<double> HubbardAtom::make_interaction_picture_destroy_op(double tau, int state_index) const {
    double Z01 = triqs::atom_diag::partition_function(this->ad, tau);
    double Z02 = triqs::atom_diag::partition_function(this->ad, -tau);

    auto const &cmat = (state_index == 0 ? this->c_dn : this->c_up);
    auto time_evol_1 = triqs::atom_diag::atomic_density_matrix(this->ad, -tau)[0];
    auto time_evol_2 = triqs::atom_diag::atomic_density_matrix(this->ad, tau)[0];
    auto ctau        = Z01 * Z02 * time_evol_1 * cmat * time_evol_2;

    return ctau;
  }

  nda::matrix<double> HubbardAtom::make_interaction_picture_create_op(double tau, int state_index) const {
    double Z01 = triqs::atom_diag::partition_function(this->ad, tau);
    double Z02 = triqs::atom_diag::partition_function(this->ad, -tau);

    auto const &cmat = (state_index == 0 ? this->cdag_dn : this->cdag_up);
    auto time_evol_1 = triqs::atom_diag::atomic_density_matrix(this->ad, -tau)[0];
    auto time_evol_2 = triqs::atom_diag::atomic_density_matrix(this->ad, tau)[0];
    auto cdagtau     = Z01 * Z02 * time_evol_1 * cmat * time_evol_2;

    return cdagtau;
  }

  double HubbardAtom::G0(cumul_args const &unprimed_args, cumul_args const &primed_args) const {
    int order = unprimed_args.size();
    if (order != primed_args.size()) { throw std::runtime_error("Error in G0: unprimed_args and primed_args must have the same size"); }

    int n_ops = 2 * order;
    std::vector<double> times;
    std::vector<int> spins;
    std::vector<int> flags;
    times.reserve(n_ops);
    spins.reserve(n_ops);
    flags.reserve(n_ops);

    for (int i = 0; i < order; ++i) {
      // Push the Creation operator first
      times.push_back(primed_args[i].first);
      spins.push_back(primed_args[i].second);
      flags.push_back(1);

      // Push the Destruction operator second
      times.push_back(unprimed_args[i].first);
      spins.push_back(unprimed_args[i].second);
      flags.push_back(0);
    }

    auto [sorted_times, sorted_spins, sorted_flags, sign] = sort_operators(times, spins, flags);
    if (!verify_consecutive_terms(sorted_spins, sorted_flags)) { return 0.0; }

    nda::matrix<double> op = this->rho0;

    for (int i = 0; i < n_ops; ++i) {
      if (sorted_flags[i] == 1) {
        op *= this->make_interaction_picture_create_op(sorted_times[i], sorted_spins[i]);
      } else {
        op *= this->make_interaction_picture_destroy_op(sorted_times[i], sorted_spins[i]);
      }
    }
    double G0_value = sign * trace(op);

    return G0_value;
  }

  bool HubbardAtom::verify_consecutive_terms_infinite_U(const std::vector<int> &sorted_spins, const std::vector<int> &sorted_flags, int n_ops) {
    int total_ops = 2 * n_ops;
    for (int i = 0; i < total_ops - 1; ++i) {
      if (sorted_flags[i] == sorted_flags[i + 1]) { return false; } //Two consecutive creates or destroys (violates alternation/no double occupancy)

      if (sorted_flags[i] == 0 && sorted_flags[i + 1] == 1) {
        if (sorted_spins[i] != sorted_spins[i + 1]) return false;
      } //Destroying a particle with different spin than the one created at the previous (later) time
    }

    // Check boundary condition for trace: first and last operators must have matching spins
    if (sorted_spins[0] != sorted_spins[total_ops - 1]) return false;

    return true;
  }

  double HubbardAtom::G0_infinite_U(cumul_args const &unprimed_args, cumul_args const &primed_args) const {

    int order = unprimed_args.size();
    if (order != primed_args.size()) { throw std::runtime_error("Error in G0: unprimed_args and primed_args must have the same size"); }

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

    if (!verify_consecutive_terms_infinite_U(sorted_spins, sorted_flags, order)) { return 0.0; }

    if (sorted_flags[0] == 1) {
      return sign * std::exp(this->beta * this->mu) / this->Z_infinite_U;
    } //First operator is create
    else {
      return sign * 1.0 / this->Z_infinite_U;
    } //First operator is destroy
  }
} // namespace sc_expansion