#pragma once
#include <array>
#include "fock_space.hpp"
#include "args.hpp"

namespace sc_expansion {

  template <typename T> struct Parameters {
    T U;
    T beta;
    T mu;
    T t            = T(0.0); // hopping inside the cluster (used by N_sites>=2 ED only)
    bool bipartite = true;
    // Shift absorbed into the perturbation: H_pert = -t sum c^dag c + delta sum n.
    // HubbardSolver still diagonalises H_0 with the `mu` stored here; the caller
    // is responsible for passing the already-shifted chemical potential (mu - delta)
    // if they want the atomic reference to match the variational choice.
    T delta = T(0.0);
  };

  // Per-cluster traits: Hilbert-space size, ED, infinite-U projection rule, and
  // an optional closed-form fast path for low-order G0n. One specialization per
  // local cluster (N_sites=1 atomic, N_sites=2 dimer, ...). HubbardSolver delegates
  // here so its body stays cluster-agnostic.
  template <int N_sites> struct HilbertTraits;

  // ----- Atomic (single site) -----
  template <> struct HilbertTraits<1> {
    static constexpr int N_SITES  = 1;
    static constexpr int N_STATES = 4;
    static constexpr int N_OPS    = 4;

    template <typename T> static void diagonalize(Parameters<T> const &params, std::array<Eigenstate<T>, N_STATES> &out);

    template <typename T> static T compute_Z_infinite_U(std::array<T, N_STATES> const &exp_beta_E);

    // Returns true and writes `out` if a closed-form G0n is available for this `args`.
    // Single-site provides one for n=2 (the analytic G01); higher orders fall through.
    template <typename T> static bool try_closed_form_G0n(Args<1, T> const &args, Parameters<T> const &params, T const &Z, T &out);
  };

  // ----- Dimer (two sites) -----
  template <> struct HilbertTraits<2> {
    static constexpr int N_SITES  = 2;
    static constexpr int N_STATES = 16;
    static constexpr int N_OPS    = 8;

    template <typename T> static void diagonalize(Parameters<T> const &params, std::array<Eigenstate<T>, N_STATES> &out);

    // Dimer has no infinite-U path; HubbardSolver still computes a value for symmetry
    // with the atomic interface, but dimer MCMC must never invoke G0n_infinite_U.
    template <typename T> static T compute_Z_infinite_U(std::array<T, N_STATES> const &exp_beta_E);

    // No closed form available.
    template <typename T> static bool try_closed_form_G0n(Args<2, T> const &args, Parameters<T> const &params, T const &Z, T &out);
  };

} // namespace sc_expansion
