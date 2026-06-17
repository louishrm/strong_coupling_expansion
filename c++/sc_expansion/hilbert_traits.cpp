#include "hilbert_traits.hpp"
#include "dual.hpp"
#include <cmath>

namespace sc_expansion {

  // ============================================================================
  // Atomic (N_sites = 1)
  // ============================================================================

  template <typename T> void HilbertTraits<1>::diagonalize(Parameters<T> const &params, std::array<Eigenstate<T>, 4> &out) {
    // Basis: |0>, |down>, |up>, |up down>
    // Bit order: bit 0 = down, bit 1 = up.
    out[0] = Eigenstate<T>{{{0, T(1.0)}}, T(0.0)};
    out[1] = Eigenstate<T>{{{1, T(1.0)}}, -params.mu};
    out[2] = Eigenstate<T>{{{2, T(1.0)}}, -params.mu};
    out[3] = Eigenstate<T>{{{3, T(1.0)}}, params.U - 2.0 * params.mu};
  }

  template <typename T> T HilbertTraits<1>::compute_Z_infinite_U(std::array<T, 4> const &exp_beta_E) {
    // Drop the doubly-occupied state (index 3).
    return exp_beta_E[0] + exp_beta_E[1] + exp_beta_E[2];
  }

  template <typename T> bool HilbertTraits<1>::try_closed_form_G0n(Args<1, T> const &args, Parameters<T> const &params, T const &Z, T &out) {
    // Closed-form one-body propagator (n = 2 only).
    // Operates on the (already time-sorted) args: taus[0] >= taus[1], delta >= 0.
    if (args.order != 2) return false;
    if (!args.operator_sequence_is_valid()) {
      out = T(0.0);
      return true;
    }

    using std::exp;
    T const &mu   = params.mu;
    T const &U    = params.U;
    T const &beta = params.beta;
    double delta  = args.taus[0] - args.taus[1];
    bool c_first  = (args.ops[0].get_action() == 0); // action 0 = annihilator

    T trace;
    if (c_first) {
      // Tr[e^{-bH} c(tau_1) c^dag(tau_1')], delta = tau_1 - tau_1' > 0
      trace = exp(mu * delta) + exp(beta * mu + (mu - U) * delta);
    } else {
      // Tr[e^{-bH} c^dag(tau_1') c(tau_1)], delta = tau_1' - tau_1 > 0
      trace = exp(beta * mu - mu * delta) + exp(-beta * (U - T(2.0) * mu) + (U - mu) * delta);
    }
    out = args.permutation_sign * trace / Z;
    return true;
  }

  // ============================================================================
  // Dimer (N_sites = 2)
  // ============================================================================

  namespace dimer_constants {
    constexpr double SQRT2_INV = 0.7071067811865475;

    template <typename T> inline T Eplus(T t, T U, T mu) {
      using std::sqrt;
      return T(0.5) * (U + sqrt(U * U + 16.0 * t * t)) - 2.0 * mu;
    }

    template <typename T> inline T Eminus(T t, T U, T mu) {
      using std::sqrt;
      return T(0.5) * (U - sqrt(U * U + 16.0 * t * t)) - 2.0 * mu;
    }
  } // namespace dimer_constants

  template <typename T> void HilbertTraits<2>::diagonalize(Parameters<T> const &params, std::array<Eigenstate<T>, 16> &out) {
    using namespace dimer_constants;
    using std::sqrt;

    T t  = params.t;
    T U  = params.U;
    T mu = params.mu;

    // N=0
    out[0] = Eigenstate<T>{{{0, T(1.0)}}, T(0.0)};

    // N=1, Sz=-1/2: (|down,0> ± |0,down>)/sqrt(2)
    out[1] = Eigenstate<T>{{{1, T(SQRT2_INV)}, {2, T(SQRT2_INV)}}, -t - mu};
    out[2] = Eigenstate<T>{{{1, T(SQRT2_INV)}, {2, -T(SQRT2_INV)}}, t - mu};

    // N=1, Sz=+1/2: (|up,0> ± |0,up>)/sqrt(2)
    out[3] = Eigenstate<T>{{{4, T(SQRT2_INV)}, {8, T(SQRT2_INV)}}, -t - mu};
    out[4] = Eigenstate<T>{{{4, T(SQRT2_INV)}, {8, -T(SQRT2_INV)}}, t - mu};

    // N=2, Sz=-1: |down,down>
    out[5] = Eigenstate<T>{{{3, T(1.0)}}, -2.0 * mu};

    // N=2, Sz=0, parity = even: 2x2 block H_even = [[U-2mu, -2t], [-2t, -2mu]].
    // Even-parity basis: |e+> = (|0101>+|1010>)/sqrt(2), |e-> = (|1001>+|0110>)/sqrt(2).
    // Eigenvector for eigenvalue E satisfies (U-2mu-E)*a = 2t*b, so b/a = (U-2mu-E)/(2t).
    T Ep = Eplus(t, U, mu);
    T Em = Eminus(t, U, mu);

    // At t=0 the block diagonalises and the ratio b/a = (U-2mu-E)/(2t) diverges,
    // so we handle the limit explicitly.
    T norm_plus, norm_minus, component_plus, component_minus;
    if (is_zero(t)) {
      // E+ = U-2mu has eigenvector |e+>, E- = -2mu has eigenvector |e->.
      norm_plus       = T(SQRT2_INV);
      component_plus  = T(0.0);
      norm_minus      = T(0.0);
      component_minus = T(SQRT2_INV);
    } else {
      T ratio_plus    = (U - T(2.0) * mu - Ep) / (T(2.0) * t);
      T ratio_minus   = (U - T(2.0) * mu - Em) / (T(2.0) * t);
      norm_plus       = T(SQRT2_INV) / sqrt(T(1.0) + ratio_plus * ratio_plus);
      norm_minus      = T(SQRT2_INV) / sqrt(T(1.0) + ratio_minus * ratio_minus);
      component_plus  = ratio_plus * norm_plus;
      component_minus = ratio_minus * norm_minus;
    }

    out[6] = Eigenstate<T>{{{5, norm_plus}, {10, norm_plus}, {9, component_plus}, {6, component_plus}}, Ep};
    out[7] = Eigenstate<T>{{{5, norm_minus}, {10, norm_minus}, {9, component_minus}, {6, component_minus}}, Em};

    // N=2, Sz=0, parity = odd
    out[8] = Eigenstate<T>{{{5, T(SQRT2_INV)}, {10, -T(SQRT2_INV)}}, U - 2.0 * mu};
    out[9] = Eigenstate<T>{{{9, T(SQRT2_INV)}, {6, -T(SQRT2_INV)}}, -2.0 * mu};

    // N=2, Sz=+1: |up,up>
    out[10] = Eigenstate<T>{{{12, T(1.0)}}, -2.0 * mu};

    // N=3, Sz=-1/2: (|down up, down> ± |down, down up>)/sqrt(2)
    out[11] = Eigenstate<T>{{{7, T(SQRT2_INV)}, {11, T(SQRT2_INV)}}, U - t - 3.0 * mu};
    out[12] = Eigenstate<T>{{{7, T(SQRT2_INV)}, {11, -T(SQRT2_INV)}}, U + t - 3.0 * mu};

    // N=3, Sz=+1/2: (|down up, up> ± |up, down up>)/sqrt(2)
    out[13] = Eigenstate<T>{{{13, T(SQRT2_INV)}, {14, T(SQRT2_INV)}}, U - t - 3.0 * mu};
    out[14] = Eigenstate<T>{{{13, T(SQRT2_INV)}, {14, -T(SQRT2_INV)}}, U + t - 3.0 * mu};

    // N=4, Sz=0: |down up, down up>
    out[15] = Eigenstate<T>{{{15, T(1.0)}}, 2.0 * U - 4.0 * mu};
  }

  template <typename T> T HilbertTraits<2>::compute_Z_infinite_U(std::array<T, 16> const &exp_beta_E) {
    // Dimer infinite-U is not used by the dimer MCMC. Return Z (full sum) so the
    // value is well-defined; dimer::Move must not invoke G0n_infinite_U regardless.
    T sum = exp_beta_E[0];
    for (int i = 1; i < 16; ++i) sum = sum + exp_beta_E[i];
    return sum;
  }

  template <typename T> bool HilbertTraits<2>::try_closed_form_G0n(Args<2, T> const & /*args*/, Parameters<T> const & /*params*/,
                                                                   T const & /*Z*/, T & /*out*/) {
    return false;
  }

  // ----- Explicit instantiations -----

  template void HilbertTraits<1>::diagonalize(Parameters<double> const &, std::array<Eigenstate<double>, 4> &);
  template void HilbertTraits<1>::diagonalize(Parameters<Dual> const &, std::array<Eigenstate<Dual>, 4> &);
  template double HilbertTraits<1>::compute_Z_infinite_U(std::array<double, 4> const &);
  template Dual HilbertTraits<1>::compute_Z_infinite_U(std::array<Dual, 4> const &);
  template bool HilbertTraits<1>::try_closed_form_G0n(Args<1, double> const &, Parameters<double> const &, double const &, double &);
  template bool HilbertTraits<1>::try_closed_form_G0n(Args<1, Dual> const &, Parameters<Dual> const &, Dual const &, Dual &);

  template void HilbertTraits<2>::diagonalize(Parameters<double> const &, std::array<Eigenstate<double>, 16> &);
  template void HilbertTraits<2>::diagonalize(Parameters<Dual> const &, std::array<Eigenstate<Dual>, 16> &);
  template double HilbertTraits<2>::compute_Z_infinite_U(std::array<double, 16> const &);
  template Dual HilbertTraits<2>::compute_Z_infinite_U(std::array<Dual, 16> const &);
  template bool HilbertTraits<2>::try_closed_form_G0n(Args<2, double> const &, Parameters<double> const &, double const &, double &);
  template bool HilbertTraits<2>::try_closed_form_G0n(Args<2, Dual> const &, Parameters<Dual> const &, Dual const &, Dual &);

} // namespace sc_expansion
