/*
 * Deterministic quadrature test for dimer expansion on the INFINITE lattice.
 *
 * Computes free-energy coefficients a_n by trapezoidal quadrature (no MC),
 * exploiting time-translation invariance: fix tau_n = 0, integrate tau_1..tau_{n-1},
 * multiply by beta.
 *
 * The cumulant has a kink at equal times, so trapezoidal convergence is O(h^2).
 * Grid sizes are chosen to give ~0.1% accuracy at each order.
 *
 * This gives reference values for the infinite lattice that can be compared
 * against MCMC (mcmc_dimer) to diagnose sign-problem bias vs computation bugs.
 *
 * Also computes the mu-derivative (density coefficient) via finite differences
 * as an independent cross-check of the Dual-number MCMC results.
 *
 * Usage:  mpirun -np 1 ./test_dimer_quadrature
 */

#include <gtest/gtest.h>
#include "sc_expansion/free_energy_calculator.hpp"
#include <mpi/mpi.hpp>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace sc_expansion;

// =====================================================================
// Helper: compute a_n by quadrature for dimer expansion on infinite lattice.
//
// Uses time-translation invariance: fix tau_n = 0, integrate tau_1..tau_{n-1}
// over [0, beta]^{n-1} with trapezoidal rule, multiply by beta.
// =====================================================================
static double dimer_coefficient_quadrature(double U, double beta, double mu, double t_intra, int order, int N_quad) {

  Parameters<double> params{U, beta, mu, t_intra, false};
  FreeEnergyCalculator<2, double> calculator(params, order);

  int n_free = order - 1; // number of free taus (one fixed at 0)
  double h   = beta / N_quad;

  // Total number of grid points = (N_quad + 1)^{n_free}
  long total_points = 1;
  for (int d = 0; d < n_free; d++) total_points *= (N_quad + 1);

  double integral = 0.0;

  for (long flat_idx = 0; flat_idx < total_points; flat_idx++) {
    std::vector<double> taus(order, 0.0); // last tau is fixed at 0
    double weight = 1.0;

    long remainder = flat_idx;
    for (int d = 0; d < n_free; d++) {
      int idx = (int)(remainder % (N_quad + 1));
      remainder /= (N_quad + 1);
      taus[d] = idx * h;
      double w = (idx == 0 || idx == N_quad) ? 0.5 : 1.0;
      weight *= w;
    }

    calculator.mark_all_dirty();
    double val = calculator.compute_sum_diagrams_dimer(taus, false);
    integral += weight * val;
  }

  integral *= std::pow(h, n_free) * beta;
  calculator.clear_all_caches();

  return integral;
}

// =====================================================================
// Order 2: cross-check against existing test
// N=5000 matches the proven 1D quadrature in test_mcmc_dimer.cpp
// =====================================================================
TEST(DimerQuadrature, Order2InfiniteLattice) {
  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order  = 2;
  int N_quad = 5000; // 1D integral, matches existing test precision

  double coeff         = dimer_coefficient_quadrature(U, beta, mu, t_intra, order, N_quad);
  double cluster_coeff = coeff / 6.0;

  std::cout << std::setprecision(12);
  std::cout << "\n=== Order 2 infinite lattice ===" << std::endl;
  std::cout << "Quadrature coeff (infinite): " << coeff << std::endl;
  std::cout << "Cluster coeff (/ fm=6):      " << cluster_coeff << std::endl;
  std::cout << "Exact cluster (Python ED):   " << -0.066467819521 << std::endl;

  EXPECT_NEAR(cluster_coeff, -0.066467819521, 1e-4);
}

// =====================================================================
// Order 3: first odd order — only the 3-cycle contributes.
// 2D integral: N=300 gives 301^2 = ~90K evaluations.
// =====================================================================
TEST(DimerQuadrature, Order3InfiniteLattice) {
  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order  = 3;
  int N_quad = 300;

  double coeff = dimer_coefficient_quadrature(U, beta, mu, t_intra, order, N_quad);

  double dmu         = 1e-4;
  double coeff_plus  = dimer_coefficient_quadrature(U, beta, mu + dmu, t_intra, order, N_quad);
  double coeff_minus = dimer_coefficient_quadrature(U, beta, mu - dmu, t_intra, order, N_quad);
  double dcoeff_dmu  = (coeff_plus - coeff_minus) / (2.0 * dmu);

  std::cout << std::setprecision(12);
  std::cout << "\n=== Order 3 infinite lattice ===" << std::endl;
  std::cout << "Quadrature coeff a_3:         " << coeff << std::endl;
  std::cout << "d(a_3)/d(mu) (finite diff):   " << dcoeff_dmu << std::endl;
  std::cout << "Density correction per site:  " << -dcoeff_dmu / 2.0 << std::endl;
  std::cout << "(Compare da_3/dmu against MCMC 'mean' in density HDF5)" << std::endl;
}

// =====================================================================
// Order 4: 3D integral
// N=50 gives 51^3 = ~133K evaluations.
// =====================================================================
TEST(DimerQuadrature, Order4InfiniteLattice) {
  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order  = 4;
  int N_quad = 50;

  double coeff = dimer_coefficient_quadrature(U, beta, mu, t_intra, order, N_quad);

  double dmu         = 1e-4;
  double coeff_plus  = dimer_coefficient_quadrature(U, beta, mu + dmu, t_intra, order, N_quad);
  double coeff_minus = dimer_coefficient_quadrature(U, beta, mu - dmu, t_intra, order, N_quad);
  double dcoeff_dmu  = (coeff_plus - coeff_minus) / (2.0 * dmu);

  std::cout << std::setprecision(12);
  std::cout << "\n=== Order 4 infinite lattice ===" << std::endl;
  std::cout << "Quadrature coeff a_4:         " << coeff << std::endl;
  std::cout << "d(a_4)/d(mu) (finite diff):   " << dcoeff_dmu << std::endl;
  std::cout << "Density correction per site:  " << -dcoeff_dmu / 2.0 << std::endl;
}

// =====================================================================
// Order 5: 4D integral
// N=15 gives 16^4 = ~65K evaluations.
// Accuracy ~1-2% (coarse grid), but sufficient to detect gross MCMC errors.
// =====================================================================
TEST(DimerQuadrature, Order5InfiniteLattice) {
  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order  = 5;
  int N_quad = 15;

  double coeff = dimer_coefficient_quadrature(U, beta, mu, t_intra, order, N_quad);

  double dmu         = 1e-4;
  double coeff_plus  = dimer_coefficient_quadrature(U, beta, mu + dmu, t_intra, order, N_quad);
  double coeff_minus = dimer_coefficient_quadrature(U, beta, mu - dmu, t_intra, order, N_quad);
  double dcoeff_dmu  = (coeff_plus - coeff_minus) / (2.0 * dmu);

  std::cout << std::setprecision(12);
  std::cout << "\n=== Order 5 infinite lattice ===" << std::endl;
  std::cout << "Quadrature coeff a_5:         " << coeff << std::endl;
  std::cout << "d(a_5)/d(mu) (finite diff):   " << dcoeff_dmu << std::endl;
  std::cout << "Density correction per site:  " << -dcoeff_dmu / 2.0 << std::endl;
  std::cout << "(~1-2%% accuracy from coarse 4D grid)" << std::endl;
}

// =====================================================================
// Mu scan at order 3: compute density correction at multiple mu values
// to compare directly against MCMC HDF5 data.
// =====================================================================
TEST(DimerQuadrature, Order3MuScan) {
  double U = 8.0, beta = 2.0, t_intra = 1.0;
  int order  = 3;
  int N_quad = 200;

  std::vector<double> mus = {2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0};
  double dmu = 1e-4;

  std::cout << std::setprecision(12);
  std::cout << "\n=== Order 3 mu scan (infinite lattice, N=" << N_quad << ") ===" << std::endl;
  std::cout << "  mu         | a_3             | da_3/dmu        | density corr /site" << std::endl;
  std::cout << "  -----------+-----------------+-----------------+-------------------" << std::endl;

  for (double mu : mus) {
    double coeff       = dimer_coefficient_quadrature(U, beta, mu, t_intra, order, N_quad);
    double coeff_plus  = dimer_coefficient_quadrature(U, beta, mu + dmu, t_intra, order, N_quad);
    double coeff_minus = dimer_coefficient_quadrature(U, beta, mu - dmu, t_intra, order, N_quad);
    double dcoeff_dmu  = (coeff_plus - coeff_minus) / (2.0 * dmu);

    std::cout << "  " << std::setw(9) << mu << "  |  " << std::setw(14) << coeff << "  |  " << std::setw(14) << dcoeff_dmu << "  |  "
              << std::setw(14) << -dcoeff_dmu / 2.0 << std::endl;
  }
}

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
