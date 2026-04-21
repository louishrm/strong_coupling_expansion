/*
 * MCMC test for AtomicConfiguration (single-site series) benchmarked against
 * the dimer free-energy coefficients from analytical/benchmark_atomic_expansion.py.
 *
 * All graphs are forced to free_multiplicity=1 via override_fm=1 so the MC
 * estimate is compared directly against the cluster coefficient (no lattice
 * embedding weight). Covers both the pure-hopping case (delta=0) and the
 * shifted expansion (mu -> mu - delta, with delta reinserted as self-loop
 * density vertices).
 *
 * Usage:  mpirun -np 4 ./test_mcmc_atom
 *         (also works with 1 rank: ./test_mcmc_atom)
 */

#include <gtest/gtest.h>
#include "sc_expansion/atomic_configuration.hpp"
#include "sc_expansion/free_energy_calculator.hpp"
#include "sc_expansion/generate_diagrams.hpp"
#include "sc_expansion/graph.hpp"
#include "sc_expansion/hubbard_solver.hpp"
#include "sc_expansion/args.hpp"
#include "sc_expansion/fock_space.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/measure.hpp"
#include "sc_expansion/dual.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <mpi/mpi.hpp>
#include <cmath>
#include <memory>

static void run_mcmc_atom_check(int order, double U, double beta, double mu, double delta, double exact_coeff, double rel_tol, int n_cycles) {

  mpi::communicator world;

  double alpha     = 0.001;
  int n_warmup     = 2000;
  int length_cycle = 1;

  bool allow_self_loops = (delta != 0.0);

  // Python ED uses V = Ht_unit + delta*Hmu = T - delta*N (since Hmu = -N).
  // C++ perturbation convention is T + params.delta*N, so set params.delta = -delta.
  // H_0 uses the shifted chemical potential mu - delta.
  sc_expansion::Parameters<double> params{U, beta, mu - delta, true, -delta};

  sc_expansion::FreeEnergyCalculator<double> calculator(params, order, /*override_fm=*/1, allow_self_loops);

  double reference_integral        = 0.0;
  double signed_reference_integral = 0.0;

  if (world.rank() == 0) {
    auto [ref_abs, ref_signed] = calculator.compute_infinite_U_coefficient();
    reference_integral         = ref_abs;
    signed_reference_integral  = ref_signed;
  }

  mpi::broadcast(reference_integral, world);
  mpi::broadcast(signed_reference_integral, world);

  // =====================================================================
  // Construct AtomicConfiguration reusing the same calculator
  // =====================================================================
  auto config = std::make_unique<AtomicConfiguration<double>>(params, order, alpha, calculator);

  // =====================================================================
  // Set up TRIQS MC loop (same pattern as apps/mcmc_atom.cpp)
  // =====================================================================
  int random_seed = 32186222 + world.rank() * 786512;
  int verbosity   = (world.rank() == 0 ? 2 : 0);

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  measure<double> meas(config.get(), reference_integral, signed_reference_integral, n_bins, block_size, mu - delta);
  mc.add_move(move<double>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "defensive_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  // =====================================================================
  // Read results from shared result struct and verify against exact coefficient
  // =====================================================================
  if (world.rank() == 0) {
    double mc_mean  = meas.result->mean;
    double mc_error = meas.result->error;

    double rel_err = std::abs(mc_mean - exact_coeff) / std::abs(exact_coeff);

    std::cout << "Exact (Python ED):  " << exact_coeff << std::endl;
    std::cout << "MC estimate:        " << mc_mean << std::endl;
    std::cout << "MC error:           " << mc_error << std::endl;
    std::cout << "Relative error:     " << rel_err << std::endl;

    EXPECT_LT(rel_err, rel_tol) << "MC estimate " << mc_mean << " deviates from exact " << exact_coeff << " by " << rel_err * 100 << "%";

    // Print vertex cache stats
    std::cout << "\n=== Vertex Cache Stats ===" << std::endl;
    auto const &vt = config->get_calculator().get_vertex_types();
    for (size_t i = 0; i < vt.size(); i++) {
      auto [lhits, lmisses] = vt[i].get_local_cache_stats();
      long ltotal           = lhits + lmisses;
      double l_hit_rate     = ltotal > 0 ? 100.0 * lhits / ltotal : 0.0;
      std::cout << "  Cumulant order " << (i + 1) << ": local hits=" << lhits << " misses=" << lmisses << " hit_rate=" << l_hit_rate << "%"
                << std::endl;
    }
  }
}

// Deterministic infinite-U reference-integral check (no MC) for the dimer.
// Uses build_dimer_corrected_graphs so graphs with self-loops get the right
// dimer multiplicity. At delta=0 this collapses to the pure-hopping benchmark.
static void run_infinite_u_check(int order, double U, double beta, double mu, double delta, double expected_signed_coeff, double tol) {
  bool allow_self_loops = (delta != 0.0);
  sc_expansion::Parameters<double> params{U, beta, mu - delta, true, -delta};

  sc_expansion::FreeEnergyCalculator<double> calculator(params, order, /*override_fm=*/1, allow_self_loops);
  auto [abs_coeff, signed_coeff] = calculator.compute_infinite_U_coefficient();

  std::cout << "Expected (Python): " << expected_signed_coeff << std::endl;
  std::cout << "C++  signed_coeff: " << signed_coeff << std::endl;
  std::cout << "abs_coeff:         " << abs_coeff << std::endl;

  EXPECT_NEAR(signed_coeff, expected_signed_coeff, tol);
}

TEST(McmcAtom, InfiniteUDimerOrder2) {
  // Pure-hopping benchmark moved from test_diagram.cpp.
  run_infinite_u_check(/*order=*/2, /*U=*/8.0, /*beta=*/2.0, /*mu=*/3.0, /*delta=*/0.0,
                       /*expected_signed_coeff=*/-0.0012363096839754632, /*tol=*/1e-12);
}
TEST(McmcAtom, InfiniteUDimerOrder4) {
  // Pure-hopping benchmark moved from test_diagram.cpp.
  run_infinite_u_check(/*order=*/4, /*U=*/8.0, /*beta=*/2.0, /*mu=*/3.0, /*delta=*/0.0,
                       /*expected_signed_coeff=*/-4.0904630472238777e-04, /*tol=*/1e-12);
}

TEST(McmcAtom, InfiniteUDimerShiftedOrder2) {
  run_infinite_u_check(/*order=*/2, /*U=*/8.0, /*beta=*/2.0, /*mu=*/3.0, /*delta=*/1.1,
                       /*expected_signed_coeff=*/-0.024175845859726222, /*tol=*/1e-10);
}
TEST(McmcAtom, InfiniteUDimerShiftedOrder3) {
  run_infinite_u_check(/*order=*/3, /*U=*/8.0, /*beta=*/2.0, /*mu=*/3.0, /*delta=*/1.1,
                       /*expected_signed_coeff=*/0.03302607921926233, /*tol=*/1e-10);
}

TEST(McmcAtom, Order4Coefficient) {
  // Pure-t expansion (delta=0), matches analytical/benchmark_atomic_expansion.py.
  run_mcmc_atom_check(/*order=*/4, /*U=*/8.0, /*beta=*/2.0, /*mu=*/3.0, /*delta=*/0.0,
                      /*exact_coeff=*/-0.019604213442906773, /*rel_tol=*/0.10, /*n_cycles=*/100000);
}

TEST(McmcAtom, ShiftedOrder2) {
  // Shifted expansion: mu -> mu - delta, with delta reinserted as self-loop vertices.
  run_mcmc_atom_check(/*order=*/2, /*U=*/8.0, /*beta=*/2.0, /*mu=*/3.0, /*delta=*/1.1,
                      /*exact_coeff=*/-0.08502669580461163, /*rel_tol=*/0.05, /*n_cycles=*/100000);
}

TEST(McmcAtom, ShiftedOrder3) {
  run_mcmc_atom_check(/*order=*/3, /*U=*/8.0, /*beta=*/2.0, /*mu=*/3.0, /*delta=*/1.1,
                      /*exact_coeff=*/0.029672838206571895, /*rel_tol=*/0.20, /*n_cycles=*/200000);
}

// Runs the defensive MCMC estimator with T=Dual so AtomicConfiguration extracts
// the .derivative of each Dual result, turning the free-energy MC into a density MC.
// compute_infinite_U_coefficient() with T=Dual likewise returns the density reference.
static void run_mcmc_density_check(int order, double U, double beta, double mu, double delta, double exact_density_coeff, double rel_tol,
                                   int n_cycles) {
  mpi::communicator world;

  double alpha     = 0.001;
  int n_warmup     = 2000;
  int length_cycle = 1;

  bool allow_self_loops = (delta != 0.0);

  sc_expansion::Parameters<Dual> params{Dual(U, 0.0), Dual(beta, 0.0), Dual(mu - delta, 1.0), true, Dual(-delta, 0.0)};

  sc_expansion::FreeEnergyCalculator<Dual> calculator(params, order, /*override_fm=*/1, allow_self_loops);

  double reference_integral        = 0.0;
  double signed_reference_integral = 0.0;

  if (world.rank() == 0) {
    auto [ref_abs, ref_signed] = calculator.compute_infinite_U_coefficient();
    reference_integral         = ref_abs;
    signed_reference_integral  = ref_signed;
  }

  mpi::broadcast(reference_integral, world);
  mpi::broadcast(signed_reference_integral, world);

  auto config = std::make_unique<AtomicConfiguration<Dual>>(params, order, alpha, calculator);

  int random_seed = 32186222 + world.rank() * 786512;
  int verbosity   = (world.rank() == 0 ? 2 : 0);

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  measure<Dual> meas(config.get(), reference_integral, signed_reference_integral, n_bins, block_size, mu - delta);
  mc.add_move(move<Dual>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "defensive_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  if (world.rank() == 0) {
    double mc_mean  = meas.result->mean;
    double mc_error = meas.result->error;
    double rel_err  = std::abs(mc_mean - exact_density_coeff) / std::abs(exact_density_coeff);

    std::cout << "Exact (ED):     " << exact_density_coeff << std::endl;
    std::cout << "MC estimate:    " << mc_mean << std::endl;
    std::cout << "MC error:       " << mc_error << std::endl;
    std::cout << "Relative error: " << rel_err << std::endl;

    EXPECT_LT(rel_err, rel_tol) << "MC density estimate " << mc_mean << " deviates from exact " << exact_density_coeff << " by " << rel_err * 100
                                << "%";
  }
}

TEST(McmcAtom, DensityOrder2MCMCDual) {
  run_mcmc_density_check(/*order=*/2, /*U=*/8.0, /*beta=*/2.0, /*mu=*/3.0, /*delta=*/0.0,
                         /*exact_density_coeff=*/0.0021180825055121203, /*rel_tol=*/0.10, /*n_cycles=*/100000);
}

TEST(McmcAtom, DensityOrder2MCMCDualShfited) {
  run_mcmc_density_check(/*order=*/2, /*U=*/8.0, /*beta=*/2.0, /*mu=*/3.0, /*delta=*/1.1,
                         /*exact_density_coeff=*/0.04456693714966627, /*rel_tol=*/0.10, /*n_cycles=*/100000);
}

// Exact first-order density correction from the Python benchmark:
//   n_1 = delta * d(n_0)/d(mu_eff)   where  mu_eff = mu - delta
// This is the leading (order-1 in delta) change in the site density.
static double exact_n1_formula(double U, double beta, double mu, double delta) {
  double mueff      = mu - delta;
  double Z          = 1.0 + 2 * std::exp(beta * mueff) + std::exp(-beta * (U - 2 * mueff));
  double dZ_dmueff  = 2 * beta * std::exp(beta * mueff) + 2 * beta * std::exp(-beta * (U - 2 * mueff));
  double n0         = 2.0 / Z * (std::exp(beta * mueff) + std::exp(-beta * (U - 2 * mueff)));
  double dn0_dmueff = -dZ_dmueff / Z * n0 + 2 * beta / Z * (std::exp(beta * mueff) + 2 * std::exp(-beta * (U - 2 * mueff)));
  return delta * dn0_dmueff;
}

// Deterministic check: computes the first-order density correction n_1 = -da_1/dmu
// using FreeEnergyCalculator<Dual> with mu differentiation (derivative seed = 1).
// At order 1 with allow_self_loops, the only diagram is a single self-loop.  Its
// integrand is constant in tau, so evaluating compute_sum_diagrams at tau=0 and
// multiplying by beta recovers the coefficient exactly (no MC noise).
static void run_density_dual_check(double U, double beta, double mu, double delta, double tol) {
  // C++ convention: params.mu = mu - delta (shifted), params.delta = -delta.
  // Set derivative seed on mu so Dual propagates d/dmu automatically.
  sc_expansion::Parameters<Dual> params{Dual(U, 0.0), Dual(beta, 0.0), Dual(mu - delta, 1.0), true, Dual(-delta, 0.0)};

  sc_expansion::FreeEnergyCalculator<Dual> calculator(params, /*order=*/1, /*override_fm=*/1, /*allow_self_loops=*/true);

  // Order-1 integrand is constant in tau; evaluate at tau=0, the same point chosen
  // by the n=1 SJT quadrature in compute_infinite_U_coefficient.
  std::vector<double> taus = {0.0};
  Dual sum_T               = calculator.compute_sum_diagrams(taus, /*infinite_U=*/false);

  // a_1 = beta * sum_T.value  =>  da_1/dmu = beta * sum_T.derivative
  // density correction  n_1 = -da_1/dmu
  double n1_cpp = -beta * sum_T.derivative;

  double n1_exact = exact_n1_formula(U, beta, mu, delta);

  std::cout << "Expected (Python n_1): " << n1_exact << std::endl;
  std::cout << "C++ n_1 (Dual):        " << n1_cpp << std::endl;
  std::cout << "Difference:            " << std::abs(n1_cpp - n1_exact) << std::endl;

  EXPECT_NEAR(n1_cpp, n1_exact, tol);
}

TEST(McmcAtom, DensityOrder1DualMuShifted) { run_density_dual_check(/*U=*/8.0, /*beta=*/2.0, /*mu=*/2.5, /*delta=*/-0.5, /*tol=*/1e-10); }

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
