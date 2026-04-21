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

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
