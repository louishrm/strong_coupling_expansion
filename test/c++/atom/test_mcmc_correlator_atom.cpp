/*
 * MCMC test for the density-density correlator series, via
 * sc_expansion::atomic::Configuration driven by SumDiagrams' rooted
 * (density-density) constructor.
 *
 * Reference values come from test_correlator_atom.cpp's exact infinite-U
 * dimer-cluster coefficients. Like that test, we pass override_lm = 1 to the
 * rooted SumDiagrams constructor so each kept diagram enters with multiplier
 * 1 instead of its actual Z² embedding count — matching the per-diagram,
 * no-lattice-sum convention the references were computed under.
 *
 * Usage:  mpirun -np 4 ./test_mcmc_correlator_atom
 *         (also works with 1 rank: ./test_mcmc_correlator_atom)
 */

#include <gtest/gtest.h>
#include "sc_expansion/atomic/configuration.hpp"
#include "sc_expansion/atomic/sum_diagrams.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/measure.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <mpi/mpi.hpp>
#include <cmath>
#include <memory>

namespace {
[[maybe_unused]] constexpr int SPIN_DOWN = 0;
constexpr int SPIN_UP                    = 1;
} // namespace

static void run_mcmc_correlator_check(int order, double U, double beta, double mu, std::vector<int> const &r, int s1, int s2,
                                      double exact_coeff, double rel_tol, int n_cycles) {

  mpi::communicator world;

  double alpha     = 0.001;
  int n_warmup     = 2000;
  int length_cycle = 1;

  sc_expansion::Parameters<double> params{U, beta, mu, 0.0, true};
  sc_expansion::atomic::SumDiagrams<double> calculator(params, order, r, s1, s2, /*override_lm=*/1);

  double reference_integral        = 0.0;
  double signed_reference_integral = 0.0;

  if (world.rank() == 0) {
    auto coeff_map      = calculator.density_density_infinite_U_coefficient();
    auto [ref_abs, ref_signed] = coeff_map.at(calculator.get_target_d_sq());
    reference_integral         = ref_abs;
    signed_reference_integral  = ref_signed;
  }

  mpi::broadcast(reference_integral, world);
  mpi::broadcast(signed_reference_integral, world);

  auto config = std::make_unique<sc_expansion::atomic::Configuration<double>>(params, order, alpha, calculator);

  int random_seed = 32186222 + world.rank() * 786512;
  int verbosity   = (world.rank() == 0 ? 2 : 0);

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  measure<double> meas(config.get(), reference_integral, signed_reference_integral, n_bins, block_size, mu);
  mc.add_move(move<double>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "defensive_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  if (world.rank() == 0) {
    double mc_mean  = meas.result->mean;
    double mc_error = meas.result->error;
    double rel_err  = std::abs(mc_mean - exact_coeff) / std::abs(exact_coeff);

    std::cout << "Exact:          " << exact_coeff << std::endl;
    std::cout << "MC estimate:    " << mc_mean << std::endl;
    std::cout << "MC error:       " << mc_error << std::endl;
    std::cout << "Relative error: " << rel_err << std::endl;

    EXPECT_LT(rel_err, rel_tol) << "MC correlator estimate " << mc_mean << " deviates from exact " << exact_coeff << " by " << rel_err * 100 << "%";
  }
}

// Sanity check: the SJT sweep over the rooted diagram list should reproduce
// test_correlator_atom.cpp's exact infinite-U coefficients to ~12 digits.
static void run_infinite_u_correlator_check(int order, double U, double beta, double mu, std::vector<int> const &r, int s1, int s2,
                                            double expected_signed_coeff, double tol) {
  sc_expansion::Parameters<double> params{U, beta, mu, 0.0, true};
  sc_expansion::atomic::SumDiagrams<double> calculator(params, order, r, s1, s2, /*override_lm=*/1);

  auto coeff_map     = calculator.density_density_infinite_U_coefficient();
  auto [abs_coeff, signed_coeff] = coeff_map.at(calculator.get_target_d_sq());

  std::cout << "Expected:      " << expected_signed_coeff << std::endl;
  std::cout << "signed_coeff:  " << signed_coeff << std::endl;
  std::cout << "abs_coeff:     " << abs_coeff << std::endl;
  EXPECT_NEAR(signed_coeff, expected_signed_coeff, tol);
}

TEST(McmcCorrelatorAtom, InfiniteUOnSiteSameSpinOrder2) {
  run_infinite_u_correlator_check(/*order=*/2, /*U=*/0.0, /*beta=*/2.0, /*mu=*/1.0,
                                  /*r=*/{0, 0}, SPIN_UP, SPIN_UP,
                                  /*expected_signed_coeff=*/-3.285401994562779e-03, /*tol=*/1e-12);
}

TEST(McmcCorrelatorAtom, InfiniteUOnSiteSameSpinOrder4) {
  run_infinite_u_correlator_check(/*order=*/4, /*U=*/0.0, /*beta=*/2.0, /*mu=*/1.0,
                                  /*r=*/{0, 0}, SPIN_UP, SPIN_UP,
                                  /*expected_signed_coeff=*/-3.002141546904583e-03, /*tol=*/1e-12);
}

// Finite-U exact references for the on-site / NN same-spin coefficients at
// U=8, β=2, μ=1. These are the exact dimer-Hubbard t^n coefficients (not the
// U→∞ limit) — used as MCMC ground truth at the same U value the sampler runs.
TEST(McmcCorrelatorAtom, OnSiteSameSpinOrder2) {
  run_mcmc_correlator_check(/*order=*/2, /*U=*/8.0, /*beta=*/2.0, /*mu=*/1.0,
                            /*r=*/{0, 0}, SPIN_UP, SPIN_UP,
                            /*exact_coeff=*/-0.002844879587836985, /*rel_tol=*/0.10, /*n_cycles=*/100000);
}

TEST(McmcCorrelatorAtom, OnSiteSameSpinOrder4) {
  // Signal ≈ 1.6e-3; 100k cycles → MC error ≈ 2e-4 (~12% rel). Bumped to 200k
  // to converge inside the 10% relative-error tolerance.
  run_mcmc_correlator_check(/*order=*/4, /*U=*/8.0, /*beta=*/2.0, /*mu=*/1.0,
                            /*r=*/{0, 0}, SPIN_UP, SPIN_UP,
                            /*exact_coeff=*/-0.001617825408717906, /*rel_tol=*/0.10, /*n_cycles=*/200000);
}

TEST(McmcCorrelatorAtom, NearestNeighborSameSpinOrder2) {
  run_mcmc_correlator_check(/*order=*/2, /*U=*/8.0, /*beta=*/2.0, /*mu=*/1.0,
                            /*r=*/{1, 0}, SPIN_UP, SPIN_UP,
                            /*exact_coeff=*/-0.05813220573604677, /*rel_tol=*/0.10, /*n_cycles=*/100000);
}

TEST(McmcCorrelatorAtom, NearestNeighborSameSpinOrder4) {
  // Signal ≈ 5e-4 — very small. 100k cycles → MC error ≈ 4e-4 (~80% rel),
  // dominated by noise. Bumped to 5M to drag rel-error under 10%; expect
  // ~1–2 min/rank. The 10% tolerance is intentionally tight to catch any
  // bias regression in the rooted-correlator path.
  run_mcmc_correlator_check(/*order=*/4, /*U=*/8.0, /*beta=*/2.0, /*mu=*/1.0,
                            /*r=*/{1, 0}, SPIN_UP, SPIN_UP,
                            /*exact_coeff=*/-0.0005279366897435804, /*rel_tol=*/0.10, /*n_cycles=*/5000000);
}

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
