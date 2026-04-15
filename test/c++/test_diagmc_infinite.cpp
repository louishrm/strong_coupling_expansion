/*
 * DiagMC tests on the infinite lattice (columnar dimer, bipartite).
 *
 * These compute the density (d/dmu of free energy), so Dual numbers are
 * used with mu as the differentiation variable.
 *
 * ED references (U=8.0, beta=2.0, mu=3.0, t_intra=1.0):
 *   Order 4 density coefficient = -0.02523099801644088
 *   Order 6 density coefficient =  0.00293783601058773
 *
 * Usage:  mpirun -np 4 ./test_diagmc_infinite
 */

#include <gtest/gtest.h>
#include "sc_expansion/diagmc_configuration.hpp"
#include "sc_expansion/diagmc_move.hpp"
#include "sc_expansion/measure_diagmc.hpp"
#include "sc_expansion/dual.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <mpi/mpi.hpp>
#include <cmath>
#include <memory>

using namespace sc_expansion;

// =====================================================================
// Order 4, infinite lattice (bipartite = true)
// 3 diagrams. Expected: -0.02523099801644088
// =====================================================================

TEST(DiagMCInfinite, Order4) {
  mpi::communicator world;

  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order    = 4;
  double alpha = 0.001;
  Parameters<Dual> params{Dual(U, 0.0), Dual(beta, 0.0), Dual(mu, 1.0), Dual(t_intra, 0.0), true};

  FreeEnergyCalculator<2, Dual> calculator(params, order);
  auto config = std::make_unique<DiagMCConfiguration<Dual>>(params, order, calculator, alpha);

  if (world.rank() == 0) {
    std::cout << "\n=== DiagMC Infinite Lattice Order 4 (density) ===" << std::endl;
    std::cout << "N_diagrams = " << config->get_n_diagrams() << std::endl;
    auto const &graphs = calculator.get_graphs();
    for (int d = 0; d < config->get_n_diagrams(); d++) {
      std::cout << "  d=" << (d + 1) << " V=" << graphs[d].get_V() << " SymFactor=" << graphs[d].get_symmetry_factor()
                << " FreeMult=" << graphs[d].get_free_multiplicity() << std::endl;
    }
  }

  int n_cycles     = 100000;
  int n_warmup     = 2000;
  int length_cycle = 1;
  int random_seed  = 32186222 + world.rank() * 786512;
  int verbosity    = (world.rank() == 0 ? 2 : 0);

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  MeasureDiagMC<Dual> meas(config.get(), n_bins, block_size, verbosity);
  mc.add_move(DiagMCMove<Dual>(config.get(), mc.get_rng()), "diagmc_move");
  mc.add_measure(meas, "diagmc_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  if (world.rank() == 0) {
    double mc_coeff = meas.result->coeff;
    double exact    = 0.02523099801644088;
    double rel_err  = std::abs(mc_coeff - exact) / std::abs(exact);

    std::cout << "\nExact:                   " << exact << std::endl;
    std::cout << "DiagMC coefficient:      " << mc_coeff << std::endl;
    std::cout << "DiagMC error:            " << meas.result->error << std::endl;
    std::cout << "Relative error:          " << rel_err << std::endl;

    EXPECT_LT(rel_err, 0.10) << "DiagMC estimate " << mc_coeff << " deviates from exact " << exact << " by " << rel_err * 100 << "%";
  }
}

// =====================================================================
// Order 6, infinite lattice (bipartite = true)
// 7 diagrams. Expected: 0.00293783601058773
// This is where DiagMC should show real benefit — the V=2 diagram
// is expensive but carries very little signal.
// =====================================================================

TEST(DiagMCInfinite, Order6) {
  mpi::communicator world;

  double U = 8.0, beta = 2.0, mu = 3.0, t_intra = 1.0;
  int order    = 6;
  double alpha = 0.001;
  Parameters<Dual> params{Dual(U, 0.0), Dual(beta, 0.0), Dual(mu, 1.0), Dual(t_intra, 0.0), true};

  FreeEnergyCalculator<2, Dual> calculator(params, order);
  auto config = std::make_unique<DiagMCConfiguration<Dual>>(params, order, calculator, alpha);

  if (world.rank() == 0) {
    std::cout << "\n=== DiagMC Infinite Lattice Order 6 (density) ===" << std::endl;
    std::cout << "N_diagrams = " << config->get_n_diagrams() << std::endl;
    auto const &graphs = calculator.get_graphs();
    for (int d = 0; d < config->get_n_diagrams(); d++) {
      std::cout << "  d=" << (d + 1) << " V=" << graphs[d].get_V() << " SymFactor=" << graphs[d].get_symmetry_factor()
                << " FreeMult=" << graphs[d].get_free_multiplicity() << std::endl;
    }
  }

  int n_cycles     = 100000;
  int n_warmup     = 10000;
  int length_cycle = 1;
  int random_seed  = 55186444 + world.rank() * 786512;
  int verbosity    = (world.rank() == 0 ? 2 : 0);

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  MeasureDiagMC<Dual> meas(config.get(), n_bins, block_size, verbosity);
  mc.add_move(DiagMCMove<Dual>(config.get(), mc.get_rng()), "diagmc_move");
  mc.add_measure(meas, "diagmc_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  if (world.rank() == 0) {
    double mc_coeff = meas.result->coeff;
    double exact    = -2.93783601e-03;
    double rel_err  = std::abs(mc_coeff - exact) / std::abs(exact);

    std::cout << "\nExact:                   " << exact << std::endl;
    std::cout << "DiagMC coefficient:      " << mc_coeff << std::endl;
    std::cout << "DiagMC error:            " << meas.result->error << std::endl;
    std::cout << "Relative error:          " << rel_err << std::endl;

    EXPECT_LT(rel_err, 0.15) << "DiagMC estimate " << mc_coeff << " deviates from exact " << exact << " by " << rel_err * 100 << "%";
  }
}

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
