/*
 * MCMC test for the staggered-dimer density-density correlator ⟨n_{r,s2}(0)
 * n_{0,s1}(0)⟩ on a finite 3-dimer cluster, via sc_expansion::dimer::Configuration
 * driven by dimer::SumDiagrams' rooted (finite-cluster) constructor.
 *
 * Same 3-dimer triangle as the free-energy MCMC (test_mcmc_dimer.cpp):
 *   A=(0,0), B=(1,0), C=(0,1) on the staggered superlattice. Physical sites:
 *   A:(0,0),(1,0)  B:(2,0),(3,0)  C:(1,1),(2,1).
 *
 * The MCMC reports the order-n coefficient with the density n(0) PINNED at a
 * single reference site (cluster_positions[0] = the pendant site (0,0)), the
 * intensive/local convention for a ⟨n(r)n(0)⟩ correlator (pin_origin=true; no
 * sweep over cells, no ÷ n_cluster_sites). Reference values are exact
 * finite-cluster ED coefficients at that same single site. Uses the
 * uniform-reference defensive estimator (measure_dimer): W = |f + alpha|.
 *
 * Orbital/spin encoding: SPIN_DOWN=0, SPIN_UP=1.
 *
 * Usage:  mpirun -np 4 ./test_mcmc_dimer_correlator
 *         (also works with 1 rank: ./test_mcmc_dimer_correlator)
 */

#include <gtest/gtest.h>
#include "sc_expansion/dimer/configuration.hpp"
#include "sc_expansion/dimer/sum_diagrams.hpp"
#include "sc_expansion/dimer/measure_dimer.hpp"
#include "sc_expansion/cumulant.hpp"
#include "sc_expansion/args.hpp"
#include "sc_expansion/move.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <mpi/mpi.hpp>
#include <cmath>
#include <iomanip>
#include <memory>
#include <utility>
#include <vector>

namespace {
  constexpr int SPIN_DOWN                = 0;
  [[maybe_unused]] constexpr int SPIN_UP = 1;
} // namespace

static void run_mcmc_dimer_correlator_check(int order, std::vector<int> const &r, int s1, int s2, double exact_coeff, double rel_tol, int n_cycles,
                                            double beta = 1.0, std::vector<std::pair<int, int>> const &cluster_positions = {{0, 0}, {1, 0}, {0, 1}}) {
  mpi::communicator world;

  double U  = 12.0;
  double mu = 3.0;
  double t  = 1.0; // intra-dimer hopping

  // 3-dimer triangular cluster on the staggered (non-bipartite) superlattice.
  // cluster_positions[0] is the pinned reference dimer (density on its site 0).
  // n_cluster_sites is the per-dimer normaliser for the sweep convention; in the
  // pinned (single-site) convention used here it is unused (no ÷n_cluster_sites).
  int n_cluster_sites = (int)cluster_positions.size();

  double alpha     = 0.001;
  int n_warmup     = 2000;
  int length_cycle = 1;

  sc_expansion::Parameters<double> params{U, beta, mu, t, /*bipartite=*/false};
  // pin_origin=true: ⟨n(r)n(0)⟩ anchored at the single pendant reference site
  // cluster_positions[0] = (0,0), matching the finite-cluster ED reference.
  sc_expansion::dimer::SumDiagrams<double> calculator(params, order, r, s1, s2, cluster_positions, n_cluster_sites, /*pin_origin=*/true);

  auto config =
     std::make_unique<sc_expansion::dimer::Configuration<double, sc_expansion::dimer::SumDiagrams<double>>>(params, order, calculator, alpha);

  int random_seed = 32186222 + world.rank() * 786512;
  int verbosity   = (world.rank() == 0 ? 2 : 0);

  triqs::mc_tools::mc_generic<double> mc("", random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  measure_dimer<double> meas(config.get(), n_bins, block_size, alpha);
  mc.add_move(move<double>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "dimer_correlator_measure");

  mc.warmup_and_accumulate(n_warmup, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  if (world.rank() == 0) {
    double mc_mean  = meas.result->coeff;
    double mc_error = meas.result->error;
    double rel_err  = std::abs(mc_mean - exact_coeff) / std::abs(exact_coeff);

    std::cout << "r = (" << r[0] << "," << r[1] << "), spins = (" << s1 << "," << s2 << "), order " << order << std::endl;
    std::cout << "Exact (ED):     " << exact_coeff << std::endl;
    std::cout << "MC estimate:    " << mc_mean << std::endl;
    std::cout << "MC error:       " << mc_error << std::endl;
    std::cout << "Relative error: " << rel_err << std::endl;

    EXPECT_LT(rel_err, rel_tol) << "MC estimate " << mc_mean << " deviates from exact " << exact_coeff << " by " << rel_err * 100 << "%";
  }
}

// ---------------------------------------------------------------------------
// Deterministic diagnostic (no MC noise). Builds the SAME catalog the MCMC uses
// (order, r, spins, 3-dimer cluster), dumps each diagram's marks/sites/spins/
// sign/spatial-configs, then integrates the diagram sum f(τ) = density_density(τ)
// over the hypercube [0,β]^order by composite Simpson quadrature — which is
// EXACTLY what the MC ratio estimator coeff = ∫ f d^order τ targets. Comparing
// this deterministic integral against both the MC value and the ED reference
// localises a discrepancy to either the diagram definition (integral ≈ MC, far
// from ED) or the MC estimator (integral ≈ ED, far from MC).
static double simpson_weight(int i, int n) { return (i == 0 || i == n) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0); }

// cluster_positions[0] is the PINNED reference dimer; the density n(0) sits on its
// within-dimer site 0. Reorder it to pin a different reference site: {{0,0},...}
// pins A's pendant site (0,0) → 2·cumB; {{1,0},{0,0},{0,1}} pins B's site (2,0),
// whose density-site carries 2 inter-bonds → 2·cumA (exercises the dirs=[0,1]
// "hop-at-density-site" config that is geometrically absent at a pendant).
static void dump_dimer_correlator_diagnostic(int order, std::vector<int> const &r, int s1, int s2, double exact_coeff, double U = 12.0,
                                             double mu = 3.0, double beta = 1.0, int n_grid = 80,
                                             std::vector<std::pair<int, int>> const &cluster_positions = {{0, 0}, {1, 0}, {0, 1}}) {
  mpi::communicator world;
  if (world.rank() != 0) return; // single-rank deterministic dump

  double t            = 1.0;
  int n_cluster_sites = (int)cluster_positions.size();

  sc_expansion::Parameters<double> params{U, beta, mu, t, /*bipartite=*/false};
  // pin_origin=true: anchor n(0) at cluster_positions[0] = the pendant site (0,0),
  // matching the single-site ED reference (no sweep over inequivalent cells, no
  // ÷n_cluster_sites). With the old sweep convention the integral over-counted by
  // ~13× (U=12) because it averaged the pendant site with the bonded interior sites.
  sc_expansion::dimer::SumDiagrams<double> calculator(params, order, r, s1, s2, cluster_positions, n_cluster_sites, /*pin_origin=*/true);
  auto const &solver = calculator.get_solver();
  int n_diag         = calculator.get_n_diagrams();

  std::cout << "\n=== DIAGNOSTIC: order " << order << ", r=(" << r[0] << "," << r[1] << "), spins=(" << s1 << "," << s2 << "), 3-dimer cluster ===\n";
  std::cout << "n_diagrams in catalog: " << n_diag << "\n";

  int di = 0;
  for (auto const &d : calculator.get_diagrams()) {
    auto const &mk = d.get_marks();
    auto const &st = d.get_sites();
    auto const &sp = d.get_mark_spins();
    auto const &sc = d.get_spatial_configurations();
    std::cout << "  diagram[" << di++ << "] V=" << d.get_graph().get_V() << " sign=" << d.get_diagram_sign() << " marks=(" << mk[0] << "," << mk[1]
              << ") sites=(" << st[0] << "," << st[1] << ") spins=(" << sp[0] << "," << sp[1] << ") n_spatial=" << sc.size()
              << " free_mult(sum w)=" << d.get_free_multiplicity() << "\n";
    for (auto const &c : sc) {
      std::cout << "      dirs=[";
      for (size_t k = 0; k < c.directions.size(); ++k) std::cout << (int)c.directions[k] << (k + 1 < c.directions.size() ? "," : "");
      std::cout << "] weight=" << c.weight << "\n";
    }
    // Decode each valid configuration's op_ids → the within-dimer SITE each
    // hopping leg attaches to (site = orbital&1, spin = orbital>>1, ACTION_BIT=4
    // ⇒ creation). Legs are ordered vertex 0, then vertex 1, ... so for V=2 the
    // last 2 entries are the MARKED vertex's legs. This makes the two physical
    // processes explicit: the marked legs sit on site 0 (= density site) in one
    // config and on site 1 (≠ density site) in the other.
    std::cout << "    valid_configs (legs ordered v0..v" << (d.get_graph().get_V() - 1) << "; density is on site " << st[0] << "):\n";
    for (auto const &vc : d.get_valid_configurations()) {
      std::cout << "      [";
      for (size_t k = 0; k < vc.config.size(); ++k) {
        uint8_t op = vc.config[k];
        bool cre   = (op & 4) != 0;
        int orb    = op & 3;
        std::cout << (cre ? "c+" : "c-") << "@site" << (orb & 1) << (((orb >> 1) & 1) ? "up" : "dn") << (k + 1 < vc.config.size() ? ", " : "");
      }
      std::cout << "]  weight=" << vc.weight << "\n";
    }
  }

  // Composite-Simpson integral over [0,beta]^order, for the full sum AND each
  // diagram separately (so we can see if one sector alone already equals ED).
  std::vector<double> taus(order, 0.0);
  double h        = beta / n_grid;
  double integral = 0.0;
  std::vector<double> per_diagram(n_diag, 0.0);
  if (order == 2) {
    for (int i = 0; i <= n_grid; ++i) {
      for (int j = 0; j <= n_grid; ++j) {
        taus[0]   = i * h;
        taus[1]   = j * h;
        double sw = simpson_weight(i, n_grid) * simpson_weight(j, n_grid);
        calculator.mark_all_dirty(); // force recompute for these taus
        integral += sw * calculator.density_density(taus);
        int k = 0;
        for (auto const &d : calculator.get_diagrams())
          per_diagram[k++] += sw * const_cast<sc_expansion::dimer::Diagram<double> &>(d).evaluate(taus, solver);
      }
    }
    double scale = (h / 3.0) * (h / 3.0);
    integral *= scale;
    for (auto &v : per_diagram) v *= scale;
  } else if (order == 1) {
    for (int i = 0; i <= n_grid; ++i) {
      taus[0]   = i * h;
      double sw = simpson_weight(i, n_grid);
      calculator.mark_all_dirty();
      integral += sw * calculator.density_density(taus);
      int k = 0;
      for (auto const &d : calculator.get_diagrams())
        per_diagram[k++] += sw * const_cast<sc_expansion::dimer::Diagram<double> &>(d).evaluate(taus, solver);
    }
    double scale = (h / 3.0);
    integral *= scale;
    for (auto &v : per_diagram) v *= scale;
  }
  for (int k = 0; k < n_diag; ++k) std::cout << "  ∫ diagram[" << k << "] = " << std::setprecision(10) << per_diagram[k] << "\n";

  // Sample value at a generic interior point, to eyeball the integrand scale.
  std::vector<double> sample(order, 0.37);
  if (order >= 2) sample[1] = 0.71;
  calculator.mark_all_dirty();
  double f_sample = calculator.density_density(sample);

  // Order-0 reference (no inter-dimer hops): the isolated-dimer connected
  // correlator κ₂(n_{site_a, σ_a}, n_{site_b, σ_b}) = ⟨n n⟩ − ⟨n⟩⟨n⟩ on the 2-site
  // Hubbard cumulant. Independent of all embedding/hopping-leg machinery, so it
  // isolates the H₀/μ/connected-definition convention. The within-dimer sites and
  // spins are read from the first catalog diagram so this matches whatever (site,
  // spin) pair the rooted generator chose: on-site r=(0,0) ⇒ sites 0/0; intra-dimer
  // r=(1,0) ⇒ sites 0/1. Compare to your ED's order-0.
  if (n_diag > 0) {
    auto const &d0 = calculator.get_diagrams().front();
    auto const &st = d0.get_sites();
    auto const &sp = d0.get_mark_spins();
    sc_expansion::Args<2, double> u({}, {}), p({}, {});
    sc_expansion::CumulantSolver<2, double> cs(u, p, calculator.get_solver(), /*infinite_U=*/false);
    cs.add_static_density(st[0] + sp[0] * 2);
    cs.add_static_density(st[1] + sp[1] * 2);
    double kappa0 = cs.compute_cumulant_decomposition();
    std::cout << "order-0 connected κ₂ (isolated dimer, sites " << st[0] << "/" << st[1] << ", spins " << sp[0] << "/" << sp[1] << "): " << kappa0
              << "\n";
  }

  std::cout << std::setprecision(10);
  std::cout << "f(sample=" << sample[0] << (order >= 2 ? ",0.71" : "") << ") = " << f_sample << "\n";
  std::cout << "Deterministic ∫f d^" << order << "τ over [0," << beta << "]^" << order << " (Simpson n=" << n_grid << "): " << integral << "\n";
  std::cout << "ED reference (order " << order << "):                                " << exact_coeff << "\n";
  std::cout << "ratio  integral/ED = " << integral / exact_coeff << "\n";
  std::cout << "=== END DIAGNOSTIC ===\n\n";
}

TEST(McmcDimerCorrelator, DiagnosticOnSiteOrder2) {
  dump_dimer_correlator_diagnostic(/*order=*/2, /*r=*/{0, 0}, SPIN_UP, SPIN_DOWN, /*exact_coeff=*/0.000675371799);
}

// Second parameter point (U=6, same beta=1, mu=3, intra-t=1) to discriminate a
// pure normalization constant (ratio stays ≈ that of U=12) from a model/value
// mismatch (ratio changes). ED reference at U=6: 0.001334786422.
TEST(McmcDimerCorrelator, DiagnosticOnSiteOrder2_U6) {
  dump_dimer_correlator_diagnostic(/*order=*/2, /*r=*/{0, 0}, SPIN_UP, SPIN_DOWN, /*exact_coeff=*/0.001334786422, /*U=*/6.0);
}

// Intra-dimer 1st-neighbour correlator at r=(1,0): origin (0,0) = pendant site 0,
// partner (1,0) = site 1 of the SAME dimer. ED at U=12, beta=2, mu=3, t_intra=1:
//   order 0:  0.041180973554   (isolated-dimer κ₂, cross-site, ↓@site0 / ↑@site1)
//   order 2: -0.001230255454
// Our convention puts s1 at the origin and s2 at r, so ED's (↑@r, ↓@0) ⇒ s1=↓, s2=↑.
TEST(McmcDimerCorrelator, DiagnosticNeighborIntraOrder2) {
  dump_dimer_correlator_diagnostic(/*order=*/2, /*r=*/{1, 0}, /*s1=*/SPIN_DOWN, /*s2=*/SPIN_UP,
                                   /*exact_coeff=*/-0.001230255454, /*U=*/12.0, /*mu=*/3.0, /*beta=*/2.0, /*n_grid=*/160);
}

// Config-A probe: on-site (r=(0,0)) correlator pinned at B's site (2,0) instead of
// A's pendant. B's density-site carries 2 inter-bonds and B's other site (3,0) is
// the pendant, so this isolates the dirs=[0,1] "hop-AT-density-site" config →
// integral should be 2·cumA (vs 2·cumB at the pendant idx0). This is the config
// that's invisible to the idx0/NN pendant tests but DOES appear (weight 3) in the
// TD limit. ED at site (2,0), on-site, order 2, U=12 beta=2 mu=3: 0.011287004262
// (order-0 κ₂ at beta=2 is -0.243162205625, exact). beta=2 to match this ED.
TEST(McmcDimerCorrelator, DiagnosticOnSiteAtBSite0Order2) {
  dump_dimer_correlator_diagnostic(/*order=*/2, /*r=*/{0, 0}, SPIN_UP, SPIN_DOWN, /*exact_coeff=*/0.011287004262, /*U=*/12.0, /*mu=*/3.0,
                                   /*beta=*/2.0, /*n_grid=*/160, /*cluster_positions=*/{{1, 0}, {0, 0}, {0, 1}});
}

// U=12, beta=1, mu=3. On-site (r=(0,0)) up-down correlator ⟨n↑ n↓⟩, order 2.
TEST(McmcDimerCorrelator, OnSiteUpDownOrder2) {
  run_mcmc_dimer_correlator_check(/*order=*/2, /*r=*/{0, 0}, SPIN_UP, SPIN_DOWN,
                                  /*exact_coeff=*/0.000675371799, /*rel_tol=*/0.12, /*n_cycles=*/1000000);
}

// U=12, beta=2, mu=3. Intra-dimer 1st-neighbour ⟨n_{(1,0),↑} n_{(0,0),↓}⟩, order 2.
// ED reference -0.001230255454 (see DiagnosticNeighborIntraOrder2). beta=2 here.
TEST(McmcDimerCorrelator, NeighborIntraUpDownOrder2) {
  run_mcmc_dimer_correlator_check(/*order=*/2, /*r=*/{1, 0}, /*s1=*/SPIN_DOWN, /*s2=*/SPIN_UP,
                                  /*exact_coeff=*/-0.001230255454, /*rel_tol=*/0.12, /*n_cycles=*/1000000, /*beta=*/2.0);
}

// U=12, beta=2, mu=3. On-site (r=(0,0)) ⟨n↑ n↓⟩ pinned at B's site (2,0), order 2.
// Isolates the dirs=[0,1] "hop-AT-density-site" config (= 2·cumA), unreachable at a
// pendant. ED reference 0.011287004262 (see DiagnosticOnSiteAtBSite0Order2).
TEST(McmcDimerCorrelator, OnSiteAtBSite0Order2) {
  run_mcmc_dimer_correlator_check(/*order=*/2, /*r=*/{0, 0}, SPIN_UP, SPIN_DOWN,
                                  /*exact_coeff=*/0.011287004262, /*rel_tol=*/0.12, /*n_cycles=*/1000000, /*beta=*/2.0,
                                  /*cluster_positions=*/{{1, 0}, {0, 0}, {0, 1}});
}

// Additional ED references to fill in as they arrive (orders 3/4, other r/spins).
// We now have a full series at site (2,0), beta=2, on-site ⟨n↑ n↓⟩ for higher-order
// validation: order 3 -0.000084368305, order 4 0.004661852767, order 5 -0.000364314888,
// order 6 0.000108298703, order 7 -0.000112168537.
// TEST(McmcDimerCorrelator, OnSiteAtBSite0Order4) {
//   run_mcmc_dimer_correlator_check(4, {0, 0}, SPIN_UP, SPIN_DOWN, 0.004661852767, 0.15, 5000000, 2.0, {{1,0},{0,0},{0,1}});
// }

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
