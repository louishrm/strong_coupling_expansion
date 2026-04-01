#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <cassert>
#include "sc_expansion/diagram2.hpp"
#include "sc_expansion/graph.hpp"

using namespace sc_expansion;

// =====================================================================
// Helper
// =====================================================================

static double single_site_free_multiplicity(std::vector<uint8_t> const &adjmat, int V, bool bipartite_only = true) {
  Graph graph(adjmat, V, bipartite_only);
  std::vector<VertexType<1, double> *> vt;
  Diagram2<1, double> diagram(graph, vt);
  return diagram.get_free_multiplicity();
}

static void check(bool condition, const char *msg) {
  if (!condition) {
    std::cerr << "FAILED: " << msg << std::endl;
    std::abort();
  }
  std::cout << "  PASS: " << msg << std::endl;
}

static void check_eq(double a, double b, const char *msg) {
  if (std::abs(a - b) > 1e-10) {
    std::cerr << "FAILED: " << msg << " (got " << a << ", expected " << b << ")" << std::endl;
    std::abort();
  }
  std::cout << "  PASS: " << msg << std::endl;
}

// =====================================================================
// Free multiplicity tests (N_sites=1, single-site expansion)
// =====================================================================

void test_free_multiplicity() {
  std::cout << "--- Free multiplicity ---" << std::endl;

  check_eq(single_site_free_multiplicity({0, 1, 1, 0}, 2), 4, "D2a fm=4");
  check_eq(single_site_free_multiplicity({0, 1, 0, 0, 0, 1, 1, 0, 0}, 3, false), 12, "D3a fm=12");
  check_eq(single_site_free_multiplicity({0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0}, 4), 36, "D4a fm=36");
  check_eq(single_site_free_multiplicity({0, 1, 1, 1, 0, 0, 1, 0, 0}, 3), 16, "D4b fm=16");
  check_eq(single_site_free_multiplicity({0, 2, 2, 0}, 2), 4, "D4c fm=4");

  // D6 variants
  check_eq(single_site_free_multiplicity({0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0}, 4), 64, "D6_1 fm=64");
  check_eq(single_site_free_multiplicity({0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0}, 4), 64, "D6_2 fm=64");
  check_eq(single_site_free_multiplicity({0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0},
                                         6),
           400, "D6a fm=400");
  check_eq(single_site_free_multiplicity({0, 3, 3, 0}, 2), 4, "D6b fm=4");
  check_eq(single_site_free_multiplicity({0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0}, 4), 64, "D6c fm=64");
  check_eq(single_site_free_multiplicity({0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0}, 5), 144, "D6d fm=144");
  check_eq(single_site_free_multiplicity({0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0}, 4), 64, "D6e fm=64");
  check_eq(single_site_free_multiplicity({0, 2, 1, 2, 0, 0, 1, 0, 0}, 3), 16, "D6f fm=16");
  check_eq(single_site_free_multiplicity({0, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0}, 4), 36, "D6g fm=36");

  // D8a (8-cycle)
  std::vector<uint8_t> D8a(64, 0);
  for (int i = 0; i < 8; ++i) {
    D8a[i * 8 + (i + 1) % 8]   = 1;
    D8a[((i + 1) % 8) * 8 + i] = 1;
  }
  check_eq(single_site_free_multiplicity(D8a, 8), 4900, "D8a fm=4900");
  check_eq(single_site_free_multiplicity({0, 4, 4, 0}, 2), 4, "D8b fm=4");
}

// =====================================================================
// Dimer spatial configuration tests (N_sites=2)
// =====================================================================

void test_dimer_spatial_configs() {
  std::cout << "--- Dimer spatial configurations ---" << std::endl;

  // D4b: 2 distinct spatial configs, weight 18 each, total 36
  {
    Graph graph({0, 1, 1, 1, 0, 0, 1, 0, 0}, 3);
    std::vector<VertexType<2, double> *> vt;
    Diagram2<2, double> diagram(graph, vt);
    auto const &spatial = diagram.get_spatial_configurations();

    check_eq(spatial.size(), 2, "D4b: 2 spatial configs");
    double total = 0;
    for (auto const &sc : spatial) total += sc.weight;
    check_eq(total, 36.0, "D4b: total weight 36");
  }

  // D2a: 1 spatial config, weight 6
  {
    Graph graph({0, 1, 1, 0}, 2);
    std::vector<VertexType<2, double> *> vt;
    Diagram2<2, double> diagram(graph, vt);
    auto const &spatial = diagram.get_spatial_configurations();

    check_eq(spatial.size(), 1, "D2a: 1 spatial config");
    check_eq(spatial[0].weight, 6.0, "D2a: weight 6");
  }

  // D4c: 1 spatial config, weight 6
  {
    Graph graph({0, 2, 2, 0}, 2);
    std::vector<VertexType<2, double> *> vt;
    Diagram2<2, double> diagram(graph, vt);
    auto const &spatial = diagram.get_spatial_configurations();

    check_eq(spatial.size(), 1, "D4c: 1 spatial config");
    check_eq(spatial[0].weight, 6.0, "D4c: weight 6");
  }

  // D6c: 2 spatial configs, weights 54 and 162
  {
    Graph graph({0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0}, 4);
    std::vector<VertexType<2, double> *> vt;
    Diagram2<2, double> diagram(graph, vt);
    auto const &spatial = diagram.get_spatial_configurations();

    check_eq(spatial.size(), 2, "D6c: 2 spatial configs");
    check_eq(spatial[0].weight, 54.0, "D6c: weight[0]=54");
    check_eq(spatial[1].weight, 162.0, "D6c: weight[1]=162");
  }
}

// =====================================================================
// Global configuration tests
// =====================================================================

void test_global_configs() {
  std::cout << "--- Global configurations ---" << std::endl;

  // D2a atom: 1 config, weight 4
  {
    Graph graph({0, 1, 1, 0}, 2);
    std::vector<VertexType<1, double> *> vt;
    Diagram2<1, double> diagram(graph, vt);
    auto const &configs = diagram.get_valid_configurations();

    check_eq(configs.size(), 1, "D2a atom: 1 config");
    check_eq(configs[0].weight, 4.0, "D2a atom: weight 4");
  }

  // D2a dimer: 1 config, weight 12
  {
    Graph graph({0, 1, 1, 0}, 2);
    std::vector<VertexType<2, double> *> vt;
    Diagram2<2, double> diagram(graph, vt);
    auto const &configs = diagram.get_valid_configurations();

    check_eq(configs.size(), 1, "D2a dimer: 1 config");
    check_eq(configs[0].weight, 12.0, "D2a dimer: weight 12");
  }

  // D4b atom: 2 configs, total weight 32
  {
    Graph graph({0, 1, 1, 1, 0, 0, 1, 0, 0}, 3);
    std::vector<VertexType<1, double> *> vt;
    Diagram2<1, double> diagram(graph, vt);
    auto const &configs = diagram.get_valid_configurations();

    check_eq(configs.size(), 2, "D4b atom: 2 configs");
    double total = 0;
    for (auto const &c : configs) total += c.weight;
    check_eq(total, 32.0, "D4b atom: total weight 32");
  }
}

// =====================================================================
// Numerical evaluation tests
// =====================================================================

void test_evaluation() {
  std::cout << "--- Numerical evaluation ---" << std::endl;

  double U = 8.0, beta = 1.0, mu = 2.0;
  Parameters<double> params{U, beta, mu, 0.0, true};
  std::vector<double> taus = {0.5, 0.0};

  // D2a atom
  {
    Graph graph({0, 1, 1, 0}, 2);
    VertexType<1, double> vt1(2);
    std::vector<VertexType<1, double> *> vt = {&vt1};
    Diagram2<1, double> diagram(graph, vt);
    HubbardSolver<1, double> solver(params);

    double val = diagram.evaluate(taus, solver, false);
    check(std::isfinite(val), "D2a atom: evaluate is finite");
    check(val != 0.0, "D2a atom: evaluate is non-zero");

    double val2 = diagram.evaluate(taus, solver, false);
    check_eq(val, val2, "D2a atom: cache gives same result");
  }

  // D2a dimer
  {
    Parameters<double> dimer_params{U, beta, mu, 1.0, true};
    Graph graph({0, 1, 1, 0}, 2);
    VertexType<2, double> vt1(2);
    std::vector<VertexType<2, double> *> vt = {&vt1};
    Diagram2<2, double> diagram(graph, vt);
    HubbardSolver<2, double> solver(dimer_params);

    double val = diagram.evaluate(taus, solver, false);
    check(std::isfinite(val), "D2a dimer: evaluate is finite");
    check(val != 0.0, "D2a dimer: evaluate is non-zero");
  }
}

// =====================================================================

int main() {
  test_free_multiplicity();
  test_dimer_spatial_configs();
  test_global_configs();
  test_evaluation();
  std::cout << "\nAll tests passed." << std::endl;
  return 0;
}
