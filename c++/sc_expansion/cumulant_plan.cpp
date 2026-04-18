#include "cumulant_plan.hpp"
#include "dual.hpp"
#include <algorithm>

namespace sc_expansion {

  // Invert an argsort: if argsort[sorted_pos] = stable_pos, then inv[stable_pos] = sorted_pos.
  static std::vector<int> invert_argsort(std::vector<int> const &argsort) {
    std::vector<int> inv(argsort.size());
    for (size_t i = 0; i < argsort.size(); ++i) inv[argsort[i]] = (int)i;
    return inv;
  }

  // Build the leaf Args in STABLE basis — ops listed in ascending stable-position order
  // as stored in the plan's leaf. The Args constructor will sort by τ internally and
  // record leaf_args.permutation_sign; G0n multiplies by that sign on the way out, so
  // the returned value is the correlator for ops listed in STABLE order. Combined with
  // the stable-basis sign bookkeeping in the plan, this gives the stable-basis cumulant.
  template <int N_sites, typename T>
  static Args<N_sites, T> build_leaf_args_stable(CumulantPlan::LeafOps const &leaf,
                                                 Args<N_sites, T> const &master_unprimed,
                                                 Args<N_sites, T> const &master_primed,
                                                 std::vector<int> const &inv_argsort_u,
                                                 std::vector<int> const &inv_argsort_p) {
    int n = (int)leaf.u_global_idx.size();
    std::vector<double> taus;
    std::vector<FermionOperator<N_sites, T>> ops;
    taus.reserve(2 * n);
    ops.reserve(2 * n);
    for (int i = 0; i < n; ++i) {
      int sp = inv_argsort_p[leaf.p_global_idx[i]]; // stable → sorted pos in eval master
      int su = inv_argsort_u[leaf.u_global_idx[i]];
      taus.push_back(master_primed.taus[sp]);
      ops.push_back(master_primed.ops[sp]);
      taus.push_back(master_unprimed.taus[su]);
      ops.push_back(master_unprimed.ops[su]);
    }
    return Args<N_sites, T>(std::move(taus), std::move(ops));
  }

  template <int N_sites, typename T>
  T evaluate_plan(CumulantPlan const &plan,
                  Args<N_sites, T> const &master_unprimed,
                  Args<N_sites, T> const &master_primed,
                  HubbardSolver<N_sites, T> const &solver,
                  bool infinite_U) {

    if (plan.root_id < 0) return T(0.0);

    std::vector<int> inv_argsort_u = invert_argsort(master_unprimed.argsort);
    std::vector<int> inv_argsort_p = invert_argsort(master_primed.argsort);

    std::vector<T> value;
    value.reserve(plan.nodes.size());

    for (size_t i = 0; i < plan.nodes.size(); ++i) {
      auto const &node = plan.nodes[i];

      Args<N_sites, T> args = build_leaf_args_stable<N_sites, T>(node.leaf, master_unprimed, master_primed, inv_argsort_u, inv_argsort_p);
      T v = infinite_U ? solver.G0n_infinite_U(args) : solver.G0n(args);

      for (auto const &term : node.subtraction_terms) {
        T prod = T((double)term.sign);
        for (int fid : term.factor_node_ids) prod = prod * value[fid];
        v = v + prod;
      }
      value.push_back(v);
    }

    // Plan works in STABLE basis. compute_cumulant_decomposition (the reference) returns
    // the τ-sorted-basis cumulant. The conversion between bases is a permutation sign
    // on each side (see vertex.cpp:59 for the same factor):
    //   C_τ_sorted = master_u.permutation_sign * master_p.permutation_sign * C_stable
    T C_stable = value[plan.root_id];
    return C_stable * T(master_unprimed.permutation_sign) * T(master_primed.permutation_sign);
  }

  // Explicit instantiations
  template double evaluate_plan<1, double>(CumulantPlan const &, Args<1, double> const &, Args<1, double> const &,
                                           HubbardSolver<1, double> const &, bool);
  template Dual evaluate_plan<1, Dual>(CumulantPlan const &, Args<1, Dual> const &, Args<1, Dual> const &,
                                       HubbardSolver<1, Dual> const &, bool);
  template double evaluate_plan<2, double>(CumulantPlan const &, Args<2, double> const &, Args<2, double> const &,
                                           HubbardSolver<2, double> const &, bool);
  template Dual evaluate_plan<2, Dual>(CumulantPlan const &, Args<2, Dual> const &, Args<2, Dual> const &,
                                       HubbardSolver<2, Dual> const &, bool);

} // namespace sc_expansion
