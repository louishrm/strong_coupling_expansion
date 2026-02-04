import numpy as np
import matplotlib.pyplot as plt
from . import robustpade
import os
import h5py
import math

# ========== utilities ==========
def orders_grid(N):
    base = int(N / 2.0)
    offsets = [-1, 0, 1]
    grid = []
    for i in offsets:
        row = []
        for j in offsets:
            m, n = base + i, base + j
            row.append((max(0, m), max(0, n)))
        grid.append(row)
    return grid  # shape 3x3


def sample_coeffs(mean, err, n_samples, rng=None):
    rng = rng or np.random.default_rng()
    return rng.normal(mean, err, size=(n_samples, len(mean)))


def compute_pade_with_consistent_orders(coeff_mean, coeff_samples, U_vals, n, m, tol=1e-5):
    """Compute Padé mean and sampled approximants keeping only those with matching (mu,nu)."""
    try:
        r_mean, a_mean, b_mean, mu_mean, nu_mean, *_ = robustpade.pade_approx(coeff_mean, [n, m], tol=tol)
    except Exception:
        return np.full_like(U_vals, np.nan), np.zeros_like(U_vals), None, None

    valid_curves = []
#    i_valid_curve = 0
    for coeff in coeff_samples:
        try:
            r, a, b, mu, nu, *_ = robustpade.pade_approx(coeff, [n, m], tol=tol)
            curve = r(U_vals).real
            #if mu == mu_mean and nu == nu_mean:
            if (np.all(abs(curve) < 2)):
                valid_curves.append(curve)
#            i_valid_curve += 1
        except Exception:
            pass

#    print(f'Acceptance of curves: {i_valid_curve/len(coeff_samples):.2f}')

    if len(valid_curves) > 0:
        valid_curves = np.array(valid_curves)
        mean_curve = np.mean(valid_curves, axis=0)
        std_curve = np.std(valid_curves, axis=0)
    else:
        mean_curve = np.full_like(U_vals, np.nan)
        std_curve = np.zeros_like(U_vals)

    return mean_curve, std_curve, mu_mean, nu_mean


# ========== main analysis ==========
def analyze_single_mu_plots(filepath, U, err_coeff, limit=8.0, n_points=400, n_bootstrap=1000, seed=0, tol=1e-5, show_plots=True):
    """Compute all Padé approximants for a given μ and produce 4x3 grid of plots."""
    rng = np.random.default_rng(seed)
    U_vals = np.linspace(0, limit, n_points)

    # Load coefficients and errors
    with h5py.File(filepath, 'r') as f:

        coeffs_raw = f["coeffs"][:-1]
        errors_raw = f["errors"][:-1]
        mu_val = f["params/mu"][()]
        alpha = f["params/alpha"][()]

    norm = alpha * 2.0

    for u in range(len(coeffs_raw)):
        coeffs_raw[u] *= norm * (-1)**u / math.factorial(u)
        errors_raw[u] *= abs(norm) * (+1)**u / math.factorial(u)

    print(mu_val)
    print(coeffs_raw)
    print(errors_raw)

    rel_err = np.abs(errors_raw) / np.abs(coeffs_raw)

    N_eff = len(coeffs_raw)
    while N_eff > 0 and rel_err[N_eff-1] >= err_coeff:
        N_eff -= 1

    print(N_eff)

    coeffs = coeffs_raw[:N_eff]
    errors = errors_raw[:N_eff]

    grid = orders_grid(N_eff)

    # Prepare figure: 4×3 grid
    fig, axes = plt.subplots(4, 3, figsize=(12, 13))
    plt.subplots_adjust(hspace=0.55, wspace=0.4)

    # Bare partial sum plot (top-left)
    ax_bare = axes[0, 0]
    orders = np.arange(N_eff)
    partial_sum = np.cumsum(coeffs * U**orders)
    partial_err = np.sqrt(np.cumsum((errors * U**orders) ** 2))

    bare_sum = np.sum(coeffs[:]*U**orders[:])
    bare_sum_err = np.sqrt(np.sum((errors[:]*U**orders[:]) ** 2))

    ax_bare.errorbar(
        orders, partial_sum, yerr=partial_err, fmt="o-", capsize=4,
        color="C3", label=rf"Bare $\sum_u c_u U^u$, $N_\mathrm{{eff}} = {N_eff}$"
    )
    #ax_bare.set_xlabel("Order $u$")
    ax_bare.set_ylabel(r"$S^{(\mathrm{bare})}(U{=}8)$")
    ax_bare.set_title("Bare series partial sums", fontsize=10)
    ax_bare.grid(True)
    ax_bare.legend(fontsize=8)

#   Relative error plot (top-right)
    ax_relerr = axes[0, 2]
    with np.errstate(divide='ignore', invalid='ignore'):
        rel_part_err = np.abs(partial_err / partial_sum)
        rel_part_err[~np.isfinite(rel_part_err)] = np.nan  # clean up infinities or NaNs

    ax_relerr.plot(orders, rel_part_err, "o-", color="C3", label="Relative error")
    #ax_relerr.set_xlabel("Order $u$")
    ax_relerr.set_ylabel(r"$\sigma(S_u) / |S_u|$")
    ax_relerr.set_title("Relative error of bare sum", fontsize=10)
    ax_relerr.set_yscale("log")  # often better for relative errors
    ax_relerr.grid(True, which="both", ls="--", alpha=0.6)
    ax_relerr.legend(fontsize=8)

    # First row → summary in center
    summary_ax = axes[0, 1]
    pade_results = []

    # ---- Loop over 3×3 Padé orders ----
    for i in range(3):
        for j in range(3):
            m, n = grid[i][j]
            ax = axes[i + 1, j]

            coeff_samples = sample_coeffs(coeffs, errors, n_bootstrap, rng)

            mean, std, mu, nu = compute_pade_with_consistent_orders(coeffs, coeff_samples, U_vals, n, m, tol)

            if mu is None or np.isnan(mean).all():
                ax.set_visible(False)
                continue

            # Value and spread at U=8
            val_U = mean[-1]
            spread_U = std[-1]
            pade_results.append({
                "i": i,
                "j": j,
                "mu": mu,
                "nu": nu,
                "val_U": val_U,
                "spread_U": spread_U,
            })

            # Plot mean ± std
            ax.plot(U_vals, mean, color="C0")
            ax.fill_between(U_vals, mean - std, mean + std, color="C0", alpha=0.3)
            #ax.set_ylim(0, 1)
            ax.set_title(f"[{mu},{nu}]  @U8={val_U:.4f}  spread={spread_U:.3g}", fontsize=9)
            ax.grid(True)

            # Add to summary
            summary_ax.plot(U_vals, mean, label=f"[{mu},{nu}]")
            summary_ax.fill_between(U_vals, mean - std, mean + std, alpha=0.15)

    # ---- Summary subplot ----
    #summary_ax.set_ylim(0, 1)
    #summary_ax.set_xlabel(r"$U$")
    summary_ax.set_ylabel(r"$S(U)$")
    summary_ax.set_title(f"μ = {mu_val:.3f} | Summary of Padé Approximants", fontsize=12)
    summary_ax.legend(fontsize=8)
    summary_ax.grid(True)

    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    # ---- Global figure title ----
    plt.suptitle(
        f"Bootstrap Padé Approximants | μ = {mu_val:.3f} | ",
        fontsize=14, y=0.93
    )
    if (show_plots):
        plt.show()

    print("\nAvailable Padé approximants:")
    for k, r in enumerate(pade_results):
        print(
            f"{k}: [{r['mu']},{r['nu']}] "
            f"U=8 → {r['val_U']:.6f} ± {r['spread_U']:.3g}"
        )

    print(len(pade_results))

    fU_values = [pade_results[i]["val_U"] for i in range(len(pade_results))]
    fU_spreads = [pade_results[i]["spread_U"] for i in range(len(pade_results))]

    # ---- Compute overall average ----
    fU_values = np.array(fU_values)
    mean_val = np.nanmean(fU_values)
    combined_err_spread = np.sqrt(np.nanmean(np.array(fU_spreads) ** 2))
    combined_err_std = np.nanstd(fU_values)
    combined_err = max(combined_err_std,combined_err_spread)

    print(f"⟨S(U=8)⟩ = {mean_val:.4f} ± {combined_err:.4f}")

    return mu_val, mean_val, combined_err, bare_sum, bare_sum_err