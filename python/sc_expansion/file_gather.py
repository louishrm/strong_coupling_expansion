import h5py
import glob
import re
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('order', type=int, help='Expansion order')
parser.add_argument('U', type=float, help='Interaction strength U')
parser.add_argument('beta', type=float, help='Inverse temperature beta')
parser.add_argument('bipartite', type=int, help='Whether the lattice is bipartite (1 for True, 0 for False)')
parser.add_argument('alpha', type=float, help='Alpha parameter for mcmc')
parser.add_argument('dual', type=int, help='Whether to use dual representation (1 for True, 0 for False)')

args = parser.parse_args()

order = args.order
U = args.U
beta = args.beta
bipartite = args.bipartite
alpha = args.alpha
dual = args.dual

# --- Find all sub-files for this order/U/beta, any delta ---
all_files = glob.glob(
    "results/full_lattice_data_order_{}_U_{:.6f}_beta_{:.6f}_delta_*_mu_*.h5".format(order, U, beta)
)

if not all_files:
    print("No subfiles found. Exiting.")
    exit()

# Extract unique delta values from filenames
delta_re = re.compile(r"_delta_([+-]?\d+\.\d+)_mu_")
deltas = sorted(set(
    float(delta_re.search(f).group(1))
    for f in all_files
    if delta_re.search(f)
))

print(f"Found {len(all_files)} subfiles across {len(deltas)} delta value(s): {deltas}")


def update_dataset(master, name, data_list):
    new_data = np.array(data_list)

    if name in master:
        dset = master[name]

        if dset.chunks is None:
            print(f"  Dataset '{name}' is not resizable. Migrating to chunked format...")
            old_data = dset[()]
            del master[name]
            master.create_dataset(name, data=old_data, maxshape=(None,), chunks=True)
            dset = master[name]

        curr_size = dset.shape[0]
        dset.resize((curr_size + len(new_data),))
        dset[curr_size:] = new_data
        print(f"  Appended {len(new_data)} items to '{name}'.")
    else:
        master.create_dataset(name, data=new_data, maxshape=(None,), chunks=True)
        print(f"  Created '{name}' with {len(new_data)} items.")


# --- Process each delta independently ---
for delta in deltas:
    delta_files = [f for f in all_files if f"_delta_{delta:.6f}_mu_" in f]
    delta_files.sort()

    print(f"\nProcessing delta={delta:.6f}: {len(delta_files)} file(s)...")

    if dual:
        extension = f"density_order_{order}_scan_mu_U={U}_beta={beta}_alpha={alpha}_delta={delta}.h5"
    else:
        extension = f"Omega_order_{order}_scan_mu_U={U}_beta={beta}_alpha={alpha}_delta={delta}.h5"

    if not bipartite:
        filename = f"./results/atom_triangular_lattice_" + extension
    else:
        filename = f"./results/atom_square_lattice_" + extension

    os.makedirs(os.path.dirname(filename), exist_ok=True)

    new_mus = []
    new_means = []
    new_errors = []
    new_infinite_U_coeffs = []

    for fname in delta_files:
        with h5py.File(fname, "r") as f:
            new_mus.append(f["mu"][()])
            new_means.append(f["mean"][()])
            new_errors.append(f["error"][()])
            if "reference_integral" in f:
                new_infinite_U_coeffs.append(f["reference_integral"][()])

    with h5py.File(filename, "a") as master:
        update_dataset(master, "mu_list", new_mus)
        update_dataset(master, "mean_list", new_means)
        update_dataset(master, "error_list", new_errors)
        if new_infinite_U_coeffs:
            update_dataset(master, "reference_integral", new_infinite_U_coeffs)

    print(f"  Cleaning up {len(delta_files)} temporary file(s)...")
    for fname in delta_files:
        try:
            os.remove(fname)
        except OSError as e:
            print(f"  Error deleting {fname}: {e}")

    print(f"  -> {filename}")

print("\nAll deltas processed.")
