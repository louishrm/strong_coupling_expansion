import h5py
import glob
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('order', type=int, help='Expansion order')
parser.add_argument('U', type=float, help='Interaction strength U')
parser.add_argument('beta', type=float, help='Inverse temperature beta')
parser.add_argument('alpha', type=float, help='Alpha parameter for mcmc')
parser.add_argument('dual', type=int, help='Whether to use dual representation (1 for True, 0 for False)')

args = parser.parse_args()

order = args.order
U = args.U
beta = args.beta
alpha = args.alpha
dual = args.dual

if dual: 
    filename = f"./results/full_lattice_density_order_{order}_scan_mu_U={U}_beta={beta}_alpha={alpha}.h5"
else:
    filename = f"./results/full_lattice_Omega_order_{order}_scan_mu_U={U}_beta={beta}_alpha={alpha}.h5"

# Creates the folder path if it's missing
os.makedirs(os.path.dirname(filename), exist_ok=True)

# 1. Find all the tiny files
files = glob.glob("results/full_lattice_data_order_{}_U_{:.6f}_beta_{:.6f}_mu_*.h5".format(order, U, beta))
files.sort() 

if not files:
    print("No subfiles found to merge. Exiting.")
    exit()

# 2. Process and Merge
with h5py.File(filename, "a") as master:
    
    new_mus = []
    new_means = []
    new_errors = []
    new_infinite_U_coeffs = []

    print(f"Reading {len(files)} subfiles...")
    for fname in files:
        with h5py.File(fname, "r") as f:
            new_mus.append(f["mu"][()])
            new_means.append(f["mean"][()])
            new_errors.append(f["error"][()])
            new_infinite_U_coeffs.append(f["reference_integral"][()])

    def update_dataset(name, data_list):
        new_data = np.array(data_list)
        
        if name in master:
            dset = master[name]
            
            # --- AUTO-REPAIR LOGIC ---
            # Check if dataset is resizable. If dset.chunks is None, it's fixed-size.
            if dset.chunks is None:
                print(f"⚠️ Dataset '{name}' is not resizable. Migrating to chunked format...")
                old_data = dset[()] # Read existing data
                del master[name]    # Delete fixed-size dataset
                # Re-create with maxshape and chunks
                master.create_dataset(
                    name, 
                    data=old_data, 
                    maxshape=(None,), 
                    chunks=True
                )
                dset = master[name] # Re-bind reference
            # -------------------------

            curr_size = dset.shape[0]
            dset.resize((curr_size + len(new_data),))
            dset[curr_size:] = new_data
            print(f"✅ Appended {len(new_data)} items to {name}.")
        else:
            master.create_dataset(
                name, 
                data=new_data, 
                maxshape=(None,), 
                chunks=True
            )
            print(f"✨ Created new resizable dataset {name} with {len(new_data)} items.")

    update_dataset("mu_list", new_mus)
    update_dataset("mean_list", new_means)
    update_dataset("error_list", new_errors)
    update_dataset("reference_integral", new_infinite_U_coeffs)

    # 3. Cleanup
    print("\nCleaning up temporary files...")
    for fname in files:
        try:
            os.remove(fname)
        except OSError as e:
            print(f"Error deleting {fname}: {e}")
            
    print("Cleanup complete.")

print(f"\nSuccessfully updated master file: {filename}")